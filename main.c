/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

#include "alloc.h"

extern int N_iter; // N iteration in each chain
extern int N_parm; // N_para of model parameters in each chain.

extern char *results_dir; //dir to save outputs

extern char *Data_file;   // user data file
extern char *Delimiter;   // delimiter of user data file
extern int ndim_data;     // dimension of user data file
int nline_data;           // lines of user data file (For data_loader)

extern double init_gp_ratio;    // init_gp_sigma to allowed_prior_range
extern unsigned init_rand_seed; // init_gp_sigma to allowed_prior_range
extern unsigned i_save_begin;   // accutually the i burning 

extern int n_iter_a_stack;       // set total N_iteration in a stack of batches
extern int n_iter_a_batch_base;  // set base N_iteration of a batch 
extern int n_iter_a_batch_rand;  // add +/- N_rand iteration to a batch 

extern int n_iter_in_tune;       // N_iteration when tuning accept rate
extern double ar_ok_lower;       // 
extern double ar_ok_upper;       // 


// to read all the input.ini parameters
int read_input_ini(char *path);
// check and mkdir for outputs
int make_dir(char *path);
// set init gaussian proposal for the sampling
int init_gaussian_proposal(double *ptr_sigma_prop, double init_gp_ratio);
// to save the sigma of gaussian proposal for each rank
int save_sigma_gauss_prop(double *ptr_sigma_prop, int i_rank);

 
// get the line numbers of a file
int mpi_get_nlines(int line_number, int my_rank, char *path, int root_rank);
// to load user data file 
void mpi_data_loader(int my_rank, int root_rank, int nline_data, int ndim_data, double *data_NlineNdim, char *datafile, char *delimiter);

// init a BetaParm arry ...
int mpi_gen_init_parm(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, int N_parm, unsigned *ptr_rolling_seed, double **transit_BetaParm_root);
// calc the init logpost
int mpi_init_calc_logllpp(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int nline_data, double *data_NlineNdim, double N_parm, double *logpost_all_ranks);
// collect sigma_prop at root
int mpi_gather_sigma_prop_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double **sigma_RanksParm_root, double *ptr_sigma_prop);



// declaration of the mpi flow function
int mpi_entire_flow(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int n_iter_a_stack, int n_iter_a_batch_base, int n_iter_a_batch_rand, unsigned i_save_begin, int nline_data, double *data_NlineNdim, double *ptr_sigma_prop, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, double *logpost_all_ranks, int N_iter, int n_iter_in_tune, double ar_ok_lower, double ar_ok_upper, double **sigma_RanksParm_root);



int main(int argc, char *argv[])
{
    // Print debug
    int debug = 1;
    
    /////////////////////////////////////////////////////
    //
    // MPI_Status: a data struct that provides info_s on the received message
    MPI_Status status;
    //
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    //
    // Get the rank of the process
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //
    // Get the number of processES.
    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
    //
    // Set the root process, can be any rank
    int root_rank = 0;
    //
    // Tag for common massges
    int rootsent_tag = 100;
    int slavereturn_tag = 200;


    /////////////////////////////////////////////////////
    // Runtime monitering start point 
    double runtime_start = MPI_Wtime();
   
 
    /////////////////////////////////////////////////////
    // Read the input parameters (every process should do that to load global parameters).
    read_input_ini("input.ini"); 
    // Check_and_make results_dir.
    if( my_rank == root_rank ) { make_dir(results_dir); }
    //
    
    /////////////////////////////////////////////////////
    // Set a rolling random seed
    unsigned rolling_seed;
    // redundant pointer to make reference clear
    unsigned *ptr_rolling_seed;
    ptr_rolling_seed = &rolling_seed;
    if (init_rand_seed <= 0)
    {
        *ptr_rolling_seed = ((unsigned int) time(NULL));
    }
    else 
    {
        *ptr_rolling_seed = init_rand_seed;
    }

    /////////////////////////////////////////////
    //
    // get the value of nline_data
    nline_data = mpi_get_nlines(nline_data, my_rank, Data_file, root_rank);
    //
    // alloc the memory to store user data 
    double * data_NlineNdim;  
    data_NlineNdim = alloc_1d_double(nline_data*ndim_data);
    // read the data using root=0 and then broadcast to other ranks.
    mpi_data_loader(my_rank, root_rank, nline_data, ndim_data, data_NlineNdim, Data_file, Delimiter);


    /////////////////////////////////////////////////////
    // Alloc a N_beta*N_parm transit array for mpi_init and mpi_swap
    double ** transit_BetaParm_root;  
    transit_BetaParm_root = alloc_2d_double(n_ranks, N_parm);
    //
    /////////////////////////////////////////////////////
    // Generate a initial random N_beta*N_parm parm array
    mpi_gen_init_parm(status, my_rank, n_ranks, root_rank, rootsent_tag, slavereturn_tag, N_parm, ptr_rolling_seed, transit_BetaParm_root);
    //
    // alloc a 1d arry to init logpost and to transit it after each batch for all ranks
    double * logpost_all_ranks;  
    logpost_all_ranks = alloc_1d_double(n_ranks);
    // calc the initial logpost_all_ranks
    mpi_init_calc_logllpp(status, my_rank, n_ranks, root_rank, rootsent_tag, slavereturn_tag, transit_BetaParm_root, nline_data, data_NlineNdim, N_parm, logpost_all_ranks);
    // 
    // 2d matrix at root to monitor, change, and distribute ptr_sigma_prop of all ranks
    double ** sigma_RanksParm_root;  
    sigma_RanksParm_root = alloc_2d_double(n_ranks, N_parm);
    // NOTE only root needs transit_BP and logpost_all_ranks, slaves do not use them
    if (my_rank != root_rank) 
    {
        free_2d_double(transit_BetaParm_root);
        transit_BetaParm_root = NULL;
        free_1d_double(logpost_all_ranks);
        logpost_all_ranks = NULL;
        free_2d_double(sigma_RanksParm_root);
        sigma_RanksParm_root = NULL;
    }
    //
    ////////////////////////////////////////
    // SYNC since we need to collect and distribute logpost to all the ranks 
    MPI_Barrier(MPI_COMM_WORLD);
    //
    /////////////////////////////////////////////////////
    // alloc a sigma_prop array for EACH RANK to store the gaussian proposals
    double * sigma_gaussian_prop;  
    // redundant poiter to make the code clear
    double * ptr_sigma_prop;
    sigma_gaussian_prop = alloc_1d_double(N_parm);
    ptr_sigma_prop = &sigma_gaussian_prop[0];
    //
    //init_gp should behind mpi_gen_init_parm that reads parms max and min 
    init_gaussian_proposal(ptr_sigma_prop, init_gp_ratio);
    // save gaussian_proposal
    save_sigma_gauss_prop(ptr_sigma_prop, my_rank);
    //
    // Collect *ptr_sigma_propS at root
    mpi_gather_sigma_prop_root(status, my_rank, n_ranks, root_rank, slavereturn_tag, sigma_RanksParm_root, ptr_sigma_prop);
 
    /////////////////////////////////////////////////////
    // accumulate i for the entire Markov chain
    unsigned i_accumul = 0; // updated in iter_batch_mh 
    // accumulate N_accept for the entire Markov chain
    unsigned i_accumul_accept = 0; // updated in iter_batch_mh 
    //
    unsigned *ptr_i_accumul;
    unsigned *ptr_i_accumul_accept;
    ptr_i_accumul = &i_accumul;
    ptr_i_accumul_accept = &i_accumul_accept;


    // Run a flow of stacks of batches. Where we tune sigma_prop between stacks and swap chains between batches
    // 
    mpi_entire_flow(status, my_rank, n_ranks, root_rank, rootsent_tag, slavereturn_tag, transit_BetaParm_root, n_iter_a_stack,  n_iter_a_batch_base, n_iter_a_batch_rand, i_save_begin, nline_data, data_NlineNdim, ptr_sigma_prop, ptr_i_accumul, ptr_i_accumul_accept, logpost_all_ranks, N_iter, n_iter_in_tune, ar_ok_lower, ar_ok_upper, sigma_RanksParm_root);
    

    // SYNC 
    MPI_Barrier(MPI_COMM_WORLD);

    if (debug)
    {
        printf("Total: my_rank,i_accumul,i_accumul_accept: %d %d %d\n", my_rank, i_accumul, i_accumul_accept);
    }

    /////////////////////////////////////////////////////
    // Runtime monitering end point
    double runtime_end = MPI_Wtime();
    double time_spent = (double)(runtime_end - runtime_start); 
    //
    if (debug)
    {    
        printf("Time spent: %f in %d rank\n", time_spent, my_rank);
    }


    /////////////////////////////////////////////////////
    // Free transit N_beta*N_parm array
    if( my_rank == root_rank ) 
    {
        free_2d_double(transit_BetaParm_root);
        transit_BetaParm_root = NULL;
        free_2d_double(sigma_RanksParm_root);
        sigma_RanksParm_root = NULL;
    }
    //
    free_1d_double(data_NlineNdim);  
    data_NlineNdim = NULL;  
    //
    free_1d_double(sigma_gaussian_prop);  
    sigma_gaussian_prop = NULL;  


    /////////////////////////////////////////////////////
    MPI_Finalize();

    return 0;
}








