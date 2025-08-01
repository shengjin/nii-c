/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#include "alloc.h"

extern double *Beta_Values;

extern int Fout_Len;   // length of filename to save chains
extern char *results_dir; //dir to save outputs
extern char *FoutPre;
extern char *FoutSuf;


int save_the_seed(unsigned seed, char *path, int i_rank);

int save_init_parm(double **transit_BetaParm, char *path, int n_ranks, int N_parm);

// init logpost(posterior)
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one);
// init log_pirior
double log_prior(double *ptr_one_chain);

// set init random parameterS at N_ITER = 0 for one chain
void init_parm_set(int seed, double *chain_parm);

// save the first chain and its logpost
int save_first_chain(double *chain_Parm, char *path, int i_rank, double logpost_first, int N_parm);


int mpi_gen_init_parm(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, int N_parm, unsigned *ptr_rolling_seed, double **transit_BetaParm_root)
{    
    int debug = 1;
    //
    if( my_rank == root_rank )
    {
    // Work for root_rank
        // Length for paring messages
        int n_to_generate;
        n_to_generate = N_parm;
        //
        // Tell slave processes N_parm and rolling_seed
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                //...tell a slave process how many parameters in a sampling model
                MPI_Send(&n_to_generate, 1 , MPI_INT, i_rank, rootsent_tag, MPI_COMM_WORLD);
                //...tell a slave process the random seed
                (*ptr_rolling_seed)++;
                //
                MPI_Send(ptr_rolling_seed, 1 , MPI_UNSIGNED, i_rank, rootsent_tag+1, MPI_COMM_WORLD);
            }
        }
        //
        // Generate a random inital parm in the segment assigned to the root process 
        // ...roll the seed
        (*ptr_rolling_seed)++;
        if (debug)
        {
            save_the_seed(*ptr_rolling_seed, results_dir, my_rank);
        }
        // ...gen the parameters
        init_parm_set(*ptr_rolling_seed, &transit_BetaParm_root[root_rank][0]);
        //
        // Collect the random initial parmES generated from the slave processes. 
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                MPI_Recv(&transit_BetaParm_root[i_rank][0], n_to_generate, MPI_DOUBLE, i_rank, slavereturn_tag, MPI_COMM_WORLD, &status);
            }
        }
        //
        if (debug)
        {
            save_init_parm(transit_BetaParm_root, results_dir, n_ranks, N_parm);
        }
        //
    } 
    else  
    {
    // Work for slave_rank
    // NOTE: they do not have to know tranist_BP array, that why their transit_BP are free-ed before.
        //
        // Receive infos
    	int n_to_generate_slave = 0;
    	unsigned rolling_seed_slave = 0;
        MPI_Recv(&n_to_generate_slave, 1, MPI_INT, root_rank, rootsent_tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&rolling_seed_slave, 1, MPI_UNSIGNED, root_rank, rootsent_tag+1, MPI_COMM_WORLD, &status);
        //
        // check rank seed
        if (debug)
        {
            save_the_seed(rolling_seed_slave, results_dir, my_rank);
        }
        //
        // Alloc memory
        double *chain_Parm_slave;
        chain_Parm_slave = alloc_1d_double(n_to_generate_slave);
        //
        // Gen the initial random parm at iter = 0
        init_parm_set(rolling_seed_slave, chain_Parm_slave);
        //
        // Send root_rank the generatED parm
        MPI_Send(&chain_Parm_slave[0], n_to_generate_slave, MPI_DOUBLE, root_rank, slavereturn_tag, MPI_COMM_WORLD); 
        //
        // Free temporary memory alloc.ed
        free_1d_double(chain_Parm_slave);
        chain_Parm_slave = NULL;
    }
    //
    return 0;
}



int mpi_init_calc_logllpp(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int nline_data, double *data_NlineNdim, double N_parm, double *logpost_all_ranks)
{
    int debug_arr_pointer = 0;
    //
    if( my_rank == root_rank )
    {
    // Work for root_rank
        //
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                //...tell a slave the init_of_batch parm 
                MPI_Send(&transit_BetaParm_root[i_rank][0], N_parm, MPI_DOUBLE, i_rank, rootsent_tag, MPI_COMM_WORLD);
            }
        }
        //
        // alloc a array to aviod the logll and lop change the values by accident
        double *chain_Parm_root;
        chain_Parm_root = alloc_1d_double(N_parm);
        for (int i=0; i<N_parm; i++)
        {
            chain_Parm_root[i] = transit_BetaParm_root[root_rank][i]; 
        }
        //
        ////////////////////////////////////////////////////////////
        // redundant debug to help understanding array and pointer
        if (debug_arr_pointer)
        {
            printf("chain[0],transit_BP[0]: %lf %lf\n", chain_Parm_root[0] , transit_BetaParm_root[root_rank][0]); 
            printf("chain[0],transit_BP[0]: %p %p\n", &chain_Parm_root[0] , &transit_BetaParm_root[root_rank][0]); 
        }
        ////////////////////////////////////////////////////////////
        //
        double logll_tempered_root;
        double logprior_root;
        double logpost_root;
        logll_tempered_root = logll_beta(&chain_Parm_root[0], nline_data, data_NlineNdim, Beta_Values[root_rank]);
        logprior_root = log_prior(&chain_Parm_root[0]);
        logpost_root = logll_tempered_root + logprior_root;
        logpost_all_ranks[root_rank] = logpost_root;
        //
        // save first chain and its logpost
        save_first_chain(chain_Parm_root, results_dir, root_rank, logpost_root, N_parm);
	//
        // Collect the logpost of other ranks the slave processes. 
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                MPI_Recv(&logpost_all_ranks[i_rank], 1, MPI_DOUBLE, i_rank, slavereturn_tag, MPI_COMM_WORLD, &status);
            }
        }
        //
        // Free temporary memory alloc.ed
        free_1d_double(chain_Parm_root);
        chain_Parm_root = NULL;
    } 
    else  
    {
    // Work for slave_rank
    // NOTE: they do not have to know tranist_BP array, that why their transit_BP are free-ed before.
        // I must be a slave process, so I must receive my array segment. storing it in a "local" array.
        // Alloc memory
        double *chain_Parm_slave;
        chain_Parm_slave = alloc_1d_double(N_parm);
        // ...receive the one parm (in tranist_B.P. array) from root
        MPI_Recv(&chain_Parm_slave[0], N_parm, MPI_DOUBLE, root_rank, rootsent_tag, MPI_COMM_WORLD, &status);
        //
        double logll_tempered_slave;
        double logprior_slave;
        double logpost_slave;
        logll_tempered_slave = logll_beta(&chain_Parm_slave[0], nline_data, data_NlineNdim, Beta_Values[my_rank]);
        logprior_slave = log_prior(&chain_Parm_slave[0]);
        logpost_slave = logll_tempered_slave + logprior_slave;
	//
        // save first chain and its logpost
        save_first_chain(chain_Parm_slave, results_dir, my_rank, logpost_slave, N_parm);
        //
        // ... Send root_rank the logpost_slave
        MPI_Send(&logpost_slave, 1, MPI_DOUBLE, root_rank, slavereturn_tag, MPI_COMM_WORLD); 
        //
        // Free temporary memory alloc.ed
        free_1d_double(chain_Parm_slave);
        chain_Parm_slave = NULL;
    }
    //
    return 0;
}






int save_the_seed(unsigned seed, char *path, int i_rank)
{
    FILE *out;
    //
    // set fnames
    char fname[Fout_Len];
    snprintf(fname, sizeof fname, "%s%s%s", path, "/",  "init.randseed");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    // save chains
    fprintf(out, "seed for rank %d is %d\n",  i_rank, seed);
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}



int save_init_parm(double **transit_BetaParm, char *path, int n_ranks, int N_parm)
{
    FILE *out;
    //
    // set fnames
    char fname[Fout_Len];
    snprintf(fname, sizeof fname, "%s%s%s", path, "/",  "init.parm");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    // save chains
    for(int i = 0; i < n_ranks; i++)
    {
        fprintf(out, "init para set %d: ", i);
        for(int j = 0; j < N_parm; j++)
        {
            fprintf(out, "%f ", transit_BetaParm[i][j]);
        }
        fprintf(out, "\n");
    }
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}




int save_first_chain(double *chain_Parm, char *path, int i_rank, double logpost_first, int N_parm)
{
    FILE *out;
    //
    // set fnames
    char fname[Fout_Len];
    snprintf(fname, sizeof fname, "%s%s%s%d%s", path, "/",  FoutPre, i_rank, FoutSuf);
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    // save chains
    for (int i=0; i<N_parm; i++)
    {
        fprintf(out, "%.12e",  chain_Parm[i]);
        fprintf(out, "  ");
    }
    //
    fprintf(out, "%lf",  logpost_first);
    fprintf(out, "  ");
    //
    fprintf(out, "%d",  0);
    fprintf(out, "  ");
    //
    fprintf(out, "%d",  0);
    fprintf(out, "  ");
    //
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}





