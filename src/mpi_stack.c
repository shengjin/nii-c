/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "alloc.h"

extern int N_parm;
extern int N_swap;
extern int Swapmode;
extern char *results_dir;

extern int Tune_Ladder;
extern int N_stopTuneLadder;

////////////////////////////////////
int mpi_adjust_ladder(MPI_Status status, int rootsent_tag, int my_rank, int n_ranks, int root_rank, double *doubleInt_Nswaps, double *doubleInt_Nprops, double *running_Beta_Values);


int max(int a, int b);
int min(int a, int b);

// uniform int between a and b
int i4_unif_ab(int a, int b);
// uniform int between 0 and a
int i4_unif_0a( int a );


// save the sequence number in a stack where swap judgement is made
// ...we use i_next rather than i_tmp since at the moment of saving the batch is done
int save_debug_stack_sequence(unsigned *ptr_i_accumul, int i_swap);
// save the value of do_swap
int save_debug_stack_doswap(int do_swap, int i_swap);
// save the Beta values
int save_debug_Beta_Values(double *running_Beta_Values, int n_ranks);

// init logpost(posterior)
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one);
// init log_pirior
double log_prior(double *ptr_one_chain);

// declaration of the mpi_batch function
int mpi_run_a_batch(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int n_iter_a_batch, unsigned i_save_begin, int nline_data, double *data_NlineNdim, double *ptr_sigma_prop, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, double *logpost_all_ranks, double *running_Beta_Values);

// judge if a swapping is needed and do it if so
int mpi_judge_and_swap(MPI_Status status, int my_rank, int root_rank, int slave_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int nline_data, double *data_NlineNdim, double *logpost_all_ranks, int i_swap, int j_swap, int Swapmode, int n_ranks, double *running_Beta_Values);

// swap the values of two chain
void swap_two_chains(double *ptr_ichain, double *ptr_jchain, int N_parm);


int mpi_static_sigma_stack(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int n_iter_a_stack, int n_iter_a_batch_base, int n_iter_a_batch_rand, unsigned i_save_begin, int nline_data, double *data_NlineNdim, double *ptr_sigma_prop, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, double *logpost_all_ranks, double *running_Beta_Values)
{    
    int save_swap_debug = 0;
    int save_Beta_debug = 1;
    //
    int i_tmp = 0; 
    int i_next = 0; 
    //
    int n_iter_a_batch_adjust;
    int n_iter_a_batch = 0;
    //
    //
    int slave_rank; // the slave_rank that needed to calc the: logpost_jchain_ibeta  
    slave_rank = 1 - min(1, root_rank);  // can be 0 or 1 depends on root_rank
    //
    //
    int i_swap = 0; // lower i_rank of pair of ranks in mpi_swapping 
    int j_swap = 0; // lower i_rank of pair of ranks in mpi_swapping 
    int do_swap = 0; 
    // 
    //
    // adjust Ladders if T_L == 1
    double * doubleInt_Nswaps;
    doubleInt_Nswaps = alloc_1d_double(n_ranks);
    double * doubleInt_Nprops;
    doubleInt_Nprops = alloc_1d_double(n_ranks);
    if ((Tune_Ladder > 0) && (my_rank == root_rank))
    {
        // alloc the array to sum swaps
        for(int i = 0; i < n_ranks; i++)
        {
            doubleInt_Nswaps[i] = 0;
            doubleInt_Nprops[i] = 0;
        }
    }
    //
    //
    // save Beta_Values for debug
    if (save_Beta_debug)
    {
        if (my_rank == root_rank)
        {
            save_debug_Beta_Values(running_Beta_Values, n_ranks);
        }
    }
    //
    //
    while ( i_next < n_iter_a_stack )
    {
        /////////////////////////////////////////////
        //
        // Calc n_iter_a_batch, i_next
        if (my_rank == root_rank)
        {
            //
            n_iter_a_batch_adjust = i4_unif_ab(-n_iter_a_batch_rand, n_iter_a_batch_rand);
            n_iter_a_batch = n_iter_a_batch_base + n_iter_a_batch_adjust;
            //
            // new i_next
            i_next = i_tmp + n_iter_a_batch;
            if ( i_next > n_iter_a_stack )
            {
                n_iter_a_batch = n_iter_a_stack - i_tmp;
                i_next = n_iter_a_stack;
            }
            //
            // Tell slaves the i_next and n_iter_a_batch
            for(int i_rank = 0; i_rank < n_ranks; i_rank++)
            {
                if (i_rank != root_rank)
                {
                    //...tell a slave process i_next
                    MPI_Send(&i_next, 1 , MPI_INT, i_rank, rootsent_tag, MPI_COMM_WORLD);
                    //...tell a slave process how many iterations in a sampling model
                    MPI_Send(&n_iter_a_batch, 1 , MPI_INT, i_rank, rootsent_tag+1, MPI_COMM_WORLD);
                }
            }
        }
        else
        {
            // Receive the i_next
            MPI_Recv(&i_next, 1, MPI_INT, root_rank, rootsent_tag, MPI_COMM_WORLD, &status);
            // Receive the n_iter_a_batch
            MPI_Recv(&n_iter_a_batch, 1, MPI_INT, root_rank, rootsent_tag+1, MPI_COMM_WORLD, &status);
        }
        //
        // 
        // sync
        MPI_Barrier(MPI_COMM_WORLD);
        //
        //
        ////////////////////////////////////////////////////
        //
        // Run a batch
        mpi_run_a_batch(status, my_rank, n_ranks, root_rank, rootsent_tag, slavereturn_tag, transit_BetaParm_root, n_iter_a_batch, i_save_begin, nline_data, data_NlineNdim, ptr_sigma_prop, ptr_i_accumul, ptr_i_accumul_accept, logpost_all_ranks, running_Beta_Values);
        //
        // sync since we need to pass logpost to the ranks in mpi_batch
        MPI_Barrier(MPI_COMM_WORLD);
        //
        //
        ///////////////////////////////////////////////
        //
        // Swap proposal. NOTE: only root has transit_BP and logpost_all_ranks
        //
        for (int i=0; i<N_swap; i++)
        {
            // let the root find i_swap
            if (my_rank == root_rank)
            {
                i_swap = i4_unif_0a(n_ranks-1);
            }
            //
            // save infos for debug
            if (save_swap_debug) 
            {
                if (my_rank == root_rank)
                {
                    // swap proposal only, should be uniform 
                    save_debug_stack_sequence(ptr_i_accumul, i_swap);
                }
            }
            //
            // Judge and make swap in case needed
            do_swap = mpi_judge_and_swap(status, my_rank, root_rank, slave_rank, rootsent_tag, slavereturn_tag, transit_BetaParm_root, nline_data, data_NlineNdim, logpost_all_ranks, i_swap, j_swap, Swapmode, n_ranks, running_Beta_Values);
            //
            //
            // adjust Ladders if T_L == 1
            if (Tune_Ladder > 0)
            {
                if (my_rank == root_rank)
                {
                    doubleInt_Nprops[i_swap] = doubleInt_Nprops[i_swap] + 1;
                    doubleInt_Nswaps[i_swap] = doubleInt_Nswaps[i_swap] + (double)(do_swap);
                }
            }
            //
        }
        //
        //
        // sync since we need to pass logpost to the ranks in mpi_batch
        MPI_Barrier(MPI_COMM_WORLD);
        //
        //
        // Update i_tmp
        i_tmp = i_next;
    }
    //
    //
    //  adust the ladder (temperature) if Tune_Ladder == 1
    if ( ( (int)*ptr_i_accumul < N_stopTuneLadder) && ( Tune_Ladder > 0) )
    {
        {
            mpi_adjust_ladder(status, rootsent_tag, my_rank, n_ranks, root_rank, doubleInt_Nswaps, doubleInt_Nprops, running_Beta_Values);
        }
    }
    //
    //
    free_1d_double(doubleInt_Nswaps);
    doubleInt_Nswaps = NULL;
    free_1d_double(doubleInt_Nprops);
    doubleInt_Nprops = NULL;
    //
    return 0;
}




int mpi_judge_and_swap(MPI_Status status, int my_rank, int root_rank, int slave_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int nline_data, double *data_NlineNdim, double *logpost_all_ranks, int i_swap, int j_swap, int Swapmode, int n_ranks, double *running_Beta_Values)
{    
    int do_swap = 0; // default not jump
    // 
    if( my_rank == root_rank )
    {
        // Work for root_rank
        //
        // save infos of swapping
        int save_swap_debug = 0;
        //
        double H;
        double rand_unif;
        //
        if (Swapmode == 0)
        {
            j_swap = i_swap +1;
        }
        if (Swapmode == 1)
        {
            j_swap = i4_unif_0a(n_ranks-1);
            if (j_swap >= i_swap)
            {
                   j_swap = j_swap + 1;
            }
        }
        //
        // redundant ptr to make code clear
        double *ptr_ichain;
        double *ptr_jchain;
        ptr_ichain = &transit_BetaParm_root[i_swap][0];
        ptr_jchain = &transit_BetaParm_root[j_swap][0];
        //
        double logpost_ichain_ibeta; // logpost of chain(i) using beta(i)
        double logpost_ichain_jbeta; // logpost of chain(i) using beta(i)
        double logpost_jchain_ibeta; // logpost of chain(i) using beta(i)
        double logpost_jchain_jbeta; // logpost of chain(i) using beta(i)
        // 1st and 2nd values
        logpost_ichain_ibeta = logpost_all_ranks[i_swap];
        logpost_jchain_jbeta = logpost_all_ranks[j_swap];
        //
        //
        ///////////////////
        // Tell the slave_rank the chain, and the Beta needed to calc: 
        //                logpost_jchain_ibeta  
        // ...tell the chain
        MPI_Send(ptr_jchain, N_parm, MPI_DOUBLE, slave_rank, rootsent_tag+1, MPI_COMM_WORLD);
        // ...tell the Beta
        MPI_Send(&i_swap, 1, MPI_INT, slave_rank, rootsent_tag+2, MPI_COMM_WORLD);
        //
        //
        // Calculate:
        //          logpost_ichain_jbeta  
        // use a redundant chain_Parm_ichain here to avoid unrealized modification of transit_BP in logll_beta.
        double *chain_Parm_ichain;
        chain_Parm_ichain = alloc_1d_double(N_parm);
        for (int i=0; i<N_parm; i++)
        {
            chain_Parm_ichain[i] = transit_BetaParm_root[i_swap][i];
        }
        //
        double logll_tempered_mix_root; 
        double logprior_root;
        // caution: <1> check i,j numbers; <2> check if logprior is needed.
        logll_tempered_mix_root = logll_beta(&chain_Parm_ichain[0], nline_data, data_NlineNdim, running_Beta_Values[j_swap]);
        logprior_root = log_prior(&chain_Parm_ichain[0]);
        // 3rd value
        logpost_ichain_jbeta = logll_tempered_mix_root + logprior_root;
        //
        // Recv the value of logpost_jchain_ibeta from slave_rank
        // 4th value
        MPI_Recv(&logpost_jchain_ibeta, 1, MPI_DOUBLE, slave_rank, slavereturn_tag, MPI_COMM_WORLD, &status); 
        //
        //
        // JUDGEMENT
        if ( (logpost_ichain_jbeta+logpost_jchain_ibeta-logpost_ichain_ibeta-logpost_jchain_jbeta) > 0)
        {
            do_swap = 1;
        }
        else
        {
            H = exp( logpost_ichain_jbeta+logpost_jchain_ibeta-logpost_ichain_ibeta-logpost_jchain_jbeta );
            // draw random uniform number
            rand_unif = drand48( );
            if (rand_unif < H)
            {
                do_swap = 1;
            }
        }

        //
        ////////////////////////////
        //
        if (do_swap)
        {
            // swap the i_swap and j_swap chain
            // NOTE: in swap_two_chains, i_accumul and i_accumul_accept DO NOT increase.
            swap_two_chains(ptr_ichain, ptr_jchain, N_parm);
            //
            // change logpost_all_ranks[i_swap] and [j_swap]
            logpost_all_ranks[i_swap] = logpost_jchain_ibeta;
            logpost_all_ranks[j_swap] = logpost_ichain_jbeta;
        }
        //
        if (save_swap_debug)
        {
            save_debug_stack_doswap(do_swap, i_swap);
        }
        //
        //
        // Free temporary memory alloc.ed
        free_1d_double(chain_Parm_ichain);
        chain_Parm_ichain = NULL;
        //
    }
    else if (my_rank == slave_rank)
    {
        // Work for the slave_rank:
        //                logpost_jchain_ibeta  
        // Receive the  the chain
        double *chain_Parm_jchain;
        chain_Parm_jchain = alloc_1d_double(N_parm);
        MPI_Recv(&chain_Parm_jchain[0], N_parm, MPI_DOUBLE, root_rank, rootsent_tag+1, MPI_COMM_WORLD, &status);
        // Receive the iBeta (it is the i_swap sent from root_rank)
        MPI_Recv(&i_swap, 1, MPI_INT, root_rank, rootsent_tag+2, MPI_COMM_WORLD, &status);
        //            
        // Calc the logpost
        double logll_tempered_mix_slave; 
        double logprior_slave;
        double logpost_jchain_ibeta_slave; 
        // caution: <1> check i,j numbers; <2> check if logprior is needed.
        logll_tempered_mix_slave = logll_beta(&chain_Parm_jchain[0], nline_data, data_NlineNdim, running_Beta_Values[i_swap]);
        logprior_slave = log_prior(&chain_Parm_jchain[0]);
        logpost_jchain_ibeta_slave = logll_tempered_mix_slave + logprior_slave;
        //            
        // Send back the value of logpost_jchain_ibeta
        MPI_Send(&logpost_jchain_ibeta_slave, 1, MPI_DOUBLE, root_rank, slavereturn_tag, MPI_COMM_WORLD); 
        //
        // Free temporary memory alloc.ed
        free_1d_double(chain_Parm_jchain);
        chain_Parm_jchain = NULL;
        //
    }
    //
    return do_swap;
}


//////////////////////////
// Find maximum between two numbers.
int max(int a, int b)
{
    return ( a > b ) ? a : b;
}


//////////////////////////
// Find minimum between two numbers.
int min(int a, int b)
{
    return ( a > b ) ? b : a;
}



//////////////////////////
// Swap two parm chains
void swap_two_chains(double *ptr_ichain, double *ptr_jchain, int N_parm)
{
    double temp;
    //
    for (int i=0; i<N_parm; i++)
    {
        temp = *(ptr_ichain+i);
        *(ptr_ichain+i) = *(ptr_jchain+i);
        *(ptr_jchain+i) = temp;
    }
}



int save_debug_stack_sequence(unsigned *ptr_i_accumul, int i_swap)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s", results_dir, "/", "swap_sequence.dat");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    fprintf(out, "%u", *ptr_i_accumul);
    fprintf(out, "   ");
    fprintf(out, "%d", i_swap);
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}
//
//
//
int save_debug_stack_doswap(int do_swap, int i_swap)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s", results_dir, "/", "swap_decision.dat");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    fprintf(out, "%d", do_swap);
    fprintf(out, "   ");
    fprintf(out, "%d", i_swap);
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}


int save_debug_Beta_Values(double *Beta_Values, int n_ranks)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s", results_dir, "/", "Beta_Values_in_stack.dat");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    //
    for(int i = 0; i < n_ranks; i++)                                                                                             
    {                                                                                                                            
        fprintf(out, "%lf ", Beta_Values[i]);                                                                         
    }                                                                                                                        
    fprintf(out, "\n");                                                                                                      
    //
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}



