/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "alloc.h"

extern char *results_dir;

extern int N_parm;


// save n_iter_a_stack, ar, rank_id, 
int save_ar_stack(int i_next_stack, int n_iter_a_stack, double accept_rate_a_stack, int my_rank);

// run a CONSTANT gaussian proposal sigma stack of many batches
int mpi_static_sigma_stack(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int n_iter_a_stack, int n_iter_a_batch_base, int n_iter_a_batch_rand, unsigned i_save_begin, int nline_data, double **data_NlineNdim, double *ptr_sigma_prop, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, double *logpost_all_ranks);

// working 
// loop A
int mpi_tune_sigma_irank(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, double *logpost_all_ranks, int nline_data, double **data_NlineNdim, double *ptr_sigma_prop, int n_iter_in_tune, int rank_in_tune, double **sigma_RanksParm_root);


// collect sigma_prop at root
int mpi_gather_sigma_prop_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double **sigma_RanksParm_root, double *ptr_sigma_prop);
//
// distribut sigma_prop from root
int mpi_distribute_sigma_prop_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, double **sigma_RanksParm_root, double *ptr_sigma_prop);

int mpi_find_ranks2tune(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double ar_ok_lower, double ar_ok_upper, double accept_rate_a_stack, int *tune_ranks); 

// to save the sigma of gaussian proposal for each rank
int save_sigma_gauss_prop(double *ptr_sigma_prop, int i_rank);


// NOTE: ptr_i_accumul, ptr_i_accumul_accept are separately updated in each rank (they have the same ptr_i_accumul, but different ptr_i_accumul_accept
int mpi_entire_flow(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int n_iter_a_stack, int n_iter_a_batch_base, int n_iter_a_batch_rand, unsigned i_save_begin, int nline_data, double **data_NlineNdim, double *ptr_sigma_prop, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, double *logpost_all_ranks, int N_iter, int n_iter_in_tune, double ar_ok_lower, double ar_ok_upper, double **sigma_RanksParm_root)
{    
    int save_debug = 1;
    int processing_debug = 1;
    //
    unsigned n_accept_old = 0;
    unsigned n_accept_now = 0;
    double accept_rate_a_stack = 0;

    ///////////////////////
    // variales needed in tune
    // for all ranks: 1 to tune, 0 not tune
    int *tune_ranks;
    tune_ranks = alloc_1d_int(n_ranks);
        
    int i_tmp_stack = 0; 
    int i_next_stack = 0; 
    //
    while ( i_next_stack < N_iter )
    {
        /////////////////////////////////////////////
        //
        // Calc i_next_stack
        i_next_stack = i_tmp_stack + n_iter_a_stack;
        if ( i_next_stack > N_iter )
        {
            n_iter_a_stack = N_iter - i_tmp_stack;
            i_next_stack = N_iter;
        }

        //distribute the sigma_RanksParm_root (can be modified in mpi tune)
        mpi_distribute_sigma_prop_root(status, my_rank, n_ranks, root_rank, rootsent_tag, sigma_RanksParm_root, ptr_sigma_prop);
        //
        MPI_Barrier(MPI_COMM_WORLD);

        //////////////////////////////////////
        // run a stack. 
        // 
        n_accept_old = *ptr_i_accumul_accept;
        //
        save_sigma_gauss_prop(ptr_sigma_prop, my_rank);
        //
        mpi_static_sigma_stack(status, my_rank, n_ranks, root_rank, rootsent_tag, slavereturn_tag, transit_BetaParm_root, n_iter_a_stack, n_iter_a_batch_base, n_iter_a_batch_rand, i_save_begin, nline_data, data_NlineNdim, ptr_sigma_prop, ptr_i_accumul, ptr_i_accumul_accept, logpost_all_ranks);
        // 
        if ( (processing_debug) && (my_rank==root_rank) )
        {
            printf("%9d out of total %d iterations have completed.\n", i_next_stack, N_iter);
        }
        //
        // calc accpet rate
        n_accept_now = *ptr_i_accumul_accept;
        accept_rate_a_stack = (double)(n_accept_now-n_accept_old)/(double)(n_iter_a_stack);
        //
        if (save_debug)
        {
            save_ar_stack(i_next_stack, n_iter_a_stack, accept_rate_a_stack, my_rank);
        }
        //
        MPI_Barrier(MPI_COMM_WORLD);
        //
        if (i_next_stack < N_iter)
        {
            //
            mpi_find_ranks2tune(status, my_rank, n_ranks, root_rank, slavereturn_tag, ar_ok_lower, ar_ok_upper, accept_rate_a_stack, tune_ranks); 
            //
            MPI_Barrier(MPI_COMM_WORLD);
            //
            //
            /////////////////////
            //  Tune sigma_prop for each i_rank needed, sigma... are modified through sigma_RanksParm_root
            //
            // LOOP A: for each i_rank
            for (int i_rank=0; i_rank<n_ranks; i_rank++)
            {
                if (tune_ranks[i_rank] == 1)
                {
                    // printf infos
                    if ( (processing_debug) && (my_rank==root_rank) )
                    {
                        printf("        tuning rank %d...\n", i_rank);
                    }
                    //
                    mpi_tune_sigma_irank(status, my_rank, n_ranks, root_rank, rootsent_tag, slavereturn_tag, transit_BetaParm_root, logpost_all_ranks, nline_data, data_NlineNdim, ptr_sigma_prop, n_iter_in_tune, i_rank, sigma_RanksParm_root);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            // end of tune LOOP A ( all n_rank needed to be tuned)
            }
        }
        // end of a stack
        //
        i_tmp_stack = i_next_stack; 
    // end of loop of a stack
    }
    //
    free_1d_int(tune_ranks);
    tune_ranks = NULL;
    //
    return 0;
}



int mpi_find_ranks2tune(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double ar_ok_lower, double ar_ok_upper, double accept_rate_a_stack, int *tune_ranks) 
{
    int debug = 0;

    if( my_rank == root_rank )
    {
        if ( (accept_rate_a_stack < ar_ok_upper) && (accept_rate_a_stack > ar_ok_lower) )
        {
            tune_ranks[root_rank] = 0;
        }
        else
        {
            tune_ranks[root_rank] = 1;
        }
        //
        if (debug)
        {
            printf("i_rank, accept_rate_a_stack, tune_ranks[irank]: %d %lf %d\n", my_rank, accept_rate_a_stack, tune_ranks[root_rank]); 
        }
        //
        // Collect the tune_ranks[i_ranks] from slaves
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                MPI_Recv(&tune_ranks[i_rank], 1, MPI_INT, i_rank, slavereturn_tag, MPI_COMM_WORLD, &status);
            }
        }
    } 
    else  
    {
    // Work for slave_rank
        int tune_irank_slave;
        //
        if ( (accept_rate_a_stack < ar_ok_upper) && (accept_rate_a_stack > ar_ok_lower) )
        {
            tune_irank_slave = 0;
        }
        else
        {
            tune_irank_slave = 1;
        }
        //
        if (debug)
        {
            printf("i_rank, accept_rate_a_stack, tune_irank_slave: %d %lf %d\n", my_rank, accept_rate_a_stack, tune_irank_slave); 
        }
        //
        // ... Send root_rank the tune_irank_slave
        MPI_Send(&tune_irank_slave, 1, MPI_INT, root_rank, slavereturn_tag, MPI_COMM_WORLD); 
    }
    //
    // Broadcast tune_ranks to all slaves
    MPI_Bcast(&tune_ranks[0], n_ranks, MPI_INT, root_rank, MPI_COMM_WORLD);
    //
    return 0;
}




int mpi_gather_sigma_prop_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double **sigma_RanksParm_root, double *ptr_sigma_prop) 
{
    int debug = 0;

    if( my_rank == root_rank )
    {
        // Copy *ptr_sigma_prop of root 
        for(int i = 0; i < N_parm; i++)
        {
            sigma_RanksParm_root[root_rank][i] = *(ptr_sigma_prop+i);
        }
        //
        // Collect the *ptr_sigma_prop of slaves
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                MPI_Recv(&sigma_RanksParm_root[i_rank][0], N_parm, MPI_DOUBLE, i_rank, slavereturn_tag, MPI_COMM_WORLD, &status);
            }
        }
        //
        // print debug
        if (debug)
        {
            for(int i_rank = 0; i_rank < n_ranks; i_rank++)
            {
                printf("i_rank: "); 
                for(int i = 0; i < N_parm; i++)
                {
                    printf("%lf ", sigma_RanksParm_root[i_rank][i]);
                }
                printf("\n"); 
            }
        }
        //
    } 
    else  
    {
    // Work for slave_rank
        // ... Send root_rank the *ptr_sigma_prop
        MPI_Send(ptr_sigma_prop, N_parm, MPI_DOUBLE, root_rank, slavereturn_tag, MPI_COMM_WORLD); 
    }
    //
    return 0;
}





int mpi_distribute_sigma_prop_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, double **sigma_RanksParm_root, double *ptr_sigma_prop) 
{
    int debug = 0;

    if( my_rank == root_rank )
    {
        // Copy *ptr_sigma_prop of root 
        for(int i = 0; i < N_parm; i++)
        {
            *(ptr_sigma_prop+i) = sigma_RanksParm_root[root_rank][i]; 
        }
        //
        // Send new values of sigma_prop to *ptr_sigma_prop of slaves
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                MPI_Send(&sigma_RanksParm_root[i_rank][0], N_parm, MPI_DOUBLE, i_rank, rootsent_tag, MPI_COMM_WORLD);
            }
        }
        //
    } 
    else  
    {
    // Work for slave_rank
        // ... Send root_rank the *ptr_sigma_prop
        MPI_Recv(ptr_sigma_prop, N_parm, MPI_DOUBLE, root_rank, rootsent_tag, MPI_COMM_WORLD, &status);
        //
        // print debug
        if (debug)
        {
            for(int i = 0; i < N_parm; i++)
            {
                printf("%lf ", *(ptr_sigma_prop+i));
            }
            printf("\n"); 
        }
        //
    }
    //
    return 0;
}




int save_ar_stack(int i_next_stack, int n_iter_a_stack, double accept_rate_a_stack, int my_rank)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s%d", results_dir, "/", "accept_rate_stacks.chain", my_rank);
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    fprintf(out, "%d", i_next_stack);
    fprintf(out, "   ");
    fprintf(out, "%d", n_iter_a_stack);
    fprintf(out, "   ");
    fprintf(out, "%lf", accept_rate_a_stack);
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}




int save_sigma_gauss_prop(double *ptr_sigma_prop, int i_rank)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s%d", results_dir, "/", "gaussian_prop.chain", i_rank);
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    for (int i = 0; i < N_parm; i++)
    {
        fprintf(out, "%lf", *(ptr_sigma_prop+i));
        fprintf(out, "   ");
    }
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}





