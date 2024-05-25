/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "alloc.h"

extern char *results_dir;

extern int Tune_Ladder;
//  NOTE  scale the tuning of Temperature
extern double scale_tune_ladder;
extern double zero_stretch;


// distribut running_Beta_Values from root
int mpi_distribute_running_Beta_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, double *running_Beta_Values);

// tune the ladder
int tune_ladder_use_3_neighb(int n_ranks, double *doubleInt_Nswaps, double *doubleInt_Nprops, double *running_Beta_Values);


/////////////////////////////////////
//int mpi_adjust_ladder(MPI_Status status, int my_rank, int n_ranks, int root_rank)
int mpi_adjust_ladder(MPI_Status status, int rootsent_tag, int my_rank, int n_ranks, int root_rank, double *doubleInt_Nswaps, double *doubleInt_Nprops, double *running_Beta_Values)
{
    int debug = 0;
    
    if (debug)
    {
        if( my_rank == root_rank )
        {
            for (int i=0; i<n_ranks-1; i++)
            {
                printf("%d, %lf.\n", i, doubleInt_Nswaps[i]/doubleInt_Nprops[i]);
            }
            // 
            for(int i_rank = 0; i_rank < n_ranks; i_rank++)
            {
                printf("the Beta of %d rank is %lf .\n", i_rank, running_Beta_Values[i_rank]);
            }
            //
        }
    } 
    //
    //
    // let root rank tune the ladder
    if( my_rank == root_rank )
    {
        tune_ladder_use_3_neighb(n_ranks, doubleInt_Nswaps, doubleInt_Nprops, running_Beta_Values);
    }
    //
    //
    // mpi distribute the ladder
    mpi_distribute_running_Beta_root(status, my_rank, n_ranks, root_rank, rootsent_tag, running_Beta_Values);
    //
    return 0;
}



int mpi_distribute_running_Beta_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, double *running_Beta_Values)
{
    int debug = 0;

    if( my_rank == root_rank )
    {
        // Send new values of running_Beta_Values to slaves
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                MPI_Send(running_Beta_Values, n_ranks, MPI_DOUBLE, i_rank, rootsent_tag, MPI_COMM_WORLD);
            }
        }
        //
    } 
    else  
    {
    // Work for slave_rank
        // ... receve from root_rank the running_Beta_Values
        MPI_Recv(running_Beta_Values, n_ranks, MPI_DOUBLE, root_rank, rootsent_tag, MPI_COMM_WORLD, &status);
        //
        // print debug
        if (debug)
        {
            for(int i = 0; i < n_ranks; i++)
            {
                printf("%lf ", *(running_Beta_Values+i));
            }
            printf("\n"); 
        }
        //
    }
    //
    return 0;
}


// tune the ladder using 3 neighbour parallel tempering chains only.
int tune_ladder_use_3_neighb(int n_ranks, double *doubleInt_Nswaps, double *doubleInt_Nprops, double *running_Beta_Values)
{
    //
    // print debug 
    int debug = 1;
    //
    ////////////////////////////////////////////////////////////////
    // temperature ladder: where T_i = 1/beta_i
    double * T_ladder;
    T_ladder = alloc_1d_double(n_ranks);
    //
    for (int i=0; i<n_ranks; i++)
    {
        T_ladder[i]  = 1/running_Beta_Values[i];
    }
    //
    ////////////////////////////////////////////////////////////////
    // swap accept rate between chain i and chain i+1
    double * swap_ar_i_ip1;
    swap_ar_i_ip1 = alloc_1d_double(n_ranks-1);
    for (int i=0; i<n_ranks-1; i++)
    {
        swap_ar_i_ip1[i] = doubleInt_Nswaps[i]/doubleInt_Nprops[i];
        if (debug)
        {
            printf("SWAP rate of i: %d,  %lf.\n", i, swap_ar_i_ip1[i]);
        }
    }
    //
    //
    //
    double devi_of_mean;
    //
    for (int i=1; i<n_ranks-1; i++)
    {
        // 
        if ( (swap_ar_i_ip1[i-1]+swap_ar_i_ip1[i]) > 0 ) 
        {
            // calc for chain i the difference between two neighbor swap acceptance rates
            // In case of a increasing beta array:
            // if devi_of_mean > 0:  then T_i and T_i_p1 too close, T_i_m1 and T_i too far, 
            //      then decrease beta, i.e., increase Ti  
            // if devi_of_mean < 0:  then T_i and T_i_m1 too close, T_i_p1 and T_i too far, 
            //      then increase beta, i.e., decrease Ti  
            devi_of_mean = (swap_ar_i_ip1[i]-swap_ar_i_ip1[i-1]) / (swap_ar_i_ip1[i]+swap_ar_i_ip1[i-1]);
            //
            if (debug)
            {
                printf("i: %d,  %lf %lf %lf.\n", i, swap_ar_i_ip1[i-1], swap_ar_i_ip1[i], devi_of_mean);
            }
            //
            ////////////////////////////////////////////////////////////////
            // Change T_ladder
            //
            printf("        Tune the ladder %d .\n", i);
            // if devi_of_mean > 0, increase Ti based on T_i_m1 (limit)
            if ( devi_of_mean > 0 )
            {
                double T_diff;
                double T_increase;
                T_diff = (1/running_Beta_Values[i-1]) - (1/running_Beta_Values[i+1]) ;  // > 0
                //
                if (devi_of_mean == 1)
                {
                    T_increase = T_diff * (1-zero_stretch);
                }
                else
                {
                    T_increase = T_diff * devi_of_mean * scale_tune_ladder;
                }
                //
                T_ladder[i] = T_ladder[i+1] + T_increase;
                //
                if (debug)
                {
                    printf("i: %d, Tdiff %lf devi_of_mean %lf  T_increase %lf T_ladder_new %lf T_old_im1 %lf T_old_i %lf.\n", i, T_diff, devi_of_mean, T_increase, T_ladder[i], (1/running_Beta_Values[i-1]),  (1/running_Beta_Values[i]) );
                }
                //
            }
            // if devi_of_mean <= 0, decrease Ti based on T_i_p1 (limit)
            else 
            {
                double T_diff;
                double T_decrease;
                T_diff = (1/running_Beta_Values[i-1]) - (1/running_Beta_Values[i+1]) ; // > 0
                //
                if (devi_of_mean == -1)
                {
                    T_decrease = T_diff *  (1-zero_stretch);
                }
                else
                {
                    T_decrease = T_diff * (-devi_of_mean) * scale_tune_ladder;
                }
                //
                T_ladder[i] = T_ladder[i-1] - T_decrease;
                //
                if (debug)
                {
                    printf("i: %d, Tdiff %lf devi_of_mean %lf  T_decrease %lf T_ladder_new %lf T_old_i %lf T_old_ip1 %lf.\n", i, T_diff, devi_of_mean, T_decrease, T_ladder[i], (1/running_Beta_Values[i]),  (1/running_Beta_Values[i+1]) );
                }
                //
            }
        }
        //
        ////////////////
        //  in case of two zero swap acceptance rates:
        else
        {
            double Total_T_multiple;
            Total_T_multiple = (1/running_Beta_Values[i-1]) / (1/running_Beta_Values[i+1]) ; 
            double each_time_T_multiple;
            each_time_T_multiple = pow(Total_T_multiple, (0.5));
            // reset the ladder
            T_ladder[i] = T_ladder[i-1] / each_time_T_multiple;
            printf("        Tune the ladder %d (double 0 swap_ARs). \n", i);
        }
        //
        // specifically tune T_ladder[1] and T_ladder[n_rank-2]
        double sub_boundary_critical;
        sub_boundary_critical = 0.01;
        if (swap_ar_i_ip1[0] < sub_boundary_critical)
        {
            T_ladder[1] = T_ladder[0] - (T_ladder[0] - T_ladder[2]) * zero_stretch;
        }
        //
        if (swap_ar_i_ip1[n_ranks-2] < sub_boundary_critical)
        {
            T_ladder[n_ranks-2] = T_ladder[n_ranks-3] - (T_ladder[n_ranks-3] - T_ladder[n_ranks-1]) * (1-zero_stretch);
        }
        //
    }
    //
    //
    ////////////////////////////////////////////////////////////////
    // Change Beta_Values
    //
    for (int i=1; i<n_ranks-1; i++)
    {
        running_Beta_Values[i] = (1/T_ladder[i]);
    }
    //
    //
    free_1d_double(swap_ar_i_ip1);
    swap_ar_i_ip1 = NULL;
    free_1d_double(T_ladder);
    T_ladder = NULL;
    //
    return 0;
}


