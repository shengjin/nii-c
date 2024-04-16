/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "alloc.h"




extern int N_parm; // N_para of model parameters in each chain.

extern double *Beta_Values;

extern double ar_ok_lower;
extern double ar_ok_upper;
extern double ar_best;
extern double ar_accept_diff;
//
extern double sigma_scale_half_ratio;
extern double *sigma_parm_min;
extern double *sigma_parm_max;
extern double sigma_jumpin_ratio;

extern double init_gp_ratio;
//
extern int Fout_Len;   // length of filename to save chains
extern char *FoutPre; //prefix of output chains_files
extern char *FoutSuf; //suffix of output chains_files
extern char *results_dir; //dir to save outputs

// log likelyhood
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, int i_rank);
// log prior
double log_prior(double *ptr_one_chain);


int gen_sigma_alltune(double **sigma_alltune_ParmNvaried, double *sigma_parm_min, double *sigma_parm_max, double sigma_jumpin_ratio, double **sigma_RanksParm_root, int rank_in_tune, int n_ranks, int N_parm); 

int check_bounceInside_sigma_boundary(double *scaled_arr, double *sigma_parm_min, double *sigma_parm_max, double sigma_jumpin_ratio, int n_ranks, int i_parm);

int decide_sigma_to_change(double **ar_ParmNvaried, int n_ranks, int *ptr_iparm_change, int *ptr_irank_change);

int calc_SD_allParm(double *std_Parm, double **ar_ParmNvaried, int n_ranks);

int argmin(double *array, int n_array);
int argmax(double *array, int n_array);

double calc_closest_ar(double **ar_ParmNvaried, int *ptr_iparm_closest, int *ptr_irank_closest, int n_ranks);
double calc_closest_ar_oneparm(double **ar_ParmNvaried, int *ptr_irank_closest_oneparm, int n_ranks, int iparm_stdmax);

int modify_sigma_prop_rankintune(int rank_in_tune, int iparm_change, int irank_change, double **sigma_RanksParm_root, double **sigma_alltune_ParmNvaried, double **ar_ParmNvaried, int n_ranks);

int race_all_parm(double *ptr_sigma_prop_rankintune, double **sigma_alltune_ParmNvaried, double **ar_ParmNvaried, int n_ranks);

int init_gaussian_proposal(double *ptr_sigma_prop, double init_gp_ratio);


// loop A
int mpi_tune_sigma_irank(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, double *logpost_all_ranks, int nline_data, double *data_NlineNdim, double *ptr_sigma_prop, int n_iter_in_tune, int rank_in_tune, double **sigma_RanksParm_root);
// save checking info
int save_tuning_sigma_ar(double **ar_ParmNvaried, double **sigma_alltune_ParmNvaried, char *path, int rank_in_tune, int n_ranks);

// create log grid scaled arr
int create_logspace_array(double base_value, double half_scale, int n_points, double *scaled_arr);

// collect sigma_prop at root
int mpi_gather_sigma_prop_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double **sigma_RanksParm_root, double *ptr_sigma_prop);
//
// distribut sigma_prop from root
int mpi_distribute_sigma_prop_root(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, double **sigma_RanksParm_root, double *ptr_sigma_prop);

// copy the sigma_prop from rank_in_tune and replicate it into a n_tune array of sigma_prop
int replicate_OneSigmaOrigin_Nvaried(double **sigma_RanksParm_root, double **sigma_tune1parm_NvariedParm, int rank_in_tune, int n_ranks); 
// tune one parm only 
int tune_oneparm_Nvaried(double **sigma_alltune_ParmNvaried, double **sigma_tune1parm_NvariedParm, int n_ranks, int i_parm); 

int mpi_find_ranks2tune(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double ar_ok_lower, double ar_ok_upper, double accept_rate_a_stack, int *tune_ranks); 

// to save the sigma of gaussian proposal for each rank
int save_sigma_gauss_prop(double *ptr_sigma_prop, int i_rank);


// loop B
int mpi_tune_sigma_iparmNrank(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double *chain_single_tune, double logpost_single_tune, int nline_data, double *data_NlineNdim, double *ptr_sigma_prop, int n_iter_in_tune, double **ar_ParmNvaried, int rank_in_tune, int j_parm);

// NOTE: in iter_batch_mh, diff rank use diff Beta that are decided by int i_rank
// NOTE: here in iter_batch_mh_tune, diff rank should use the same Beta that are decided by rank_in_tune
double iter_batch_mh_tune(double **chain_IterParm, double *ptr_sigma_prop, int n_iter_in_tune, int nline_data, double *data_NlineNdim, int rank_in_tune, unsigned *ptr_accumul_accept_tune, double logpost_old, int my_rank);

// ...save the batch_in tune: save the tune chain
int save_the_batch_tune(double **chain_IterParm, int n_iter_a_batch, char *path, int my_rank, int rank_in_tune, double *logpost);
//
// ...save the details of the all the sampled chain(including abandoned ones), logll, logpri.
int save_log_posterior(double *ptr_one_chain, double logll_tempered, double logprior,char *path, int i_rank);
//
// ...to check and bounce the parms in the allowed range
int para_boundary(double *ptr_one_chain_new);

// ...make a gaussian_pro
double do_gaussian_propose(double *ptr_one_chain, double *ptr_one_chain_new, double *ptr_sigma_prop);

// save debug in do_gaussian_propose
int save_debug_gaussian_proposal(double *ptr_one_chain_new, double *ptr_sigma_prop);

// returns a unit pseudonormal R8 useing the Box Muller method.
double r8_normal_01();




// tune 1 rank only
int mpi_tune_sigma_irank(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, double *logpost_all_ranks, int nline_data, double *data_NlineNdim, double *ptr_sigma_prop, int n_iter_in_tune, int rank_in_tune, double **sigma_RanksParm_root)
{
// LOOP A is a wrapper of this function (i.e., for all ranks needed to be tune 
    //
    // save debugging 
    int save_tuning_debug = 0;
    //
    //////////////////////////////////////////////////////
    // Every slave have these two varibles
    double logpost_single_tune = 0;
    double *chain_single_tune;
    chain_single_tune = alloc_1d_double(N_parm);
    //
    //////////////////////////////////////////////////////
    // Alloc a N_parm*Nranks matrix to store the ar of all the tuning runs at root rank only
    double **ar_ParmNvaried;
    ar_ParmNvaried = alloc_2d_double(N_parm, n_ranks);
    // Alloc a N_parm*Nranks matrix to store the varying sigma (each parm varied n_ranks times)
    //   it contains all values needed to be tuned. In loop B, cut a slice of it and distribute
    //   the slice to Nparmvaired parallel runs that have different variations of one parm.
    //   Then, after Parm(Nparm) batches of parallel runs, we would have tuned all the Parm and 
    //   calculated the ar_ParmNvaried that contains the ar corresponding to all tuned runs
    double **sigma_alltune_ParmNvaried;
    sigma_alltune_ParmNvaried = alloc_2d_double(N_parm, n_ranks);
    //
    if (my_rank != root_rank)
    {
        free_2d_double(ar_ParmNvaried);
        ar_ParmNvaried = NULL;
        free_2d_double(sigma_alltune_ParmNvaried);
        sigma_alltune_ParmNvaried = NULL;
    }
    //
    //
    //////////////////////////////////////////////////////
    // Calc the N_parm*Nranks matrix 
    if (my_rank == root_rank)
    {
        // been checked
        gen_sigma_alltune(sigma_alltune_ParmNvaried, sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, sigma_RanksParm_root, rank_in_tune, n_ranks, N_parm); 
    }
    //
    //////////////////////////////////////////////////////
    // Set starting values
    if (my_rank == root_rank)
    {
        /////////////
        // Copy the start point values at root rank 
        // logpost 
        logpost_single_tune = logpost_all_ranks[rank_in_tune];
        // parm chain
        for (int j=0; j<N_parm; j++)
        {
            chain_single_tune[j] =  transit_BetaParm_root[rank_in_tune][j];
        }
        //
    }
    // Broadcast start point values (same for N_parm*Nranks runs) to all ranks
    MPI_Bcast(&chain_single_tune[0], N_parm, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    MPI_Bcast(&logpost_single_tune, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //
    //////////////////////////////////////////////////////
    // LOOP B: for each parm, tune 1 parm for a group of n_ranks parallel runs
    for (int j_parm=0; j_parm<N_parm; j_parm++)
    {
        //////////////////////
        // Ntune of different parm chains, root rank only 
        //(sigma_tune1parm_NvariedParm);
        double **sigma_tune1parm_NvariedParm;
        // reverse the order of n_ranks and N_parm to easy the mpi msg passing 
        sigma_tune1parm_NvariedParm = alloc_2d_double(n_ranks, N_parm);
        if (my_rank != root_rank)
        {
            free_2d_double(sigma_tune1parm_NvariedParm);
            sigma_tune1parm_NvariedParm = NULL;
        }
        //
        // tune one parm
        if (my_rank == root_rank)
        {
        // been checked
            // Replicate a tune matrix sigma_tune1parm_NvariedParm(n_ranks x ptr_sigma_prop)
            replicate_OneSigmaOrigin_Nvaried(sigma_RanksParm_root, sigma_tune1parm_NvariedParm, rank_in_tune, n_ranks); 
            // Modify the tune matrix sigma_tune1parm_NvariedParm
            tune_oneparm_Nvaried(sigma_alltune_ParmNvaried, sigma_tune1parm_NvariedParm, n_ranks, j_parm); 
        }
        // NOTE: distributon the tuning sigma tmp matrix, not the MCMC_sigma_matrix(sigma_RanksParm_root)
        mpi_distribute_sigma_prop_root(status, my_rank, n_ranks, root_rank, rootsent_tag, sigma_tune1parm_NvariedParm, ptr_sigma_prop);
        MPI_Barrier(MPI_COMM_WORLD);
        //
        //////////////////////
        // ar[][] filling has been checked
        // MPI TUNE fill a slice of ar_all_ranks[:][one_parm]
        mpi_tune_sigma_iparmNrank(status, my_rank, n_ranks, root_rank, slavereturn_tag, chain_single_tune, logpost_single_tune, nline_data, data_NlineNdim, ptr_sigma_prop, n_iter_in_tune, ar_ParmNvaried, rank_in_tune, j_parm);
        MPI_Barrier(MPI_COMM_WORLD);
        //
        //////////////////////
        if (my_rank == root_rank)
        {
            free_2d_double(sigma_tune1parm_NvariedParm);
            sigma_tune1parm_NvariedParm = NULL;
        }
        // end of tune a parm LOOP B
    }
    //
    //////////////////////////////////////////////////////
    free_1d_double(chain_single_tune);
    chain_single_tune = NULL;
    //
    //////////////////////////////////////////////////////
    //
    if (my_rank == root_rank)
    {
        int iparm_change;
        int irank_change;
        int *ptr_iparm_change;
        int *ptr_irank_change;
        ptr_iparm_change = &iparm_change;
        ptr_irank_change = &irank_change;
        // choose the parm(s) to change for the rank_in_tune using ar_ParmNvaried
        decide_sigma_to_change(ar_ParmNvaried, n_ranks, ptr_iparm_change, ptr_irank_change);
        //
        // update new sigma_prop for rankintune
        modify_sigma_prop_rankintune(rank_in_tune, iparm_change, irank_change, sigma_RanksParm_root, sigma_alltune_ParmNvaried, ar_ParmNvaried, n_ranks);
        //
    }
    //
    //////////////////////////////////////////////////////
    //
    if (my_rank == root_rank)
    {
        if (save_tuning_debug)
        {
            save_tuning_sigma_ar(ar_ParmNvaried, sigma_alltune_ParmNvaried, results_dir, rank_in_tune, n_ranks);
        }
        // free the 2d array of all tuned runs
        free_2d_double(ar_ParmNvaried);
        ar_ParmNvaried = NULL;
        // 
        free_2d_double(sigma_alltune_ParmNvaried);
        sigma_alltune_ParmNvaried = NULL;
    }
    //
    MPI_Barrier(MPI_COMM_WORLD);
    // end of tune a rank (N_parm done) inside LOOP A
    //
    return 0;
}





// loop B
int mpi_tune_sigma_iparmNrank(MPI_Status status, int my_rank, int n_ranks, int root_rank, int slavereturn_tag, double *chain_single_tune, double logpost_single_tune, int nline_data, double *data_NlineNdim, double *ptr_sigma_prop, int n_iter_in_tune, double **ar_ParmNvaried, int rank_in_tune, int j_parm)
{    
    //
    if( my_rank == root_rank )
    {
    // Work for root_rank
        //
        // Run the batch for root
        double **chain_IterParm;
        chain_IterParm = alloc_2d_double(n_iter_in_tune, N_parm);
        // ...copy the first parm_set from chain_single_tune
        for(int i = 0; i < N_parm; i++)
	{
            chain_IterParm[0][i] = chain_single_tune[i];
        }
        // ...copy the logpost_old for root
        double logpost_old;
        logpost_old = logpost_single_tune;
        //
        // accumulate N_accept for the tuning part
        unsigned i_accumul_accept_tune; // updated in iter_batch_mh_tune
        i_accumul_accept_tune = 0; 
        unsigned *ptr_accumul_accept_tune;
        ptr_accumul_accept_tune = &i_accumul_accept_tune;
        //
        //
        // Iterate array : chain, n_iter_a_batch, beta
        //           chain_IterParm filled
        ar_ParmNvaried[j_parm][root_rank] = iter_batch_mh_tune(chain_IterParm, ptr_sigma_prop, n_iter_in_tune, nline_data, data_NlineNdim, rank_in_tune, ptr_accumul_accept_tune, logpost_old, my_rank);
	//
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                // Collect logpost_all_ranks[i_rank]s from slaves
                MPI_Recv(&ar_ParmNvaried[j_parm][i_rank], 1, MPI_DOUBLE, i_rank, slavereturn_tag, MPI_COMM_WORLD, &status);
            }
        }
        //
        // Free temporary memory alloc.ed
        free_2d_double(chain_IterParm);
        chain_IterParm = NULL;
        //
    } 
    else  
    {
    // Work for slave_rank
        //
        double **chain_IterParm;
        chain_IterParm = alloc_2d_double(n_iter_in_tune, N_parm);
        // ...copy the first parm_set from chain_single_tune
        for(int i = 0; i < N_parm; i++)
	{
            chain_IterParm[0][i] = chain_single_tune[i];
        }
        // ...copy the logpost_old for root
        double logpost_old;
        logpost_old = logpost_single_tune;
        //
        // accumulate N_accept for the tuning part
        unsigned i_accumul_accept_tune; // updated in iter_batch_mh_tune
        i_accumul_accept_tune = 0; 
        unsigned *ptr_accumul_accept_tune;
        ptr_accumul_accept_tune = &i_accumul_accept_tune;
        //
        double ar;
        //
        // Iter the chains
        ar = iter_batch_mh_tune(chain_IterParm, ptr_sigma_prop,  n_iter_in_tune, nline_data, data_NlineNdim, rank_in_tune, ptr_accumul_accept_tune, logpost_old, my_rank);
        //
        // Send root_rank the final parm set
        MPI_Send(&ar, 1, MPI_DOUBLE, root_rank, slavereturn_tag, MPI_COMM_WORLD); 
        //
        // Free temporary memory alloc.ed
        free_2d_double(chain_IterParm);
        chain_IterParm = NULL;
    }
    //
    return 0;
}



// Here we use rank_in_tune rather than i_rank because in tuning process, all ranks using the same Beta Value 
double iter_batch_mh_tune(double **chain_IterParm, double *ptr_sigma_prop, int n_iter_in_tune, int nline_data, double *data_NlineNdim, int rank_in_tune, unsigned *ptr_accumul_accept_tune, double logpost_old, int my_rank)
{
    double ar;
    //
    int save_debug = 0;

    //save_sigma_gauss_prop(ptr_sigma_prop, my_rank);
 
    double *logpost;
    logpost = alloc_1d_double(n_iter_in_tune);
    //
    double logll_tempered_new;
    double logprior_new;
    double logpost_new;
    //
    double q_factor;
    double H;
    double rand_unif;
    //
    // redundant poiterS to make the code clear
    double *ptr_one_chain_old;
    double *ptr_one_chain_new;
    double *ptr_one_chain_toset;
    // the old chain for i=0, here the parms are directly passed by mpi_ (transit_2d)
    ptr_one_chain_old = &chain_IterParm[0][0];
    // chain new is a seperate memory alloc, free_ed in the ends
    ptr_one_chain_new = alloc_1d_double(N_parm);
    //
    //debug int used to turn on int save_log_posterior
    int save_allch_ll = 0; 
    //
    // set the initial first logpost 
    *(logpost+0) = logpost_old;
    //
    int i = 0;
    //
    while  (i<n_iter_in_tune)
    { 
        // check if copy old chain to ptr_one_chain is needed
        int copy_old = 1;
        //
        ////////////////////////////
        if (i>0)
        {
            // set the old logpost
            logpost_old = *(logpost+i-1);
            // set the address of of ptr_one_chain to the (i-1)_th set of chains
            ptr_one_chain_old = &chain_IterParm[i-1][0];
        }
        ptr_one_chain_toset = &chain_IterParm[i][0];
        //
        //
        /////////////////////////
        // ...update parm chain using a gaussian proposal
        q_factor = do_gaussian_propose(ptr_one_chain_old, ptr_one_chain_new, ptr_sigma_prop);
	// ...check and bounce the parms in range
        para_boundary(ptr_one_chain_new);
        //
        /////////////////////////
        // ...calc logll, log_p, and the logpost for the proposed chain 
        logll_tempered_new = logll_beta(ptr_one_chain_new, nline_data, data_NlineNdim, rank_in_tune);
        logprior_new = log_prior(ptr_one_chain_new);
        logpost_new = logll_tempered_new + logprior_new;
        //
        /////////// MH-algorithm ///////////
        // NOTE: not compute Hastings ratio in the 1st_criterion to avoid overflow of exp function
        //
        //
        if ( (logpost_new - logpost_old) > log(1/q_factor) )
        {
            for (int j=0; j<N_parm; j++)
            {
                *(ptr_one_chain_toset + j) = *(ptr_one_chain_new+j); // 1st setting, total 3
            }
            //
            *(logpost+i) = logpost_new; // 1st setting, total 3
            (*ptr_accumul_accept_tune)++;
            copy_old = 0;
        }
        else
        {
            H = exp(logpost_new - logpost_old) * q_factor;
            rand_unif = drand48( );
            if (rand_unif < H)
            {
                for (int j=0; j<N_parm; j++)
                {
                    *(ptr_one_chain_toset + j) = *(ptr_one_chain_new+j); // 2nd setting, total 3
                }
                //
                *(logpost+i) = logpost_new; // 2nd setting, total 3
                (*ptr_accumul_accept_tune)++;
                copy_old = 0;
            }
        }
        //
        //// copy old values ot ptr_one_chain if no jump made, 0th direct passed thus no need
        if ( (copy_old) && (i>0) )
        {
            for (int j=0; j<N_parm; j++)
            {
                *(ptr_one_chain_toset + j) = chain_IterParm[i-1][j]; // 3rd setting, total 3
            }
            //
            // set i_th logpost to old value if not updated in the MH part
            *(logpost+i) = logpost_old;  // 3rd setting, total 3
        }
        // 
        //
        ///////////////////////////////
        // ...save all the proposaed chains (including abandaned)
        if (save_allch_ll)
        {
            // manually add 100 to avoid conflict of filenames
            save_log_posterior(ptr_one_chain_new, logll_tempered_new, logprior_new, results_dir, rank_in_tune+100);
        }
        /////////////////////////////
        //
        i++;
        //
    }
    //
    //////// save chain stuff or not ////
    if (save_debug)
    {
        save_the_batch_tune(chain_IterParm, n_iter_in_tune, results_dir, my_rank, rank_in_tune, logpost);
    }
    //
    ar = (double)(*ptr_accumul_accept_tune)/(double)(n_iter_in_tune);
    //
    //////////////////
    free_1d_double(ptr_one_chain_new);
    ptr_one_chain_new = NULL;
    //
    free_1d_double(logpost);
    logpost = NULL;
    //
    //////////////////
    return ar;
}





int replicate_OneSigmaOrigin_Nvaried(double **sigma_RanksParm_root, double **sigma_tune1parm_NvariedParm, int rank_in_tune, int n_ranks) 
{
    int N_tune = 0;
    N_tune = n_ranks;
    //
    // Copy the first line to other lines
    for(int i = 0; i < N_tune; i++)
    {
        for(int j = 0; j < N_parm; j++)
        {
            sigma_tune1parm_NvariedParm[i][j] = sigma_RanksParm_root[rank_in_tune][j];
        }
    }
    //
    return 0;
}



int tune_oneparm_Nvaried(double **sigma_alltune_ParmNvaried, double **sigma_tune1parm_NvariedParm, int n_ranks, int i_parm) 
{ 
    // change the Nvaried sigma by varing i_parm Parm only
    for(int i = 0; i < n_ranks; i++)
    {
        // NOTE sigma_alltune... and sigma_tune1parm... have opposite dimensions
        // been checked
        sigma_tune1parm_NvariedParm[i][i_parm] = sigma_alltune_ParmNvaried[i_parm][i];
    }
    //
    return 0;
}


int create_logspace_array(double base_value, double half_scale, int n_points, double *scaled_arr)
{ 
    //printf("INPUT %lf %lf %d \n", base_value, half_scale, n_points);
    double log_interval, log_min, log_max;
    log_min = log(base_value/half_scale);
    log_max = log(base_value*half_scale);
    log_interval = (log_max-log_min) / (double)(n_points-1);
   
    //printf("log min, max, interval, %lf %lf %lf\n", log_min, log_max, log_interval);
    //
    for(int i = 0; i < n_points; i++)
    {
        scaled_arr[i] = exp(log_min+log_interval*(double)i);
        //printf("%lf ", scaled_arr[i]);
    }
    //
    return 0;
}


int check_bounceInside_sigma_boundary(double *scaled_arr, double *sigma_parm_min, double *sigma_parm_max, double sigma_jumpin_ratio, int n_ranks, int i_parm)
{
    if (scaled_arr[0] < sigma_parm_min[i_parm])
    {
        //printf("lowerBOUND, %lf, %lf \n", scaled_arr[0], sigma_parm_min[i_parm]);
        double half_scale;
        half_scale = sqrt(sigma_jumpin_ratio);
        create_logspace_array(sigma_parm_min[i_parm]*half_scale, half_scale, n_ranks, scaled_arr);
    }
    //
    if (scaled_arr[n_ranks-1] > sigma_parm_max[i_parm])
    {
        //printf("upperBOUND, %lf, %lf \n", scaled_arr[n_ranks-1], sigma_parm_max[i_parm]);
        double half_scale;
        half_scale = sqrt(sigma_jumpin_ratio);
        create_logspace_array(sigma_parm_max[i_parm]/half_scale, half_scale, n_ranks, scaled_arr);
    }
    //
    return 0;
}



int gen_sigma_alltune(double **sigma_alltune_ParmNvaried, double *sigma_parm_min, double *sigma_parm_max, double sigma_jumpin_ratio, double **sigma_RanksParm_root, int rank_in_tune, int n_ranks, int N_parm) 
{
    double *scaled_arr;
    scaled_arr = alloc_1d_double(n_ranks);
    //
    for(int i=0; i<N_parm; i++)
    {
        create_logspace_array(sigma_RanksParm_root[rank_in_tune][i], sigma_scale_half_ratio, n_ranks, scaled_arr);
        //
        check_bounceInside_sigma_boundary(scaled_arr, sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, n_ranks, i);
        //
        for(int j=0; j<n_ranks; j++)
        {
            sigma_alltune_ParmNvaried[i][j] = scaled_arr[j];
        }
        //
    }
    //
    free_1d_double(scaled_arr);
    scaled_arr = NULL;
    //   
    return 0;
}


int modify_sigma_prop_rankintune(int rank_in_tune, int iparm_change, int irank_change, double **sigma_RanksParm_root, double **sigma_alltune_ParmNvaried, double **ar_ParmNvaried, int n_ranks)
{
    //
    if (iparm_change >= 0)
    {
        sigma_RanksParm_root[rank_in_tune][iparm_change] = sigma_alltune_ParmNvaried[iparm_change][irank_change];
    }
    else
    {
        /////////////
        // gen values
        //
        double *sigma_gaussian_prop_rankintune;
        sigma_gaussian_prop_rankintune = alloc_1d_double(N_parm);
        double * ptr_sigma_prop_rankintune;
        ptr_sigma_prop_rankintune = &sigma_gaussian_prop_rankintune[0];
        //
        if (irank_change >= 0)
        {
            init_gaussian_proposal(ptr_sigma_prop_rankintune, init_gp_ratio);
        }
        else
        {
            race_all_parm(ptr_sigma_prop_rankintune, sigma_alltune_ParmNvaried, ar_ParmNvaried, n_ranks);
        }
        //
        /////////////
        // copy values
        //
        for (int i=0; i<N_parm; i++)
        {
            sigma_RanksParm_root[rank_in_tune][i] = sigma_gaussian_prop_rankintune[i];
        }
        //
        /////////////
        //
        free_1d_double(sigma_gaussian_prop_rankintune);
        sigma_gaussian_prop_rankintune = NULL;
    }
    //
    return 0;
}



int race_all_parm(double *ptr_sigma_prop_rankintune, double **sigma_alltune_ParmNvaried, double **ar_ParmNvaried, int n_ranks)
{
    double *oneparm_arerr_Nvaried;
    oneparm_arerr_Nvaried = alloc_1d_double(n_ranks);
    //
    int j_min_oneparm;
    //
    for (int i=0; i<N_parm; i++)
    {
        // calc oneparm_arerr_Nvaried corresponding to parm i.
        for (int j=0; j<n_ranks; j++)
        {
            oneparm_arerr_Nvaried[j] = fabs(ar_ParmNvaried[i][j] - ar_best);
        }
        // find the closest j_Nvaried to use
        j_min_oneparm = argmin(oneparm_arerr_Nvaried, n_ranks);
        // modify the parm i
        *(ptr_sigma_prop_rankintune + i) = sigma_alltune_ParmNvaried[i][j_min_oneparm];
    }
    //
    free_1d_double(oneparm_arerr_Nvaried);
    oneparm_arerr_Nvaried = NULL;
    //
    return 0;
}


int decide_sigma_to_change(double **ar_ParmNvaried, int n_ranks, int *ptr_iparm_change, int *ptr_irank_change)
{
    double *std_Parm;
    std_Parm = alloc_1d_double(N_parm);
    //
    int iparm_stdmax;
    int irank_closest_oneparm;
    int *ptr_irank_closest_oneparm;
    double v_oneparm_closest;
    int iparm_change1;
    int irank_change1;
    //
    int iparm_closest;
    int irank_closest;
    int *ptr_iparm_closest;
    int *ptr_irank_closest;
    double v_2d_closest;
    int iparm_change2;
    int irank_change2;
    //
    double rand_unif;
    //
    calc_SD_allParm(std_Parm, ar_ParmNvaried, n_ranks);
    //
    /////////////
    //
    irank_closest_oneparm = 0;
    ptr_irank_closest_oneparm = &irank_closest_oneparm;
    //
    iparm_stdmax = argmax(std_Parm, N_parm);
    iparm_change1 = iparm_stdmax;
    //
    v_oneparm_closest = calc_closest_ar_oneparm(ar_ParmNvaried, ptr_irank_closest_oneparm, n_ranks, iparm_stdmax);
    irank_change1 = irank_closest_oneparm;
    //
    /////////////
    //
    iparm_closest = 0;
    irank_closest = 0;
    ptr_iparm_closest = &iparm_closest;
    ptr_irank_closest = &irank_closest;
    //
    v_2d_closest = calc_closest_ar(ar_ParmNvaried, ptr_iparm_closest, ptr_irank_closest, n_ranks);
    iparm_change2 = iparm_closest;
    irank_change2 = irank_closest;
    //
    ////////////////
    // judgement
    //
    if (v_oneparm_closest < ar_accept_diff)
    {
        if (v_2d_closest < ar_accept_diff)
        {
            rand_unif = drand48( );
            if (rand_unif > 0.5)
            {
                *ptr_iparm_change = iparm_change1;
                *ptr_irank_change = irank_change1;
            }
            else
            {
                *ptr_iparm_change = iparm_change2;
                *ptr_irank_change = irank_change2;
            }
        }
        else
        {
            *ptr_iparm_change = iparm_change1;
            *ptr_irank_change = irank_change1;
        }
    }
    else
    {
        if (v_2d_closest < ar_accept_diff)
        {
            *ptr_iparm_change = iparm_change2;
            *ptr_irank_change = irank_change2;
        }
        else
        {
            rand_unif = drand48( );
            if (rand_unif > 0.5)
            {
                *ptr_iparm_change = -1;
                *ptr_irank_change = 0;
            }
            else
            {
                *ptr_iparm_change = -1;
                *ptr_irank_change = -1;
            }
        }
    }
    //
    free_1d_double(std_Parm);
    std_Parm = NULL;
    //
    return 0;
}


double calc_closest_ar_oneparm(double **ar_ParmNvaried, int *ptr_irank_closest_oneparm, int n_ranks, int iparm_stdmax)
{
    double closest_value_oneparm;
    //
    closest_value_oneparm = fabs(ar_ParmNvaried[iparm_stdmax][0] - ar_best);
    //
    for (int j=0; j<n_ranks; j++)
    {
        if ( fabs(ar_ParmNvaried[iparm_stdmax][j] - ar_best) < closest_value_oneparm )
        {
            closest_value_oneparm = fabs(ar_ParmNvaried[iparm_stdmax][j] - ar_best);
            //
            *ptr_irank_closest_oneparm = j;
        }
    }
    //
    return closest_value_oneparm;
}



double calc_closest_ar(double **ar_ParmNvaried, int *ptr_iparm_closest, int *ptr_irank_closest, int n_ranks)
{
    double closest_value;
    //
    closest_value = fabs(ar_ParmNvaried[0][0] - ar_best);
    //
    for (int i=0; i<N_parm; i++)
    {
        for (int j=0; j<n_ranks; j++)
        {
            if ( fabs(ar_ParmNvaried[i][j] - ar_best) < closest_value )
            {
                closest_value = fabs(ar_ParmNvaried[i][j] - ar_best);
                //
                *ptr_iparm_closest = i;
                *ptr_irank_closest = j;
            }
        }
    }
    //
    return closest_value;
}



int argmin(double *array, int n_array)
{
    int i_min;
    double value_min;
    //
    i_min = 0;
    value_min = array[0];
    //
    for (int i=1; i<n_array; i++)
    {
        if (array[i] < value_min)
        {
            i_min = i;
            value_min = array[i];
        }
    }
    //
    return i_min;
}


int argmax(double *array, int n_array)
{
    int i_max;
    double value_max;
    //
    i_max = 0;
    value_max = array[0];
    //
    for (int i=1; i<n_array; i++)
    {
        if (array[i] > value_max)
        {
            i_max = i;
            value_max = array[i];
        }
    }
    //
    return i_max;
}




int calc_SD_allParm(double *std_Parm, double **ar_ParmNvaried, int n_ranks)
{
    //
    double sum_oneparm;
    double sum_variance;
    double mean_oneparm;
    //
    for (int i=0; i<N_parm; i++)
    {
        sum_oneparm = 0;
        for (int j=0; j<n_ranks; j++)
        {
            sum_oneparm += ar_ParmNvaried[i][j];
        }
        mean_oneparm = sum_oneparm/(double)n_ranks;
        //
        sum_variance = 0;
        for (int j=0; j<n_ranks; j++)
        {
            sum_variance += pow((ar_ParmNvaried[i][j]-mean_oneparm), 2);
        }
        //
        std_Parm[i] = sqrt(sum_variance/(double)n_ranks);
    }
    //
    return 0;
}



int save_the_batch_tune(double **chain_IterParm, int n_iter_a_batch, char *path, int my_rank, int rank_in_tune, double *logpost)
{
    FILE *out;
    //
    // set fnames
    char fname[Fout_Len];
    snprintf(fname, sizeof fname, "%s%s%s%d%s%d%s", path, "/tune.",  FoutPre, rank_in_tune, ".", my_rank, FoutSuf);
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
    for (int i=0; i<n_iter_a_batch; i++)
    {
        for (int j=0; j<N_parm; j++)
        {
            fprintf(out, "%.12e",  chain_IterParm[i][j]);
            fprintf(out, "  ");
        }
        //
        fprintf(out, "%lf",  logpost[i]);
        //
        fprintf(out, "\n");
    }
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}



int save_tuning_sigma_ar(double **ar_ParmNvaried, double **sigma_alltune_ParmNvaried, char *path, int rank_in_tune, int n_ranks)
{
    FILE *out;
    //
    // set fnames
    char fname[Fout_Len];
    snprintf(fname, sizeof fname, "%s%s%d%s", path, "/sigma_ar_intune.rank",  rank_in_tune, ".dat");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    fprintf(out, "################################################\n");
    fprintf(out, "ALLTUNE sigma variations:\n");
    for (int j=0; j<N_parm; j++)
    {
        for (int k=0; k<n_ranks; k++)
        {
            fprintf(out, "%lf ", sigma_alltune_ParmNvaried[j][k]);
        }
        fprintf(out, "\n");
    }
    //
    fprintf(out, "ARs of them:\n");
    for (int j=0; j<N_parm; j++)
    {
        for (int k=0; k<n_ranks; k++)
        {
            fprintf(out, "%lf ", ar_ParmNvaried[j][k]);
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




