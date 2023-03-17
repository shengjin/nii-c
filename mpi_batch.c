/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "alloc.h"



extern int N_parm; // N_para of model parameters in each chain.


extern int Fout_Len;   // length of filename to save chains
extern char *FoutPre; //prefix of output chains_files
extern char *FoutSuf; //suffix of output chains_files
extern char *results_dir; //dir to save outputs


// ...iter the batch, main work done here
double iter_batch_mh(double **chain_IterParm, double *ptr_sigma_prop, int n_iter_a_batch, int nline_data, double **data_NlineNdim, int i_rank, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, unsigned i_save_begin, double logpost_old);
//
// ...save the batch: save the chain sellected
int save_the_batch(double **chain_IterParm, int n_iter_a_batch, char *path, int i_rank, double *logpost, int *accumul, int *accumul_accept);
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
//
// returns a unit pseudonormal R8 useing the Box Muller method.
double r8_normal_01();


// log likelyhood
double logll_beta(double *ptr_one_chain, int nline_data, double **data_NlineNdim, int i_rank);
// log prior
double log_prior(double *ptr_one_chain);



int mpi_run_a_batch(MPI_Status status, int my_rank, int n_ranks, int root_rank, int rootsent_tag, int slavereturn_tag, double **transit_BetaParm_root, int n_iter_a_batch, unsigned i_save_begin, int nline_data, double **data_NlineNdim, double *ptr_sigma_prop, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, double *logpost_all_ranks)
{    
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
                MPI_Send(&transit_BetaParm_root[i_rank][0], N_parm, MPI_DOUBLE, i_rank, rootsent_tag+2, MPI_COMM_WORLD);
                //...tell a slave the logpost
                MPI_Send(&logpost_all_ranks[i_rank], 1, MPI_DOUBLE, i_rank, rootsent_tag+3, MPI_COMM_WORLD);
            }
        }
        //
        // Run the batch for root
        double **chain_IterParm_root;
        chain_IterParm_root = alloc_2d_double(n_iter_a_batch, N_parm);
        // ...copy the first parm_set from tranist_B.P. array
        for(int i = 0; i < N_parm; i++)
	{
            chain_IterParm_root[0][i] = transit_BetaParm_root[root_rank][i];
        }
        // ...copy the logpost_old for root
        double logpost_old_root;
        logpost_old_root = logpost_all_ranks[root_rank];
        //
        // Iterate array : chain, n_iter_a_batch, beta
        //           chain_IterParm filled, logpost_all_ranks updated         
        logpost_all_ranks[root_rank] = iter_batch_mh(chain_IterParm_root, ptr_sigma_prop, n_iter_a_batch, nline_data, data_NlineNdim, root_rank, ptr_i_accumul, ptr_i_accumul_accept, i_save_begin, logpost_old_root);
        // ...copy the final parm_set to tranist_B.P. array
        for(int i = 0; i < N_parm; i++)
	{
 	    transit_BetaParm_root[root_rank][i] = chain_IterParm_root[n_iter_a_batch-1][i];
	}
	//
        for(int i_rank = 0; i_rank < n_ranks; i_rank++)
        {
            if (i_rank != root_rank)
            {
                // Collect the final iteration parmES from the slave processes. 
                MPI_Recv(&transit_BetaParm_root[i_rank][0], N_parm, MPI_DOUBLE, i_rank, slavereturn_tag, MPI_COMM_WORLD, &status);
                // Collect logpost_all_ranks[i_rank]s from slaves
                MPI_Recv(&logpost_all_ranks[i_rank], 1, MPI_DOUBLE, i_rank, slavereturn_tag+1, MPI_COMM_WORLD, &status);
            }
        }
        //
        // Free temporary memory alloc.ed
        free_2d_double(chain_IterParm_root);
        chain_IterParm_root = NULL;
        //
    } 
    else  
    {
    // Work for slave_rank
        // I must be a slave process, so I must receive my array segment. storing it in a "local" array.
        //
        // Receive the first parm_set (in tranist_B.P. array) from root
        double **chain_IterParm_slave;
        chain_IterParm_slave = alloc_2d_double(n_iter_a_batch, N_parm);
        MPI_Recv(&chain_IterParm_slave[0][0], N_parm, MPI_DOUBLE, root_rank, rootsent_tag+2, MPI_COMM_WORLD, &status);
        // Receive the logpost_old
        double logpost_old_slave;
        MPI_Recv(&logpost_old_slave, 1, MPI_DOUBLE, root_rank, rootsent_tag+3, MPI_COMM_WORLD, &status);
        //
        // Iter the chains
        double logpost_final_slave;
        logpost_final_slave = iter_batch_mh(chain_IterParm_slave, ptr_sigma_prop,  n_iter_a_batch, nline_data, data_NlineNdim, my_rank, ptr_i_accumul, ptr_i_accumul_accept, i_save_begin, logpost_old_slave);
        //
        // Send root_rank the final parm set
        MPI_Send(&chain_IterParm_slave[n_iter_a_batch-1][0], N_parm, MPI_DOUBLE, root_rank, slavereturn_tag, MPI_COMM_WORLD); 
        MPI_Send(&logpost_final_slave, 1, MPI_DOUBLE, root_rank, slavereturn_tag+1, MPI_COMM_WORLD); 
        //
        //
        // Free temporary memory alloc.ed
        free_2d_double(chain_IterParm_slave);
        chain_IterParm_slave = NULL;
    }
    //
    return 0;
}






double iter_batch_mh(double **chain_IterParm, double *ptr_sigma_prop, int n_iter_a_batch, int nline_data, double **data_NlineNdim, int i_rank, unsigned *ptr_i_accumul, unsigned *ptr_i_accumul_accept, unsigned i_save_begin, double logpost_old)
{
    //
    int *accumul_accept;
    int *accumul;
    accumul_accept = alloc_1d_int(n_iter_a_batch);
    accumul = alloc_1d_int(n_iter_a_batch);
    //
    double *logpost;
    logpost = alloc_1d_double(n_iter_a_batch);
    //
    // the logpost of the final chain to pass out
    double logpost_final;
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
    while  (i<n_iter_a_batch)
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
        logll_tempered_new = logll_beta(ptr_one_chain_new, nline_data, data_NlineNdim, i_rank);
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
            (*ptr_i_accumul_accept)++;
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
                (*ptr_i_accumul_accept)++;
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
        ///////////////////////////////
        // ...save all the proposaed chains (including abandaned)
        if (save_allch_ll)
        {
            save_log_posterior(ptr_one_chain_new, logll_tempered_new, logprior_new, results_dir, i_rank);
        }
        /////////////////////////////
        //
        ////// update i_*s ////
        (*ptr_i_accumul)++;
        *(accumul+i) = *ptr_i_accumul;
        //
        *(accumul_accept+i) = *ptr_i_accumul_accept;
        //
        i++;
        //
    }
    //
    //////// save chain stuff or not ////
    if (*ptr_i_accumul >= i_save_begin)
    {
        save_the_batch(chain_IterParm, n_iter_a_batch, results_dir, i_rank, logpost, accumul, accumul_accept);
    }
    //
    ///////////////////
    //
    // pass out logpost final
    logpost_final = logpost[n_iter_a_batch-1];
    //
    //////////////////
    free_1d_double(ptr_one_chain_new);
    ptr_one_chain_new = NULL;
    //
    free_1d_double(logpost);
    logpost = NULL;
    //
    free_1d_int(accumul_accept);
    accumul_accept = NULL;
    //
    free_1d_int(accumul);
    accumul_accept = NULL;
    //
    //////////////////
    return logpost_final;
}




double do_gaussian_propose(double *ptr_one_chain, double *ptr_one_chain_new, double *ptr_sigma_prop)
{
/*
    Gaussian proposal distribution.
    This proposal is nearly as simple as it gets, mathematically it is:
        q(x∗∣xi)=Normal(xi,σ2),
    that is, a Gaussian centered on the current position xi with variance given by a standard deviation parameter σ.
    //
    Input:
    ptr_one_chain: poiter to the part of Parameter array in 2d chain_iterparm.
    ptr_one_chain_new: poiter to a Parameter array of parms.
    ptr_sigma_prop: poiter to the array of the Standard deviation of Gaussian distribution.
    //
    Main Work: propose parameter set of in ptr_one_chain_new.
    //
    Return: 
    q_factor: ratio of proposal densities.
*/
    //
    double x;
    double dx;
    double sigma_prop_x;
    //
    double q_factor;
    //
    int debug = 0;
    //
    for (int i = 0; i < N_parm; i++)
    {
        x = *(ptr_one_chain+i);
        //
        sigma_prop_x = *(ptr_sigma_prop+i);
        dx = r8_normal_01()*sigma_prop_x;
        // 
        // prop new parm chain
        *(ptr_one_chain_new+i) = x + dx;
        //
    }
    //
    if (debug)
    {
        save_debug_gaussian_proposal(ptr_one_chain_new, ptr_sigma_prop);
    }
        //
    //
    //proposal ratio factor is 1 since jump is symmetric
    q_factor = 1;
    return q_factor;
}




int save_debug_gaussian_proposal(double *ptr_one_chain, double *ptr_sigma_prop)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s", results_dir, "/", "debug_gaussian_prop");
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
        fprintf(out, "%.12e", *(ptr_one_chain+i));
        fprintf(out, "   ");
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





int save_log_posterior(double *ptr_one_chain, double logll_tempered_new, double logprior_new, char *path, int i_rank)
{
    FILE *out;
    //
    // set fnames
    char fname[Fout_Len];
    snprintf(fname, sizeof fname, "%s%s%s%d%s%s", path, "/",  FoutPre, i_rank, FoutSuf, ".all.ll");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    // save the parm and ll whether the chain jumps or not
    for (int i=0; i<N_parm; i++)
    {
        fprintf(out, "%.12e", *(ptr_one_chain+i));
        fprintf(out, "  ");
    }
    fprintf(out, "%lf", logll_tempered_new);
    fprintf(out, "  ");
    // logprior_new inf debuging
    fprintf(out, "%lf", logprior_new);
    fprintf(out, "  ");
    fprintf(out, "%lf", logprior_new + logll_tempered_new);
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}




int save_the_batch(double **chain_IterParm, int n_iter_a_batch, char *path, int i_rank, double *logpost, int *accumul,  int *accumul_accept)
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
    for (int i=0; i<n_iter_a_batch; i++)
    {
        for (int j=0; j<N_parm; j++)
        {
            fprintf(out, "%.12e",  chain_IterParm[i][j]);
            fprintf(out, "  ");
        }
        //
        fprintf(out, "%lf",  logpost[i]);
        fprintf(out, "  ");
        //
        fprintf(out, "%d",  accumul[i]);
        fprintf(out, "  ");
        //
        fprintf(out, "%d",  accumul_accept[i]);
        fprintf(out, "  ");
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



