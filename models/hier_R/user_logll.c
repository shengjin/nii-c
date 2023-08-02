////////////
#include <stdio.h>
#include <math.h>

extern double *Beta_Values;

//////////////////////////////
// NOTE: change the declaration of logll_beta in mpi_batch/init if func_prototype is changed.


double logllhd_Vdevi_V(double Vdevi, double V, double lambda);

double logllhd_Vdevi_V_plusgamma(double Vdevi, double V, double lambda, double gamma_lambda_Y);

double logllhd_oneMeasure(double gamma_power, double logC, double Yobs, double Ytr, double Xobs, double Xtr, double lambda_Yobs, double lambda_Y, double lambda_Xobs, double gamma_lambda_Y);

double logll_one_planet(int N_measures, int ndim_data, double Xmu, double Ymu, int i_x, int i_y, double *ptr_df_pl_once, double *ptr_one_chain);

double logll_beta(double *ptr_one_chain, int nline_data, double **data_NlineNdim, int i_rank);

double calc_mean(double *numeric_array, int N_array, int N_dim, int i_dim);




////////////////////////////
// X,Y relation 
//
/*
// in case of F: Ftr = Ftr/FJupiter
double log_muY(double logC, double gamma_power, Xtr)
{
    double log_value;
    log_value = logC + gamma_power*log(Xtr);
    return log_value;
}
*/
//
// has been checked
// V is mean V, Vdevi is a measured V with a scattering of lambda
double logllhd_Vdevi_V(double Vdevi, double V, double lambda)
{
    double lambda_square;
    lambda_square = lambda*lambda;
    static const double half_log_2pi = 0.9189385332046727;
    double half_log_lambSQ;
    half_log_lambSQ = 0.5*log(lambda_square);
    double err;
    err = (Vdevi-V);
    double log_value;
    log_value = -half_log_2pi 
                -half_log_lambSQ 
                -err*err/2/lambda_square; 
    return log_value;
}
//
double logllhd_Vdevi_V_plusgamma(double Vdevi, double V, double lambda, double gamma_lambda_Y)
{
    double lambda_scaled;
    lambda_scaled = lambda*pow(V, gamma_lambda_Y);
    double lambda_square;
    lambda_square = lambda_scaled*lambda_scaled;
    static const double half_log_2pi = 0.9189385332046727;
    double half_log_lambSQ;
    half_log_lambSQ = 0.5*log(lambda_square);
    double err;
    err = (Vdevi-V);
    double log_value;
    log_value = -half_log_2pi 
                -half_log_lambSQ 
                -err*err/2/lambda_square; 
    return log_value;
}
//
// passing value has been checked.
double logllhd_oneMeasure(double gamma_power, double logC, double Yobs, double Ytr, double Xobs, double Xtr, double lambda_Yobs, double lambda_Y, double lambda_Xobs, double gamma_lambda_Y)
{
    //////////
    double ll_Yobs_Ytr;
    ll_Yobs_Ytr = logllhd_Vdevi_V(Yobs, Ytr, lambda_Yobs);
    //////////
    // Ymu (mean Y) at a specific X
    double Ymu;
    Ymu = exp(logC)*pow(Xtr, gamma_power);
    //////////
    // Xtr entering through Ymu
    double ll_Ytr_YmuXtr;
    ll_Ytr_YmuXtr = logllhd_Vdevi_V_plusgamma(Ytr, Ymu, lambda_Y, gamma_lambda_Y);
    //////////
    double ll_Xobs_Xtr;
    ll_Xobs_Xtr = logllhd_Vdevi_V(Xobs, Xtr, lambda_Xobs);
    //
    double logll;
    logll = ll_Yobs_Ytr + ll_Ytr_YmuXtr + ll_Xobs_Xtr;
    //
    return logll;
}
//
//
//
// it has been checked
double logll_one_planet(int N_measures, int ndim_data, double Xmu, double Ymu, int i_x, int i_y, double *ptr_df_pl_once, double *ptr_one_chain)
{
    double gamma_power;
    double logC;
    double lambda_Y;
    double lambda_Yobs;
    double lambda_Xobs;
    // doing
    double gamma_lambda_Y;
    //
    gamma_power = *ptr_one_chain;
    logC = *(ptr_one_chain+1);
    // to be checked exp right
    lambda_Y = exp(*(ptr_one_chain+2));
    lambda_Yobs = exp(*(ptr_one_chain+3));
    lambda_Xobs = exp(*(ptr_one_chain+4));
    gamma_lambda_Y = *(ptr_one_chain+5);
    //
    double Xtr;
    double Ytr;
    double Xobs;
    double Yobs;
    //
    double logll;
    //
    if (N_measures == 1)
    {
        Xtr = Xmu;
        Ytr = Ymu;
        Xobs = Xtr;
        Yobs = *(ptr_df_pl_once+i_y);  // tbd
        logll = logllhd_oneMeasure(gamma_power, logC, Yobs, Ytr, Xobs, Xtr, lambda_Yobs, lambda_Y, lambda_Xobs, gamma_lambda_Y);
        return logll;
    }
    else
    {
        logll = 0;
        double logll_once;
        for (int i=0; i<N_measures; i++)
        {
            Xtr = Xmu;
            Ytr = Ymu;
            Xobs = *(ptr_df_pl_once+i*ndim_data+i_x);  // tbd
            Yobs = *(ptr_df_pl_once+i*ndim_data+i_y);  // tbd
            logll_once = logllhd_oneMeasure(gamma_power, logC, Yobs, Ytr, Xobs, Xtr, lambda_Yobs, lambda_Y, lambda_Xobs, gamma_lambda_Y);
            logll = logll + logll_once;
        }
        return logll;
    }
}
//
//
//
// it has been checked
// N_parm not explicitly declare here because we will use each of parms in very detail, thus for sure we know the number.
double logll_beta(double *ptr_one_chain, int nline_data, double **data_NlineNdim, int i_rank)
{
    ////////////////////////////
    // planet index
    int i;
    //
    // number of measures
    int N_measures;
    //
    // the i_line of the actual data of the planet
    int dealing_line;
    //
    // local ndim_data, we donot pass it in for simplicity
    int ndim_data = 2;
    // index of the columns for x and y
    int i_x = 0;
    int i_y = 1;
    // mean values
    double Xmu;
    double Ymu;
    //
    // ptr to the starting data of a planet in a 2D array
    double *ptr_df_pl_once; 
    //
    ///////////////////////////////////////////////////////
    //
    double logll;
    double logll_one;
    logll = 0;
    //
    //////////////////////////
    // 0st planet
    i = 0;
    //
    N_measures = (int)data_NlineNdim[i][1];
    dealing_line = (int)data_NlineNdim[i][0];
    ptr_df_pl_once = &data_NlineNdim[dealing_line][0];
    Xmu = calc_mean(ptr_df_pl_once, N_measures, ndim_data, i_x);
    Ymu = calc_mean(ptr_df_pl_once, N_measures, ndim_data, i_y);
    //
    // logll of one planet
    logll_one = logll_one_planet(N_measures, ndim_data, Xmu, Ymu, i_x, i_y, ptr_df_pl_once, ptr_one_chain);
    logll = logll + logll_one;
    //
    ////////////////////////  
    // loop for other planets
    i = 1;
    //
    while (   ( ((int)data_NlineNdim[i][0]+(int)data_NlineNdim[i][1]) <= nline_data ) 
            &&            ( (int)data_NlineNdim[i][0] > (int)data_NlineNdim[i-1][0] )  )
    {
        //////////
        N_measures = (int)data_NlineNdim[i][1];
        dealing_line = (int)data_NlineNdim[i][0];
        ptr_df_pl_once = &data_NlineNdim[dealing_line][0];
        Xmu = calc_mean(ptr_df_pl_once, N_measures, ndim_data, i_x);
        Ymu = calc_mean(ptr_df_pl_once, N_measures, ndim_data, i_y);
        //
        // logll of one planet
        logll_one = logll_one_planet(N_measures, ndim_data, Xmu, Ymu, i_x, i_y, ptr_df_pl_once, ptr_one_chain);
        logll = logll + logll_one;
        //
        i++;
    }
    //
    // final tempering 
    logll = logll*Beta_Values[i_rank];
    //
    return logll;
}




////////////////////////
//
double calc_mean(double *numeric_array, int N_array, int N_dim, int i_dim)
{
    double sum_all = 0;
    for (int i=0; i<N_array; i++)
    {
        sum_all = sum_all + *(numeric_array+i*N_dim+i_dim);
    }
    return sum_all/(double)N_array;
}




