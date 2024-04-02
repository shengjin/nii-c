////////////
#include <stdio.h>
#include <math.h>

extern double *Beta_Values;

//////////////////////////////
// NOTE: change logll_beta in mpi_batch/init if func_prototype is changed.
double logll_beta(double *ptr_one_chain, int nline_data, double **data_NlineNdim, int i_rank);

double value_calc(double a, double b);

// N_parm not explicitly declare here because we will use each of parms invery detail, thus for sure we know the number.
double logll_beta(double *ptr_one_chain, int nline_data, double **data_NlineNdim, int i_rank)
{
   /* model description:
    * 
    * 
    */
    //
    double a;
    double b;
    a = *ptr_one_chain;
    b = *(ptr_one_chain+1);
    // 
    double value;
    double logll;

    value = value_calc(a, b);
    logll = value*Beta_Values[i_rank];
    //    
    return logll;
}



double value_calc(double a, double b)
{
   return sin(b/2.0)*cos(a/2.0)+2;
}



