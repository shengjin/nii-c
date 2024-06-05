////////////
#include <stdio.h>
#include <math.h>


extern int ndim_data; 

//////////////////////////////
// NOTE: change logll_beta in mpi_batch/init if func_prototype is changed.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one);

double value_calc(double a, double b);

// N_parm not explicitly declare here because we will use each of parms invery detail, thus for sure we know the number.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one)
{
   /* model description:
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
    logll = value*beta_one;
    //    
    return logll;
}




double value_calc(double a, double b)
{
    //
    double mx1 = 10.0;
    double my1 = 10.0;
    double r1 = 1.0;
    double gauss1;
    //
    double mx2 = -10.0;
    double my2 = -10.0;
    double r2 = 1.0;
    double gauss2;
    //
    gauss1 = exp( - ( pow((a-mx1),2.0)/2/(pow(r1,2.0)) + pow((b-my1),2.0)/2/(pow(r1,2.0)) ) );
    gauss2 = exp( - ( pow((a-mx2),2.0)/2/(pow(r2,2.0)) + pow((b-my2),2.0)/2/(pow(r2,2.0)) ) );
    //
    return gauss1 + gauss2;
}



