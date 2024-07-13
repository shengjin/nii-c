////////////
#include <stdio.h>
#include <math.h>

//////////////////////////////
// NOTE: change logll_beta in mpi_batch/init if func_prototype is changed.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one);

double value_calc(double a, double b);

// N_parm not explicitly declare here because we will use each of parms invery detail, thus for sure we know the number.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one)
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
    logll = value*beta_one;
    //    
    return logll;
}



double value_calc(double a, double b)
{
   double ll;
   ll = exp(-pow(a,2) - pow((9+4*pow(a,2) + 9*b),2)) + 0.5 * exp(-8*pow(a,2)-8*pow((b-2),2));
   return log(ll);
   //return exp(-pow(a,2) - pow((9+4*pow(a,2) + 9*b),2)) + 0.5 * exp(-8*pow(a,2)-8*pow((b-2),2));
}



    
        
        
        
    
    






