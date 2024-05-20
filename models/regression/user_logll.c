////////////
#include <stdio.h>
#include <math.h>


extern int ndim_data; 

//////////////////////////////
// NOTE: change logll_beta in mpi_batch/init if func_prototype is changed.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one);

// N_parm not explicitly declare here because we will use each of parms invery detail, thus for sure we know the number.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one)
{
   /* model description:
    * 
    * x = np.random.rand(1000)*10
    * noise = np.random.randn(1000)*2
    * y = x*3.5+4 +noise
    * 
    * i.e., y = ax+b+N(0,d) 
    * 
    */
    //
    const double PI = 3.141592653589793;
    double a;
    double b;
    double d;
    a = *ptr_one_chain;
    b = *(ptr_one_chain+1);
    d = *(ptr_one_chain+2);
    // 
    double sig_power;
    double logll;
    double logll_one;
    logll = 0;
    for (int i=0; i<nline_data; i++)
    {
        sig_power = pow(d, 2.0);
        logll_one = -0.5*log(2*PI) -0.5*log(sig_power) - pow((data_NlineNdim[i*ndim_data+0]*a+b-data_NlineNdim[i*ndim_data+1]), 2.0)/2.0/sig_power;
        logll = logll + logll_one;
    }
    //
    logll = logll*beta_one;
    //    
    return logll;
}





