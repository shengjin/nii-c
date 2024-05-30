////////////
#include <stdio.h>
#include <math.h>


//////////////////////////////
// NOTE: change logll_beta in mpi_batch/init if func_prototype is changed.
double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one);

double value_calc(double a, double b, double c, double d, double f);

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
    double c;
    double d;
    double f;
    a = *ptr_one_chain;
    b = *(ptr_one_chain+1);
    c = *(ptr_one_chain+2);
    d = *(ptr_one_chain+3);
    f = *(ptr_one_chain+4);
    // 
//    printf("a b d beta %lf %lf %lf %lf \n", a, b, d, beta_one);
    double value;
    double logll;


    value = value_calc(a, b, c, d, f);
    logll = value*beta_one;
    //printf("Beta of one rank is %lf\n", beta_one);
    //    
    return logll;
}



double value_calc(double a, double b, double c, double d, double e)
{
   const double PI = 3.141592653589793;
   //
   double ll;
   ll = pow(a,2.0)-10*cos(2*PI*a) + pow(b,2.0)-10*cos(2*PI*b) + pow(c,2.0)-10*cos(2*PI*c) + pow(d,2.0)-10*cos(2*PI*d) + pow(e,2.0)-10*cos(2*PI*e)+50 ;
   return log(ll);
}



    
        
        

