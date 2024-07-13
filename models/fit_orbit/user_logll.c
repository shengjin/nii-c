////////////
#include <stdio.h>
#include <math.h>
//#include <stdlib.h>
#include <string.h>

#include "alloc.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//   NOTE this user_logll.c only works for the Nii-C v 1.0.0.
//   We haven't update it to Nii-C v 1.1.0 yet.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


extern double *Beta_Values;

int get_nlines_of_file(char *path);
double r8_normal_ab ( double a, double b );
double r8_normal_01 ( );

double period_to_au(double period,double ms_Msun);
double calc_osi(double mp_Mearth, double ms_Msun, double a_AU, double d_pc);
double newton_solver(double EE, double e);
int func_as(double *time_con, int Nline_time, double ecc, double osi, double cos_inc, double OmegaO, double M0, double omega, double per, double **Ntime_radec);

double log_likelihood(int Nline_time, double **Ntime_radec, double **data_radec, double var_uk);

//////////////////////////////
// NOTE: change the declaration of logll_beta in mpi_batch/init if func_prototype is changed.
double logll_beta(double *ptr_one_chain, int nline_data, double **data_NlineNdim, int i_rank);



//
// N_parm not explicitly declare here because we will use each of parms in very detail, thus for sure we know the number.
double logll_beta(double *ptr_one_chain, int nline_data, double **data_NlineNdim, int i_rank)
{
    ////////////////////////////
    //
    /////////////////////////
    // system parameter
    double ms_Msun = 1.0;
    double d_pc = 3.0;
    //
    double a_AU1;
    double osi1;
    //
    double a_AU2;
    double osi2;
    //
    // data dim, ra and dec
    int radec_int2 = 2;
    //
    //
    // local ndim_data, we donot pass this parameter in, for simplicity
    int ndim_data = 2;
    // 
    double cos_inc1;
    double ecc1;
    double an_Omg1;
    double p_omg1 ;
    double M0_1 ;
    double pl_m1;
    //
    ////////////////////////////////////////////////////////
    // NOTE: here we set "var_uk^2" to "sigma_neta^2 + epsilon_x^2" 
    double var_uk;
    ////////////////////////////////////////////////////////
    //
    double period1;
    double cos_inc2;
    double ecc2;
    double an_Omg2;
    double p_omg2 ;
    double M0_2 ;
    double pl_m2;
    double period2;
    //
    cos_inc1 = *ptr_one_chain;
    ecc1     = *(ptr_one_chain+1);
    an_Omg1  = *(ptr_one_chain+2);
    p_omg1   = *(ptr_one_chain+3);
    M0_1     = *(ptr_one_chain+4);
    pl_m1    = *(ptr_one_chain+5);
    var_uk   = *(ptr_one_chain+6);
    period1  = *(ptr_one_chain+7);
    cos_inc2 = *(ptr_one_chain+8);
    ecc2     = *(ptr_one_chain+9);
    an_Omg2  = *(ptr_one_chain+10);
    p_omg2   = *(ptr_one_chain+11);
    M0_2     = *(ptr_one_chain+12);
    pl_m2    = *(ptr_one_chain+13);
    period2  = *(ptr_one_chain+14);
    //
    //
    //
    double logll;
    //
    logll = 0;
    // 
    //
    /////////////////////////////////////////////
    //
    int Nline_time;
    Nline_time = nline_data/2;
    //printf("Nline of time_seq: %d \n", Nline_time);
    //
    double ** Ntime_radec1;
    double ** Ntime_radec2;
    double ** Ntime_radec;
    double ** data_radec;
    Ntime_radec1 = alloc_2d_double(Nline_time, radec_int2);
    Ntime_radec2 = alloc_2d_double(Nline_time, radec_int2);
    Ntime_radec = alloc_2d_double(Nline_time, radec_int2);
    data_radec = alloc_2d_double(Nline_time, radec_int2);
    //
    //
    double * time_con;
    time_con = alloc_1d_double(Nline_time);
    for (int i = 0; i < Nline_time; i++)
    {
        time_con[i] = data_NlineNdim[i*ndim_data][0];
        data_radec[i][0] = data_NlineNdim[i*ndim_data][1];
        data_radec[i][1] = data_NlineNdim[i*ndim_data+1][1];
    }
    //
    //
    a_AU1 = period_to_au(period1, ms_Msun);
    osi1 = calc_osi(pl_m1, ms_Msun, a_AU1, d_pc);
    //
    a_AU2 = period_to_au(period2, ms_Msun);
    osi2 = calc_osi(pl_m2, ms_Msun, a_AU2, d_pc);
    //
    //
    func_as(time_con, Nline_time, ecc1, osi1, cos_inc1, an_Omg1, M0_1, p_omg1, period1, Ntime_radec1);
    func_as(time_con, Nline_time, ecc2, osi2, cos_inc2, an_Omg2, M0_2, p_omg2, period2, Ntime_radec2);
    //
    for (int i = 0; i < Nline_time; i++)
    {
        Ntime_radec[i][0] = Ntime_radec1[i][0] + Ntime_radec2[i][0];
        Ntime_radec[i][1] = Ntime_radec1[i][1] + Ntime_radec2[i][1];
    }
    //
    //
    // calc log likelihood
    logll = log_likelihood(Nline_time, Ntime_radec, data_radec, var_uk); 
    //
    // final tempering 
    logll = logll*Beta_Values[i_rank];
    //
    free_2d_double(Ntime_radec1);
    free_2d_double(Ntime_radec2);
    free_2d_double(Ntime_radec);
    free_2d_double(data_radec);
    Ntime_radec1=NULL;
    Ntime_radec2=NULL;
    Ntime_radec=NULL;
    data_radec=NULL;
    //
    free_1d_double(time_con);
    time_con=NULL;
    //
    return logll;
}




////////////////////
// 
double period_to_au(double period_days,double ms_Msun)
{
    double PI;
    double GG;
    double Msun;
    double AU2cm;
    PI = 3.14159265358;
    GG = 6.67259e-8;
    Msun = 1.9891e33;
    AU2cm = 1.4959787e13;
    //
    double a_AU;
    double period_second;
    period_second = period_days*24.0*3600.0;
    double a = period_second/2.0/PI;
    double b = GG*Msun*ms_Msun;
    double r = pow(a, 2.0)*b;
    a_AU = pow(r, 1.0/3.0)/AU2cm;
    return a_AU;
}


////////////////////////////////////////////////////////////
// python function: calc_osi_etc(mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega) 1
double calc_osi(double mp_Mearth, double ms_Msun, double a_AU, double d_pc)
{
    //osi_alpha = 3.0*(mp_Mearth)*(ms_Msun**-1.0)*(a_AU)*(d_pc**-1) # mu as
    //#osi_alpha = 5.0*(mplanet_JupiterMass)*(msun_SolarMass**-1.0)*(a_5AU)*(d_pc**-1) # mas
    double osi_alpha; 
    osi_alpha = 3.0*(mp_Mearth)*(pow(ms_Msun, -1.0))*(a_AU)*(pow(d_pc, -1.0));
    return osi_alpha;
}


////////////////////
double newton_solver(double EE, double e)
{
    double tolerance = 1e-8;
    double eps = 1;
    double M;
    M = EE;
    while(fabs(eps)>tolerance)
    {
        double E1;
        E1=EE-(EE-e*sin(EE)-M)/(1-e*cos(EE));
        eps=E1-EE;
        EE=E1;
    }
    //
    return EE; 
}


int func_as(double *time_con, int Nline_time, double ecc, double osi, double cosi, double OmegaO, double M0, double omega, double per, double **Ntime_radec)
{
//  # alpha, e*sin(omega), e*cos(omega), i, Omega, M0, period
    //
    double deg2rad = 0.0174532925;
    double PI = 3.14159265358;
    //
    omega = omega*deg2rad;
    M0 = M0*deg2rad;
    OmegaO = OmegaO*deg2rad;
    double coso;
    double sino;
    double cosOg;
    double sinOg;
    coso=cos(omega);
    sino=sin(omega);
    cosOg=cos(OmegaO);
    sinOg=sin(OmegaO);
    //cosi=cos_inc
    //
    double A;
    double B;
    double F;
    double G;
    A=osi*(coso*cosOg-sino*sinOg*cosi);
    B=osi*(coso*sinOg+sino*cosOg*cosi);
    F=osi*(-sino*cosOg-coso*sinOg*cosi);
    G=osi*(-sino*sinOg+coso*cosOg*cosi);
    //
    double X;
    double Y;
    double EE;
    for (int i = 0; i < Nline_time; i++)
    {
        EE = newton_solver(((2*PI)/per*time_con[i]-M0), ecc);
        X = cos(EE)-ecc;
        Y = sqrt(1-pow(ecc,2))*sin(EE);
        Ntime_radec[i][0] = (B*X+G*Y);
        Ntime_radec[i][1] = (A*X+F*Y);
        //
        //printf("%f %f %f\n", time_con[i], Ntime_radec[i][0], Ntime_radec[i][1]);
    }
    //
    return 0;
}


double log_likelihood(int Nline_time, double **Ntime_radec, double **data_radec, double var_uk)
{
    //
    //
    double PI2;
    PI2 = 6.28318530716;
    //
    double log_llhd = 0;
    //
    double ra_once;
    double dec_once;
    //
    double sig_power;
    sig_power = var_uk*var_uk;
    //
    double AC_twice;
    AC_twice = 2.0*log(pow(PI2, -0.5)) + 2.0*log(pow(sig_power, -0.5));
    //
    for (int i = 0; i < Nline_time; i++)
    {
        ra_once = -pow((data_radec[i][0]-Ntime_radec[i][0]), 2.0)/2.0/sig_power;
        dec_once = -pow((data_radec[i][1]-Ntime_radec[i][1]), 2.0)/2.0/sig_power;
        //printf("time[i]: %d %lf\n", i, log_llhd);
        log_llhd = log_llhd + ra_once + dec_once + AC_twice;
    }
    //
    return log_llhd;
}





