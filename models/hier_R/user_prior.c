///////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEST_SIYE 100

///////////////////////////////////////
double gamma_mean;
double gamma_sd;
double gamma_max;
double gamma_min;
double logC_min;
double logC_max;
//lambda_Y_min = e**log_lambda_Y_min
//lambda_Y_max = e**log_lambda_Y_max
double log_lambda_Y_min;
double log_lambda_Y_max;
//lambda_Yobs_min = e**log_lambda_Yobs_min
//lambda_Yobs_max = e**log_lambda_Yobs_max
double log_lambda_Yobs_min;
double log_lambda_Yobs_max;
//lambda_Xobs_min = e**log_lambda_Xobs_min
//lambda_Xobs_max = e**log_lambda_Xobs_max
double log_lambda_Xobs_min;
double log_lambda_Xobs_max;
///////
extern int N_parm;
///////////////////////////////////////
double *sigma_parm_min;
double *sigma_parm_max;
extern double sigma_scale_min;
extern double sigma_scale_max;
///////////////////////////////////////
extern char *results_dir; //dir to save outputs
///////////////////////////////////////
extern double init_gp_ratio;

//////////////////
// new added
double gamma_lambda_Y;
double gamma_lambda_Y_min;
double gamma_lambda_Y_max;

/////////////////////////////////////
////////////////////
//NOTE: change log_prior in mpi_batch/init if func_prototype is changed
double log_prior(double *ptr_one_chain);
////////////////////
//calc sigma scale boundary (function in user_prior.c)
int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max);
////////////////////
double r8_normal_ab ( double a, double b );
double r8_unif_ab ( double a, double b );
double r8_logunif_ab ( double a, double b );
double normal_pdf(double x, double mean, double sd);
///////////////////
void read_parm_range(char *path);
///////////////////
int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max);


////////////////////////////////////////////
////////////////////////////////////////////
///// read prior range
////////////////////////////////////////////
////////////////////////////////////////////

void read_parm_range(char *path)
{
    char *para_line; 
    char *para_name;
    //
    char *read_onepara(char *path, char *para_name);
    //
    char dummy[DEST_SIYE] = {0};
    //
    para_name = "gamma_mean";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &gamma_mean);
    //
    para_name = "gamma_sd";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &gamma_sd);
    //
    para_name = "gamma_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &gamma_max);
    //
    para_name = "gamma_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &gamma_min);
    //
    para_name = "logC_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &logC_min);
    //
    para_name = "logC_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &logC_max);
    //
    para_name = "gamma_lambda_Y_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &gamma_lambda_Y_min);
    //
    para_name = "gamma_lambda_Y_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &gamma_lambda_Y_max);
    //
    para_name = "log_lambda_Y_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &log_lambda_Y_min);
    //
    para_name = "log_lambda_Y_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &log_lambda_Y_max);
    //
    para_name = "log_lambda_Yobs_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &log_lambda_Yobs_min);
    //
    para_name = "log_lambda_Yobs_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &log_lambda_Yobs_max);
    //
    para_name = "log_lambda_Xobs_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &log_lambda_Xobs_min);
    //
    para_name = "log_lambda_Xobs_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &log_lambda_Xobs_max);
    //
    free(para_line);
    para_line = NULL;
}






////////////////////////////////////////////
////////////////////////////////////////////
///// initialize prior
////////////////////////////////////////////
////////////////////////////////////////////

double gamma_power_init(double gamma_mean, double gamma_sd)
{
    double gamma_power;
    gamma_power = r8_normal_ab ( gamma_mean, gamma_sd );
    return gamma_power;
}

double logC_init(double logC_min, double logC_max)
{
    double logC;
    logC = r8_unif_ab( logC_min, logC_max );
    return logC;
}

double log_lambda_Y_init(double log_lambda_Y_min, double log_lambda_Y_max)
{
    double log_lambda_Y;
    log_lambda_Y = r8_unif_ab( log_lambda_Y_min, log_lambda_Y_max );
    return log_lambda_Y;
}

double gamma_lambda_Y_init(double gamma_lambda_Y_min, double gamma_lambda_Y_max)
{
    double gamma_lambda_Y;
    gamma_lambda_Y = r8_unif_ab( gamma_lambda_Y_min, gamma_lambda_Y_max );
    return gamma_lambda_Y;
}

double log_lambda_Yobs_init(double log_lambda_Yobs_min, double log_lambda_Yobs_max)
{
    double log_lambda_Yobs;
    log_lambda_Yobs = r8_unif_ab( log_lambda_Yobs_min, log_lambda_Yobs_max );
    return log_lambda_Yobs;
}

double log_lambda_Xobs_init(double log_lambda_Xobs_min, double log_lambda_Xobs_max)
{
    double log_lambda_Xobs;
    log_lambda_Xobs = r8_unif_ab( log_lambda_Xobs_min, log_lambda_Xobs_max );
    return log_lambda_Xobs;
}


// set init parameter at N_ITER = 0
void init_parm_set(int seed, double* chain_parm)
{
    // rand seed
    srand48(seed);
    //
    char *input_file; 
    input_file = "input.ini";
    // read the range
    read_parm_range(input_file);
    //
    //
    chain_parm[0] = gamma_power_init(gamma_mean, gamma_sd);
    chain_parm[1] = logC_init(logC_min, logC_max);
    chain_parm[2] = log_lambda_Y_init(log_lambda_Y_min, log_lambda_Y_max);
    chain_parm[3] = log_lambda_Yobs_init(log_lambda_Yobs_min, log_lambda_Yobs_max);
    chain_parm[4] = log_lambda_Xobs_init(log_lambda_Xobs_min, log_lambda_Xobs_max);
    chain_parm[5] = gamma_lambda_Y_init(gamma_lambda_Y_min, gamma_lambda_Y_max);
    //
    ////// set sigma tuning range of parms
    //
    sigma_parm_min = (double *) malloc(sizeof(double)*N_parm);
    sigma_parm_max = (double *) malloc(sizeof(double)*N_parm);
    calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max);
}


////////////////////////////////////////////
////////////////////////////////////////////
///// init sigma ...
////////////////////////////////////////////
////////////////////////////////////////////

int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max)
{
    //
    int sigma_bound_ok = 0;

    sigma_parm_min[0] = sigma_scale_min * (gamma_sd);
    sigma_parm_min[1] = sigma_scale_min * (logC_max - logC_min);
    sigma_parm_min[2] = sigma_scale_min * (log_lambda_Y_max - log_lambda_Y_min);
    sigma_parm_min[3] = sigma_scale_min * (log_lambda_Yobs_max - log_lambda_Yobs_min);
    sigma_parm_min[4] = sigma_scale_min * (log_lambda_Xobs_max - log_lambda_Xobs_min);
    sigma_parm_min[5] = sigma_scale_min * (gamma_lambda_Y_max - gamma_lambda_Y_min);
    //
    sigma_parm_max[0] = sigma_scale_max * (gamma_sd);
    sigma_parm_max[1] = sigma_scale_max * (logC_max - logC_min);
    sigma_parm_max[2] = sigma_scale_max * (log_lambda_Y_max - log_lambda_Y_min);
    sigma_parm_max[3] = sigma_scale_max * (log_lambda_Yobs_max - log_lambda_Yobs_min);
    sigma_parm_max[4] = sigma_scale_max * (log_lambda_Xobs_max - log_lambda_Xobs_min);
    sigma_parm_max[5] = sigma_scale_max * (gamma_lambda_Y_max - gamma_lambda_Y_min);
    //
    //
    for (int i=0; i<N_parm; i++)
    {
        if (sigma_parm_min[i] < 0)
        {
            sigma_bound_ok++;
        }
        //
        if (sigma_parm_max[i] < 0)
        {
            sigma_bound_ok++;
        }
        //
        if ( (sigma_parm_max[i]-sigma_parm_min[i]) < 0)
        {
            sigma_bound_ok++;
        }
    }
    //
    ///////
    //
    if (sigma_bound_ok > 0)
    {
        printf("ERR: bad sigma_parm_min/max %d TIMES!\n", sigma_bound_ok);
        return 1;
    }
    else
    {
        return 0;
    }
}

// set init gaussian proposal for the sampling
int init_gaussian_proposal(double *ptr_sigma_prop, double init_gp_ratio)
{
    //
    double gamma_power_gp;
    double logC_gp;
    double log_lambda_Y_gp;
    double log_lambda_Yobs_gp;
    double log_lambda_Xobs_gp;
    double gamma_lambda_Y_gp;
    //
    gamma_power_gp = (gamma_max - gamma_min) * init_gp_ratio; 
    logC_gp = (logC_max - logC_min) * init_gp_ratio; 
    log_lambda_Y_gp = (log_lambda_Y_max - log_lambda_Y_min) * init_gp_ratio; 
    log_lambda_Yobs_gp = (log_lambda_Yobs_max - log_lambda_Yobs_min) * init_gp_ratio; 
    log_lambda_Xobs_gp = (log_lambda_Xobs_max - log_lambda_Xobs_min) * init_gp_ratio; 
    gamma_lambda_Y_gp = (gamma_lambda_Y_max - gamma_lambda_Y_min) * init_gp_ratio; 
    //
    *(ptr_sigma_prop) = gamma_power_gp;
    *(ptr_sigma_prop+1) = logC_gp;
    *(ptr_sigma_prop+2) = log_lambda_Y_gp;
    *(ptr_sigma_prop+3) = log_lambda_Yobs_gp;
    *(ptr_sigma_prop+4) = log_lambda_Xobs_gp;
    *(ptr_sigma_prop+5) = gamma_lambda_Y_gp;
    //
    return 0; 
    //
}



////////////////////////////////////////////
////////////////////////////////////////////
///// check boundary
////////////////////////////////////////////
////////////////////////////////////////////


int para_boundary(double *ptr_one_chain_new)
{
    //
    // PingPong outilers inside the allowed prior range
    //
    double gamma_power;
    double logC;
    double log_lambda_Y;
    double log_lambda_Yobs;
    double log_lambda_Xobs;
    double gamma_lambda_Y;
    gamma_power = *ptr_one_chain_new;
    logC = *(ptr_one_chain_new+1);
    log_lambda_Y = *(ptr_one_chain_new+2);
    log_lambda_Yobs = *(ptr_one_chain_new+3);
    log_lambda_Xobs = *(ptr_one_chain_new+4);
    gamma_lambda_Y = *(ptr_one_chain_new+5);
    //
    int debug = 0; // output debug infos 
    //
    //////////////////
    if (gamma_power > gamma_max)
    {
        gamma_power = gamma_max - (gamma_power - gamma_max);
        gamma_power = fmax(gamma_power,gamma_min);
        *ptr_one_chain_new = gamma_power;
    }
    //
    if (gamma_power<gamma_min)
    {
        gamma_power = gamma_min + (gamma_min-gamma_power);
        gamma_power = fmin(gamma_power,gamma_max);
        *ptr_one_chain_new = gamma_power;
    }
    //
    ///////////////////
    if (logC>logC_max)
    {
        logC = logC_max - (logC-logC_max);
        logC = fmax(logC,logC_min);
        *(ptr_one_chain_new+1) = logC;
    }
    //
    if (logC<logC_min)
    {
        logC = logC_min + (logC_min-logC);
        logC = fmin(logC,logC_max);
        *(ptr_one_chain_new+1) = logC;
    }
    //
    ///////////////////
    if (log_lambda_Y>log_lambda_Y_max)
    {
        log_lambda_Y = log_lambda_Y_max - (log_lambda_Y-log_lambda_Y_max);
        log_lambda_Y = fmax(log_lambda_Y,log_lambda_Y_min);
        *(ptr_one_chain_new+2) = log_lambda_Y;
    }
    //
    if (log_lambda_Y<log_lambda_Y_min)
    {
        log_lambda_Y = log_lambda_Y_min + (log_lambda_Y_min-log_lambda_Y);
        log_lambda_Y = fmin(log_lambda_Y,log_lambda_Y_max);
        *(ptr_one_chain_new+2) = log_lambda_Y;
    }
    //
    ///////////////////
    if (log_lambda_Yobs>log_lambda_Yobs_max)
    {
        log_lambda_Yobs = log_lambda_Yobs_max - (log_lambda_Yobs-log_lambda_Yobs_max);
        log_lambda_Yobs = fmax(log_lambda_Yobs,log_lambda_Yobs_min);
        *(ptr_one_chain_new+3) = log_lambda_Yobs;
    }
    //
    if (log_lambda_Yobs<log_lambda_Yobs_min)
    {
        log_lambda_Yobs = log_lambda_Yobs_min + (log_lambda_Yobs_min-log_lambda_Yobs);
        log_lambda_Yobs = fmin(log_lambda_Yobs,log_lambda_Yobs_max);
        *(ptr_one_chain_new+3) = log_lambda_Yobs;
    }
    //
    ///////////////////
    if (log_lambda_Xobs>log_lambda_Xobs_max)
    {
        log_lambda_Xobs = log_lambda_Xobs_max - (log_lambda_Xobs-log_lambda_Xobs_max);
        log_lambda_Xobs = fmax(log_lambda_Xobs,log_lambda_Xobs_min);
        *(ptr_one_chain_new+4) = log_lambda_Xobs;
    }
    //
    if (log_lambda_Xobs<log_lambda_Xobs_min)
    {
        log_lambda_Xobs = log_lambda_Xobs_min + (log_lambda_Xobs_min-log_lambda_Xobs);
        log_lambda_Xobs = fmin(log_lambda_Xobs,log_lambda_Xobs_max);
        *(ptr_one_chain_new+4) = log_lambda_Xobs;
    }
    //
    ///////////////////
    if (gamma_lambda_Y>gamma_lambda_Y_max)
    {
        gamma_lambda_Y = gamma_lambda_Y_max - (gamma_lambda_Y-gamma_lambda_Y_max);
        gamma_lambda_Y = fmax(gamma_lambda_Y,gamma_lambda_Y_min);
        *(ptr_one_chain_new+5) = gamma_lambda_Y;
    }
    //
    if (gamma_lambda_Y<gamma_lambda_Y_min)
    {
        gamma_lambda_Y = gamma_lambda_Y_min + (gamma_lambda_Y_min-gamma_lambda_Y);
        gamma_lambda_Y = fmin(gamma_lambda_Y,gamma_lambda_Y_max);
        *(ptr_one_chain_new+5) = gamma_lambda_Y;
    }
    //
    //
    if (debug)
    {
        //printf("gm,lC,ldY,ldYo,ldXo: %lf  %lf  %lf  %lf  %lf\n", gamma_power, logC, log_lambda_Y, log_lambda_Yobs, log_lambda_Xobs);
        save_debug_para_boundary("gm", gamma_power, gamma_min, gamma_max);
        save_debug_para_boundary("lC", logC, logC_min, logC_max);
        save_debug_para_boundary("ldY", log_lambda_Y, log_lambda_Y_min, log_lambda_Y_max);
        save_debug_para_boundary("ldYo", log_lambda_Yobs, log_lambda_Yobs_min, log_lambda_Yobs_max);
        save_debug_para_boundary("ldXo", log_lambda_Xobs, log_lambda_Xobs_min, log_lambda_Xobs_max);
        save_debug_para_boundary("gmldY", gamma_lambda_Y, gamma_lambda_Y_min, gamma_lambda_Y_max);
    }
    //
    return 0;
}

////////////////////////////////////////////
////////////////////////////////////////////
///// log prior
////////////////////////////////////////////
////////////////////////////////////////////
//
// para_boundary ensures: ((a<=a_max) && (a>=a_min))
//
double prior_gamma_power(double gamma_power, double gamma_mean, double gamma_sd)
{
    return normal_pdf(gamma_power, gamma_mean, gamma_sd);
}
//
double prior_logC(void)
{
    return 1/(logC_max-logC_min);
}
//
double prior_lambda_Y(double lambda_Y, double log_lambda_Y_min, double log_lambda_Y_max)
{
    double lambda_Y_min = exp(log_lambda_Y_min);
    double lambda_Y_max = exp(log_lambda_Y_max);
    return 1/(sqrt(fabs(lambda_Y)) * log(exp(lambda_Y_max)/exp(lambda_Y_min)) );
}
//
double prior_lambda_Yobs(double lambda_Yobs, double log_lambda_Yobs_min, double log_lambda_Yobs_max)
{
    double lambda_Yobs_min = exp(log_lambda_Yobs_min);
    double lambda_Yobs_max = exp(log_lambda_Yobs_max);
    return 1/(sqrt(fabs(lambda_Yobs)) * log(exp(lambda_Yobs_max)/exp(lambda_Yobs_min)) );
}
//
double prior_lambda_Xobs(double lambda_Xobs, double log_lambda_Xobs_min, double log_lambda_Xobs_max)
{
    double lambda_Xobs_min = exp(log_lambda_Xobs_min);
    double lambda_Xobs_max = exp(log_lambda_Xobs_max);
    return 1/(sqrt(fabs(lambda_Xobs)) * log(exp(lambda_Xobs_max)/exp(lambda_Xobs_min)) );
}
//
double prior_gamma_lambda_Y(void)
{
    return 1/(gamma_lambda_Y_max-gamma_lambda_Y_min);
}
//
double log_prior(double *ptr_one_chain)
{
    //
    double gamma_power;
    double logC;
    double log_lambda_Y;
    double log_lambda_Yobs;
    double log_lambda_Xobs;
    //
    double lambda_Y;
    double lambda_Yobs;
    double lambda_Xobs;
    //
    double gamma_lambda_Y;
    //
    double log_prior;
    //
    int debug = 0; // output debug infos 
    //
    gamma_power = *ptr_one_chain;
    logC = *(ptr_one_chain+1);
    log_lambda_Y = *(ptr_one_chain+2);
    log_lambda_Yobs = *(ptr_one_chain+3);
    log_lambda_Xobs = *(ptr_one_chain+4);
    gamma_lambda_Y = *(ptr_one_chain+5);
    //
    lambda_Y = exp(log_lambda_Y);
    lambda_Yobs = exp(log_lambda_Yobs);
    lambda_Xobs = exp(log_lambda_Xobs);

/// to be checked in MODEL logC double log?
    log_prior = log(prior_gamma_power(gamma_power, gamma_mean, gamma_sd)) + 
                log(prior_logC()) + 
                log(prior_lambda_Y(lambda_Y, log_lambda_Y_min, log_lambda_Y_max)) + 
                log(prior_lambda_Yobs(lambda_Yobs, log_lambda_Yobs_min, log_lambda_Yobs_max)) + 
                log(prior_lambda_Xobs(lambda_Xobs, log_lambda_Xobs_min, log_lambda_Xobs_max)) +
                log(prior_gamma_lambda_Y());
    //
    if (debug)
    {
      printf("gm,lC,ldY,ldYo,ldXo,gmldY: %lf  %lf  %lf  %lf  %lf  %lf\n", gamma_power, logC, lambda_Y, lambda_Yobs, lambda_Xobs, gamma_lambda_Y);
      printf("their prior and sum: %lf  %.12e  %lf  %lf  %lf  %lf  %lf\n", log(prior_gamma_power(gamma_power, gamma_mean, gamma_sd)) , log(prior_logC()) , log(prior_lambda_Y(lambda_Y, log_lambda_Y_min, log_lambda_Y_max)) , log(prior_lambda_Yobs(lambda_Yobs, log_lambda_Yobs_min, log_lambda_Yobs_max)) , log(prior_lambda_Xobs(lambda_Xobs, log_lambda_Xobs_min, log_lambda_Xobs_max)), log(prior_gamma_lambda_Y()), log_prior);
    }
    //
    return log_prior;
    // 
}


////////////////////////////////////////////
////////////////////////////////////////////
///// para debug
////////////////////////////////////////////
////////////////////////////////////////////
//
//
int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s%s", results_dir, "/", p_s, ".debug_para_boundary");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    fprintf(out, "%f", p);
    fprintf(out, "  ");
    fprintf(out, "%f", p_min);
    fprintf(out, "  ");
    fprintf(out, "%f", p_max);
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}





