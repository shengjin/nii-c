/////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double a_min;
double a_max;
double b_min;
double b_max;
double c_min;
double c_max;
double d_min;
double d_max;
double f_min;
double f_max;

double *sigma_parm_min;
double *sigma_parm_max;

extern double sigma_scale_min;
extern double sigma_scale_max;

extern int N_parm;

extern char *results_dir; //dir to save outputs

extern double init_gp_ratio;

/////////////////////
// NOTE: change log_prior in mpi_batch/init if func_prototype is changed 
double log_prior(double *ptr_one_chain);


// calc sigma scale boundary (function in user_prior.c)
int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max);

double r8_logunif_ab ( double a, double b ); 

double r8_unif_ab( double a, double b );

void read_parm_range(char *path);

int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max);


#define DEST_SIZE 100

// read parm range
void read_parm_range(char *path)
{
    char *para_line; 
    char *para_name;
    //
    char *read_onepara(char *path, char *para_name);
    //
    //
    char dummy[DEST_SIZE] = {0};
    //
    para_name = "a_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &a_min);
    //
    para_name = "a_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &a_max);
    //
    para_name = "b_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &b_min);
    //
    para_name = "b_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &b_max);
    //
    para_name = "c_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &c_min);
    //
    para_name = "c_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &c_max);
    //
    para_name = "d_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &d_min);
    //
    para_name = "d_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &d_max);
    //
    para_name = "f_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &f_min);
    //
    para_name = "f_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &f_max);
    //
    free(para_line);
    para_line = NULL;
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
    // a 
    chain_parm[0] = r8_unif_ab(a_min, a_max);
    // b 
    chain_parm[1] = r8_unif_ab(b_min, b_max);
    // c 
    chain_parm[2] = r8_unif_ab(c_min, c_max);
    // d 
    chain_parm[3] = r8_unif_ab(d_min, d_max);
    // e 
    chain_parm[4] = r8_unif_ab(f_min, f_max);
    //
    ////// set sigma tuning range of parms
    //
    sigma_parm_min = (double *) malloc(sizeof(double)*N_parm);
    sigma_parm_max = (double *) malloc(sizeof(double)*N_parm);
    calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max);
}


// set init gaussian proposal for the sampling
int init_gaussian_proposal(double *ptr_sigma_prop, double init_gp_ratio)
{
    //
    double a_gp;
    double b_gp;
    double c_gp;
    double d_gp;
    double f_gp;
    //
    a_gp = (a_max-a_min) * init_gp_ratio; 
    b_gp = (b_max-b_min) * init_gp_ratio; 
    c_gp = (c_max-c_min) * init_gp_ratio; 
    d_gp = (d_max-d_min) * init_gp_ratio; 
    f_gp = (f_max-f_min) * init_gp_ratio; 
    //
    *(ptr_sigma_prop) = a_gp;
    *(ptr_sigma_prop+1) = b_gp;
    *(ptr_sigma_prop+2) = c_gp;
    *(ptr_sigma_prop+3) = d_gp;
    *(ptr_sigma_prop+4) = f_gp;
    //
    return 0; 
    //
}


//////////////////////////////////
//
// para_boundary ensures: ((a<=a_max) && (a>=a_min))
double prior_a(void)
{
    return 1/(a_max-a_min);
}

// para_boundary ensures: ((b<=b_max) && (b>=b_min))
double prior_b(void)
{
    return 1/(b_max-b_min);
}
//
double prior_c(void)
{
    return 1/(c_max-c_min);
}
//
double prior_d(void)
{
    return 1/(d_max-d_min);
}
//
double prior_f(void)
{
    return 1/(f_max-f_min);
}


double log_prior(double *ptr_one_chain)
{
    //
    double a;
    double log_a;
    double b;
    double log_b;
    double c;
    double log_c;
    double d;
    double log_d;
    double f;
    double log_f;
    //
    double log_prior;
    //
    int debug = 0; // output debug infos 
    //
    a = *ptr_one_chain;
    b = *(ptr_one_chain+1);
    c = *(ptr_one_chain+2);
    d = *(ptr_one_chain+3);
    f = *(ptr_one_chain+4);
    log_a = log(prior_a());
    log_b = log(prior_b());
    log_c = log(prior_c());
    log_d = log(prior_d());
    log_f = log(prior_f());
    log_prior = log_a + log_b +log_c + log_d +log_f;
    //
    if (debug)
    {
        printf("a,b,d,la,lb,ld: %lf  %lf  %lf  %lf  %lf  %lf\n", a, b, d, log_a, log_b, log_d);
        printf("a,b,d,la,lb,ld: %lf  %lf  %lf  %lf  %lf  %lf\n", a, c, f, log_c, log_f, log_d);
    }
    return log_prior;
    // 
}

//////////////////////////////////


////
int para_boundary(double *ptr_one_chain_new)
{
    //
    // PingPong outilers inside the allowed prior range
    //
    double a;
    double b;
    double c;
    double d;
    double f;
    a = *ptr_one_chain_new;
    b = *(ptr_one_chain_new+1);
    c = *(ptr_one_chain_new+2);
    d = *(ptr_one_chain_new+3);
    f = *(ptr_one_chain_new+4);
    //
    int debug = 0; // output debug infos 
    //
    //////////////////
    if (a>a_max)
    {
        a = a_max - (a-a_max);
        a = fmax(a,a_min);
        *ptr_one_chain_new = a;
    }
    //
    if (a<a_min)
    {
        a = a_min + (a_min-a);
        a = fmin(a,a_max);
        *ptr_one_chain_new = a;
    }
    //
    ///////////////////
    if (b>b_max)
    {
        b = b_max - (b-b_max);
        b = fmax(b,b_min);
        *(ptr_one_chain_new+1) = b;
    }
    //
    if (b<b_min)
    {
        b = b_min + (b_min-b);
        b = fmin(b,b_max);
        *(ptr_one_chain_new+1) = b;
    }
    //
    ///////////////////
    if (c>c_max)
    {
        c = c_max - (c-c_max);
        c = fmax(c,c_min);
        *(ptr_one_chain_new+2) = c;
    }
    //
    if (c<c_min)
    {
        c = c_min + (c_min-c);
        c = fmin(c,c_max);
        *(ptr_one_chain_new+2) = c;
    }
    //
    ///////////////////
    if (d>d_max)
    {
        d = d_max - (d-d_max);
        d = fmax(d,d_min);
        *(ptr_one_chain_new+3) = d;
    }
    //
    if (d<d_min)
    {
        d = d_min + (d_min-d);
        d = fmin(d,d_max);
        *(ptr_one_chain_new+3) = d;
    }
    //
    ///////////////////
    if (f>f_max)
    {
        f = f_max - (f-f_max);
        f = fmax(f,f_min);
        *(ptr_one_chain_new+4) = f;
    }
    //
    if (f<f_min)
    {
        f = f_min + (f_min-f);
        f = fmin(f,f_max);
        *(ptr_one_chain_new+4) = f;
    }
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    if (debug)
    {
        printf("a,b,d: %lf  %lf  %lf\n", a, b, d);
        save_debug_para_boundary("a", a, a_min, a_max);
        save_debug_para_boundary("b", b, b_min, b_max);
        save_debug_para_boundary("d", d, d_min, d_max);
    }
    //
    return 0;
}


/////////////////////////////////
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






int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max)
{
    //
    int sigma_bound_ok = 0;

    sigma_parm_min[0] = sigma_scale_min * (a_max-a_min);
    sigma_parm_min[1] = sigma_scale_min * (b_max-b_min);
    sigma_parm_min[2] = sigma_scale_min * (c_max-c_min);
    sigma_parm_min[3] = sigma_scale_min * (d_max-d_min);
    sigma_parm_min[4] = sigma_scale_min * (f_max-f_min);
    //
    sigma_parm_max[0] = sigma_scale_max * (a_max-a_min);
    sigma_parm_max[1] = sigma_scale_max * (b_max-b_min);
    sigma_parm_max[2] = sigma_scale_max * (c_max-c_min);
    sigma_parm_max[3] = sigma_scale_max * (d_max-d_min);
    sigma_parm_max[4] = sigma_scale_max * (f_max-f_min);
    //
    //////// check
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




