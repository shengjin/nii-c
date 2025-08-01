#/usr/bin/python 


# define the number of parameters
n = 3 # int(input("Set the number of parameters : "))



####################### Part 0 ##########
# the header of user_prior.c 
c = """
//////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEST_SIZE 100

// NOTE: the description of all parameters starts from line 714

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////// Part 0
////// General part, no need to change in most cases
//////////////////////////////////////////////////////////////////////////////
//
// Number of parameters
extern int N_parm;
//
// Gaussian proposal array
double *sigma_parm_min;
double *sigma_parm_max;
//
// Range to scale gaussian proposal
extern double sigma_scale_min;
extern double sigma_scale_max;
//
// initial gaussian proposal ratio to the paramter range
extern double init_gp_ratio;
//
// Define results directory
extern char *results_dir; //dir to save outputs
//
/////////////////////////////////////
//NOTE: change log_prior in mpi_batch/init if func_prototype is changed
double log_prior(double *ptr_one_chain);
//
////////////////////
//calc sigma scale boundary (function in user_prior.c)
int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max);
//
//
// useful random number genetors
double r8_normal_ab ( double a, double b );
double r8_unif_ab ( double a, double b );
double r8_logunif_ab ( double a, double b );
double normal_pdf(double x, double mean, double sd);
//
//
///////////////////
// prototype no need to chaneg, contents need change
void read_parm_range(char *path);
//
//
///////////////////////////////////////////////
// check if modification is needed or not
///////////////////
int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max);
//
//

double bounce_inside(double para, double para_min, double para_max)
{
//
    if (para > para_max)
    {
        para = para_max - (para - para_max);
        para = fmax(para, para_min);
    }
    if (para < para_min)
    {
        para = para_min + (para_min - para);
        para = fmin(para, para_max);
    }
    //
    return para;
}

"""

# writer the header
with open("user_prior.c", "w", encoding="utf-8") as f:
    f.write(c)


################################## Part 1 #### 
c = """
//////////////////////////////////////////////////////////////////////////////
////// Part 1
////// define how many parameters in the model
////// their values are read from input files by read_parm_range(char *path)
//
"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c)
    for i in range(n):
        s = "%s%d%s\n" % ("double para", i, "_min;")
        f.write(s)
    f.write("//\n")
    for i in range(n):
        s = "%s%d%s\n" % ("double para", i, "_max;")
        f.write(s)
    f.write("//")




################################## Part 2 #### 
c0 = """
/////// Part 2: read prior range
///////   NOTE: Input Min and Max in default.
///////   NOTE: Add your own control para (mean, std, etc) in case of needed.
////////////////////////////////////////////
//
void read_parm_range(char *path)
{
    char *para_line;
    char *para_name;
    //
    char *read_onepara(char *path, char *para_name);
    //
    char dummy[DEST_SIZE] = {0};
    //
    //
"""
#
c1 = """
    //
    free(para_line);
    para_line = NULL;
}
"""
#
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    for i in range(n):
        s0 = "%s%d%s\n" % ("    para_name = \"para", i, "_max\";")
        s1 = "    para_line = read_onepara(path, para_name);\n" 
        s2 = "%s%d%s\n" % ("    sscanf(para_line, \"%[^:]:%lf\", dummy, &para", i, "_max);")
        s3 = "%s%d%s\n" % ("    para_name = \"para", i, "_min\";")
        s4 = "    para_line = read_onepara(path, para_name);\n" 
        s5 = "%s%d%s\n" % ("    sscanf(para_line, \"%[^:]:%lf\", dummy, &para", i, "_min);")
        f.write(s0)
        f.write(s1)
        f.write(s2)
        f.write(s3)
        f.write(s4)
        f.write(s5)
        f.write("    //\n")
    f.write(c1)
    
        
        
################################## Part 3 #### 
c0 = """
///////////////////////////////////////////////////////////////////////////////
///////   Part 3: initialize prior, part 1 
///////   NOTE: uniform initialization in default.
///////   NOTE: Specify your own init function in case of needed.
////////////////////////////////////////////
//
// set init parameter at N_ITER = 0
"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s\n" % ("double para",i, "_init(double para" ,i,"_min,double para" ,i,"_max)")
        s1 = "{\n" 
        s2 = "%s%d%s\n" % ("    double para",i,";")
        s3 = "%s%d%s%d%s%d%s\n" % ("    para",i," = r8_unif_ab(para",i,"_min, para",i,"_max);")
        s4 = "%s%d%s\n" %("    return para",i,";") 
        s5 = "}\n"
        f.write(s0)
        f.write(s1)
        f.write(s2)
        f.write(s3)
        f.write(s4)
        f.write(s5)
        f.write("//\n")
   
        
 

################################## Part 4 #### 
c0 = """
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
///////   Part 4 initialize, part 2
///////   NOTE: uniform initialization in default.
///////   NOTE: Specify your own init function in case of needed.
//
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
"""
c1 = """    ////// set sigma tuning range of parms
    //
    sigma_parm_min = (double *) malloc(sizeof(double)*N_parm);
    sigma_parm_max = (double *) malloc(sizeof(double)*N_parm);
    //
    // initilize
    for (int i=0; i<N_parm; i++)
    {
       sigma_parm_min[i] = 0;
       sigma_parm_max[i] = 0;
    }
    //
    calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max);
}
"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s%d%s\n" % ("    chain_parm[",i,"] = para",i,"_init(para",i,"_min, para",i,"_max);")        
        f.write(s0)    
    f.write("    //\n")
    f.write(c1)



################################## Part 5 #### 
c0="""
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 5: calc the scaling boundary of each para
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max)
{
    //
    int sigma_bound_ok = 0;
    
"""
c1 = """
    //
    for (int i=0; i<N_parm; i++)
    {
        // each of the following case is a error, sum how many bad sigma
        if (sigma_parm_min[i] < 0)
        {
            sigma_bound_ok++;
        }
        //
        if (sigma_parm_max[i] < 0)
        {
            sigma_bound_ok++;
        }
        if ((sigma_parm_max[i] - sigma_parm_min[i]) < 0)
        {
            sigma_bound_ok++;
        }
    }
    //
    ///////
    //
    if (sigma_bound_ok > 0)
    {
        printf("ERR: bad sigma_parm_min/max TIMES!\\n");
        return 1;
    }
    else
    {
        return 0;
    }
}

"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s\n" % ("    sigma_parm_min[",i,"]= sigma_scale_min * (para",i,"_max - para",i,"_min);")
        
        f.write(s0)
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s\n" % ("    sigma_parm_max[",i,"]= sigma_scale_max * (para",i,"_max - para",i,"_min);")
        f.write(s0)
    f.write(c1)
 
 
 
 ################################## Part 6.1:init gaussian proposal ##########    
c0 ="""
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 6: init gaussian proposal
///////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
//
// set init gaussian proposal for the sampling
int init_gaussian_proposal(double *ptr_sigma_prop, double init_gp_ratio)
{
    //
    // PingPong outilers inside the allowed prior range
    """
c3 = """
    //
    return 0; 
    //
}
   """
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    f.write("    //\n")
    for i in range(n):
        s0 = "%s%d%s\n" % ("    double para",i,"_gp;")
        f.write(s0)
    f.write("    //\n")
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s\n" % ("    para",i,"_gp = (para",i,"_max - para",i,"_min) * init_gp_ratio;")
        f.write(s0)
    f.write("    //\n")
    for i in range(n):
        s0 = "%s%d%s%d%s\n" % ("    *(ptr_sigma_prop+",i,") = para",i,"_gp;")
        f.write(s0)
    f.write(c3)



################################## Part 6.2:heck if a proposed point is within the r   #################################
c0 ="""
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 6: check if a proposed point is within the range
///////
////////////////////////////////////////////
////////////////////////////////////////////
//
// TODO: add a new function for this ugly long boring function.
//
int para_boundary(double *ptr_one_chain_new)
{
    //
    // PingPong outilers inside the allowed prior range
    //
"""
c3 ="""
    //
    return 0;
}
"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    for i in range(n):
        s0 = "%s%d%s\n" % ("    double para",i,";")
        f.write(s0)
    f.write("    //\n")
    for i in range(n):
        s0 = "%s%d%s%d%s\n" % ("    para",i," = *(ptr_one_chain_new+",i,");")
        f.write(s0)
    f.write("    //\n")
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s%d%s\n" % ("    *(ptr_one_chain_new+",i,") = bounce_inside(para",i,", para",i,"_min, para",i,"_max);")
        f.write(s0)
        f.write("    //\n")
    f.write(c3)




################################## Part 7   #################################  
c0 = """       
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 7: log prior
///////   NOTE: Write your own prior function of all parameters
///////
////////////////////////////////////////////
////////////////////////////////////////////
//
// NOTE: para_boundary ensures: ((a<=a_max) && (a>=a_min))
//
////////////////////////////////////////////
//

/////////////////////////
"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s\n" % ("double prior_para",i,"(double para",i,"_min, double para",i,"_max)")
        s1 = "{\n" 
        s2 = "%s%d%s%d%s\n" % ("    return 1/(para",i,"_max-para",i,"_min);")
        s3 = "}\n"
        f.write(s0)
        f.write(s1)
        f.write(s2)
        f.write(s3)
        f.write("//\n")



################################## Part 8   ################################# 
        
c0="""
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////// Part 8: Combine all prior distributions 
///////
////////////////////////////////////////////
//
double log_prior(double *ptr_one_chain)
{
    //
"""
c1="""
    //
    double log_prior;
    //
"""
c4="""
    //
    return log_prior;
    // 
}
"""
c5="""
    log_prior = 
"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)
    for i in range(n):
        s0 = "%s%d%s\n" % ("    double para",i,";")
        f.write(s0)
    f.write("    //\n")
    f.write(c1)
    for i in range(n):
        s0 = "%s%d%s%d%s\n" % ("    para",i," = *(ptr_one_chain+",i,");")
        f.write(s0)
    f.write("    //\n")
    f.write(c5)  
    for i in range(n):
        s0 = "%s%d%s%d%s%d%s\n" % ("                log(prior_para",i,"(para",i,"_min, para",i,"_max)) +")
        f.write(s0)
    f.write("                +0;")
    f.write("    //\n")
    f.write(c4)




c0="""
////////////////////////////////////////////
////////////////////////////////////////////
/////  Save para debug
/////  NO need to change
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
    //printf("%s\\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\\n");
        exit(3);
    }
    //
    fprintf(out, "%f", p);
    fprintf(out, "  ");
    fprintf(out, "%f", p_min);
    fprintf(out, "  ");
    fprintf(out, "%f", p_max);
    fprintf(out, "\\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\\n");
    //
    return 0;
}
"""
with open("user_prior.c", "a", encoding="utf-8") as f:
    f.write(c0)



print("DONE.")

