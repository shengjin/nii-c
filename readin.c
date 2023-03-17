////////////
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>


#define DEST_SIZE 100


int N_parm; // N parameters of a sampling model
int N_iter; // N iterations of each chain
int N_beta; // N beta for N_beta parallel tempering chains
double *Beta_Values;
//
int n_iter_a_stack;      // 
int n_iter_a_batch_base; // 
int n_iter_a_batch_rand; //
int N_swap;              // 
//
int n_iter_in_tune;      // 
//
double ar_ok_lower; 
double ar_ok_upper; 
double ar_best; 
double ar_accept_diff;
//
double sigma_scale_min;
double sigma_scale_max;
double sigma_scale_half_ratio;
//
double sigma_jumpin_ratio;
//
unsigned i_save_begin; // burning

double init_gp_ratio; // init_gp_sigma / prior_range
unsigned init_rand_seed; // init_random seed

int Fout_Len;
char *FoutPre; //prefix of output chains_files
char *FoutSuf; //suffix of output chains_files
char *results_dir;  //dir to save outputs

// user data file
char *Data_file;  
char *Delimiter; 
int ndim_data; 

// read the line containing a specific para_name
char *read_onepara(char *path, char *para_name);

// wrapper: to read all the input.ini parameters
int read_input_ini(char *path);

// read grid info of model chains from input file
void read_chains_grid(char *path);

// read beta values
int read_beta_values(char *path);

// setting prefix and suffix of the output chains_dat
void read_chain_out_name(char *path);

// read describtion of user's datafile
void read_data_desc(char *path);

// read ini and control parms of bayes sampling 
void read_sampling_para(char *path);

// check and mkdir for outputs
int make_dir(char *path);



char *read_onepara(char *path, char *para_name)
{
    //getline() is initially called with no buffer allocated. When calling it, getline() allocates a buffer, reads the first line and places the line's contents in the new buffer. On subsequent calls, getline() updates the same buffer and only reallocates the buffer when it is no longer large enough to fit the whole line. The temporary buffer is then freed when we are done with the file.
    char *line_buf = NULL;
    size_t line_buf_size = 0;
    //
    char *para_line; 
    //
    FILE *fp = fopen(path, "r");
    //
    int found = 0;
    //
    para_line = (char*)malloc(DEST_SIZE);
    //
    if (fp == NULL)
    {
	fprintf(stderr, "Error opening file '%s'\n", path);
    }
    //
    // Loop through until we are done with the file.
    while ((getline(&line_buf, &line_buf_size, fp)) >= 0)
    {
	if (strstr(line_buf, para_name))
	{
	    strcpy(para_line, line_buf);
	    found++;
	}
    }
    //
    if (found > 1)
    {
        fprintf(stderr, "Error(readin.c): duplicate defination of '%s'!\n", para_name);
    	para_line = NULL;
    } 
    else if (found == 0)
    {
        fprintf(stderr, "Error: '%s' not found!\n", para_name);
    	para_line = NULL;
    }
    //
    // Close the file now that we are done with it.
    fclose(fp);
    //
    // Free the allocated line buffer (required by strstr())
    free(line_buf);
    // Release the dangling pointer
    line_buf = NULL;
    // 
    return para_line;
}


/*******************************************/

void read_chains_grid(char *path)
{
    char *para_line; 
    //
    char *read_onepara(char *path, char *para_name);
    //
    char *para_name;
    //
    char dummy[DEST_SIZE] = {0};
    //
    para_name = "N_parm";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &N_parm);
    //
    para_name = "N_iter";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &N_iter);
    //
    para_name = "N_beta";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &N_beta);
    //
    free(para_line);
    para_line = NULL;
}

void read_chain_out_name(char *path)
{
    char *para_line; 
    //
    char *read_onepara(char *path, char *para_name);
    //
    char *para_name;
    //
    char dummy[DEST_SIZE] = {0};
    //
    FoutPre = (char*)malloc(DEST_SIZE);
    FoutSuf = (char*)malloc(DEST_SIZE);
    results_dir = (char*)malloc(DEST_SIZE);
    //
    para_name = "Fout_Len";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &Fout_Len);
    //
    para_name = "FoutPre";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%s", dummy, FoutPre);
    //
    para_name = "FoutSuf";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%s", dummy, FoutSuf);
    //
    para_name = "results_dir";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%s", dummy, results_dir);
    //
    free(para_line);
    para_line = NULL;
}



/*******************************************/

void read_data_desc(char *path)
{
    char *para_line; 
    //
    char *read_onepara(char *path, char *para_name);
    //
    char *para_name;
    //
    char dummy[DEST_SIZE] = {0};
    //
    Data_file = (char*)malloc(DEST_SIZE);
    Delimiter = (char*)malloc(DEST_SIZE);
    //
    para_name = "Data_file";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%s", dummy, Data_file);
    //
    para_name = "ndim_data";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &ndim_data);
    //
    para_name = "Delimiter";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%s", dummy, Delimiter);
    if (strcmp(Delimiter, "blank")==0)
    {
        Delimiter = " ";
    }
    //
    free(para_line);
    para_line = NULL;
}



int read_beta_values(char *path)
{
    char *para_line; 
    //
    char *read_onepara(char *path, char *para_name);
    //
    char *para_name;
    //
    char* token;
    //
    int found = 0;
    //
    Beta_Values = (double *) malloc(sizeof(double)*N_beta);
    //
    para_name = "Beta_Values";
    para_line = read_onepara(path, para_name);
    //
    token = strtok(para_line, ":");
 
    while ((token = strtok(NULL, ",")))
    {
	Beta_Values[found] = atof(token);
	found++;
    }
    //
    free(para_line);
    para_line = NULL;
    //
    // check number of found and number of N_beta
    if (found == N_beta)
    {
        return EXIT_SUCCESS; 
    } 
    else 
    {
	Beta_Values = NULL;
        return EXIT_FAILURE; 
    }
}


void read_sampling_para(char *path)
{
    char *para_line; 
    //
    char *read_onepara(char *path, char *para_name);
    //
    char *para_name;
    //
    char dummy[DEST_SIZE] = {0};
    //
    para_name = "init_gp_ratio";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &init_gp_ratio);
    //
    para_name = "init_rand_seed";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%u", dummy, &init_rand_seed);
    //
    para_name = "i_save_begin";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%u", dummy, &i_save_begin);
    //
    para_name = "n_iter_a_stack";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &n_iter_a_stack);
    //
    para_name = "n_iter_a_batch_base";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &n_iter_a_batch_base);
    //
    para_name = "n_iter_a_batch_rand";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &n_iter_a_batch_rand);
    //
    para_name = "n_iter_in_tune";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &n_iter_in_tune);
    //
    para_name = "ar_ok_lower";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &ar_ok_lower);
    //
    para_name = "ar_ok_upper";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &ar_ok_upper);
    //
    para_name = "ar_best";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &ar_best);
    //
    para_name = "ar_accept_diff";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &ar_accept_diff);
    //
    para_name = "sigma_scale_half_ratio";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &sigma_scale_half_ratio);
    //
    para_name = "sigma_scale_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &sigma_scale_min);
    //
    para_name = "sigma_scale_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &sigma_scale_max);
    //
    para_name = "sigma_jumpin_ratio";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &sigma_jumpin_ratio);
    //
    para_name = "N_swap";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%d", dummy, &N_swap);
    //
    free(para_line);
    para_line = NULL;
}



int  make_dir(char *path)
{
    if (mkdir(path, 0755) == -1)
    {
        printf("\nWARNING: dir: %s exists!\n", path);
        printf("-------  may overwrite previous results!\n");
        printf("-------  MSG from readin.c.\n");
        printf("\n");
        //
        return 1;
        //
    }
    //
    return 0;
}


// wraper: read all the input parameters
int read_input_ini(char *path)
{
    read_chains_grid(path);
    read_chain_out_name(path);
    // make sure read_chains_grid is called before read_beta_values.
    read_beta_values(path);
    // read  describtion of user Data file
    read_data_desc(path);
    // read control parms of bayes sampling
    read_sampling_para(path);
    //
    return 0;
}


