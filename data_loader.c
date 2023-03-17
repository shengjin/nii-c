/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int get_nlines_of_file(char *path);

int read_2d_data(char *path, double **data_NlineNdim, int nline_data, int ndim_data, char *delimiter);

// load data 
void mpi_data_loader(int my_rank, int root_rank, int nline_data, int ndim_data, double **data_NlineNdim, char *path, char *delimiter)
{
    //
    if( my_rank == root_rank )
    {
        // read data
        read_2d_data(path, data_NlineNdim, nline_data, ndim_data, delimiter);
    } 
    // 
    // MPI_Bcast the user_data
    MPI_Bcast(&data_NlineNdim[0][0], nline_data*ndim_data, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
}



int read_2d_data(char *path, double **data_NlineNdim, int nline_data, int ndim_data, char *delimiter)
{
    FILE *fp = fopen(path, "r");
    //
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file '%s'\n", path);
    }
    // 
    // Read lines using POSIX function getline
    char *line = NULL;
    size_t len = 0;
    // 
    int iline_local;
    iline_local = 0;
    //
    while(getline(&line, &len, fp) != -1)
    {
        char *token = strtok(line, delimiter);
        // read the first element
        data_NlineNdim[iline_local][0] = atof(token);
        // read the other elements
        for(int j=1; j<ndim_data; j++)
        {
            token = strtok(NULL, delimiter);
            data_NlineNdim[iline_local][j] = atof(token);
        }
        iline_local++;
    }
    //
    fclose(fp);
    free(line);     
    //
    if (nline_data == iline_local) 
    { 
        return 0;
    } 
    else
    {
        return 1;
    }
    //
}






// read the number of line by root_rank and broadcast it
int mpi_get_nlines(int nline_data, int my_rank, char *path, int root_rank)
{
    //
    if( my_rank == root_rank )
    {
        // get the number of line of user_data
        nline_data = get_nlines_of_file(path);
    } 
    // 
    // MPI_Bcast the line_number
    MPI_Bcast(&nline_data, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
    //
    return nline_data;
}




int get_nlines_of_file(char *path)
{
    //
    FILE *fp = fopen(path, "r");
    //
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file '%s'\n", path);
    }
    // 
    int line_number = 0;
    //
    // Extract characters from file and store in character c
    char c;
    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n') // Increment line_number if this character is newline
            line_number++;
    // 
    // Close the file now that we are done with it.
    fclose(fp);
    //
    return line_number; 
}


