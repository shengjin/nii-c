/* ---- memory allocation ---- */
#include <stdlib.h>
#include <stdio.h>


#define alloc_error_check(p) { \
    if ((p) == NULL) { \
        puts("Memory allocation Failed!\n"); \
        exit(1); \
    } \
}


/* 1d int */
int *alloc_1d_int(int n1)
{
    int *i;
    
    i = (int *) malloc(sizeof(int) * n1);
    alloc_error_check(i);
    return i;
}

void free_1d_int(int *i)
{
    free(i);
}



/* 1d dob */
double *alloc_1d_double(int n1)
{
    double *d;
    
    d = (double *) malloc(sizeof(double) * n1);
    alloc_error_check(d);
    return d;
}

void free_1d_double(double *d)
{
    free(d);
}



/* 2d int */
int **alloc_2d_int(int n1, int n2)
{
    int **ii, *i;
    int j;
    
    ii = (int **) malloc(sizeof(int *) * n1);
    // ii == &ii[0] // ii address
    // ii[0] // poiter to int
    // ii: ptr2int_ptr2int_ptr2int... (n1) 
    alloc_error_check(ii);

    i = (int *) malloc(sizeof(int) * n1 * n2);
    // i == &i[0] // i address
    // i[0] // int
    // i: int_int_int ... (n1*n2)
    alloc_error_check(i);

    ii[0] = i; 
    //ii[0](1st_ptr) point to i_address
    for (j = 1; j < n1; j++) {
        ii[j] = ii[j - 1] + n2;
        //ii[j](jst_ptr)  ( n2*int jump)
        // since ii[0] point to i, this means n2*int jump in i
    }
    return ii;
}

void free_2d_int(int **ii)
{
    //You can free only pointer returned from malloc, which points on the first byte of allocated array, and as well, it will free all the allocated array of bytes. 
    free(ii[0]); // free i
    free(ii);
}



/* 2d dob */
double **alloc_2d_double(int n1, int n2)
{
    double **dd, *d;
    int j;
    
    dd = (double **) malloc(sizeof(double *) * n1);
    alloc_error_check(dd);
    d = (double *) malloc(sizeof(double) * n1 * n2);
    alloc_error_check(d);
    dd[0] = d;
    for (j = 1; j < n1; j++) {
        dd[j] = dd[j - 1] + n2;
    }
    return dd;
}
/* allternative 2d dob
double **alloc_2d_double(int n1, int n2) 
{
    double *d = (double *)malloc(n1*n2*sizeof(double));
    double **dd= (double **)malloc(n1*sizeof(double*));
    for (int i=0; i<n1; i++)
        dd[i] = &(d[n2*i]);
    return dd;
}
*/

void free_2d_double(double **dd)
{
    free(dd[0]);
    free(dd);
}



/* 3d int */
int ***alloc_3d_int(int n1, int n2, int n3)
{
    int ***iii, **ii, *i;
    int j;
    
    iii = (int ***) malloc(sizeof(int **) * n1);
    alloc_error_check(iii);
    ii = (int **) malloc(sizeof(int *) * n1 * n2);
    alloc_error_check(ii);
    iii[0] = ii;
    for (j = 1; j < n1; j++) {
        iii[j] = iii[j - 1] + n2;
    }
    i = (int *) malloc(sizeof(int) * n1 * n2 * n3);
    alloc_error_check(i);
    ii[0] = i;
    for (j = 1; j < n1 * n2; j++) {
        ii[j] = ii[j - 1] + n3;
    }
    return iii;
}

void free_3d_int(int ***iii)
{
    free(iii[0][0]);
    free(iii[0]);
    free(iii);
}



/* 3d dob */
double ***alloc_3d_double(int n1, int n2, int n3)
{
    double ***ddd, **dd, *d;
    int j;
    
    ddd = (double ***) malloc(sizeof(double **) * n1);
    alloc_error_check(ddd);
    dd = (double **) malloc(sizeof(double *) * n1 * n2);
    alloc_error_check(dd);
    ddd[0] = dd;
    for (j = 1; j < n1; j++) {
        ddd[j] = ddd[j - 1] + n2;
    }
    d = (double *) malloc(sizeof(double) * n1 * n2 * n3);
    alloc_error_check(d);
    dd[0] = d;
    for (j = 1; j < n1 * n2; j++) {
        dd[j] = dd[j - 1] + n3;
    }
    return ddd;
}

void free_3d_double(double ***ddd)
{
    free(ddd[0][0]);
    free(ddd[0]);
    free(ddd);
}


