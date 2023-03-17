#include <math.h>
#include <stdlib.h>


double r8_normal_01 ( );

double r8_normal_ab ( double a, double b );

double r8_logunif_ab ( double a, double b ); 

double r8_unif_ab( double a, double b );

int i4_unif_ab(int a, int b);

double normal_pdf(double x, double mean, double sd);

/******************************************************************************/
/*
  Purpose:

    r8_normal_01() returns a unit pseudonormal R8.

    The standard normal probability distribution function (PDF) has 
    mean 0 and standard deviation 1.

    Because this routine uses the Box Muller method, it requires pairs
    of uniform random values to generate a pair of normal random values.
    This means that on every other call, the code can use the second
    value that it calculated.
*/

double r8_normal_01 ( )
{
  double r1;
  double r2;
  const double r8_pi = 3.141592653589793;
  double x;

  r1 = drand48( );
  r2 = drand48( );
  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r8_pi * r2 );

  return x;
}




/******************************************************************************/
/*
  Purpose:

    r8_normal_ab() returns a scaled pseudonormal R8.

    The normal probability distribution function (PDF) is sampled,
    with mean A and standard deviation B.

  Input:

    double A, the mean of the PDF.

    double B, the standard deviation of the PDF.

  Output:

    double R8_NORMAL_AB, a sample of the normal PDF.
*/

double r8_normal_ab ( double a, double b )
{
  double value;

  value = a + b * r8_normal_01 ( );

  return value;
}




/******************************************************************************/
/*
  Purpose:

    r8_logunif_ab() returns a log uniform pseudo random number R8.
    
    here we use the Inverse_transform_sampling method

  Input:

    double A, the lower limit of the PDF.

    double B, the upper limit of the PDF.

  Output:

    double R8_logunif_ab, a sample of the log uniform PDF.
*/

double r8_logunif_ab( double a, double b )
{
    double value;
    double r8_unif;
    //
    // generate an uniform random number 
    r8_unif = drand48( );
    //
    // calling the inverse cdf
    value = a*exp(r8_unif*log(b/a));
    //
    return value;
}


/******************************************************************************/
/*
  Purpose:

    r8_unif_ab() returns a uniform pseudo random number R8 between a and b.
    
  Input:

    double A, the lower limit of the PDF.

    double B, the upper limit of the PDF.

  Output:

    double R8_unif_ab, a sample of the log uniform PDF.
*/

double r8_unif_ab( double a, double b )
{
    double value;
    double r8_unif;
    //
    // generate an uniform random number 
    r8_unif = drand48( );
    //
    value = r8_unif*(b-a)+a;
    //
    return value;
}



/******************************************************************************/
/*
  Purpose:
    i4_unif_ab() returns a uniform pseudo random number i4 between a and b.
    
  Input:
    int a, the lower limit of the PDF.
    int b, the lower limit of the PDF.

  Output:
    int i4_unif_ab
*/
int i4_unif_ab( int a, int b )
{
    double r8_unif;
    int value;
    //
    // generate an uniform random number 
    r8_unif = drand48( );
    //
    value = round(r8_unif*(double)(b-a))+a;
    //
    return value;
}



/******************************************************************************/
/*
  Purpose:
    i4_unif_0a() returns a rand int between i=0 to i=a
    
  Input:
    int a, the upper limit of the PDF.

  Output:
    int i4_unif_0a
*/
int i4_unif_0a( int a )
{
    double r8_unif;
    int value;
    //
    // generate an uniform random number 
    r8_unif = drand48( );
    //
    value = floor(r8_unif*(double)(a));
    //
    return value;
}



double normal_pdf(double x, double mean, double sd)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - mean) / sd;

    return inv_sqrt_2pi / sd * exp(-0.5 * a * a);
}


