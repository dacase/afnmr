#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "sff.h"

#define	EQ(a,b)	(strcmp((a),(b))==0)
#define	NE(a,b)	(strcmp((a),(b))!=0)

double x0[10], x[10], b[3000], a[3000 * 10];
double sig[3000];
double ystar[10], xall[10], ybar[10], ybar2[10];
int fix[10], ij[3000], ib[3000];
int i, iarg, ier, ijack, iboot, j, n, m, used;
double denom, f, chi, chi0, calc, calc0, sig0, chii, chi0i, ave, std;
double rmsgrad, dfpred;
int maxiter;

double lsqfunc( double x[], double g[], int *niter ) {

   double f, calc, z, hz2p1, rho, drho;
   int ivar, iobs, iobsu;

   f = 0.000000E+00;used = 0;
   for( ivar = 1;ivar <= n; ++ ivar )g[ivar - 1] = 0.000000E+00;

   for( iobs = 1;iobs <= m;iobs ++  ){
      if( ijack > 0 && ij[iobs - 1] )continue;
      if( iboot > 0 )iobsu = ib[iobs - 1]; else iobsu = iobs;
      if( !(  iobsu > 0 && iobsu <= m ) ){
        fprintf( stderr, "lsq.nab:59 assert 'iobsu > 0 && iobsu <= m' failed.\n" );
        fprintf( stderr, "lsq.nab:59        'iobsu' = %d\n",  iobsu );
        fprintf( stderr, "lsq.nab:59        'iobsu' = %d\n",  iobsu );
        fprintf( stderr, "lsq.nab:59        'm' = %d\n",  m );
        exit( 1 );
      }

      calc = 0.000000E+00;
      for( ivar = 1;ivar <= n;ivar ++  ){
         calc += a[iobsu - 1 + 3000 * ( ivar - 1 )] * x[ivar - 1];
      }
      z = ( b[iobsu - 1] - calc ) / sig[iobsu - 1];
      hz2p1 = 5.000000E-01 * z * z + 1.000000E+00;
      rho = log( hz2p1 );
      drho = z / hz2p1;
      for( ivar = 1;ivar <= n;ivar ++  ){
         if( fix[ivar - 1] )continue;
         g[ivar - 1] -= drho * a[iobsu - 1 + 3000 * ( ivar - 1 )] / sig[iobsu - 1];
      }
      f += rho;
      used ++ ;
   }

   return( ( f ) );
}

int main( int argc, char **argv ) {

   if( argc < 2 || ( EQ( argv[2 - 1], "-help" ) ) ){
      fprintf( stderr, "Usage: lsq n sig x1 x2 .... < input-stream\n" );
      fprintf( stderr, "\nn is number of parameters; sig an error parameter\n" );
      fprintf( stderr, "x1, x2...are initial values of the paramters\n" );
      fprintf( stderr, "   If an initial guess is exactly 1.0, that parameter will be fixed\n" );
      fprintf( stderr, "\nEach line of the input-stream contains the observed value plus\n" );
      fprintf( stderr, "one row of the coefficient matrix\n" );
      fprintf( stderr, "\nstderr will contain fits and statistics;\n" );
      fprintf( stderr, "stdout will have observed and fitted values\n" );
      exit( 1 );
   }

   n = atoi( argv[2 - 1] );
   if( !( n <= 10 ) ){
     fprintf( stderr, "lsq.nab:96 assert 'n <= 10' failed.\n" );
     fprintf( stderr, "lsq.nab:96        'n' = %d\n", n );
     exit( 1 );
   }

   sig0 = atof( argv[3 - 1] );
   for( iarg = 4;iarg <= argc;iarg = iarg + 1 ){
      x0[iarg - 3 - 1] = atof( argv[iarg - 1] );
   }
   ijack = iboot = 0;

   m = 0;
   for( i = 1;;i = i + 1 ){
      if( ( scanf( "%lf",  &b[i - 1] ) ) == EOF )break;
      m = m + 1;
      if( !( m <= 3000 ) ){
        fprintf( stderr, "lsq.nab:109 assert 'm <= 3000' failed.\n" );
        fprintf( stderr, "lsq.nab:109        'm' = %d\n", m );
        exit( 1 );
      }

      for( j = 1;j <= n;j = j + 1 ){
         scanf( "%lf",  &a[i - 1 + 3000 * ( j - 1 )] );
      }
      sig[i - 1] = sig0;
   }
   fprintf( stderr, "%d data points, %d variables\n", m, n );

   for( j = 1;j <= n;j = j + 1 )x[j - 1] = x0[j - 1];

   for( i = 1;i <= n;i ++  ){
      if( x[i - 1] == 1.000000E+00 )fix[i - 1] = 1;
         else fix[i - 1] = 0;
   }

   rmsgrad = 1.000000E-04;
   dfpred = 1.000000E-01;
   maxiter = 2000;
   ier = conjgrad( x,  &n,  &f, lsqfunc,  &rmsgrad,  &dfpred,  &maxiter );
   if( ier < 0 ){
      fprintf( stderr, "conjgrad returns %d\n", ier );
      exit( 1 );
   }

   fprintf( stderr, "Cauchy fit of %d parameters with sigma = %8.3f\n\n", n, sig0 );
   fprintf( stderr, "Parameters:  fit         orig\n" );
   for( i = 1;i <= n;i = i + 1 )
      fprintf( stderr, "        %12.5f %12.5f\n", x[i - 1], x0[i - 1] );

   chi = 0.000000E+00;
   chi0 = 0.000000E+00;
   for( i = 1;i <= m;i = i + 1 ){
      calc = 0.000000E+00;calc0 = 0.000000E+00;
      for( j = 1;j <= n;j = j + 1 ){
         calc = calc + a[i - 1 + 3000 * ( j - 1 )] * x[j - 1];
         calc0 = calc0 + a[i - 1 + 3000 * ( j - 1 )] * x0[j - 1];
      }
      chii = log( 5.000000E-01 * ( ( b[i - 1] - calc ) / sig[i - 1] ) * ( ( b[i - 1] - calc ) / sig[i - 1] ) + 1.000000E+00 );
      chi0i = log( 5.000000E-01 * ( ( b[i - 1] - calc0 ) / sig[i - 1] ) * ( ( b[i - 1] - calc0 ) / sig[i - 1] ) + 1.000000E+00 );

      printf( "%12.5f %12.5f\n", b[i - 1], calc );
      chi += chii;chi0 += chi0i;
   }

   fprintf( stderr, "\nError fn:    fit         orig\n        %12.5f %12.5f\n", chi, chi0 );

   ijack = 0;
   for( i = 1;i <= n;i = i + 1 ){
      ybar[i - 1] = 0.000000E+00;ybar2[i - 1] = 0.000000E+00;
   }

   for( iboot = 1;iboot <= 50;iboot = iboot + 1 ){

      for( i = 1;i <= m;i = i + 1 )ib[i - 1] = m * ( rand2(  ) ) + 1;
      for( i = 1;i <= n;i ++  )x[i - 1] = xall[i - 1];
      ier = conjgrad( x,  &n,  &f, lsqfunc,  &rmsgrad,  &dfpred,  &maxiter );
      if( ier < 0 )fprintf( stderr, "conjgrad returns %d\n", ier );

      for( i = 1;i <= n;i ++  ){
         ybar[i - 1] += x[i - 1];
         ybar2[i - 1] += x[i - 1] * x[i - 1];
      }
   }

   denom = ( 50 - 1 ) * 50;
   fprintf( stderr, "\nBootstrap:    av.          std.  \n" );
   for( i = 1;i <= n;i ++  ){
      ave = ybar[i - 1] / 50;
      ybar2[i - 1] = ybar2[i - 1] / 50;
      std = sqrt( ybar2[i - 1] - ave * ave );
      fprintf( stderr, "        %12.5f %12.5f\n", ave, std );
   }
   exit( 0 );
}
