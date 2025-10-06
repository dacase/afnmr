#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sff.h"
FILE* nabout;

//  Very primitive short minimizer, to illustrate the C-api; no argument-
//     checking, etc., etc.
//
//  Usage:  txmin  <parm-file>  <input-coord-file> <output-coord-file>
//       output to stdout

int main( int argc, char *argv[] )
{

   PARMSTRUCT_T *prm;   //  struct to hold info from a prmtop file
   XMIN_OPT_T xo;       //  options for the minimizer
   int natm, iter;
   double *xyz,  *grad, *xyz_ref;
   double energy, grms;
   double start_time = 0.0;   // dummy, since this is minimization, not md

   nabout = stdout;    // change to redirect output (historical kludge)

//   options for the minimizer:

   xmin_opt_init( &xo );  // sets the default parameters;

   xo.maxiter = 2;                  // non-default minimization options:
   xo.grms_tol = 0.0005;
   xo.ls_maxatmov = 0.05;
   xo.print_level = 1;
   xo.method = 2;

//   read in the prmtop file and the coordinates:

   prm = rdparm( argv[1] );    // reads the prmtop file
   natm = prm->Natom;
   xyz = malloc( 3 * natm * (sizeof(double)) );
   xyz_ref = malloc( 3 * natm * (sizeof(double)) );
   grad = malloc( 3 * natm * (sizeof(double)) );
   getxv( argv[2], natm, start_time, xyz, grad );  // reads a restart file
   getxv( argv[2], natm, start_time, xyz_ref, grad );  // reads a restart file

//   setup the force field parameters, and get an initial energy:

   mm_options( "ntpr=1, gb=8, kappa=0.10395, rgbmax=9., cut=9.0, wcons=0. " );

   // solvent frozen; constrain non-hydrogens:
   int* frozen = parseMaskString( ":WAT,Na+,Cl-", prm, xyz, 2 );
   int* constrained = parseMaskString( "!@H*", prm, xyz, 2 );

   mme_init_sff( prm, frozen, constrained, xyz_ref, NULL );
   iter = -1;   // historical flag to give more verbose output
   energy = mme( xyz, grad, &iter );
   energy = mme_rattle( xyz, grad, &iter );

//   run the minimization:

   char title[] = "minimization";
   energy = xmin( mme_rattle,  &natm, xyz, grad,  &energy,  &grms,  &xo );
   energy = mme_rattle( xyz, grad, &iter );
   energy = mme( xyz, grad, &iter );
   putxv( argv[3], title, natm, start_time, xyz, xyz );

}
