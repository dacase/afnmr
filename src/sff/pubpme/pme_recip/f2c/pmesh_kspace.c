/* pmesh_kspace.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h>
#include <stdio.h>
#define max(a,b) ((a) >= (b) ? (a) : (b))

/* Subroutine */ int pmesh_kspace_get_sizes_(int *nfft1, int *nfft2, 
	int *nfft3, int *numatoms, int *order, int *sizfftab, 
	int *sizffwrk, int *siztheta, int *siz_q__, int *
	sizheap, int *sizstack)
{
    static int nfftable;
    extern /* Subroutine */ int get_fftdims__(int *, int *, int *,
	     int *, int *, int *, int *, int *, int *,
	     int *);
    static int nffwork, nfftdim1, nfftdim2, nfftdim3;

/* INPUT */
/*      nfft1,nfft2,nfft3,numatoms,order */
/*      nfft1,nfft2,nfft3 are the dimensions of the charge grid array */
/*      numatoms is number of atoms */
/*      order is the order of B-spline interpolation */
/* OUTPUT */
/*      sizfftab,sizffwrk,siztheta,siz_Q */
/*      sizfftab is permanent 3d fft table storage */
/*      sizffwrk is temporary 3d fft work storage */
/*      siztheta is size of arrays theta1-3 dtheta1-3 */
/*      sizheap is total size of permanent storage */
/*      sizstack is total size of temporary storage */
/* This routine computes the above output parameters needed for */
/* heap or stack allocation. */
    get_fftdims__(nfft1, nfft2, nfft3, &nfftdim1, &nfftdim2, &nfftdim3, &
	    nfftable, &nffwork, sizfftab, sizffwrk);
    *siztheta = *numatoms * *order;
    *siz_q__ = (nfftdim1 << 1) * nfftdim2 * nfftdim3;
    *sizheap = *nfft1 + *nfft2 + *nfft3 + *sizfftab;
    *sizstack = *siz_q__ + *siztheta * 6 + *sizffwrk + *numatoms * 3;
/*     write(6,*)'total HEAP storage needed = ',sizheap */
/*     write(6,*)'total STACK storage needed = ',sizstack */
    return 0;
} /* pmesh_kspace_get_sizes_ */

/* ---------------------------------------------------- */
/* Subroutine */ int pmesh_kspace_setup_(double *bsp_mod1__, double *
	bsp_mod2__, double *bsp_mod3__, double *fftable, double *
	ffwork, int *nfft1, int *nfft2, int *nfft3, int *
	order, int *sizfftab, int *sizffwrk)
{
    static int nfftable;
    extern /* Subroutine */ int fft_setup__(double *, double *, 
	    double *, int *, int *, int *, int *, int 
	    *, int *, int *, int *), get_fftdims__(int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *);
    static int sfft, sffw;
    static double dummy;
    static int nffwork;
    extern /* Subroutine */ int load_bsp_moduli__(double *, double *, 
	    double *, int *, int *, int *, int *);
    static int nfftdim1, nfftdim2, nfftdim3;

/*  see DO_PMESH_KSPACE for explanation of arguments */
    /* Parameter adjustments */
    --bsp_mod1__;
    --bsp_mod2__;
    --bsp_mod3__;
    --fftable;
    --ffwork;

    /* Function Body */
    get_fftdims__(nfft1, nfft2, nfft3, &nfftdim1, &nfftdim2, &nfftdim3, &
	    nfftable, &nffwork, &sfft, &sffw);
    load_bsp_moduli__(&bsp_mod1__[1], &bsp_mod2__[1], &bsp_mod3__[1], nfft1, 
	    nfft2, nfft3, order);
    fft_setup__(&dummy, &fftable[1], &ffwork[1], nfft1, nfft2, nfft3, &
	    nfftdim1, &nfftdim2, &nfftdim3, &nfftable, &nffwork);
    return 0;
} /* pmesh_kspace_setup_ */

/* ---------------------------------------------------- */
/* Subroutine */ int do_pmesh_kspace_(int *numatoms, double *x, 
	double *y, double *z__, double *charge, double *recip,
	 double *volume, double *ewald_coeff__, int *order, 
	int *nfft1, int *nfft2, int *nfft3, double *eer, 
	double *dx, double *dy, double *dz, double *virial, 
	int *sizfftab, int *sizffwrk, int *siztheta, int *
	siz_q__, double *bsp_mod1__, double *bsp_mod2__, double *
	bsp_mod3__, double *fftable, double *q, double *ffwork, 
	double *theta1, double *theta2, double *theta3, 
	double *dtheta1, double *dtheta2, double *dtheta3, 
	double *fr1, double *fr2, double *fr3, double *time)
{
    extern /* Subroutine */ int fft_back__(double *, double *, 
	    double *, int *, int *, int *, int *, int 
	    *, int *, int *, int *);
    static int nfftable;
    extern /* Subroutine */ int grad_sum__(int *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, int *, int *, int *, int *, int 
	    *, int *, int *, double *), fill_charge_grid__(
	    int *, double *, double *, double *, double *,
	     double *, double *, double *, int *, int *, 
	    int *, int *, int *, int *, int *, double 
	    *), scalar_sum__(double *, double *, double *, 
	    double *, double *, double *, double *, int *,
	     int *, int *, int *, int *, int *, 
	    double *, double *), get_bspline_coeffs__(int *, 
	    double *, double *, double *, int *, double *,
	     double *, double *, double *, double *, 
	    double *), get_fftdims__(int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *), fft_forward__(double *, double *, double *
	    , int *, int *, int *, int *, int *, int *
	    , int *, int *);
    static float tim1, tim2;
    static int sfft, sffw;
    extern /* Subroutine */ int second_(float *), get_scaled_fractionals__(
	    int *, double *, double *, double *, double *,
	     int *, int *, int *, double *, double *, 
	    double *);
    static int nffwork, nfftdim1, nfftdim2, nfftdim3;

/* INPUT */
/*       numatoms:  number of atoms */
/*       x,y,z:   atomic coords */
/*       charge  atomic charges */
/*       recip: 3x3 array of reciprocal unit cell vectors (stored as columns) */
/*       volume: the volume of the unit cell */
/*       ewald_coeff:   ewald convergence parameter */
/*       order: the order of Bspline interpolation. E.g. cubic is order 4 */
/*          fifth degree is order 6 etc. The order must be an even number */
/*          and at least 4. */
/*       nfft1,nfft2,nfft3: the dimensions of the charge grid array */
/* OUTPUT */
/*       eer:  ewald reciprocal or k-space  energy */
/*       dx,dy,dz: forces incremented by k-space sum */
/*       virial:  virial due to k-space sum (valid for atomic scaling; */
/*                rigid molecule virial needs a correction term not */
/*                computed here */
/*       time: used to profile the different component routines */
/* SIZES of some arrays */
/* HEAP STORAGE:  These arrays need to be preserved throughout simulation */
/* STACK STORAGE: These arrays can be tossed after leaving this routine */
/*  get some int array dimensions */
    /* Parameter adjustments */
    --fr3;
    --fr2;
    --fr1;
    --dz;
    --dy;
    --dx;
    --charge;
    --z__;
    --y;
    --x;
    recip -= 4;
    --bsp_mod1__;
    --bsp_mod2__;
    --bsp_mod3__;
    --virial;
    --fftable;
    --ffwork;
    --dtheta3;
    --dtheta2;
    --dtheta1;
    --theta3;
    --theta2;
    --theta1;
    --q;
    --time;

    /* Function Body */
    get_fftdims__(nfft1, nfft2, nfft3, &nfftdim1, &nfftdim2, &nfftdim3, &
	    nfftable, &nffwork, &sfft, &sffw);
    get_scaled_fractionals__(numatoms, &x[1], &y[1], &z__[1], &recip[4], 
	    nfft1, nfft2, nfft3, &fr1[1], &fr2[1], &fr3[1]);
    second_(&tim1);
    get_bspline_coeffs__(numatoms, &fr1[1], &fr2[1], &fr3[1], order, &theta1[
	    1], &theta2[1], &theta3[1], &dtheta1[1], &dtheta2[1], &dtheta3[1])
	    ;
    second_(&tim2);
    time[1] = time[1] + tim2 - tim1;
    tim1 = tim2;
    fill_charge_grid__(numatoms, &charge[1], &theta1[1], &theta2[1], &theta3[
	    1], &fr1[1], &fr2[1], &fr3[1], order, nfft1, nfft2, nfft3, &
	    nfftdim1, &nfftdim2, &nfftdim3, &q[1]);
    second_(&tim2);
    time[2] = time[2] + tim2 - tim1;
    tim1 = tim2;
    fft_back__(&q[1], &fftable[1], &ffwork[1], nfft1, nfft2, nfft3, &nfftdim1,
	     &nfftdim2, &nfftdim3, &nfftable, &nffwork);
    second_(&tim2);
    time[3] = time[3] + tim2 - tim1;
    tim1 = tim2;
    scalar_sum__(&q[1], ewald_coeff__, volume, &recip[4], &bsp_mod1__[1], &
	    bsp_mod2__[1], &bsp_mod3__[1], nfft1, nfft2, nfft3, &nfftdim1, &
	    nfftdim2, &nfftdim3, eer, &virial[1]);
    second_(&tim2);
    time[4] = time[4] + tim2 - tim1;
    tim1 = tim2;
    fft_forward__(&q[1], &fftable[1], &ffwork[1], nfft1, nfft2, nfft3, &
	    nfftdim1, &nfftdim2, &nfftdim3, &nfftable, &nffwork);
    second_(&tim2);
    time[3] = time[3] + tim2 - tim1;
    tim1 = tim2;
    grad_sum__(numatoms, &charge[1], &recip[4], &theta1[1], &theta2[1], &
	    theta3[1], &dtheta1[1], &dtheta2[1], &dtheta3[1], &dx[1], &dy[1], 
	    &dz[1], &fr1[1], &fr2[1], &fr3[1], order, nfft1, nfft2, nfft3, &
	    nfftdim1, &nfftdim2, &nfftdim3, &q[1]);
    second_(&tim2);
    time[5] = time[5] + tim2 - tim1;
    return 0;
} /* do_pmesh_kspace_ */

/* ---------------------------------------------------------------------- */
/* Subroutine */ int get_scaled_fractionals__(int *numatoms, double *
	x, double *y, double *z__, double *recip, int *nfft1, 
	int *nfft2, int *nfft3, double *fr1, double *fr2, 
	double *fr3)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int n;
    static double w;

/* INPUT: */
/*      numatoms: number of atoms */
/*      x,y,z: arrays of cartesian coords */
/*      recip: the 3x3 array of reciprocal vectors stored as columns */
/* OUTPUT: */
/*     fr1,fr2,fr3 the scaled and shifted fractional coords */
    /* Parameter adjustments */
    --fr3;
    --fr2;
    --fr1;
    --z__;
    --y;
    --x;
    recip -= 4;

    /* Function Body */
    i__1 = *numatoms;
    for (n = 1; n <= i__1; ++n) {
	w = x[n] * recip[4] + y[n] * recip[5] + z__[n] * recip[6];
	fr1[n] = *nfft1 * (w - round(w) + .5);
	w = x[n] * recip[7] + y[n] * recip[8] + z__[n] * recip[9];
	fr2[n] = *nfft2 * (w - round(w) + .5);
	w = x[n] * recip[10] + y[n] * recip[11] + z__[n] * recip[12];
	fr3[n] = *nfft3 * (w - round(w) + .5);
/* L100: */
    }
    return 0;
} /* get_scaled_fractionals__ */

/* --------------------------------------------------------------- */
/* Subroutine */ int load_bsp_moduli__(double *bsp_mod1__, double *
	bsp_mod2__, double *bsp_mod3__, int *nfft1, int *nfft2, 
	int *nfft3, int *order)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;
    static double w;
    static int maxn;
    extern /* Subroutine */ int fill_bspline__(double *, int *, 
	    double *, double *);
    static double array[25];
    extern /* Subroutine */ int dftmod_(double *, double *, int *)
	    ;
    static double darray[25], bsp_arr__[1000];

/* this routine loads the moduli of the inverse DFT of the B splines */
/* bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions, */
/* Order is the order of the B spline approx. */
    /* Parameter adjustments */
    --bsp_mod1__;
    --bsp_mod2__;
    --bsp_mod3__;

    /* Function Body */
    if (*order > 25) {
/*      write(6,*)'order too large! check on MAXORDER' */
		fprintf( stderr, "order too large! check on MAXORDER\n" );
		exit(1);
    }
/* Computing MAX */
    i__1 = max(*nfft1,*nfft2);
    maxn = max(i__1,*nfft3);
    if (maxn > 1000) {
/*      write(6,*)'nfft1-3 too large! check on MAXNFFT' */
		fprintf( stderr, "nfft1-3 too large! check on MAXNFFT\n" );
		exit(1);
    }
    w = 0.;
    fill_bspline__(&w, order, array, darray);
    i__1 = maxn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bsp_arr__[i__ - 1] = 0.;
/* L100: */
    }
    i__1 = *order + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	bsp_arr__[i__ - 1] = array[i__ - 2];
/* L150: */
    }
    dftmod_(&bsp_mod1__[1], bsp_arr__, nfft1);
    dftmod_(&bsp_mod2__[1], bsp_arr__, nfft2);
    dftmod_(&bsp_mod3__[1], bsp_arr__, nfft3);
    return 0;
} /* load_bsp_moduli__ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int dftmod_(double *bsp_mod__, double *bsp_arr__, 
	int *nfft)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    double cos(double), sin(double);

    /* Local variables */
    static int j, k;
    static double arg, sum1, sum2, tiny, twopi;

/* Computes the modulus of the discrete fourier transform of bsp_arr, */
/*  storing it into bsp_mod */
    /* Parameter adjustments */
    --bsp_arr__;
    --bsp_mod__;

    /* Function Body */
    twopi = 6.2831853071795862;
    tiny = 1e-7;
    i__1 = *nfft;
    for (k = 1; k <= i__1; ++k) {
	sum1 = 0.;
	sum2 = 0.;
	i__2 = *nfft;
	for (j = 1; j <= i__2; ++j) {
	    arg = twopi * (k - 1) * (j - 1) / *nfft;
	    sum1 += bsp_arr__[j] * cos(arg);
	    sum2 += bsp_arr__[j] * sin(arg);
/* L250: */
	}
/* Computing 2nd power */
	d__1 = sum1;
/* Computing 2nd power */
	d__2 = sum2;
	bsp_mod__[k] = d__1 * d__1 + d__2 * d__2;
/* L300: */
    }
    i__1 = *nfft;
    for (k = 1; k <= i__1; ++k) {
	if (bsp_mod__[k] < tiny) {
	    bsp_mod__[k] = (bsp_mod__[k - 1] + bsp_mod__[k + 1]) * .5;
	}
/* L400: */
    }
    return 0;
} /* dftmod_ */

