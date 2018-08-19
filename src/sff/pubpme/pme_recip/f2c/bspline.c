/* bspline.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* --------------------------------------------------------------------- */
/* Subroutine */ int get_bspline_coeffs__(int *numatoms, double *fr1, 
	double *fr2, double *fr3, int *order, double *theta1, 
	double *theta2, double *theta3, double *dtheta1, 
	double *dtheta2, double *dtheta3)
{
    /* System generated locals */
    int theta1_dim1, theta1_offset, theta2_dim1, theta2_offset, 
	    theta3_dim1, theta3_offset, dtheta1_dim1, dtheta1_offset, 
	    dtheta2_dim1, dtheta2_offset, dtheta3_dim1, dtheta3_offset, i__1;

    /* Local variables */
    static int n;
    static double w;
    extern /* Subroutine */ int fill_bspline__(double *, int *, 
	    double *, double *);

/* --------------------------------------------------------------------- */
/* INPUT: */
/*      numatoms: number of atoms */
/*      fr1,fr2,fr3 the scaled and shifted fractional coords */
/*      order: the order of spline interpolation */
/* OUTPUT */
/*      theta1,theta2,theta3: the spline coeff arrays */
/*      dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays */
/* --------------------------------------------------------------------- */
    /* Parameter adjustments */
    --fr3;
    --fr2;
    --fr1;
    dtheta3_dim1 = *order;
    dtheta3_offset = 1 + dtheta3_dim1;
    dtheta3 -= dtheta3_offset;
    dtheta2_dim1 = *order;
    dtheta2_offset = 1 + dtheta2_dim1;
    dtheta2 -= dtheta2_offset;
    dtheta1_dim1 = *order;
    dtheta1_offset = 1 + dtheta1_dim1;
    dtheta1 -= dtheta1_offset;
    theta3_dim1 = *order;
    theta3_offset = 1 + theta3_dim1;
    theta3 -= theta3_offset;
    theta2_dim1 = *order;
    theta2_offset = 1 + theta2_dim1;
    theta2 -= theta2_offset;
    theta1_dim1 = *order;
    theta1_offset = 1 + theta1_dim1;
    theta1 -= theta1_offset;

    /* Function Body */
    i__1 = *numatoms;
    for (n = 1; n <= i__1; ++n) {
	w = fr1[n] - (int) fr1[n];
	fill_bspline__(&w, order, &theta1[n * theta1_dim1 + 1], &dtheta1[n * 
		dtheta1_dim1 + 1]);
	w = fr2[n] - (int) fr2[n];
	fill_bspline__(&w, order, &theta2[n * theta2_dim1 + 1], &dtheta2[n * 
		dtheta2_dim1 + 1]);
	w = fr3[n] - (int) fr3[n];
	fill_bspline__(&w, order, &theta3[n * theta3_dim1 + 1], &dtheta3[n * 
		dtheta3_dim1 + 1]);
/* L100: */
    }
    return 0;
} /* get_bspline_coeffs__ */

/* --------------------------------------------------- */
/* Subroutine */ int fill_bspline__(double *w, int *order, double 
	*array, double *darray)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    extern /* Subroutine */ int one_pass__(double *, double *, 
	    int *);
    static int k;
    extern /* Subroutine */ int diff_(double *, double *, int *), 
	    init_(double *, double *, int *);

/* ---------- use standard B-spline recursions: see doc file */
/* do linear case */
    /* Parameter adjustments */
    --darray;
    --array;

    /* Function Body */
    init_(&array[1], w, order);
/* compute standard b-spline recursion */
    i__1 = *order - 1;
    for (k = 3; k <= i__1; ++k) {
	one_pass__(&array[1], w, &k);
/* L10: */
    }
/* perform standard b-spline differentiation */
    diff_(&array[1], &darray[1], order);
/* one more recursion */
    one_pass__(&array[1], w, order);
    return 0;
} /* fill_bspline__ */

/* --------------------------------------------------- */
/* Subroutine */ int init_(double *c__, double *x, int *order)
{
    /* Parameter adjustments */
    --c__;

    /* Function Body */
    c__[*order] = 0.;
    c__[2] = *x;
    c__[1] = 1. - *x;
    return 0;
} /* init_ */

/* ------------------------------------- */
/* Subroutine */ int one_pass__(double *c__, double *x, int *k)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int j;
    static double div;

    /* Parameter adjustments */
    --c__;

    /* Function Body */
    div = 1. / (*k - 1);
    c__[*k] = div * *x * c__[*k - 1];
    i__1 = *k - 2;
    for (j = 1; j <= i__1; ++j) {
	c__[*k - j] = div * ((*x + j) * c__[*k - j - 1] + (*k - j - *x) * c__[
		*k - j]);
/* L100: */
    }
    c__[1] = div * (1 - *x) * c__[1];
    return 0;
} /* one_pass__ */

/* ------------------------------------- */
/* Subroutine */ int diff_(double *c__, double *d__, int *order)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int j;

    /* Parameter adjustments */
    --d__;
    --c__;

    /* Function Body */
    d__[1] = -c__[1];
    i__1 = *order;
    for (j = 2; j <= i__1; ++j) {
	d__[j] = c__[j - 1] - c__[j];
/* L10: */
    }
    return 0;
} /* diff_ */

