/* charge_grid.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#define i_sign(a,b) (*(b) >= 0 ? abs(*(a)) : -abs(*(a)))

/* Subroutine */ int fill_charge_grid__(int *numatoms, double *charge,
	 double *theta1, double *theta2, double *theta3, 
	double *fr1, double *fr2, double *fr3, int *order, 
	int *nfft1, int *nfft2, int *nfft3, int *nfftdim1, 
	int *nfftdim2, int *nfftdim3, double *q)
{
    /* System generated locals */
    int theta1_dim1, theta1_offset, theta2_dim1, theta2_offset, 
	    theta3_dim1, theta3_offset, q_dim2, q_dim3, q_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static int i__, j, k, n, i0, j0, k0, ith1, ith2, ith3;
    static double prod;
    static int ntot;
    extern /* Subroutine */ int clearq_(double *, int *);

/* --------------------------------------------------------------------- */
/* INPUT: */
/*      numatoms:  number of atoms */
/*      charge: the array of atomic charges */
/*      theta1,theta2,theta3: the spline coeff arrays */
/*      fr1,fr2,fr3 the scaled and shifted fractional coords */
/*      nfft1,nfft2,nfft3: the charge grid dimensions */
/*      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims */
/*      order: the order of spline interpolation */
/* OUTPUT: */
/*      Q the charge grid */
/* --------------------------------------------------------------------- */
    /* Parameter adjustments */
    --fr3;
    --fr2;
    --fr1;
    --charge;
    theta3_dim1 = *order;
    theta3_offset = 1 + theta3_dim1;
    theta3 -= theta3_offset;
    theta2_dim1 = *order;
    theta2_offset = 1 + theta2_dim1;
    theta2 -= theta2_offset;
    theta1_dim1 = *order;
    theta1_offset = 1 + theta1_dim1;
    theta1 -= theta1_offset;
    q_dim2 = *nfftdim1;
    q_dim3 = *nfftdim2;
    q_offset = 1 + 2 * (1 + q_dim2 * (1 + q_dim3));
    q -= q_offset;

    /* Function Body */
    ntot = (*nfftdim1 << 1) * *nfftdim2 * *nfftdim3;
    clearq_(&q[q_offset], &ntot);
    i__1 = *numatoms;
    for (n = 1; n <= i__1; ++n) {
	k0 = (int) fr3[n] - *order;
	i__2 = *order;
	for (ith3 = 1; ith3 <= i__2; ++ith3) {
	    ++k0;
	    k = k0 + 1 + (*nfft3 - i_sign(nfft3, &k0)) / 2;
	    j0 = (int) fr2[n] - *order;
	    i__3 = *order;
	    for (ith2 = 1; ith2 <= i__3; ++ith2) {
		++j0;
		j = j0 + 1 + (*nfft2 - i_sign(nfft2, &j0)) / 2;
		prod = theta2[ith2 + n * theta2_dim1] * theta3[ith3 + n * 
			theta3_dim1] * charge[n];
		i0 = (int) fr1[n] - *order;
		i__4 = *order;
		for (ith1 = 1; ith1 <= i__4; ++ith1) {
		    ++i0;
		    i__ = i0 + 1 + (*nfft1 - i_sign(nfft1, &i0)) / 2;
		    q[(i__ + (j + k * q_dim3) * q_dim2 << 1) + 1] += theta1[
			    ith1 + n * theta1_dim1] * prod;
/* L100: */
		}
/* L150: */
	    }
/* L200: */
	}
/* L300: */
    }
    return 0;
} /* fill_charge_grid__ */

/* ----------------------------------------------------------- */
/* Subroutine */ int clearq_(double *q, int *ntot)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

    /* Parameter adjustments */
    --q;

    /* Function Body */
    i__1 = *ntot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = 0.;
/* L10: */
    }
    return 0;
} /* clearq_ */

/* ----------------------------------------------------------- */
/* Subroutine */ int scalar_sum__(double *q, double *ewaldcof, 
	double *volume, double *recip, double *bsp_mod1__, 
	double *bsp_mod2__, double *bsp_mod3__, int *nfft1, 
	int *nfft2, int *nfft3, int *nfftdim1, int *nfftdim2, 
	int *nfftdim3, double *eer, double *vir)
{
    /* System generated locals */
    int q_dim2, q_dim3, q_offset, i__1;
    double d__1, d__2;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static int k, k1, k2, k3, m1, m2, m3;
    static double pi;
    static int nf1, nf2, nf3;
    static double fac;
    static int nff, ind, jnd;
    static double msq, mhat1, mhat2, mhat3, denom, eterm, vterm, struc2, 
	    energy;
    static int indtop;

    /* Parameter adjustments */
    recip -= 4;
    --bsp_mod1__;
    --bsp_mod2__;
    --bsp_mod3__;
    q_dim2 = *nfftdim1;
    q_dim3 = *nfftdim2;
    q_offset = 1 + 2 * (1 + q_dim2 * (1 + q_dim3));
    q -= q_offset;
    --vir;

    /* Function Body */
    indtop = *nfft1 * *nfft2 * *nfft3;
    pi = 3.14159265358979323846f;
/* Computing 2nd power */
    d__1 = pi;
/* Computing 2nd power */
    d__2 = *ewaldcof;
    fac = d__1 * d__1 / (d__2 * d__2);
    nff = *nfft1 * *nfft2;
    nf1 = *nfft1 / 2;
    if (nf1 << 1 < *nfft1) {
	++nf1;
    }
    nf2 = *nfft2 / 2;
    if (nf2 << 1 < *nfft2) {
	++nf2;
    }
    nf3 = *nfft3 / 2;
    if (nf3 << 1 < *nfft3) {
	++nf3;
    }
    energy = 0.;
    for (k = 1; k <= 6; ++k) {
	vir[k] = 0.;
/* L10: */
    }
    i__1 = indtop - 1;
    for (ind = 1; ind <= i__1; ++ind) {
/* get k1,k2,k3 from the relationship */
/*           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1 */
	k3 = ind / nff + 1;
	jnd = ind - (k3 - 1) * nff;
	k2 = jnd / *nfft1 + 1;
	k1 = jnd - (k2 - 1) * *nfft1 + 1;
	m1 = k1 - 1;
	if (k1 > nf1) {
	    m1 = k1 - 1 - *nfft1;
	}
	m2 = k2 - 1;
	if (k2 > nf2) {
	    m2 = k2 - 1 - *nfft2;
	}
	m3 = k3 - 1;
	if (k3 > nf3) {
	    m3 = k3 - 1 - *nfft3;
	}
	mhat1 = recip[4] * m1 + recip[7] * m2 + recip[10] * m3;
	mhat2 = recip[5] * m1 + recip[8] * m2 + recip[11] * m3;
	mhat3 = recip[6] * m1 + recip[9] * m2 + recip[12] * m3;
	msq = mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3;
	denom = pi * *volume * bsp_mod1__[k1] * bsp_mod2__[k2] * bsp_mod3__[
		k3] * msq;
	eterm = exp(-fac * msq) / denom;
	vterm = (fac * msq + 1.) * 2. / msq;
/* Computing 2nd power */
	d__1 = q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 1];
/* Computing 2nd power */
	d__2 = q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 2];
	struc2 = d__1 * d__1 + d__2 * d__2;
	energy += eterm * struc2;
	vir[1] += eterm * struc2 * (vterm * mhat1 * mhat1 - 1.);
	vir[2] += eterm * struc2 * (vterm * mhat1 * mhat2);
	vir[3] += eterm * struc2 * (vterm * mhat1 * mhat3);
	vir[4] += eterm * struc2 * (vterm * mhat2 * mhat2 - 1.);
	vir[5] += eterm * struc2 * (vterm * mhat2 * mhat3);
	vir[6] += eterm * struc2 * (vterm * mhat3 * mhat3 - 1.);
	q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 1] = eterm * q[(k1 + (k2 
		+ k3 * q_dim3) * q_dim2 << 1) + 1];
	q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 2] = eterm * q[(k1 + (k2 
		+ k3 * q_dim3) * q_dim2 << 1) + 2];
/* L100: */
    }
    *eer = energy * .5;
    for (k = 1; k <= 6; ++k) {
	vir[k] *= .5;
/* L150: */
    }
    return 0;
} /* scalar_sum__ */

/* ----------------------------------------------------------- */
/* Subroutine */ int grad_sum__(int *numatoms, double *charge, 
	double *recip, double *theta1, double *theta2, double 
	*theta3, double *dtheta1, double *dtheta2, double *
	dtheta3, double *fx, double *fy, double *fz, double *
	fr1, double *fr2, double *fr3, int *order, int *nfft1,
	 int *nfft2, int *nfft3, int *nfftdim1, int *nfftdim2,
	 int *nfftdim3, double *q)
{
    /* System generated locals */
    int theta1_dim1, theta1_offset, theta2_dim1, theta2_offset, 
	    theta3_dim1, theta3_offset, dtheta1_dim1, dtheta1_offset, 
	    dtheta2_dim1, dtheta2_offset, dtheta3_dim1, dtheta3_offset, 
	    q_dim2, q_dim3, q_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int i__, j, k, n;
    static double f1, f2;
    static int i0, j0, k0;
    static double f3;
    static int ith1, ith2, ith3;
    static double term;

/* $DOACROSS LOCAL(f1,f2,f3,k0,k,j0,j,i0,i,term,n,ith1,ith2,ith3), */
/* $&  SHARE(numatoms,fr1,fr2,fr3,charge,Q,fx,fy,fz,recip,order, */
/* $&   nfft1,nfft2,nfft3,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3) */
    /* Parameter adjustments */
    --fr3;
    --fr2;
    --fr1;
    --fz;
    --fy;
    --fx;
    --charge;
    recip -= 4;
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
    q_dim2 = *nfftdim1;
    q_dim3 = *nfftdim2;
    q_offset = 1 + 2 * (1 + q_dim2 * (1 + q_dim3));
    q -= q_offset;

    /* Function Body */
    i__1 = *numatoms;
    for (n = 1; n <= i__1; ++n) {
	f1 = 0.;
	f2 = 0.;
	f3 = 0.;
	k0 = (int) fr3[n] - *order;
	i__2 = *order;
	for (ith3 = 1; ith3 <= i__2; ++ith3) {
	    ++k0;
	    k = k0 + 1 + (*nfft3 - i_sign(nfft3, &k0)) / 2;
	    j0 = (int) fr2[n] - *order;
	    i__3 = *order;
	    for (ith2 = 1; ith2 <= i__3; ++ith2) {
		++j0;
		j = j0 + 1 + (*nfft2 - i_sign(nfft2, &j0)) / 2;
		i0 = (int) fr1[n] - *order;
		i__4 = *order;
		for (ith1 = 1; ith1 <= i__4; ++ith1) {
		    ++i0;
		    i__ = i0 + 1 + (*nfft1 - i_sign(nfft1, &i0)) / 2;
		    term = charge[n] * q[(i__ + (j + k * q_dim3) * q_dim2 << 
			    1) + 1];
/* AARON: here is the change needed--- note you need to define an array to */
/* hold the electrostatic potentials */
/* assume esp(n) is initted to zero like f1-f3 */
/* note I am putting new code in comments-- no time to make it legal */
/*          esp(n) = esp(n) + */
/*            theta1(ith1,n)*theta2(ith2,n)*theta3(ith3,n)*Q(1,i,j,k) */

/* end of new code */
/* force is negative of grad */
		    f1 -= *nfft1 * term * dtheta1[ith1 + n * dtheta1_dim1] * 
			    theta2[ith2 + n * theta2_dim1] * theta3[ith3 + n *
			     theta3_dim1];
		    f2 -= *nfft2 * term * theta1[ith1 + n * theta1_dim1] * 
			    dtheta2[ith2 + n * dtheta2_dim1] * theta3[ith3 + 
			    n * theta3_dim1];
		    f3 -= *nfft3 * term * theta1[ith1 + n * theta1_dim1] * 
			    theta2[ith2 + n * theta2_dim1] * dtheta3[ith3 + n 
			    * dtheta3_dim1];
/* L100: */
		}
/* L150: */
	    }
/* L200: */
	}
	fx[n] = fx[n] + recip[4] * f1 + recip[7] * f2 + recip[10] * f3;
	fy[n] = fy[n] + recip[5] * f1 + recip[8] * f2 + recip[11] * f3;
	fz[n] = fz[n] + recip[6] * f1 + recip[9] * f2 + recip[12] * f3;
/* L400: */
    }
    return 0;
} /* grad_sum__ */

