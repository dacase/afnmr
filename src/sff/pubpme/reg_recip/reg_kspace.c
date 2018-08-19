/* reg_kspace.f -- translated by f2c (version 20061008).
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

/* Subroutine */ int reg_kspace_(int *numwats, double *x, double 
	*y, double *z__, double *cg, double *ewaldcof, double 
	*ene, double *force, double *frac, double *recip, 
	double *maxexp, int *mlimit, double *volume, double *
	box, double *xx, int *navail, double *virial)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int numatoms, i__;
    extern /* Subroutine */ int recip_reg__(int *, double *, 
	    double *, double *, double *, double *, int *,
	     int *, double *, double *, double *, double *
	    , double *, double *);
    static double vir, tvir[9]	/* was [3][3] */;

    /* Parameter adjustments */
    --x;
    --y;
    --z__;
    --cg;
    force -= 4;
    frac -= 4;
    recip -= 4;
    --mlimit;
    --xx;
    --virial;

    /* Function Body */
    numatoms = *numwats * 3;
    *ene = 0.;
    i__1 = numatoms;
    for (i__ = 1; i__ <= i__1; ++i__) {
	force[i__ * 3 + 1] = 0.;
	force[i__ * 3 + 2] = 0.;
	force[i__ * 3 + 3] = 0.;
	frac[i__ * 3 + 1] = x[i__] / *box;
	frac[i__ * 3 + 2] = y[i__] / *box;
	frac[i__ * 3 + 3] = z__[i__] / *box;
/* L50: */
    }
    recip_reg__(&numatoms, &cg[1], ene, &vir, tvir, &xx[1], navail, &mlimit[1]
	    , volume, &recip[4], &force[4], &frac[4], ewaldcof, maxexp);
    virial[1] = tvir[0];
    virial[2] = tvir[3];
    virial[3] = tvir[6];
    virial[4] = tvir[4];
    virial[5] = tvir[7];
    virial[6] = tvir[8];
    return 0;
} /* reg_kspace__ */

/* -------------------------------------------------------- */
/* Subroutine */ int recip_reg__(int *numatoms, double *charge, 
	double *eer, double *vir, double *tvir, double *x, 
	int *navail, int *mlimit, double *volume, double *
	recip, double *force, double *fraction, double *ewaldcof, 
	double *maxexp)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, j, need, lckc, lelc, lemc, lenc, lfrc, lclm, lcks, 
	    lels, lems, lens, lslm, mlimax;
    extern /* Subroutine */ int ew_ccp5__(int *, int *, int *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *);

    /* Parameter adjustments */
    fraction -= 4;
    force -= 4;
    recip -= 4;
    --mlimit;
    --x;
    tvir -= 4;
    --charge;

    /* Function Body */
    *vir = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    tvir[i__ + j * 3] = 0.;
/* L20: */
	}
/* L10: */
    }
/* .. check if the the array space in X is large eneough */
/*   to supply space for the regular ewald summations. */
/* Computing MAX */
    i__1 = max(mlimit[1],mlimit[2]);
    mlimax = max(i__1,mlimit[3]) + 1;
    need = (mlimax * 7 + 6) * *numatoms;
    if (*navail < need) {
/*        write(6,*)'Memory space in the FFT-arrays not large enough ' */
/*        write(6,*)'Provided space: ',navail,' Needed space: ',need */
/*        stop */
		fprintf( stderr, "Memory space in the FFT-arrays not large enough\n" );			exit(1);
    }
    lckc = 1;
    lcks = lckc + *numatoms;
    lclm = lcks + *numatoms;
    lslm = lclm + *numatoms;
    lelc = lslm + mlimax * *numatoms;
    lemc = lelc + mlimax * *numatoms;
    lenc = lemc + mlimax * *numatoms;
    lels = lenc + mlimax * *numatoms;
    lems = lels + mlimax * *numatoms;
    lens = lems + mlimax * *numatoms;
    lfrc = lens + mlimax * *numatoms;
    ew_ccp5__(numatoms, &mlimax, &mlimit[1], volume, &recip[4], ewaldcof, 
	    maxexp, &charge[1], &fraction[4], &force[4], &x[lckc], &x[lcks], &
	    x[lclm], &x[lslm], &x[lelc], &x[lemc], &x[lenc], &x[lels], &x[
	    lems], &x[lens], eer, vir, &tvir[4]);
    return 0;
} /* recip_reg__ */

/* Subroutine */ int ew_ccp5__(int *maxsit, int *kmax, int *
	mlimit, double *dvol1, double *hi, double *alphad, 
	double *rksmax, double *schg, double *fraction, 
	double *ffxyz, double *ckc, double *cks, double *clm, 
	double *slm, double *elc, double *emc, double *enc, 
	double *els, double *ems, double *ens, double *pe, 
	double *vir, double *tvr)
{
    /* Initialized data */

    static double cl = 1.;

    /* System generated locals */
    int elc_dim1, elc_offset, emc_dim1, emc_offset, enc_dim1, enc_offset, 
	    els_dim1, els_offset, ems_dim1, ems_offset, ens_dim1, ens_offset, 
	    i__1, i__2, i__3, i__4;
    double d__1, d__2, d__3;

    /* Builtin functions */
    double atan(double), sqrt(double), cos(double), sin(
	    double), exp(double);

    /* Local variables */
    static double rksmaxsq;
    static int i__, l, m, n;
    static double ak;
    static int ll, mm, nn;
    static double rl, rm, rn, fac, akv, rcl, qpe, omg[9];
    static int mmm, nnn;
    static double qvf, rkx1, rkx2, rky2, rkz2, rky1, rkz1, rkx3, rky3, 
	    rkz3, ckcs;
    static int mmin, nmin;
    static double ckss, sqpi, rksq, rvol, ralph;
    static int klimx, klimy, klimz;
    static double twopi;
    static int klim2y, klim2z;
    static double qforce;


/* ----------------------------------------------------------------------- */
/*     Arguments: */
/*                 maxsit   = number of sites ( atoms ) */
/*                 kmax     = maximum number of k-vectors */
/*                 dvol1    = volume of md-cell */
/*                 hi       = reciprocal cell-vectors */
/*                 alphad   = Ewald - convergence parameter */
/*                 schg     = fractional coordinates */
/*                 fraction = dimensionless ( fractional ) coordinates */
/*                 ffxyz    = forces on sites */
/*                 ckc      = work-array to store cosines ( size : maxsit ) */
/*                 cks      = work-array to store sines   ( size : maxsit ) */
/*                 clm      = work-array for cosine-products */
/*                 slm      = work-array for sine-products */
/*                 elc      = work-array ( size : (kmax+1)*maxsit ) */
/*                 emc      = work-array ( size : (kmax+1)*maxsit ) */
/*                 enc      = work-array ( size : (kmax+1)*maxsit ) */
/*                 els      = work-array ( size : (kmax+1)*maxsit ) */
/*                 ems      = work-array ( size : (kmax+1)*maxsit ) */
/*                 ens      = work-array ( size : (kmax+1)*maxsit ) */
/*                 pe       = potential energy */
/*                 vir      = virial */
/*                 tvr      = virial-tensor */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --slm;
    --clm;
    --cks;
    --ckc;
    ffxyz -= 4;
    fraction -= 4;
    --schg;
    ens_dim1 = *kmax;
    ens_offset = 1 + ens_dim1;
    ens -= ens_offset;
    ems_dim1 = *kmax;
    ems_offset = 1 + ems_dim1;
    ems -= ems_offset;
    els_dim1 = *kmax;
    els_offset = 1 + els_dim1;
    els -= els_offset;
    enc_dim1 = *kmax;
    enc_offset = 1 + enc_dim1;
    enc -= enc_offset;
    emc_dim1 = *kmax;
    emc_offset = 1 + emc_dim1;
    emc -= emc_offset;
    elc_dim1 = *kmax;
    elc_offset = 1 + elc_dim1;
    elc -= elc_offset;
    --mlimit;
    --hi;
    --tvr;

    /* Function Body */

    twopi = 8. * atan(1.);
    sqpi = sqrt(atan(1.) * 4.);
    klimx = mlimit[1] + 1;
    klimy = mlimit[2] + 1;
    klimz = mlimit[3] + 1;
    klim2y = (mlimit[2] << 1) + 1;
    klim2z = (mlimit[3] << 1) + 1;
    rcl = twopi / cl;
    rvol = twopi / *dvol1;
/* Computing 2nd power */
    d__1 = *alphad;
    ralph = -1. / (cl * 4 * cl * (d__1 * d__1));
/* Computing 2nd power */
    d__1 = *rksmax;
    rksmaxsq = twopi * twopi * (d__1 * d__1);
    fac = 2.;

/*     INITIALISE TENSORS */
    for (i__ = 1; i__ <= 9; ++i__) {
	tvr[i__] = 0.;
	omg[i__ - 1] = 0.;
/* L10: */
    }

/*     INITIALISE ACCUMULATORS */
    *pe = 0.;
    qpe = 0.;
    *vir = 0.;

/*     CALCULATE AND STORE EXPONENTIAL FACTORS */
    i__1 = *maxsit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	elc[i__ * elc_dim1 + 1] = 1.;
	emc[i__ * emc_dim1 + 1] = 1.;
	enc[i__ * enc_dim1 + 1] = 1.;
	els[i__ * els_dim1 + 1] = 0.;
	ems[i__ * ems_dim1 + 1] = 0.;
	ens[i__ * ens_dim1 + 1] = 0.;
	elc[i__ * elc_dim1 + 2] = cos(rcl * fraction[i__ * 3 + 1]);
	emc[i__ * emc_dim1 + 2] = cos(rcl * fraction[i__ * 3 + 2]);
	enc[i__ * enc_dim1 + 2] = cos(rcl * fraction[i__ * 3 + 3]);
	els[i__ * els_dim1 + 2] = sin(rcl * fraction[i__ * 3 + 1]);
	ems[i__ * ems_dim1 + 2] = sin(rcl * fraction[i__ * 3 + 2]);
	ens[i__ * ens_dim1 + 2] = sin(rcl * fraction[i__ * 3 + 3]);
/* L100: */
    }
    i__1 = *kmax;
    for (l = 3; l <= i__1; ++l) {
	i__2 = *maxsit;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    elc[l + i__ * elc_dim1] = elc[l - 1 + i__ * elc_dim1] * elc[i__ * 
		    elc_dim1 + 2] - els[l - 1 + i__ * els_dim1] * els[i__ * 
		    els_dim1 + 2];
	    emc[l + i__ * emc_dim1] = emc[l - 1 + i__ * emc_dim1] * emc[i__ * 
		    emc_dim1 + 2] - ems[l - 1 + i__ * ems_dim1] * ems[i__ * 
		    ems_dim1 + 2];
	    enc[l + i__ * enc_dim1] = enc[l - 1 + i__ * enc_dim1] * enc[i__ * 
		    enc_dim1 + 2] - ens[l - 1 + i__ * ens_dim1] * ens[i__ * 
		    ens_dim1 + 2];
	    els[l + i__ * els_dim1] = els[l - 1 + i__ * els_dim1] * elc[i__ * 
		    elc_dim1 + 2] + elc[l - 1 + i__ * elc_dim1] * els[i__ * 
		    els_dim1 + 2];
	    ems[l + i__ * ems_dim1] = ems[l - 1 + i__ * ems_dim1] * emc[i__ * 
		    emc_dim1 + 2] + emc[l - 1 + i__ * emc_dim1] * ems[i__ * 
		    ems_dim1 + 2];
	    ens[l + i__ * ens_dim1] = ens[l - 1 + i__ * ens_dim1] * enc[i__ * 
		    enc_dim1 + 2] + enc[l - 1 + i__ * enc_dim1] * ens[i__ * 
		    ens_dim1 + 2];
/* L130: */
	}
/* L140: */
    }

/*     LOOP OVER ALL K VECTORS  K=2PI(LL/CL,MM/CL,NN/CL) */
    mmin = klimy;
    nmin = klimz + 1;
    i__1 = klimx;
    for (l = 1; l <= i__1; ++l) {
	ll = l - 1;
	rl = rcl * (double) ll;
	rkx1 = rl * hi[1];
	rky1 = rl * hi[4];
	rkz1 = rl * hi[7];
	i__2 = klim2y;
	for (mmm = mmin; mmm <= i__2; ++mmm) {
	    mm = mmm - klimy;
	    m = abs(mm) + 1;
	    rm = rcl * (double) mm;
	    rkx2 = rkx1 + rm * hi[2];
	    rky2 = rky1 + rm * hi[5];
	    rkz2 = rkz1 + rm * hi[8];
/*     SET TEMPORARY PRODUCTS OF EXPONENTIAL TERMS */
	    if (mm >= 0) {
		i__3 = *maxsit;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    clm[i__] = elc[l + i__ * elc_dim1] * emc[m + i__ * 
			    emc_dim1] - els[l + i__ * els_dim1] * ems[m + i__ 
			    * ems_dim1];
		    slm[i__] = els[l + i__ * els_dim1] * emc[m + i__ * 
			    emc_dim1] + ems[m + i__ * ems_dim1] * elc[l + i__ 
			    * elc_dim1];
/* L150: */
		}
	    } else {
		i__3 = *maxsit;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    clm[i__] = elc[l + i__ * elc_dim1] * emc[m + i__ * 
			    emc_dim1] + els[l + i__ * els_dim1] * ems[m + i__ 
			    * ems_dim1];
		    slm[i__] = els[l + i__ * els_dim1] * emc[m + i__ * 
			    emc_dim1] - ems[m + i__ * ems_dim1] * elc[l + i__ 
			    * elc_dim1];
/* L160: */
		}
	    }
	    i__3 = klim2z;
	    for (nnn = nmin; nnn <= i__3; ++nnn) {
		nn = nnn - klimz;
		n = abs(nn) + 1;
		rn = rcl * (double) nn;
		rkx3 = rkx2 + rn * hi[3];
		rky3 = rky2 + rn * hi[6];
		rkz3 = rkz2 + rn * hi[9];
/*     CALCULATE AK COEFFICIENTS */
/* Computing 2nd power */
		d__1 = rkx3;
/* Computing 2nd power */
		d__2 = rky3;
/* Computing 2nd power */
		d__3 = rkz3;
		rksq = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;

/*     BYPASS K VECTORS OUTSIDE CUTOFF */
		if (rksq <= rksmaxsq) {
		    ak = exp(ralph * rksq) / rksq;
		    akv = ak * 2. * (1. / rksq - ralph);
/*     CALCULATE EXP(IKR) TERMS AND PRODUCT WITH SITE CHARGES */
		    if (nn >= 0) {
			i__4 = *maxsit;
			for (i__ = 1; i__ <= i__4; ++i__) {
			    ckc[i__] = schg[i__] * (clm[i__] * enc[n + i__ * 
				    enc_dim1] - slm[i__] * ens[n + i__ * 
				    ens_dim1]);
			    cks[i__] = schg[i__] * (slm[i__] * enc[n + i__ * 
				    enc_dim1] + clm[i__] * ens[n + i__ * 
				    ens_dim1]);
/* L170: */
			}
		    } else {
			i__4 = *maxsit;
			for (i__ = 1; i__ <= i__4; ++i__) {
			    ckc[i__] = schg[i__] * (clm[i__] * enc[n + i__ * 
				    enc_dim1] + slm[i__] * ens[n + i__ * 
				    ens_dim1]);
			    cks[i__] = schg[i__] * (slm[i__] * enc[n + i__ * 
				    enc_dim1] - clm[i__] * ens[n + i__ * 
				    ens_dim1]);
/* L180: */
			}
		    }
/*     CALCULATE VECTOR SUMS */
		    ckcs = 0.;
		    ckss = 0.;
		    i__4 = *maxsit;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			ckcs += ckc[i__];
			ckss += cks[i__];
/* L190: */
		    }
/*     ACCUMULATE POTENTIAL ENERGY AND VIRIAL */
		    qpe += ak * (ckcs * ckcs + ckss * ckss);
		    qvf = akv * (ckcs * ckcs + ckss * ckss);
		    omg[0] -= qvf * rkx3 * rkx3;
		    omg[1] -= qvf * rkx3 * rky3;
		    omg[2] -= qvf * rkx3 * rkz3;
		    omg[3] -= qvf * rky3 * rkx3;
		    omg[4] -= qvf * rky3 * rky3;
		    omg[5] -= qvf * rky3 * rkz3;
		    omg[6] -= qvf * rkz3 * rkx3;
		    omg[7] -= qvf * rkz3 * rky3;
		    omg[8] -= qvf * rkz3 * rkz3;
/*     CALCULATE FORCE ON EACH SITE */
		    i__4 = *maxsit;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			qforce = ak * (cks[i__] * ckcs - ckc[i__] * ckss);
			ffxyz[i__ * 3 + 1] += rl * qforce;
			ffxyz[i__ * 3 + 2] += rm * qforce;
			ffxyz[i__ * 3 + 3] += rn * qforce;
/* L200: */
		    }

/*     END VECTOR LOOP */
		}
/* L210: */
	    }
	    nmin = 1;
/* L220: */
	}
	mmin = 1;
/* L230: */
    }

/*     CALCULATE EWALD FORCE ARRAYS */
    i__1 = *maxsit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ffxyz[i__ * 3 + 1] = fac * 2. * rvol * (hi[1] * ffxyz[i__ * 3 + 1] + 
		hi[2] * ffxyz[i__ * 3 + 2] + hi[3] * ffxyz[i__ * 3 + 3]);
	ffxyz[i__ * 3 + 2] = fac * 2. * rvol * (hi[4] * ffxyz[i__ * 3 + 1] + 
		hi[5] * ffxyz[i__ * 3 + 2] + hi[6] * ffxyz[i__ * 3 + 3]);
	ffxyz[i__ * 3 + 3] = fac * 2. * rvol * (hi[7] * ffxyz[i__ * 3 + 1] + 
		hi[8] * ffxyz[i__ * 3 + 2] + hi[9] * ffxyz[i__ * 3 + 3]);
/* L240: */
    }

/*     CALCULATE FINAL CORRECTED POTENTIAL */
    *pe = rvol * qpe * fac;

/*     CALCULATE FINAL VIRIAL TENSOR */
    tvr[1] = -rvol * fac * (qpe + omg[0]);
    tvr[2] = -rvol * fac * omg[1];
    tvr[3] = -rvol * fac * omg[2];
    tvr[4] = -rvol * fac * omg[3];
    tvr[5] = -rvol * fac * (qpe + omg[4]);
    tvr[6] = -rvol * fac * omg[5];
    tvr[7] = -rvol * fac * omg[6];
    tvr[8] = -rvol * fac * omg[7];
    tvr[9] = -rvol * fac * (qpe + omg[8]);

    *vir = tvr[1] + tvr[5] + tvr[9];

    return 0;
} /* ew_ccp5__ */

