/* fftcalls.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#define max(a,b) ((a) >= (b) ? (a) : (b))

/* Subroutine */ int get_fftdims__(int *nfft1, int *nfft2, int *
	nfft3, int *nfftdim1, int *nfftdim2, int *nfftdim3, 
	int *nfftable, int *nffwork, int *sizfftab, int *
	sizffwrk)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int n, nfftmax;

/* Computing MAX */
    i__1 = max(*nfft1,*nfft2);
    nfftmax = max(i__1,*nfft3);
    *nfftdim1 = *nfft1;
    n = *nfft1 / 2;
    if (*nfft1 == n << 1) {
	*nfftdim1 = *nfft1 + 1;
    }
    *nfftdim2 = *nfft2;
    n = *nfft2 / 2;
    if (*nfft2 == n << 1) {
	*nfftdim2 = *nfft2 + 1;
    }
    *nfftdim3 = *nfft3;
    n = *nfft3 / 2;
    if (*nfft3 == n << 1) {
	*nfftdim3 = *nfft3 + 1;
    }
    *nfftable = (nfftmax << 2) + 15;
    *nffwork = nfftmax;
    *sizfftab = *nfftable * 3;
    *sizffwrk = nfftmax << 1;
    return 0;
} /* get_fftdims__ */

/* --------------------------------------------------------------- */
/* Subroutine */ int fft_setup__(double *array, double *fftable, 
	double *ffwork, int *nfft1, int *nfft2, int *nfft3, 
	int *nfftdim1, int *nfftdim2, int *nfftdim3, int *
	nfftable, int *nffwork)
{
    extern /* Subroutine */ int pubz3di_(int *, int *, int *, 
	    double *, int *);

/*     write(6,*)'using public domain fft code' */
    /* Parameter adjustments */
    --ffwork;
    --fftable;
    --array;

    /* Function Body */
    pubz3di_(nfft1, nfft2, nfft3, &fftable[1], nfftable);
    return 0;
} /* fft_setup__ */

/* ----------------------------------------------------------- */
/* Subroutine */ int fft_forward__(double *array, double *fftable, 
	double *ffwork, int *nfft1, int *nfft2, int *nfft3, 
	int *nfftdim1, int *nfftdim2, int *nfftdim3, int *
	nfftable, int *nffwork)
{
    static int isign;
    extern /* Subroutine */ int pubz3d_(int *, int *, int *, 
	    int *, double *, int *, int *, double *, 
	    int *, double *, int *);

    /* Parameter adjustments */
    --ffwork;
    --fftable;
    --array;

    /* Function Body */
    isign = 1;
    pubz3d_(&isign, nfft1, nfft2, nfft3, &array[1], nfftdim1, nfftdim2, &
	    fftable[1], nfftable, &ffwork[1], nffwork);
    return 0;
} /* fft_forward__ */

/* ----------------------------------------------------------- */
/* Subroutine */ int fft_back__(double *array, double *fftable, 
	double *ffwork, int *nfft1, int *nfft2, int *nfft3, 
	int *nfftdim1, int *nfftdim2, int *nfftdim3, int *
	nfftable, int *nffwork)
{
    static int isign;
    extern /* Subroutine */ int pubz3d_(int *, int *, int *, 
	    int *, double *, int *, int *, double *, 
	    int *, double *, int *);

    /* Parameter adjustments */
    --ffwork;
    --fftable;
    --array;

    /* Function Body */
    isign = -1;
    pubz3d_(&isign, nfft1, nfft2, nfft3, &array[1], nfftdim1, nfftdim2, &
	    fftable[1], nfftable, &ffwork[1], nffwork);
    return 0;
} /* fft_back__ */

