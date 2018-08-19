/* pub3dfft.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/


typedef struct { double r, i; } doublecomplex;

/* **************************************************************************** */

/* 	3D (slow) Fourier Transform */
/*   this 1d->3d code is brute force approach */
/*   the 1d code is a double precision version of fftpack from netlib */
/*   due to Paul N Swartztrauber at NCAR Boulder Coloraso */

/* **************************************************************************** */
/* Subroutine */ int pubz3di_(int *n1, int *n2, int *n3, 
	double *table, int *ntable)
{
    /* System generated locals */
    int table_dim1, table_offset;

    /* Local variables */
    extern /* Subroutine */ int cffti_(int *, double *);

/* ntable should be 4*max(n1,n2,n3) +15 */
    /* Parameter adjustments */
    table_dim1 = *ntable;
    table_offset = 1 + table_dim1;
    table -= table_offset;

    /* Function Body */
    cffti_(n1, &table[table_dim1 + 1]);
    cffti_(n2, &table[(table_dim1 << 1) + 1]);
    cffti_(n3, &table[table_dim1 * 3 + 1]);
    return 0;
} /* pubz3di_ */

/* **************************************************************************** */
/* Subroutine */ int pubz3d_(int *isign, int *n1, int *n2, 
	int *n3, doublecomplex *w, int *ld1, int *ld2, double 
	*table, int *ntable, doublecomplex *work, int *nwork)
{
    /* System generated locals */
    int w_dim1, w_dim2, w_offset, table_dim1, table_offset, i__1, i__2, 
	    i__3, i__4, i__5;

    /* Local variables */
    static int i__, j, k;
    extern /* Subroutine */ int cfftb_(int *, doublecomplex *, double 
	    *), cfftf_(int *, doublecomplex *, double *);

/* ntable should be 4*max(n1,n2,n3) +15 */
/* nwork should be max(n1,n2,n3) */

/*   transform along X  first ... */

    /* Parameter adjustments */
    w_dim1 = *ld1;
    w_dim2 = *ld2;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    table_dim1 = *ntable;
    table_offset = 1 + table_dim1;
    table -= table_offset;
    --work;

    /* Function Body */
    i__1 = *n3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *n2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *n1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = i__;
		i__5 = i__ + (j + k * w_dim2) * w_dim1;
		work[i__4].r = w[i__5].r, work[i__4].i = w[i__5].i;
/* L70: */
	    }
	    if (*isign == -1) {
		cfftf_(n1, &work[1], &table[table_dim1 + 1]);
	    }
	    if (*isign == 1) {
		cfftb_(n1, &work[1], &table[table_dim1 + 1]);
	    }
	    i__3 = *n1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = i__ + (j + k * w_dim2) * w_dim1;
		i__5 = i__;
		w[i__4].r = work[i__5].r, w[i__4].i = work[i__5].i;
/* L80: */
	    }
/* L90: */
	}
/* L100: */
    }

/*   transform along Y then ... */

    i__1 = *n3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *n1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *n2;
	    for (j = 1; j <= i__3; ++j) {
		i__4 = j;
		i__5 = i__ + (j + k * w_dim2) * w_dim1;
		work[i__4].r = w[i__5].r, work[i__4].i = w[i__5].i;
/* L170: */
	    }
	    if (*isign == -1) {
		cfftf_(n2, &work[1], &table[(table_dim1 << 1) + 1]);
	    }
	    if (*isign == 1) {
		cfftb_(n2, &work[1], &table[(table_dim1 << 1) + 1]);
	    }
	    i__3 = *n2;
	    for (j = 1; j <= i__3; ++j) {
		i__4 = i__ + (j + k * w_dim2) * w_dim1;
		i__5 = j;
		w[i__4].r = work[i__5].r, w[i__4].i = work[i__5].i;
/* L180: */
	    }
/* L190: */
	}
/* L200: */
    }

/*   transform along Z finally ... */

    i__1 = *n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *n3;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = k;
		i__5 = i__ + (j + k * w_dim2) * w_dim1;
		work[i__4].r = w[i__5].r, work[i__4].i = w[i__5].i;
/* L270: */
	    }
	    if (*isign == -1) {
		cfftf_(n3, &work[1], &table[table_dim1 * 3 + 1]);
	    }
	    if (*isign == 1) {
		cfftb_(n3, &work[1], &table[table_dim1 * 3 + 1]);
	    }
	    i__3 = *n3;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = i__ + (j + k * w_dim2) * w_dim1;
		i__5 = k;
		w[i__4].r = work[i__5].r, w[i__4].i = work[i__5].i;
/* L280: */
	    }
/* L290: */
	}
/* L300: */
    }
    return 0;
} /* pubz3d_ */

