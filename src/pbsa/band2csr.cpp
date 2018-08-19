/*
 * Convert banded matrix, typically from PB fortran code, into CSR format to use with CUDA C conjugate gradient code.
 * Ruxi Qi, UC Irvine
 *
 * 1st version
 * Normal matrix
 * Jun 27-30, 2016
 *
 * 2nd version
 * PBC matrix
 * Oct 3-12 2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "band2csr.h"

/********************
 * Two parts:
 * Part 1 - non-PBC
 * Part 2 - PBC
 ********************/

/*************** Part 1. non-PBC *****************/
void band2csr(int *I, int *J, float *val, int N, int nz, float **band, int bandwidth, int xm, int ym, int zm)
{
	// Check input
	assert(N == xm * ym * zm);
	int xmym = xm * ym;
	int upperside = (int)((bandwidth + 1) / 2);

	for (int i = 1; i < 4; i++) {
		for (int j = 0; j < N; j++) {
			// From Prof. Luo, band data do not contain signs!
			band[i][j] = - band[i][j];
		}
	}

	/* Now converting banded matrix  into CSR format - I, J and val are nonzero row, column and value arraies in CSR */
	// Separate i into 7 ranges, handle each at a time. Generate I[] first, so that it can be used later for generating J[] and val[]
	// Use offset[] array for lower corner val[] assignment (i.e. extracting band[1][] ~ band[3][] element correctly)
	int offset[3] = {1, xm, xmym}; // For b1, b2 & b3 respectively

	// I
	I[0] = 0;
	I[1] = upperside;
	for (int i = 2; i < N; i++) {
		// Head
		if (i <= xm) {
			I[i] = I[i-1] + upperside + 1;
		}

		if (i > xm && i <= xmym) {
			I[i] = I[i-1] + upperside + 2;
		}

		// Body
		if (i > xmym && i <= N - xmym) {
			I[i] = I[i-1] + upperside + 3;
		}

		// Tail
		if (i > N - xmym && i <= N - xm) {
			I[i] = I[i-1] + upperside + 2;
		}

		if (i > N - xm && i <= N - 1) {
			I[i] = I[i-1] + upperside + 1;
		}
	}
	I[N] = I[N-1] + upperside;

	// Check with input
	assert(I[N] == nz);

	// J & val
	for (int i = 0; i < N; i++) {
		// Head
		if (i == 0) {
			// J
			J[I[i]	  ] = i;
			J[I[i] + 1] = i + 1;
			J[I[i] + 2] = i + xm;
			J[I[i] + 3] = i + xmym;

			// val
			for (int j = 0; j < upperside; j++) {
				val[j] = band[j][0];
			}
		}
		if (i >= 1 && i < xm) {
			// J
			J[I[i]	  ] = i - 1;
			J[I[i] + 1] = i;
			J[I[i] + 2] = i + 1;
			J[I[i] + 3] = i + xm;
			J[I[i] + 4] = i + xmym;

			// val
			for (int j = 0; j < bandwidth - 2; j++) {
				val[I[i] + j] = (j < upperside - 3) ? band[upperside - 3 - j][i - offset[upperside - 3 - j - 1]] : band[j - upperside + 3][i];
			}
		}
		if (i >= xm && i < xmym) {
			// J
			J[I[i]	  ] = i - xm;
			J[I[i] + 1] = i - 1;
			J[I[i] + 2] = i;
			J[I[i] + 3] = i + 1;
			J[I[i] + 4] = i + xm;
			J[I[i] + 5] = i + xmym;

			// val
			for (int j = 0; j < bandwidth - 1; j++) {
				val[I[i] + j] = (j < upperside - 2) ? band[upperside - 2 - j][i - offset[upperside - 2 - j - 1]] : band[j - upperside + 2][i];
			}
		}

		// Body
		if (i >= xmym && i <= N - xmym - 1) {
			// J
			J[I[i]    ] = i - xmym;
			J[I[i] + 1] = i - xm;
			J[I[i] + 2] = i - 1;
			J[I[i] + 3] = i;
			J[I[i] + 4] = i + 1;
			J[I[i] + 5] = i + xm;
			J[I[i] + 6] = i + xmym;

			// val
			for (int j = 0; j < bandwidth; j++) {
				val[I[i] + j] = (j < upperside - 1) ? band[upperside - 1 - j][i - offset[upperside - 1 - j - 1]] : band[j - upperside + 1][i];
			}
		}

		// Tail
		if (i > N - xmym - 1 && i <= N - xm - 1) {
			// J
			J[I[i]	  ] = i - xmym;
			J[I[i] + 1] = i - xm;
			J[I[i] + 2] = i - 1;
			J[I[i] + 3] = i;
			J[I[i] + 4] = i + 1;
			J[I[i] + 5] = i + xm;

			// val
			for (int j = 0; j < bandwidth - 1; j++) {
				val[I[i] + j] = (j < upperside - 1) ? band[upperside - 1 - j][i - offset[upperside - 1 - j - 1]] : band[j - upperside + 1][i];
			}
		}
		if (i > N - xm - 1 && i <= N - 2) {
			// J
			J[I[i]	  ] = i - xmym;
			J[I[i] + 1] = i - xm;
			J[I[i] + 2] = i - 1;
			J[I[i] + 3] = i;
			J[I[i] + 4] = i + 1;

			// val
			for (int j = 0; j < bandwidth - 2; j++) {
				val[I[i] + j] = (j < upperside - 1) ? band[upperside - 1 - j][i - offset[upperside - 1 - j - 1]] : band[j - upperside + 1][i];
			}
		}
		if (i == N - 1) {
			// J
			J[I[i]	  ] = i - xmym;
			J[I[i] + 1] = i - xm;
			J[I[i] + 2] = i - 1;
			J[I[i] + 3] = i;

			// val
			for (int j = 0; j < upperside; j++) {
				val[I[i] + j] = (j < upperside - 1) ? band[upperside - 1 - j][i - offset[upperside - 1 - j - 1]] : band[j - upperside + 1][i];
			}
		}
	}
}

/*************** Part 2. PBC *****************/

void band2csr_pbc(int *I, int *J, float *val, int N, int nz, float **band, int bandwidth, int xm, int ym, int zm)
{
	// Check input
	assert(N == xm * ym * zm);
	int xmym = xm * ym;
	int ii, jj;

	for (int i = 1; i < 7; i++) {
		for (int j = 0; j < N; j++) {
			// From Prof. Luo, band data do not contain signs!
			band[i][j] = - band[i][j];
		}
	}

	// Algorithm:
	//     New PBC version with band 0 ~ 6. The algorithm used here is quite different from its non-PBC band2csr counterpart. I
	// treat the matrix as constant 7-banded, put new surface mirror points into their original sites as in old 7 bands.
	// Also I handle [ 6 (corresponding to six surfaces of grid cube) + 1 ] original bands, check each bands with/without surface
	// points, and feed them J[] and val[] array values accordingly. For easy-coding, I use non-ordered row storage in CSR.
	// The new algorithm is more neat and bypasses the complexity. Ask Ruxi for detail, if he still remembers this.

	/* Now converting banded matrix into CSR format - I, J and val are nonzero row, column and value arraies in CSR */

	// Loop over three indices - keep in mind the big grid cube and Left-hand coordination system
	for (int k = 0; k < zm; k++) {
		for (int j = 0; j < ym; j++) {
			for (int i = 0; i < xm; i++) {
				// Main grid site index
				ii = i + j * xm + k * xmym;

				// Assign array I
				I[ii] = bandwidth * ii;

				// First band, top z surface
				if (k == 0) {
					jj = ii + (zm -1) * xmym;
					J[I[ii]] = jj;
					val[I[ii]] = band[6][ii];
				} else {
					J[I[ii]] = ii - xmym;
					val[I[ii]] = band[3][ii - xmym];
				}

				// Second band, left y surface
				if (j == 0) {
					jj = ii + (ym -1) * xm;
					J[I[ii] + 1] = jj;
					val[I[ii] + 1] = band[5][ii];
				} else {
					J[I[ii] + 1] = ii - xm;
					val[I[ii] + 1] = band[2][ii - xm];
				}

				// Third band, left x surface
				if (i == 0) {
					jj = ii + xm - 1;
					J[I[ii] + 2] = jj;
					val[I[ii] + 2] = band[4][ii];
				} else {
					J[I[ii] + 2] = ii - 1;
					val[I[ii] + 2] = band[1][ii - 1];
				}

				// Fourth band, non-surface
				J[I[ii] + 3] = ii;
				val[I[ii] + 3] = band[0][ii];

				// Fifth band, right x surface
				if (i == xm - 1) {
					jj = ii - xm + 1;
					J[I[ii] + 4] = jj;
					val[I[ii] + 4] = band[4][jj];
				} else {
					J[I[ii] + 4] = ii + 1;
					val[I[ii] + 4] = band[1][ii];
				}

				// Sixth band, right y surface
				//if (i == ym - 1) { Bug located! Heaven! Alarm here for three days!
				if (j == ym - 1) {
					jj = ii - (ym - 1) * xm;
					J[I[ii] + 5] = jj;
					val[I[ii] + 5] = band[5][jj];
				} else {
					J[I[ii] + 5] = ii + xm;
					val[I[ii] + 5] = band[2][ii];
				}

				// Seventh band, bottom z surface
				if (k == zm - 1) {
					jj = ii - (zm - 1) * xmym;
					J[I[ii] + 6] = jj;
					val[I[ii] + 6] = band[6][jj];
				} else {
					J[I[ii] + 6] = ii + xmym;
					val[I[ii] + 6] = band[3][ii];
				}
			}
		}
	}
	I[N] = bandwidth * N;
}
