#ifndef BAND2CSR_H_
#define BAND2CSR_H_

/* This function converts banded matrix into CSR format */
void band2csr(int *I, int *J, float *val, int N, int nz, float **band, int bandwidth, int xm, int ym, int zm);
void band2csr_pbc(int *I, int *J, float *val, int N, int nz, float **band, int bandwidth, int xm, int ym, int zm);

#endif
