/*
 * CUDA MG solver with Unified Memory
 * Ruxi Qi @ UC Irvine, 2017
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
// For gdb
#include <signal.h>
// For timing
#include <sys/time.h>
// CUDA Runtime
#include <cuda_runtime.h>
// CUDA BLAS Library
#include "cublas_v2.h"
// For error handling and device pickup
#include "helper_cuda.h"
#include "cuda_mg_wrapper.h"

// Global vairables
int l_xm;
int l_ym;
int l_zm;
int l_xmym;
int l_xmymzm;
int l_maxitn;
// For initialization
int l_l;
int l_m;
int l_n;
__device__ __managed__ int l_bcopt;
float l_accept;
//int mg_nlevel;
int ncyc_before;
int ncyc_after;
float l_pbkappa;
float l_epsout;
float l_h;
float l_wsor;
int l_itn;
float l_inorm;
float l_norm;

int threshold;

#define MG_NLEVEL 4
int mg_index[MG_NLEVEL + 1];
int mg_index_ext[MG_NLEVEL + 1];
int mg_x_idx[MG_NLEVEL + 1];
int mg_size[MG_NLEVEL][3];
float mg_onorm[MG_NLEVEL];

float *l_zv;
float *l_ad;
float *l_bv;
float *l_rv;
float *l_iv;
float *l_bz;
float *l_am1;
float *l_am2;
float *l_am3;
float *l_xv;

int devThreadsPerBlock = 8;

extern "C"
void init_param_c_(int *nx, int *ny, int *nz, int *p_maxitn, int *p_bcopt, float *p_accept, float *p_pbkappa, float *p_epsout, float *p_h, float *p_wsor) {
    l_xm = *nx;
    l_ym = *ny;
    l_zm = *nz;
    l_xmym = *nx * *ny;
    l_xmymzm = *nx * *ny * *nz;
    l_maxitn = *p_maxitn;
    l_bcopt = *p_bcopt;
    l_accept = *p_accept;
    ncyc_before = 10;
    ncyc_after = 10;
    l_pbkappa = *p_pbkappa;
    l_epsout = *p_epsout;
    l_h = *p_h;
    l_wsor = *p_wsor;

    threshold = 2;
}

extern "C"
void allocate_array_cuda_(int *solvopt) {
    if (!(*solvopt == 2 || *solvopt == 4)) {
        printf("Error: Only MG/SOR is supported now.\n");
        exit(2);
    }

    int m, l, n;

    // set indices for the finest level for all solvers
    mg_index_ext[0] = 0;
    mg_index[0] = 0;
    mg_x_idx[0] = 0;
    mg_size[0][0] = l_xm;
    mg_size[0][1] = l_ym;
    mg_size[0][2] = l_zm;
    m = l_xmymzm;
    l = m + l_xmym;
    n = l + l_xmym;

    // set indices for all other levels for MG only
    if (*solvopt == 2) {
        for (int i = 1; i < MG_NLEVEL; i++) {
            mg_index_ext[i] = l;
            mg_index[i] = m;
            mg_x_idx[i] = n;

            //l_bcopt != 10 for now
            for (int j = 0; j < 3; j++) {
                mg_size[i][j] = mg_size[i - 1][j] / 2;
            }
            m += mg_size[i][0] * mg_size[i][1] * mg_size[i][2];
            l += mg_size[i][0] * mg_size[i][1] * mg_size[i][2] + mg_size[i][0] * mg_size[i][1];
            n += mg_size[i][0] * mg_size[i][1] * mg_size[i][2] + 2 * mg_size[i][0] * mg_size[i][1];
        }

        mg_index_ext[MG_NLEVEL] = l;
        mg_index[MG_NLEVEL] = m;
        mg_x_idx[MG_NLEVEL] = n;
    }

    // Now for all arrays
    // Try __managed__ declaration later for performance tuning
    //__device__ __managed__ l_xv[n];
    // Note in Fortran these arrays index from 1, not 1-xmym etc.
    cudaErrorCheck(cudaMallocManaged(&l_zv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_ad, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_bv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_rv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_iv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_bz, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_am1, sizeof(float) * l));
    cudaErrorCheck(cudaMallocManaged(&l_am2, sizeof(float) * l));
    cudaErrorCheck(cudaMallocManaged(&l_am3, sizeof(float) * l));
    cudaErrorCheck(cudaMallocManaged(&l_xv, sizeof(float) * n));

    l_l = l;
    l_m = m;
    l_n = n;
}

extern "C"
void deallocate_array_cuda_() {
    cudaFree(l_zv);
    cudaFree(l_ad);
    cudaFree(l_bv);
    cudaFree(l_rv);
    cudaFree(l_iv);
    cudaFree(l_bz);
    cudaFree(l_am1);
    cudaFree(l_am2);
    cudaFree(l_am3);
    cudaFree(l_xv);
    cudaDeviceReset();
}

extern "C"
void init_array_cuda_(int *solvopt, float *epsx, float *epsy, float *epsz, float *p_bv, float *p_iv, float *p_xs) {
    if (!(*solvopt == 2 || *solvopt == 4) ) {
        printf("Error: Only MG/SOR is supported now.\n");
        exit(2);
    }

    // Initialize arrays l_ad, l_am*, l_*v to 0 on device. Use 1D thread block.
    int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;

    // m
    int nblocks = (l_m - 1) / blocksize + 1;
    init_vector_kernel<<<nblocks, blocksize>>>(l_zv, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_ad, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_bv, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_rv, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_iv, l_m);
    cudaLaunchErrorCheck();

    // l
    nblocks = (l_l - 1) / blocksize + 1;
    init_vector_kernel<<<nblocks, blocksize>>>(l_am1, l_l);
    init_vector_kernel<<<nblocks, blocksize>>>(l_am2, l_l);
    init_vector_kernel<<<nblocks, blocksize>>>(l_am3, l_l);
    cudaLaunchErrorCheck();

    // n
    nblocks = (l_n - 1) / blocksize + 1;
    init_vector_kernel<<<nblocks, blocksize>>>(l_xv, l_n);
    cudaLaunchErrorCheck();

    // Ready to array assignment: do the easiest arrays first, directly copying from the caller
    // starting index passed from F90: l_xv(1:n) which covers 2*xmym buffer; p_xs(1:l_xmymzm) just core part
    // so l_xv(1+xmym:l_xmymzm+xmym), p_xs(1:l_xmymzm)
    // Eidt: copy to device memory to limit CPU page faults
    float *lp_bv, *lp_iv, *lp_xs;
    cudaErrorCheck(cudaMalloc(&lp_bv, sizeof(float) * l_xmymzm));
    cudaErrorCheck(cudaMalloc(&lp_iv, sizeof(float) * l_xmymzm));
    cudaErrorCheck(cudaMalloc(&lp_xs, sizeof(float) * l_xmymzm));
    cudaMemcpy(lp_bv, p_bv, l_xmymzm * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lp_iv, p_iv, l_xmymzm * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lp_xs, p_xs, l_xmymzm * sizeof(float), cudaMemcpyHostToDevice);

    nblocks = (l_xmymzm  - 1) / blocksize + 1;
    copy_vector_kernel<<<nblocks, blocksize>>>(l_xv + l_xmym, lp_xs, l_xmymzm);
    cudaLaunchErrorCheck();

    copy_vector_kernel<<<nblocks, blocksize>>>(l_bv, lp_bv, l_xmymzm);
    cudaLaunchErrorCheck();

    // Set up local eps arrays for data assignment on the kernel
    int m = 0;
    for (int i = 0; i < MG_NLEVEL; i++) {
        m += (mg_size[i][0] + 1) * (mg_size[i][1] + 1) * (mg_size[i][2] + 1);
    }

    float *lepsx, *lepsy, *lepsz;
    cudaErrorCheck(cudaMallocManaged(&lepsx, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&lepsy, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&lepsz, sizeof(float) * m));
    float *epsx_f, *epsy_f, *epsz_f;
    cudaErrorCheck(cudaMalloc(&epsx_f, (l_xmymzm + l_ym * l_zm) * sizeof(float)));
    cudaErrorCheck(cudaMalloc(&epsy_f, (l_xmymzm + l_xm * l_zm) * sizeof(float)));
    cudaErrorCheck(cudaMalloc(&epsz_f, (l_xmymzm + l_xm * l_ym) * sizeof(float)));

    // Copy passed array to UM
    // Wanring: This can make managed epsx_f[i] accessable from CPU, but not from kernel, which will cause
    // CUDA_EXCEPTION_1/14 errors in the kernel.
    //epsx_f = epsx;
    //epsy_f = epsy;
    //epsz_f = epsz;
    cudaMemcpy(epsx_f, epsx, (l_xmymzm + l_ym * l_zm) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(epsy_f, epsy, (l_xmymzm + l_xm * l_zm) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(epsz_f, epsz, (l_xmymzm + l_xm * l_ym) * sizeof(float), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(devThreadsPerBlock, devThreadsPerBlock, devThreadsPerBlock);
    dim3 blocks((l_xm - 1)/devThreadsPerBlock + 1, (l_ym - 1) / devThreadsPerBlock + 1, (l_zm - 1) / devThreadsPerBlock + 1);
    feedepsintoam_kernel<<<blocks, threadsPerBlock>>>(l_xm, l_ym, l_zm, lepsx, lepsy, lepsz, epsx_f, epsy_f, epsz_f);
    cudaLaunchErrorCheck();

    // Next the salt term
    copy_vector_kernel<<<nblocks, blocksize>>>(l_iv, lp_iv, l_xmymzm);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();

    float lfactor = l_epsout * (l_h * l_pbkappa) * (l_h * l_pbkappa);

    // Finally we are ready to set up the A matrix
    // set up am/ad arrays at the finest level for all solvers
    // so only 1_xmymzm elements of leps* are initialized
    int j = 0;
    m = mg_index[j]; // m == 0 here
    int n = mg_index_ext[j];
    int lxmym = mg_size[j][0] * mg_size[j][1];
    // 1-D grid
    int lxmymzm = mg_size[j][0] * mg_size[j][1] * mg_size[j][2];
    dim3 h_threadsPerBlock (512);
    dim3 h_blocks((lxmymzm - 1) / 512 + 1);
    set_am_ad_kernel_head<<<h_blocks, h_threadsPerBlock>>>(lepsx + m, lepsy + m, lepsz + m, l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, lxmymzm);
    cudaLaunchErrorCheck();
    // 3-D grid
    set_am_ad_kernel_body<<<blocks, threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, l_ad + m, l_bz + m, l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2], lfactor, l_epsout);
    cudaLaunchErrorCheck();
    // 2-D grid
    dim3 t_threadsPerBlock(16, 16);
    dim3 t_blocks((max(mg_size[j][0], mg_size[j][1]) - 1) / 16 + 1, (max(mg_size[j][1], mg_size[j][2]) - 1) / 16 + 1);
    set_am_ad_kernel_tail<<<t_blocks, t_threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
    cudaLaunchErrorCheck();

    cudaDeviceSynchronize();

    if (*solvopt == 2) {
        for (j = 1; j < MG_NLEVEL; j++) {
            int l = mg_index[j-1];
            m = mg_index[j];
            n = mg_index_ext[j];

            lfactor *= 4;
            lxmym = mg_size[j][0] * mg_size[j][1];
            lxmymzm = mg_size[j][0] * mg_size[j][1] * mg_size[j][2];

            if (j < threshold) {
                // Resize
                dim3 blocks((mg_size[j][0] - 1)/devThreadsPerBlock + 1, (mg_size[j][1] - 1) / devThreadsPerBlock + 1, (mg_size[j][2] - 1) / devThreadsPerBlock + 1);
                restrict_eps_map_kernel<<<blocks, threadsPerBlock>>>(lepsx + l, lepsy + l, lepsz + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], lepsx + m, lepsy + m, lepsz+ m, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
                cudaLaunchErrorCheck();
                restrict_v_kernel<<<blocks, threadsPerBlock>>>(64.0, l_iv + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2]); //iv
                cudaLaunchErrorCheck();

                // 1-D grid
                dim3 h_blocks ((lxmymzm - 1) / 512 + 1);
                set_am_ad_kernel_head<<<h_blocks, h_threadsPerBlock>>>(lepsx + m, lepsy + m, lepsz + m, l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, lxmymzm);
                cudaLaunchErrorCheck();
                // 3-D grid
                set_am_ad_kernel_body<<<blocks, threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, l_ad + m, l_bz + m, l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2], lfactor, l_epsout);
                cudaLaunchErrorCheck();
                // 2-D grid
                dim3 t_blocks((max(mg_size[j][0], mg_size[j][1]) - 1) / 16 + 1, (max(mg_size[j][1], mg_size[j][2]) - 1) / 16 + 1);
                set_am_ad_kernel_tail<<<t_blocks, t_threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
                cudaLaunchErrorCheck();

                cudaDeviceSynchronize();

            } else {
                restrict_eps_map(lepsx + l, lepsy + l, lepsz + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], lepsx + m, lepsy + m, lepsz + m, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
                restrict_v(64.0, l_iv + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2]); //iv
                set_am_ad(lepsx + m, lepsy + m, lepsz + m, l_iv + m, l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, l_ad + m, l_bz + m, mg_size[j][0], mg_size[j][1], mg_size[j][2], lfactor, l_epsout);
            }
        }
    }
    cudaFree(lp_bv);
    cudaFree(lp_iv);
    cudaFree(lp_xs);
    cudaFree(lepsx);
    cudaFree(lepsy);
    cudaFree(lepsz);
    cudaFree(epsx_f);
    cudaFree(epsy_f);
    cudaFree(epsz_f);
}

__global__
void init_vector_kernel(float *vec, int m) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) vec[i] = 0.0;
}

__global__
void copy_vector_kernel(float *vec, float *vec_f, int m) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) vec[i] = vec_f[i];
}

__global__
void inv_vector_kernel(float *vec, float *inv, int m) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) inv[i] = 1.0/vec[i];
}

__host__
void restrict_eps_map(float *epsxf, float *epsyf, float *epszf, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr) {
    if (l_bcopt == 10) {
        printf("Not yet for PBC.\n");
        //Do nothing;
    } else {
        for(int k = 0; k < znr; k++) {
            int k2 = 2 * k + 1;
            for(int j = 0; j < ynr; j++) {
                int j2 = 2 * j + 1;
                for(int i = 0; i < xnr; i++) {
                    int i2 = 2 * i + 1;
                    int flatid = i + xnr * j + xnr * ynr * k;
                    // eps*r causes CUDA Exception 15 error. Solved
                    epsxr[flatid] = r_map_exp_x(epsxf, i2, j2, k2, xn, yn);
                    epsyr[flatid] = r_map_exp_y(epsyf, i2, j2, k2, xn, yn);
                    epszr[flatid] = r_map_exp_z(epszf, i2, j2, k2, xn, yn);
                }
            }
        }
    }
}

__host__ __device__
float r_map_exp_x(float *epsxmp, int i2, int j2, int k2, int xn, int yn) {
    float exp_x = hmav(epsxmp[f_id(i2  , j2  , k2  , xn, yn)], epsxmp[f_id(i2+1, j2  , k2  , xn, yn)]) / 4.0 +
                 (hmav(epsxmp[f_id(i2  , j2-1, k2  , xn, yn)], epsxmp[f_id(i2+1, j2-1, k2  , xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2+1, k2  , xn, yn)], epsxmp[f_id(i2+1, j2+1, k2  , xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2  , k2-1, xn, yn)], epsxmp[f_id(i2+1, j2  , k2-1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2  , k2+1, xn, yn)], epsxmp[f_id(i2+1, j2  , k2+1, xn, yn)])) / 8.0 +
                 (hmav(epsxmp[f_id(i2  , j2-1, k2-1, xn, yn)], epsxmp[f_id(i2+1, j2-1, k2-1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2+1, k2-1, xn, yn)], epsxmp[f_id(i2+1, j2+1, k2-1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2-1, k2+1, xn, yn)], epsxmp[f_id(i2+1, j2-1, k2+1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2+1, k2+1, xn, yn)], epsxmp[f_id(i2+1, j2+1, k2+1, xn, yn)])) / 16.0;

    // Debug
    //raise(SIGINT);
    return exp_x;
}

__host__ __device__
float r_map_exp_y(float *epsymp, int i2, int j2, int k2, int xn, int yn) {
    float exp_y = hmav(epsymp[f_id(i2  , j2  , k2  , xn, yn)], epsymp[f_id(i2  , j2+1, k2  , xn, yn)]) / 4.0 +
                 (hmav(epsymp[f_id(i2-1, j2  , k2  , xn, yn)], epsymp[f_id(i2-1, j2+1, k2  , xn, yn)]) +
                  hmav(epsymp[f_id(i2+1, j2  , k2  , xn, yn)], epsymp[f_id(i2+1, j2+1, k2  , xn, yn)]) +
                  hmav(epsymp[f_id(i2  , j2  , k2-1, xn, yn)], epsymp[f_id(i2  , j2+1, k2-1, xn, yn)]) +
                  hmav(epsymp[f_id(i2  , j2  , k2+1, xn, yn)], epsymp[f_id(i2  , j2+1, k2+1, xn, yn)])) / 8.0 +
                 (hmav(epsymp[f_id(i2-1, j2  , k2-1, xn, yn)], epsymp[f_id(i2-1, j2+1, k2-1, xn, yn)]) +
                  hmav(epsymp[f_id(i2+1, j2  , k2-1, xn, yn)], epsymp[f_id(i2+1, j2+1, k2-1, xn, yn)]) +
                  hmav(epsymp[f_id(i2-1, j2  , k2+1, xn, yn)], epsymp[f_id(i2-1, j2+1, k2+1, xn, yn)]) +
                  hmav(epsymp[f_id(i2+1, j2  , k2+1, xn, yn)], epsymp[f_id(i2+1, j2+1, k2+1, xn, yn)])) / 16.0;
    return exp_y;
}

__host__ __device__
float r_map_exp_z(float *epszmp, int i2, int j2, int k2, int xn, int yn) {
    float exp_z = hmav(epszmp[f_id(i2  , j2  , k2  , xn, yn)], epszmp[f_id(i2  , j2  , k2+1, xn, yn)]) / 4.0 +
                 (hmav(epszmp[f_id(i2  , j2-1, k2  , xn, yn)], epszmp[f_id(i2  , j2-1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2  , j2+1, k2  , xn, yn)], epszmp[f_id(i2  , j2+1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2-1, j2  , k2  , xn, yn)], epszmp[f_id(i2-1, j2  , k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2+1, j2  , k2  , xn, yn)], epszmp[f_id(i2+1, j2  , k2+1, xn, yn)])) / 8.0 +
                 (hmav(epszmp[f_id(i2-1, j2-1, k2  , xn, yn)], epszmp[f_id(i2-1, j2-1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2-1, j2+1, k2  , xn, yn)], epszmp[f_id(i2-1, j2+1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2+1, j2-1, k2  , xn, yn)], epszmp[f_id(i2+1, j2-1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2+1, j2+1, k2  , xn, yn)], epszmp[f_id(i2+1, j2+1, k2+1, xn, yn)])) / 16.0;
    return exp_z;
}

__host__ __device__
float hmav(float a, float b) {
    return 2.0 * a * b / (a + b);
}

__global__
void restrict_eps_map_kernel(float *epsxf, float *epsyf, float *epszf, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int i2 = 2 * i + 1;
    int j2 = 2 * j + 1;
    int k2 = 2 * k + 1;
    if (i < xnr  && j < ynr  && k < znr){
        int flatid = i + xnr * j + xnr * ynr * k;
        epsxr[flatid] = r_map_exp_x(epsxf, i2, j2, k2, xn, yn);
        epsyr[flatid] = r_map_exp_y(epsyf, i2, j2, k2, xn, yn);
        epszr[flatid] = r_map_exp_z(epszf, i2, j2, k2, xn, yn);
    }
}

// Assign on device
__global__
void feedepsintoam_kernel(int lxm, int lym, int lzm, float *am1, float *am2, float *am3, float *eps1, float *eps2, float *eps3) {
    // eps1 has (0:xm, ym, zm) dimension (0:ym/0:zm for eps2/eps3), which contains am1 (xm, ym, zm), so need to recheck the passed array index.
    // Edit: Solved, using 3-D mapping
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < lxm && j < lym && k < lzm) {
        int flatid = i + lxm * j + lxm * lym * k;
        // Warning of bug, x/y dimensions are incorrect, that's why only lepsz is good
        //int flatid_x_plus1 = (i + 1) + lxm * j + lxm * lym * k;
        //int flatid_y_plus1 = i + lxm * (j + 1) + lxm * lym * k;
        int flatid_x_plus1 = (i + 1) + (lxm + 1) * j + (lxm + 1) * lym * k;
        int flatid_y_plus1 = i + lxm * (j + 1) + lxm * (lym + 1) * k;
        int flatid_z_plus1 = i + lxm * j + lxm * lym * (k + 1);

        // eps* causes CUDA Exception 14 error. Solved
        am1[flatid] = eps1[flatid_x_plus1];
        am2[flatid] = eps2[flatid_y_plus1];
        am3[flatid] = eps3[flatid_z_plus1];
    }
}

__host__
void set_am_ad(float *epsx, float *epsy, float *epsz, float *iv, float *lam1, float *lam2, float *lam3, float *lad, float *lbz, int xn, int yn, int zn, float lfactor, float epsout) {
    for (int i = 0; i < xn * yn * zn; i++) {
        lam1[i] = epsx[i];
        lam2[i] = epsy[i];
        lam3[i] = epsz[i];
    }

    for (int k = 0; k < zn; k++) {
        for (int j = 0; j < yn; j++) {
            for (int i = 0; i < xn; i++) {
                int flatid = i + xn * j + xn * yn * k;
                int flatid_x_minus1 = (i - 1) + xn * j + xn * yn * k;
                int flatid_y_minus1 = i + xn * (j - 1) + xn * yn * k;
                int flatid_z_minus1 = i + xn * j + xn * yn * (k - 1);

                lad[flatid] = lam1[flatid] + lam2[flatid] + lam3[flatid];
                if (i == 0) lad[flatid] += epsout; else lad[flatid] += lam1[flatid_x_minus1];
                if (j == 0) lad[flatid] += epsout; else lad[flatid] += lam2[flatid_y_minus1];
                if (k == 0) lad[flatid] += epsout; else lad[flatid] += lam3[flatid_z_minus1];

                lbz[flatid] = lfactor * iv[flatid];
                lad[flatid] += lbz[flatid];
            }
        }
    }

    if (l_bcopt != 10) {
        for (int k = 0; k < zn; k++) {
            for (int j = 0; j < yn; j++) {
                for (int i = 0; i < xn; i++) {
                    int flatid = i + xn * j + xn * yn * k;
                    if (i == xn - 1) lam1[flatid] = 0;
                    if (j == yn - 1) lam2[flatid] = 0;
                    if (k == zn - 1) lam3[flatid] = 0;
                }
            }
        }
    }
}

__global__
void set_am_ad_kernel_head(float *epsx, float *epsy, float *epsz, float *lam1, float *lam2, float *lam3, int xnynzn) {
    // 1-D, xnynzn
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < xnynzn) {
        lam1[i] = epsx[i];
        lam2[i] = epsy[i];
        lam3[i] = epsz[i];
    }
}

__global__
void set_am_ad_kernel_body(float *lam1, float *lam2, float *lam3, float *lad, float *lbz, float *iv, int xn, int yn, int zn, float lfactor, float epsout) {
    // 3-D, (xn, yn, zn)
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < xn && j < yn && k < zn) {
        int flatid = i + xn * j + xn * yn * k;
        int flatid_x_minus1 = (i - 1) + xn * j + xn * yn * k;
        int flatid_y_minus1 = i + xn * (j - 1) + xn * yn * k;
        int flatid_z_minus1 = i + xn * j + xn * yn * (k - 1);

        lad[flatid] = lam1[flatid] + lam2[flatid] + lam3[flatid];
        if (i == 0) lad[flatid] += epsout; else lad[flatid] += lam1[flatid_x_minus1];
        if (j == 0) lad[flatid] += epsout; else lad[flatid] += lam2[flatid_y_minus1];
        if (k == 0) lad[flatid] += epsout; else lad[flatid] += lam3[flatid_z_minus1];

        lbz[flatid] = lfactor * iv[flatid];
        lad[flatid] += lbz[flatid];
    }
}

__global__
void set_am_ad_kernel_tail(float *lam1, float *lam2, float *lam3, int xn, int yn, int zn) {
    // 2-D, (i, j) -> (max_xy, max_yz)
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < yn && j < zn) {
        lam1[xn - 1 + xn * i + xn * yn * j] = 0.0;
    }

    if (i < xn && j < zn) {
        lam2[i + xn * (yn - 1) + xn * yn * j] = 0.0;
    }

    if (i < xn && j < yn) {
        lam2[i + xn * j + xn * yn * (zn - 1)] = 0.0;
    }
}

__host__
void restrict_v(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr) {
    for(int k = 0; k < nzr; k++) {
        int k2 = 2 * k + 1;
        for(int j = 0; j < nyr; j++) {
            int j2 = 2 * j + 1;
            for(int i = 0; i < nxr; i++) {
                int i2 = 2 * i + 1;
                int flatid = i + nxr * j + nxr * nyr * k;
                bvr[flatid] = ( bvf[f_id(i2-1, j2-1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2-1, nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2  , k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2-1, nx, ny)] + bvf[f_id(i2+1, j2  , k2-1, nx, ny)] ) +
                              ( bvf[f_id(i2-1, j2+1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2-1, nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2-1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2  , nx, ny)] + bvf[f_id(i2+1, j2-1, k2  , nx, ny)] ) +
                          4 * ( bvf[f_id(i2-1, j2  , k2  , nx, ny)] + 2 * bvf[f_id(i2, j2  , k2  , nx, ny)] + bvf[f_id(i2+1, j2  , k2  , nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2+1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2  , nx, ny)] + bvf[f_id(i2+1, j2+1, k2  , nx, ny)] ) +
                              ( bvf[f_id(i2-1, j2-1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2+1, nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2  , k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2+1, nx, ny)] + bvf[f_id(i2+1, j2  , k2+1, nx, ny)] ) +
                              ( bvf[f_id(i2-1, j2+1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2+1, nx, ny)] );
                bvr[flatid] /= divider;
            }
        }
    }
}

__global__
void restrict_v_kernel(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int i2 = 2 * i + 1;
    int j2 = 2 * j + 1;
    int k2 = 2 * k + 1;
    if (i < nxr  && j < nyr  && k < nzr){
        int b_id = i + nxr * j + nxr * nyr * k;

        bvr[b_id] = ( bvf[f_id(i2-1, j2-1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2-1, nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2  , k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2-1, nx, ny)] + bvf[f_id(i2+1, j2  , k2-1, nx, ny)] ) +
                    ( bvf[f_id(i2-1, j2+1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2-1, nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2-1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2  , nx, ny)] + bvf[f_id(i2+1, j2-1, k2  , nx, ny)] ) +
                4 * ( bvf[f_id(i2-1, j2  , k2  , nx, ny)] + 2 * bvf[f_id(i2, j2  , k2  , nx, ny)] + bvf[f_id(i2+1, j2  , k2  , nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2+1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2  , nx, ny)] + bvf[f_id(i2+1, j2+1, k2  , nx, ny)] ) +
                    ( bvf[f_id(i2-1, j2-1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2+1, nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2  , k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2+1, nx, ny)] + bvf[f_id(i2+1, j2  , k2+1, nx, ny)] ) +
                    ( bvf[f_id(i2-1, j2+1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2+1, nx, ny)] );
        bvr[b_id] /= divider;
    }
}

__host__ __device__
int f_id(int i, int j, int k, int nx, int ny) {
    return i + nx * j + nx * ny * k;
}

__host__
void restrict_cuda(int level) {
    float div = 16.0; //restrict_bv

    int nx = mg_size[level][0];
    int ny = mg_size[level][1];
    int nz = mg_size[level][2];
    int nxr = mg_size[level + 1][0];
    int nyr = mg_size[level + 1][1];
    int nzr = mg_size[level + 1][2];

    dim3 threadsPerBlock(devThreadsPerBlock, devThreadsPerBlock, devThreadsPerBlock);
    dim3 blocks((nx - 1)/devThreadsPerBlock + 1, (ny - 1) / devThreadsPerBlock + 1, (nz - 1) / devThreadsPerBlock +1);
    restrict_v_kernel<<<blocks, threadsPerBlock>>>(div, l_rv + mg_index[level], nx, ny, nz, l_bv + mg_index[level + 1], nxr, nyr, nzr);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();
}

__host__
void interpolate(int level) {
    // Derive all vectors
    int xn = mg_size[level + 1][0];
    int yn = mg_size[level + 1][1];
    int zn = mg_size[level + 1][2];
    int xni = mg_size[level][0];
    int yni = mg_size[level][1];
    int zni = mg_size[level][2];
    int xniyni = xni * yni;
    int xniynizni = xni * yni * zni;

    int p1 = mg_x_idx[level + 1] + mg_size[level + 1][0] * mg_size[level + 1][1];
    int p2 = mg_x_idx[level] + mg_size[level][0] * mg_size[level][1];

    float *v = &l_xv[p1];
    float *vi = &l_xv[p2];
    float *lam1 = &l_am1[mg_index_ext[level]]; // 1-xniyni:~
    float *lam2 = &l_am2[mg_index_ext[level]];
    float *lam3 = &l_am3[mg_index_ext[level]];
    float *lbz = &l_bz[mg_index[level]];
    float epsout = l_epsout;

    if (xn * 2 + 1 != xni || yn * 2 + 1 != yni || zn * 2 + 1 != zni) {
        printf("Interpolation failed because of incorrect dimension (interpolate_host)\n");
        printf("xn %d, yn %d, zn %d\n", xn, yn, zn);
        printf("xni %d, yni %d, zni %d\n", xni, yni, zni);
        exit(2);
    }

    //raise(SIGINT);

    /*
    for(int k = 0; k < zni; k++) {
        for(int j = 0; j < yni; j++) {
            for(int i = 0; i < xni; i++) {
                int flatid_c = i + xni * j + xni * yni * k; // Coarse grid
                // Later handel lam value of epsout at points on surfaces (i-1/j-1/k-1). TBD
                // Use 3-D index to manually add epsout value with outter surface
                // Check first what lam values were used inside each function
                // The original code wasted too much with 1-xniyni:0 range. Need to use 3-D surface to computate the
                // contribution at each point.
                // For now, temporarily use original index, add any offset needed here to fetch the correct position
                // Since am1/2/3 is used as 1-D array only in ipl_comp* functions, add offset there will do the trick.
                // Edit: the index of am* is linear to l, so let l += offset, and finially flatid_f += offset.
                // Edit: Don't add to flatid_f, which will cause error to vi; add to am* index in ipl_com* func only.

                //if (i == 0) lam1[flatid_c] += ;

                // Surfaces around vertices
                // Mapping: Thread block --> Coarse grid block --> Fine grid block
                if (i == xni) lam1[flatid_c] = epsout;
                if (j == yni) lam2[flatid_c] = epsout;
                if (k == zni) lam3[flatid_c] = epsout;

            }
        }
    }*/

    // Initialize to epsout, should be before ipl_chain function calls
    // 1-xniyni:0, outside 3-D loop
    for ( int i = 0; i < xniyni; i++) {
        lam1[i] = epsout;
        lam2[i] = epsout;
        lam3[i] = epsout;
    }

    // Three surfaces
    // Need adding index offset. TBF. Fixed
    for (int k = 0; k < zni; k++) {
        for (int j = 0; j < yni; j++) {
            // [(xni - 1) + j * xni + k * xniyni] + xniyni
            lam1[-1 + (j + 1) * xni + (k + 1) * xniyni] = epsout;
        }
    }
    for (int k = 0; k < zni; k++) {
        for (int i = 0; i < xni; i++) {
            // [i + (yni - 1) * xni + k * xniyni] + xniyni
            lam2[i - xni + (k + 2) * xniyni] = epsout;
        }
    }
    for (int j = 0; j < yni; j++) {
        for (int i = 0; i < xni; i++) {
            lam3[i + j * xni + xniynizni] = epsout;
        }
    }

    for(int k = 0; k < zn; k++) {
        for(int j = 0; j < yn; j++) {
            for(int i = 0; i < xn; i++) {
                // Caution with the indexes! Starting from 0
                int flatid_c = i + xn * j + xn * yn * k; // Coarse grid
                int flatid_f = (2 * i + 1) + (2 * j + 1) * xni + (2 * k + 1) * xniyni; // Fine grid

                vi[flatid_f] += v[flatid_c];

                // offset for offsetting the index
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      -1, lam2, xni, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      +1, lam2, xni, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    -xni, lam1,   1, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    +xni, lam1,   1, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, -xniyni, lam2, xni, lam1,      1, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, +xniyni, lam2, xni, lam1,      1, xn, yn, zn);
            }
        }
    }

    // 1-xniyni:0, outside 3-D loop
    for ( int i = 0; i < xniyni; i++) {
        lam1[i] = 0.0;
        lam2[i] = 0.0;
        lam3[i] = 0.0;
    }

    // Three surfaces
    for (int k = 0; k < zni; k++) {
        for (int j = 0; j < yni; j++) {
            lam1[-1 + (j + 1) * xni + (k + 1) * xniyni] = 0.0;
        }
    }
    for (int k = 0; k < zni; k++) {
        for (int i = 0; i < xni; i++) {
            lam2[i - xni + (k + 2) * xniyni] = 0.0;
        }
    }
    for (int j = 0; j < yni; j++) {
        for (int i = 0; i < xni; i++) {
            lam3[i + j * xni + xniynizni] = 0.0;
        }
    }

    // A better way to handel and fix all index issues, using index 1 ~ xnynzn, not
    // 1-nxny:nxnynz+nxny.
    // i.e. grid[level]->xs[0:lxlylz-1] v.s. (arrays + index dereference),
    // but at the expense of struct indrection overhead. Figure out to what degree.
    // May define a general Grid structure with all vector pointers inside, then create an array
    // Grid *grid[MG_NLEVEL],
    // then initialize all matrix vectors via grid[level]->xs = (float *) malloc();
}

__global__
void interpolate_kernel_head(int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float epsout) {
    int xniyni = xni * yni;
    // 2-D, [i, j] -> [max(xni, yni), max(yni, zni)]
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // Padding offset <xniyni> included
    if (i < yni && j < zni) {
        lam1[-1 + (i + 1) * xni + (j + 1) * xniyni] = epsout;
    }

    if (i < xni && j < zni) {
        lam2[i - xni + (j + 2) * xniyni] = epsout;
    }

    if (i < xni && j < yni) {
        lam3[i + j * xni + xni * yni * zni] = epsout;
        // 0:xniyni-1
        lam1[i + j * xni] = epsout;
        lam2[i + j * xni] = epsout;
        lam3[i + j * xni] = epsout;
    }
}

__global__
void interpolate_kernel_body(float *v, int xn, int yn, int zn, float *vi, int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float *lbz) {
    int xniyni = xni * yni;
    int xniynizni = xni * yni * zni;

    // 3-D, (xn, yn, zn)
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < xn  && j < yn  && k < zn){
        int flatid_c = i + xn * j + xn * yn * k; // Coarse grid
        int flatid_f = (2 * i + 1) + (2 * j + 1) * xni + (2 * k + 1) * xniyni; // Fine grid

        vi[flatid_f] += v[flatid_c];

        // offset for offsetting the index
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      -1, lam2, xni, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      +1, lam2, xni, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    -xni, lam1,   1, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    +xni, lam1,   1, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, -xniyni, lam2, xni, lam1,      1, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, +xniyni, lam2, xni, lam1,      1, xn, yn, zn);
    }

}

// Refine the four types point value assignment using directly thread block, layer by layer,
// pass single array value distributed on single thread. TBD
__host__
void ipl_chain_h(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    // Here v is value of coarse grid; vi is array of fine grid
    float v1 = ipl_comp1(v, l, lbz, am_1, xnyn, xnynzn, shift_1);
    int l1 = l + shift_1;
    vi[l1] += v1;
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2, -shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2,  shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3, -shift_3, am_2, shift_2, xn, yn, zn);
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3,  shift_3, am_2, shift_2, xn, yn, zn);
}

__device__
void ipl_chain_d(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    // Here v is value of coarse grid; vi is array of fine grid
    float v1 = ipl_comp1(v, l, lbz, am_1, xnyn, xnynzn, shift_1);
    int l1 = l + shift_1;
    atomicAdd(&vi[l1], v1);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2, -shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2,  shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3, -shift_3, am_2, shift_2, xn, yn, zn);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3,  shift_3, am_2, shift_2, xn, yn, zn);
}

__host__
void ipl_chain2_h(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    float v2 = ipl_comp2(v1, l1, lbz, am_1, am_2, xnyn, xnynzn, shift_1, shift_2);
    int l2 = l1 + shift_2;
    vi[l2] += v2;
    vi[l2 - shift_3] += ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, -shift_3);
    vi[l2 + shift_3] += ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, +shift_3);
}

__device__
void ipl_chain2_d(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    float v2 = ipl_comp2(v1, l1, lbz, am_1, am_2, xnyn, xnynzn, shift_1, shift_2);
    int l2 = l1 + shift_2;
    atomicAdd(&vi[l2], v2);
    float v3 = ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, -shift_3);
    atomicAdd(&vi[l2 - shift_3], v3);
    float v4 = ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, +shift_3);
    atomicAdd(&vi[l2 + shift_3], v4);
}

__host__ __device__
float ipl_comp1(float v, int l, float *lbz, float *am_1, int xnyn, int xnynzn, int shift_1) {
    float ipl_comp1_v;
    // For offsetting
    int bz_l = l;
    l += xnyn;

    if (shift_1 < 0)
        ipl_comp1_v = v * am_1[l + shift_1] / ( lbz[bz_l + shift_1] + am_1[l + 2 * shift_1] + am_1[l + shift_1] );
    else
        ipl_comp1_v = v * am_1[l] / ( lbz[bz_l + shift_1] + am_1[l] + am_1[l + shift_1] );
    return ipl_comp1_v;
}

__host__ __device__
float ipl_comp2(float v, int l, float *lbz, float *am_1, float *am_2, int xnyn, int xnynzn, int shift_1, int shift_2) {
    // For offsetting
    int bz_l = l;
    l += xnyn;

    float lad = am_1[l + shift_2] + am_1[l + shift_2 - abs(shift_1)] + lbz[bz_l + shift_2];
    float ipl_comp2_v;

    if (shift_2 < 0)
        ipl_comp2_v = v * am_2[l + shift_2] / ( am_2[l + 2 * shift_2] + am_2[l + shift_2] + lad);
    else
        ipl_comp2_v = v * am_2[l] / ( am_2[l] + am_2[l + shift_2] + lad );
    return ipl_comp2_v;
}

__host__ __device__
float ipl_comp3(float v, int l, float *lbz, float *am_1, float *am_2, float *am_3, int xnyn, int xnynzn, int shift_1, int shift_2, int shift_3) {
    // For offsetting
    int bz_l = l;
    l += xnyn;

    float lad = am_1[l + shift_3] + am_1[l + shift_3 - abs(shift_1)] + am_2[l + shift_3] + am_2[l + shift_3 - abs(shift_2)] + lbz[bz_l + shift_3];
    float ipl_comp3_v;

    if (shift_3 < 0)
        ipl_comp3_v = v * am_3[l + shift_3] / ( am_3[l + 2 * shift_3] + am_3[l + shift_3] + lad);
    else
        ipl_comp3_v = v * am_3[l] / ( am_3[l] + am_3[l + shift_3] + lad );
    return ipl_comp3_v;
}

__host__
void interpolate_cuda(int level) {
    int xn = mg_size[level + 1][0];
    int yn = mg_size[level + 1][1];
    int zn = mg_size[level + 1][2];
    int xni = mg_size[level][0];
    int yni = mg_size[level][1];
    int zni = mg_size[level][2];

    int p1 = mg_x_idx[level + 1] + mg_size[level + 1][0] * mg_size[level + 1][1];
    int p2 = mg_x_idx[level] + mg_size[level][0] * mg_size[level][1];
    int p3 = mg_index_ext[level];
    int p4 = mg_index[level];

    if (xn * 2 + 1 != xni || yn * 2 + 1 != yni || zn * 2 + 1 != zni) {
        printf("Interpolation failed because of incorrect dimension (interpolate_cuda)\n");
        printf("xn %d, yn %d, zn %d\n", xn, yn, zn);
        printf("xni %d, yni %d, zni %d\n", xni, yni, zni);
        exit(2);
    }

    /*
     * warning: Cuda API error detected: cudaLaunch returned (0x7).
     *
     * From /usr/local/cuda-7.5/include/driver_types.h:
     *      cudaErrorLaunchOutOfResources         =      7,
     *      This indicates that a launch did not occur because it did not have
     *      appropriate resources. Although this error is similar to
     *      ::cudaErrorInvalidConfiguration, this error usually indicates that the
     *      user has attempted to pass too many arguments to the device kernel, or the
     *      kernel launch specifies too many threads for the kernel's register count.
     *
     * Fixed temporarily. Use nvcc --maxrregcount=128
     * In release code this should be specified inside code, best using launch bounds:
     * http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#launch-bounds
     */

    // 2-D grid
    dim3 h_threadsPerBlock(16, 16);
    dim3 h_blocks((max(xni, yni) - 1) / 16 + 1, (max(yni, zni) - 1) / 16 + 1);
    interpolate_kernel_head<<<h_blocks, h_threadsPerBlock>>>(xni, yni, zni, l_am1 + p3, l_am2 + p3, l_am3 + p3, l_epsout);
    cudaLaunchErrorCheck();

    // Using (8, 8, 4) grid solved the 'CUDA launch failed: too many resources requested' issue.
    // Optimize the register use, block, grid size later. TBD
    dim3 threadsPerBlock(devThreadsPerBlock, devThreadsPerBlock, devThreadsPerBlock/2);
    dim3 blocks((xn - 1) / devThreadsPerBlock + 1, (yn - 1) / devThreadsPerBlock + 1, 2 * (zn - 1) / devThreadsPerBlock + 1);
    interpolate_kernel_body<<<blocks, threadsPerBlock>>>(l_xv + p1, xn, yn, zn, l_xv + p2, xni, yni, zni, l_am1 + p3, l_am2 + p3, l_am3 + p3, l_bz + p4);
    cudaLaunchErrorCheck();

    interpolate_kernel_head<<<h_blocks, h_threadsPerBlock>>>(xni, yni, zni, l_am1 + p3, l_am2 + p3, l_am3 + p3, 0.0);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();
}

__host__
void relax(int level, int ncyc) {
    int nx = mg_size[level][0];
    int ny = mg_size[level][1];
    int nz = mg_size[level][2];
    int nxny = nx * ny;
    int nxnynz = nxny * nz;
    float *xs = l_xv + mg_x_idx[level];
    float *lam1 = l_am1 + mg_index_ext[level];
    float *lam2 = l_am2 + mg_index_ext[level];
    float *lam3 = l_am3 + mg_index_ext[level];
    float *lzv = l_zv + mg_index[level];
    float *lad = l_ad + mg_index[level];
    float *lbv = l_bv + mg_index[level];
    float *lrv = l_rv + mg_index[level];
    float accept = l_accept;

    float onorm = mg_onorm[level];

    int itn_checknorm;
    float wsor;//, wsor1;
    float linorm = 0.0;
    float lnorm;
    int itmax = 1000; // Should move to upper layer

    if (ncyc > 0) {
        itn_checknorm = ncyc;
        wsor = 1.0; //1.0; Debug: 1.0 is Gauss-Seidel
    } else {
        itn_checknorm = 10;
        wsor = 1.9; // 1.0; Debug of omega - use GS to solve on coarsest level
    }

    for (int i = 0; i < nxnynz; i++) {
        linorm += abs(lbv[i]);
        lzv[i] = 1.0 / lad[i];
    }

    bool converged = false;
    int litn = 0;
    while (!converged) {
        int i;
        for (i = nxny;  i < nxnynz+nxny; i++) {
            xs[i] -= wsor * (xs[i] - (lam1[i - 1   ] * xs[i - 1   ] + lam1[i         ] * xs[i + 1   ] +
                                      lam2[i - nx  ] * xs[i - nx  ] + lam2[i         ] * xs[i + nx  ] +
                                      lam3[i - nxny] * xs[i - nxny] + lam3[i         ] * xs[i + nxny] + lbv[i - nxny]) * lzv[i - nxny]);
        }

        litn++;

        // Check convergence
        if (litn % itn_checknorm == 0) {
            // residue
            for (int i = nxny; i < nxnynz+nxny; i++) {
                lrv[i - nxny] = lam1[i - 1   ] * xs[i - 1   ] + lam1[i       ] * xs[i + 1   ] +
                                lam2[i - nx  ] * xs[i - nx  ] + lam2[i       ] * xs[i + nx  ] +
                                lam3[i - nxny] * xs[i - nxny] + lam3[i       ] * xs[i + nxny] + lbv[i - nxny] - lad[i - nxny] * xs[i];
                }

            // norm
            lnorm = 0.0;
            for (int i = 0; i < nxnynz; i++) {
                lnorm += abs(lrv[i]);
            }

            if (litn >= itmax || (ncyc > 0 && (litn >= ncyc && lnorm < onorm)) || lnorm <= accept * linorm) {
                converged = true;
                if (ncyc > 0 && litn >= ncyc && lnorm > onorm) {
                    printf("PB_MG FAILED: ncyc %d\t, itn %d\t, norm %e\t, onorm %e\n", ncyc, litn, lnorm, onorm);
                    exit(2);
                }

                if (ncyc > 0) mg_onorm[level] = lnorm; // Update global array
                if (litn >= itmax) printf("PB_MG WARNING: SOR maxitn exceeded (relax_host)!\n");
            }

        }
    } // while
}

__host__
void relax_cuda(int level, int ncyc) {
    int threadsPerBlock = 512; // This shouldn't be final, RL

    int nx = mg_size[level][0];
    int ny = mg_size[level][1];
    int nz = mg_size[level][2];
    int nxny = nx * ny;
    int nxnynz = nxny * nz;

    float *xs = l_xv + mg_x_idx[level];
    float *lam1 = l_am1 + mg_index_ext[level];
    float *lam2 = l_am2 + mg_index_ext[level];
    float *lam3 = l_am3 + mg_index_ext[level];
    float *lzv = l_zv + mg_index[level];
    float *lad = l_ad + mg_index[level];
    float *lbv = l_bv + mg_index[level];
    float *lrv = l_rv + mg_index[level];
    float accept = l_accept;
    float onorm = mg_onorm[level];

    int itn_checknorm = 10;
    float wsor;//, wsor1;
    int itmax = 1000;

    if (ncyc > 0) {
        itn_checknorm = ncyc;
        wsor = 1.0;// 1.0 Test
    } else {
        itn_checknorm = 10;
        wsor = 1.9; // 1.9 Test
    }

    ncyc = 10; // This shouldn't be final, RL

    // inverse AD for fast processing later
    int blocks = (nxnynz - 1) / threadsPerBlock + 1;
    inv_vector_kernel<<<blocks, threadsPerBlock>>>(lad, lzv, nxnynz);
    cudaLaunchErrorCheck();

    // initial norm. Create cuBlAS context
    float linorm = 0.0;
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    cublasErrorCheck(status);
    status = cublasSasum(handle, nxnynz, lbv, 1, &linorm);
    cublasErrorCheck(status);

    bool converged = false;
    int litn = 0;
    int sblocks = ((int)(nxnynz / 2) + (nxnynz & 1) - 1) / threadsPerBlock + 1;
    while (!converged) {
        // non-periodic
        solver_red_kernel<<<sblocks, threadsPerBlock>>>(xs, lam1, lam2, lam3, lzv, lbv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        solver_black_kernel<<<sblocks, threadsPerBlock>>>(xs, lam1, lam2, lam3, lzv, lbv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        //cudaDeviceSynchronize(); // Warning of bug: CUDA_EXCEPTION_15

        litn++;

        // Check convergence
        if (litn % itn_checknorm == 0) {
            residue_kernel<<<blocks, threadsPerBlock>>>(xs, lam1, lam2, lam3, lad, lbv, nx, ny, nz, lrv);
            cudaLaunchErrorCheck();

            // Create cuBlAS context
            float lnorm = 0.0;
            status = cublasSasum(handle, nxnynz, lrv, 1, &lnorm);
            cublasErrorCheck(status);

            if (litn >= itmax || (ncyc > 0 && (litn >= ncyc && lnorm < onorm)) || lnorm <= accept * linorm) {
                converged = true;
                if (ncyc > 0 && litn >= ncyc && lnorm > onorm) {
                    printf("PB_MG FAILED: ncyc %d\t, itn %d\t, norm %e\t, onorm %e\n", ncyc, litn, lnorm, onorm);
                    exit(2);
                }

                if (ncyc > 0) mg_onorm[level] = lnorm; // Update global array
                if (litn >= itmax) printf("PB_MG WARNING: SOR maxitn exceeded (relax_kernel)!\n");
            }
        } // if
    } // while

    // Destroy context
    cublasDestroy(handle);
}

__global__
void solver_black_kernel(float *phi, float *epsi, float *epsj, float *epsk,float *repsc, float *rho, float wsor, int xm, int ym, int zm) {

    int xmym = xm * ym;
    int i = 2 * (blockIdx.x * blockDim.x + threadIdx.x) + xmym;

    if (i < xm * ym * zm + xmym) {
        phi[i] -= wsor * (phi[i] - (epsi[i - 1   ] * phi[i - 1   ] + epsi[i      ] * phi[i + 1   ] +
                                    epsj[i - xm  ] * phi[i - xm  ] + epsj[i      ] * phi[i + xm  ] +
                                    epsk[i - xmym] * phi[i - xmym] + epsk[i      ] * phi[i + xmym] + rho[i - xmym]) * repsc[i - xmym]);
    }

}

__global__
void solver_red_kernel(float *phi, float *epsi,float *epsj, float *epsk, float *repsc, float *rho, float wsor, int xm, int ym, int zm) {

    int xmym = xm * ym;
    int i = 2 * (blockIdx.x * blockDim.x + threadIdx.x) + 1 + xmym;

    if (i < xm * ym * zm + xmym) {
        phi[i] -= wsor * (phi[i] - (epsi[i - 1   ] * phi[i - 1   ] + epsi[i      ] * phi[i + 1   ] +
                                    epsj[i - xm  ] * phi[i - xm  ] + epsj[i      ] * phi[i + xm  ] +
                                    epsk[i - xmym] * phi[i - xmym] + epsk[i      ] * phi[i + xmym] + rho[i - xmym]) * repsc[i - xmym]);
    }

}

__global__
void residue_kernel(float *phi, float *epsi,float *epsj, float *epsk, float *epsc, float *rho, int xm, int ym, int zm, float* res) {

    int xmym = xm * ym;
    int i = (blockIdx.x * blockDim.x + threadIdx.x) + xmym;

    if (i < xm * ym * zm + xmym) {
        res[i - xmym] = epsi[i - 1   ] * phi[i - 1   ] + epsi[i       ] * phi[i + 1   ] +
                        epsj[i - xm  ] * phi[i - xm  ] + epsj[i       ] * phi[i + xm  ] +
                        epsk[i - xmym] * phi[i - xmym] + epsk[i       ] * phi[i + xmym] + rho[i - xmym] - epsc[i - xmym] * phi[i];
    }

}

//*****************************************************
// Recursive V-Cycle
__host__
void VCycle(int level) {

    if (level == MG_NLEVEL - 1) {
        // Solve on coarsest grid
        relax(level, -1);
    } else {
        // Relax & restrict
        if (level < threshold) {
            // On CUDA
            relax_cuda(level, ncyc_before);
            restrict_cuda(level);

        } else {
            // On CPU
            relax(level, ncyc_before);

            restrict_v(16.0, l_rv + mg_index[level], mg_size[level][0], mg_size[level][1], mg_size[level][2], l_bv + mg_index[level + 1], mg_size[level + 1][0], mg_size[level + 1][1], mg_size[level + 1][2]); //bv

        }

        // Reinitialize l_xv on level+1 to zero
        int vnx = mg_size[level + 1][0];
        int vny = mg_size[level + 1][1];
        int vnz = mg_size[level + 1][2];
        int vnxny = vnx * vny;
        int vnxnynz = vnxny * vnz;
        int pa = mg_x_idx[level + 1] + vnxny;
        if (level < threshold - 1) {
            int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;
            int nblocks = (vnxnynz - 1) / blocksize + 1;
            init_vector_kernel<<<nblocks, blocksize>>>(l_xv + pa, vnxnynz);
            cudaLaunchErrorCheck();
            cudaDeviceSynchronize();
        } else {
            int pb = pa + vnxnynz;
            for (int i = pa; i < pb; i++) {
                l_xv[i] = 0.0;
            }
        }

        // Recursive call
        VCycle(level + 1);

        // Interpolate & relax
        if (level < threshold) {
            // On CUDA
            interpolate_cuda(level);
            relax_cuda(level, ncyc_after);

            /*@@ Debug of omega
            if (level == 0) {
                relax_cuda(level, 0); // GS
            } else { // level 1, 2
                relax_cuda(level, 2); // SOR, upgoing
            }
            */
        } else {
            // On CPU
            interpolate(level);
            relax(level, ncyc_after);
        }
    }
}

// PB/MG CUDA driver
extern "C"
void pb_mg_cuda_(float *phi_f, float *xs_f) {
    l_itn = 0;

    // Initial norm
    l_inorm = 0.0;
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    cublasErrorCheck(status);
    status = cublasSasum(handle, l_xmymzm, l_bv, 1, &l_inorm);
    cublasErrorCheck(status);

    bool mgconverged = false;
    while (!mgconverged) {
        for (int i = 0; i < MG_NLEVEL; i++) {
            mg_onorm[i] = 9.9E99;
        }

        l_itn++;

        VCycle(0);

        // Norm
        l_norm = 0.0;
        status = cublasSasum(handle, l_xmymzm, l_rv, 1, &l_norm);
        cublasErrorCheck(status);

        if (l_itn >= l_maxitn || l_norm <= l_inorm * l_accept ) {
            mgconverged = true;

            if(l_itn >= l_maxitn) {
                printf("PB_MG WARNING: maxitn exceeded (pb_mg_cuda)!\n");
            }
        }
    } // while

    // Destroy context
    cublasDestroy(handle);

    cudaMemcpy(xs_f, l_xv + l_xmym, l_xmymzm * sizeof(float), cudaMemcpyDeviceToHost);

    if (l_bcopt != 10 || l_pbkappa != 0) memcpy(phi_f, xs_f, l_xmymzm * sizeof(float));
}

// PB/SOR CUDA driver
extern "C"
void pb_sor_cuda_(float *phi_f, float *xs_f) {
    int threadsPerBlock = 512; // This shouldn't be final, RL

    int nx = l_xm; //mg_size[0][0];
    int ny = l_ym; //mg_size[0][1];
    int nz = l_zm; //mg_size[0][2];
    int nxny = nx * ny;
    int nxnynz = nxny * nz;

    int itn_checknorm = 10;
    float wsor = 1.9;//, wsor1;

    // inverse AD for fast processing later
    int blocks = (nxnynz - 1) / threadsPerBlock + 1;
    inv_vector_kernel<<<blocks, threadsPerBlock>>>(l_ad, l_zv, nxnynz);
    cudaLaunchErrorCheck();

    // initial norm
    float linorm = 0.0;
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    cublasErrorCheck(status);
    status = cublasSasum(handle, nxnynz, l_bv, 1, &linorm);
    cublasErrorCheck(status);

    bool converged = false;
    int litn = 0;
    float lnorm = 0.0;
    int sblocks = ((int)(nxnynz / 2) + (nxnynz & 1) - 1) / threadsPerBlock + 1;
    while (!converged) {
        solver_red_kernel<<<sblocks, threadsPerBlock>>>(l_xv, l_am1, l_am2, l_am3, l_zv, l_bv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        solver_black_kernel<<<sblocks, threadsPerBlock>>>(l_xv, l_am1, l_am2, l_am3, l_zv, l_bv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        //cudaDeviceSynchronize(); // Warning of bug: CUDA_EXCEPTION_15

        litn++;

        // Check convergence
        if (litn % itn_checknorm == 0) {
            // residue
            residue_kernel<<<blocks, threadsPerBlock>>>(l_xv, l_am1, l_am2, l_am3, l_ad, l_bv, nx, ny, nz, l_rv);
            cudaLaunchErrorCheck();

            // norm
            status = cublasSasum(handle, nxnynz, l_rv, 1, &lnorm);
            cublasErrorCheck(status);

            if (litn >= l_maxitn || lnorm <= l_accept * linorm) {
                converged = true;
                if (litn >= l_maxitn) printf("PB_SOR WARNING: maxitn exceeded (kernel)!\n");
            } // if
        }
    } // while

    // Destroy context
    cublasDestroy(handle);

    l_itn = litn;
    l_inorm = linorm;
    l_norm = lnorm;

    cudaMemcpy(xs_f, l_xv + l_xmym, l_xmymzm * sizeof(float), cudaMemcpyDeviceToHost);
    if (l_bcopt != 10 || l_pbkappa != 0) memcpy(phi_f, xs_f, l_xmymzm * sizeof(float));
}

// Return values
extern "C"
int get_itn_() {
    return l_itn;
}

extern "C"
float get_inorm_() {
    return l_inorm;
}

extern "C"
float get_norm_() {
    return l_norm;
}

