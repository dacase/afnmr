/*
 * Conjugate gradient method with cuda CUSPARSE/CUSP to solve sparse matrix in PB program,
 * with CG and Preconditioned CG - ILU0, IC02, Jacobi & Smooth. Tested matrix formats - CSR,
 * DIA, COO, ELL, HYB.
 *
 * Coded by Ruxi Qi @ UC Irvine, Jul 2016
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
// For timing
#include <sys/time.h>
// CUDA Runtime
#include <cuda_runtime.h>
// For converting banded matrix into CSR
#include "band2csr.h"
// Using updated (v2) interfaces for CUBLAS and CUSPARSE
#include <cusparse_v2.h>
#include <cublas_v2.h>
// For error handling and device pickup
#include "helper_cuda.h"

/********************
 * Two parts:
 * Part 1 - non-PBC
 * Part 2 - PBC
 ********************/

/*************** Part 1. non-PBC (lines 31 ~ 791) ******************/

#ifdef CUSPARSE

	/*** CG without preconditioner. Called from pb_lsolver.F90 program ***/
	#ifdef CG
		extern "C" void cuda_cg_wrapper_(float *x, float *b0, float *b1, float *b2, float *b3, float *rhs, int *bwidth, int *xm, int *ym, int *zm, int *maxitn, float *acpt, int *itn, float *residual)
		{
		    const int maxiter = *maxitn;
			// CSR matrix parameters
			int xmym = *xm * *ym;
			int N = *xm * *ym * *zm;
			int nz = N + 2 * (N - 1 + N - *xm + N - xmym);
			int *I = NULL, *J = NULL;
			float *val = NULL;
			float *band[4];

			band[0] = b0;
			band[1] = b1;
			band[2] = b2;
			band[3] = b3;

		    int k, *d_col, *d_row; // On device
		    const float tol = *acpt; // Tolerance
		    float r0, r1, alpha, beta; // CG initial parameters
			float mod_b2; // For comparing with residual r1
		    float *d_val, *d_x;
		    float *d_r, *d_p, *d_omega;
			float dot, nalpha;
		    const float constONE = 1.0;
		    const float constZERO = 0.0;

			// Select and use best GPU
			setBestGPU();

			/* Generate CSR matrix A and vector rhs (b) */
			I = (int *)malloc(sizeof(int) * (N + 1));
			J = (int *)malloc(sizeof(int) * nz);
			val = (float *)malloc(sizeof(float) * nz);
			// Caution! DO NOT allocate memory to x anymore! It was already fed by Fortran!

			// Initial approximation of solution
			for (int i = 0; i < N; i++) {
		        x[i] = 0.0;
		    }

		/*	// Check band2csr timing
			struct timeval b2cBegin, b2cEnd;
			gettimeofday(&b2cBegin, 0);
		*/
			band2csr(I, J, val, N, nz, band, *bwidth, *xm, *ym, *zm);

		/*	gettimeofday(&b2cEnd, 0);
			float b2cTime = (b2cEnd.tv_sec - b2cBegin.tv_sec) * 1000.0 + (b2cEnd.tv_usec - b2cBegin.tv_usec) / 1000.0;
			FILE *b2c_pt = fopen("band2csrTiming", "a");
			fprintf(b2c_pt, "Time spent on band2csr: %f ms\n", b2cTime);
			fclose(b2c_pt);
		*/

		/* -- Testing band2csr generated values in AX=b --
		 *
		 * Test passed. Jun 24, 28 2016
		 *
			// x from pb
			FILE *fort5 = fopen("fort.106", "r");
			float bandx[N];
			if(fort5 == NULL) {
				printf("Error reading fort.106");
				exit(0);
			}

			for (int i = 0; i < N; i++) {
				fscanf(fort5, "%f", &bandx[i]);
				//printf("bandx %E\n", bandx[i]);
				//printf("rhs %E\n", rhs[i]);
				printf("%d : ", i);
				for ( int j = I[i]; j < I[i+1]; j++) {
					printf("%E ", val[j]);
				}
				printf("\n");

				//printf("J %d\n", J[i]);
				//printf("val %E\n", val[i]);
			}
			fclose(fort5);

			// A from my band2csr, calculate A * x into b_test
			FILE *fp_bt = fopen("b.test", "w");
			float b_test[N], r_xsum;
			for (int i = 0; i < N; i++) {
				r_xsum = 0.0;
				for (int j = I[i]; j < I[i + 1]; j++) {
					r_xsum += val[j] * bandx[J[j]];
				}
				b_test[i] = r_xsum;
				fprintf(fp_bt, "b_test %E\n", b_test[i]);
			}
			fclose(fp_bt);
			exit(0);
		 //---------------------------------------------
		 */

			/* Timing begin */
			//struct timeval begin, end;
			//gettimeofday(&begin, 0);

		    /* Create CUBLAS context */
		    cublasHandle_t cublasHandle = 0;
		    cublasStatus_t cublasStatus;

			/*// Timing test
			struct timeval begin, end;
			gettimeofday(&begin, 0);
			*/

			// The culprit is here!!! This will cost 10 seconds! Nov 7, 2016
		    cublasStatus = cublasCreate(&cublasHandle);

			/*// --Timing end. Turn off for efficiency
		 	gettimeofday(&end, 0);
			float cgtime = (end.tv_sec - begin.tv_sec) * 1000.0 + (end.tv_usec - begin.tv_usec) / 1000.0;
			printf("\nTime elapse: %f ms.\n", cgtime);
			*/

		    cublasErrorCheck(cublasStatus);

		    /* Create CUSPARSE context */
		    cusparseHandle_t cusparseHandle = 0;
		    cusparseStatus_t cusparseStatus;
		    cusparseStatus = cusparseCreate(&cusparseHandle);

		    cusparseErrorCheck(cusparseStatus);

		    /* Description of the A matrix*/
		    cusparseMatDescr_t descr = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descr);

		    cusparseErrorCheck(cusparseStatus);

		    /* Define the properties of the matrix */
		    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
		    //cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_SYMMETRIC); // Avoid using this, 10x slower - need extra transpose
		    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

		    /* Allocate required memory */
		    cudaErrorCheck(cudaMalloc((void **)&d_col, nz*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_val, nz*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_r, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_p, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_omega, N*sizeof(float)));

		    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_r, rhs, N*sizeof(float), cudaMemcpyHostToDevice);

		    /* CG Algorithm. Reference. Golub and Van Loan, <Matrix Computations> */

		//	clock_t start = clock();
		    k = 0;
		    r0 = 0.0;
		    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &mod_b2); // Module b square
			r1 = mod_b2;

		    while (r1 > tol * tol * mod_b2 && k <= maxiter) {
		        k++;
		        if (k == 1) {
		            cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
		        }
		        else {
		            beta = r1/r0;
		            cublasSscal(cublasHandle, N, &beta, d_p, 1);
		            cublasSaxpy(cublasHandle, N, &constONE, d_r, 1, d_p, 1) ;
		        }

		        cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &constONE, descr, d_val, d_row, d_col, d_p, &constZERO, d_omega);
		        cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &dot);
		        alpha = r1/dot;
		        cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		        nalpha = -alpha;
		        cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		        r0 = r1;
		        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		    }

		/*	clock_t end = clock();
			FILE *fp = fopen("time-cuda.dat", "w");
			fprintf(fp, "cuda timing: %f\n", (float) (end - start) * 1000 / CLOCKS_PER_SEC);
			fclose(fp);
		*/

		    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

			*itn = k;
			*residual = r1;

		    /* Destroy contexts */
		    cusparseDestroy(cusparseHandle);
		    cublasDestroy(cublasHandle);

		    /* Free device memory */
		    free(I);
		    free(J);
		    free(val);
			// Avoid dupilcate memory deallocation with pb_lsolver.F90
		    //free(x);
		    //free(rhs);
		    cudaFree(d_col);
		    cudaFree(d_row);
		    cudaFree(d_val);
		    cudaFree(d_x);
		    cudaFree(d_r);
		    cudaFree(d_p);
		    cudaFree(d_omega);

		    // clean up all state, flush all profile data
		    cudaDeviceReset();

		}
	#endif // CG

	#ifdef PCG
		extern "C" void cuda_cg_wrapper_(float *x, float *b0, float *b1, float *b2, float *b3, float *rhs, int *bwidth, int *xm, int *ym, int *zm, int *maxitn, float *acpt, int *itn, float *residual)
		{
		    const int maxiter = *maxitn;
			// CSR matrix parameters
			int xmym = *xm * *ym;
			int N = *xm * *ym * *zm;
			int nz = N + 2 * (N - 1 + N - *xm + N - xmym);
			int *I = NULL, *J = NULL;
			float *val = NULL;
			float *band[4];

			band[0] = b0;
			band[1] = b1;
			band[2] = b2;
			band[3] = b3;

		    int k, *d_col, *d_row; // On device
		    const float tol = *acpt; // Tolerance
		    float r1, alpha, beta; // CG initial parameters
			float mod_b2; // For comparing with residual r1, Ruxi
		    float *d_val, *d_x, *d_valsIncomp;
		    float *d_z1, *d_z2, *d_rm2;
		    float *d_r, *d_p, *d_omega, *d_y;
		    float numerator, denominator, nalpha;
		    const float constONE = 1.0;
		    const float constZERO = 0.0;

			// Select and use best GPU
			setBestGPU();

			/* Generate CSR matrix A and vector rhs (b) */
			I = (int *)malloc(sizeof(int) * (N + 1));
			J = (int *)malloc(sizeof(int) * nz);
			val = (float *)malloc(sizeof(float) * nz);
			// Caution! DO NOT allocate memory to x anymore! It was already fed by Fortran!

			// Initial approximation of solution
			for (int i = 0; i < N; i++) {
		        x[i] = 0.0;
		    }
		/*
			// Check band2csr timing
			struct timeval b2cBegin, b2cEnd;
			gettimeofday(&b2cBegin, 0);
		*/
			band2csr(I, J, val, N, nz, band, *bwidth, *xm, *ym, *zm);
		/*
			gettimeofday(&b2cEnd, 0);
			float b2cTime = (b2cEnd.tv_sec - b2cBegin.tv_sec) * 1000.0 + (b2cEnd.tv_usec - b2cBegin.tv_usec) / 1000.0;
			printf("Time spent on band2csr: %f ms\n\n", b2cTime);
		*/

			/* Timing begin */
			//struct timeval begin, end;
			//gettimeofday(&begin, 0);

		    /* Create CUBLAS context */
		    cublasHandle_t cublasHandle = 0;
		    cublasStatus_t cublasStatus;
		    cublasStatus = cublasCreate(&cublasHandle);

		    cublasErrorCheck(cublasStatus);

		    /* Create CUSPARSE context */
		    cusparseHandle_t cusparseHandle = 0;
		    cusparseStatus_t cusparseStatus;
		    cusparseStatus = cusparseCreate(&cusparseHandle);

		    cusparseErrorCheck(cusparseStatus);

		    /* Allocate required memory */
		    cudaErrorCheck(cudaMalloc((void **)&d_col, nz*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_val, nz*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_r, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_valsIncomp, nz*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_z1, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_z2, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_rm2, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_p, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_omega, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_y, N*sizeof(float)));

		    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_r, rhs, N*sizeof(float), cudaMemcpyHostToDevice);

		#ifdef IC02
		    /*
			 * Preconditioned CG using IC02.
			 * Reference. NVIDIA cuSPARSE documentation, code example under csric02.
			 */

			cusparseMatDescr_t descr_A = 0;
			cusparseMatDescr_t descr_L = 0;
			csric02Info_t info_A  = 0;
			csrsv2Info_t  info_L  = 0;
			csrsv2Info_t  info_Lt = 0;
			int pBufferSize_A;
			int pBufferSize_L;
			int pBufferSize_Lt;
			int pBufferSize;
			void *pBuffer = 0;
			int structural_zero;
			int numerical_zero;
			//const float constONE = 1.;
			const cusparseSolvePolicy_t policy_A  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
			const cusparseSolvePolicy_t policy_L  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
			const cusparseSolvePolicy_t policy_Lt = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
			const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
			const cusparseOperation_t trans_Lt = CUSPARSE_OPERATION_TRANSPOSE;

			// step 1: create a descriptor which contains
			// - matrix A is base-0
			// - matrix L is base-0
			// - matrix L is lower triangular
			// - matrix L has non-unit diagonal
			cusparseCreateMatDescr(&descr_A);
			cusparseSetMatIndexBase(descr_A, CUSPARSE_INDEX_BASE_ZERO);
			cusparseSetMatType(descr_A, CUSPARSE_MATRIX_TYPE_GENERAL);

			cusparseCreateMatDescr(&descr_L);
			cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
			cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
			cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
			cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_NON_UNIT);

			// step 2: create a empty info structure
			// one info for csric02 and two info's for csrsv2
			cusparseCreateCsric02Info(&info_A);
			cusparseCreateCsrsv2Info(&info_L);
			cusparseCreateCsrsv2Info(&info_Lt);

			// step 3: query how much memory used in csric02 and csrsv2, and allocate the buffer
			cusparseScsric02_bufferSize(cusparseHandle, N, nz, descr_A, d_val, d_row, d_col, info_A, &pBufferSize_A);
			cusparseScsrsv2_bufferSize(cusparseHandle, trans_L, N, nz, descr_L, d_val, d_row, d_col, info_L, &pBufferSize_L);
			cusparseScsrsv2_bufferSize(cusparseHandle, trans_Lt, N, nz, descr_L, d_val, d_row, d_col, info_Lt,&pBufferSize_Lt);

			pBufferSize = std::max(pBufferSize_A, std::max(pBufferSize_L, pBufferSize_Lt));
			// pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
			cudaMalloc((void**)&pBuffer, pBufferSize);

			// Copy A data to IC02 vals as input - Ruxi
			cudaMemcpy(d_valsIncomp, d_val, nz*sizeof(float), cudaMemcpyDeviceToDevice);

			// step 4: perform analysis of incomplete Cholesky on A
			//         perform analysis of triangular solve on L
			//         perform analysis of triangular solve on L'
			// The lower triangular part of A has the same sparsity pattern as L, so
			// we can do analysis of csric02 and csrsv2 simultaneously.

			cusparseScsric02_analysis(cusparseHandle, N, nz, descr_A, d_val, d_row, d_col, info_A, policy_A, pBuffer);
			cusparseStatus = cusparseXcsric02_zeroPivot(cusparseHandle, info_A, &structural_zero);
			if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus) {
			   printf("A(%d,%d) is missing\n", structural_zero, structural_zero);
			}

			cusparseErrorCheck(cusparseStatus);

			cusparseScsrsv2_analysis(cusparseHandle, trans_L, N, nz, descr_L, d_val, d_row, d_col, info_L, policy_L, pBuffer);

			cusparseScsrsv2_analysis(cusparseHandle, trans_Lt, N, nz, descr_L, d_val, d_row, d_col, info_Lt, policy_Lt, pBuffer);

			// step 5: A ~= L * L'
			cusparseScsric02(cusparseHandle, N, nz, descr_A, d_valsIncomp, d_row, d_col, info_A, policy_A, pBuffer);

			cusparseStatus = cusparseXcsric02_zeroPivot(cusparseHandle, info_A, &numerical_zero);
			if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus) {
			   printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
			}

			cusparseErrorCheck(cusparseStatus);

		// =======IC02 iteration Starts=======
		    k = 0;
		    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &mod_b2); // Module b square
		//    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
			r1 = mod_b2;

		    while (r1 > tol * tol * mod_b2 && k <= maxiter) {
				// step 6: solve L * y = r
				cusparseScsrsv2_solve(cusparseHandle, trans_L, N, nz, &constONE, descr_L, d_valsIncomp, d_row, d_col, info_L, d_r, d_y, policy_L, pBuffer);

				// step 7: solve L' * z1 = y
				cusparseScsrsv2_solve(cusparseHandle, trans_Lt, N, nz, &constONE, descr_L, d_valsIncomp, d_row, d_col, info_Lt, d_y, d_z1, policy_Lt, pBuffer);

		        k++;

		        if (k == 1) {
		            cublasScopy(cublasHandle, N, d_z1, 1, d_p, 1);
		        }
		        else {
		            cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		            cublasSdot(cublasHandle, N, d_rm2, 1, d_z2, 1, &denominator);
		            beta = numerator/denominator;
		            cublasSscal(cublasHandle, N, &beta, d_p, 1);
		            cublasSaxpy(cublasHandle, N, &constONE, d_z1, 1, d_p, 1) ;
		        }

		        cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &constONE, descr_A, d_val, d_row, d_col, d_p, &constZERO, d_omega);
		        cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		        cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
		        alpha = numerator / denominator;
		        cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		        cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
		        cublasScopy(cublasHandle, N, d_z1, 1, d_z2, 1);
		        nalpha = -alpha;
		        cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		    }

		// ======= Over =======

		    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

			*itn = k;
			*residual = r1;

		/*
		// Testing 14dd template - Aug 4, 2016
			int num = N;
			float *rx = (float *)malloc(sizeof(float) * num);
			//int *rx = (int *)malloc(sizeof(int) * num);
		    cudaMemcpy(rx, d_z1, num * sizeof(float), cudaMemcpyDeviceToHost);
		    //cudaMemcpy(rx, d_row, num * sizeof(float), cudaMemcpyDeviceToHost);
			FILE *qi = fopen("tmp.dat", "w");
			for (int i = 0; i < num; i++) {
				fprintf(qi, "data: %e\n", rx[i]);
			}
			fclose(qi);
			free(rx);
		// TESTing over
		*/
			// step 8: free resources
			cudaFree(pBuffer);
			cusparseDestroyMatDescr(descr_A);
			cusparseDestroyMatDescr(descr_L);
			cusparseDestroyCsric02Info(info_A);
			cusparseDestroyCsrsv2Info(info_L);
			cusparseDestroyCsrsv2Info(info_Lt);

		#endif // IC02

		#ifdef ILU0
			/*
			 * Preconditioned CG using ILU.
			 * CG Algorithm. Reference. Golub and Van Loan, <Matrix Computations>
			 */

		    /* Description of the A matrix*/
		    cusparseMatDescr_t descr = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descr);

		    cusparseErrorCheck(cusparseStatus);

		    /* Define the properties of the matrix */
		    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
		    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

		    /* create the analysis info object for the A matrix */
		    cusparseSolveAnalysisInfo_t infoA = 0;
		    cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);

		    cusparseErrorCheck(cusparseStatus);

		    /* Perform the analysis for the Non-Transpose case */
		    cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descr, d_val, d_row, d_col, infoA);

		    cusparseErrorCheck(cusparseStatus);

		    /* Copy A data to ILU0 vals as input*/
		    cudaMemcpy(d_valsIncomp, d_val, nz*sizeof(float), cudaMemcpyDeviceToDevice);

		    /* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
		    cusparseStatus = cusparseScsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, d_valsIncomp, d_row, d_col, infoA);

		    cusparseErrorCheck(cusparseStatus);


		    /* Create info objects for the ILU0 preconditioner */
		    cusparseSolveAnalysisInfo_t info_u;
		    cusparseCreateSolveAnalysisInfo(&info_u);

		    cusparseMatDescr_t descrL = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descrL);
			cusparseErrorCheck(cusparseStatus);
		    cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL);
		    cusparseSetMatIndexBase(descrL,CUSPARSE_INDEX_BASE_ZERO);
		    cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
		    cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

		    cusparseMatDescr_t descrU = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descrU);
			cusparseErrorCheck(cusparseStatus);
		    cusparseSetMatType(descrU,CUSPARSE_MATRIX_TYPE_GENERAL);
		    cusparseSetMatIndexBase(descrU,CUSPARSE_INDEX_BASE_ZERO);
		    cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
		    cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
		    cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val, d_row, d_col, info_u);
			cusparseErrorCheck(cusparseStatus);

			// ILU0 Iteration starts
		    k = 0;
		    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &mod_b2); // Module b square
			//cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
			r1 = mod_b2;

		    while (r1 > tol * tol * mod_b2 && k <= maxiter) {
		        // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
		        cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &constONE, descrL, d_valsIncomp, d_row, d_col, infoA, d_r, d_y);

		        cusparseErrorCheck(cusparseStatus);

		        // Back Substitution
		        cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &constONE, descrU, d_valsIncomp, d_row, d_col, info_u, d_y, d_z1);

		        cusparseErrorCheck(cusparseStatus);

		        k++;

		        if (k == 1) {
		            cublasScopy(cublasHandle, N, d_z1, 1, d_p, 1);
		        }
		        else {
		            cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		            cublasSdot(cublasHandle, N, d_rm2, 1, d_z2, 1, &denominator);
		            beta = numerator/denominator;
		            cublasSscal(cublasHandle, N, &beta, d_p, 1);
		            cublasSaxpy(cublasHandle, N, &constONE, d_z1, 1, d_p, 1) ;
		        }

		        cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &constONE, descrU, d_val, d_row, d_col, d_p, &constZERO, d_omega);
		        cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		        cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
		        alpha = numerator / denominator;
		        cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		        cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
		        cublasScopy(cublasHandle, N, d_z1, 1, d_z2, 1);
		        nalpha = -alpha;
		        cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		    }

		    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

			*itn = k;
			*residual = r1;

		    /* Destroy parameters */
		    cusparseDestroySolveAnalysisInfo(infoA);
		    cusparseDestroySolveAnalysisInfo(info_u);

		#endif // ILU0

		    /* Destroy contexts */
		    cusparseDestroy(cusparseHandle);
		    cublasDestroy(cublasHandle);

		    /* Free host memory */
		    free(I);
		    free(J);
		    free(val);
			// Avoid dupilcate memory deallocation with pb_lsolver.F90
		    //free(x);
		    //free(rhs);

			/* Free device memory */
		    cudaFree(d_col);
		    cudaFree(d_row);
		    cudaFree(d_val);
		    cudaFree(d_x);
		    cudaFree(d_y);
		    cudaFree(d_r);
		    cudaFree(d_p);
		    cudaFree(d_omega);
		    cudaFree(d_valsIncomp);
		    cudaFree(d_z1);
		    cudaFree(d_z2);
		    cudaFree(d_rm2);

		    // clean up all state, flush all profile data
		    cudaDeviceReset();

			// Timing end
		/*	Turn off for efficiency
		 	gettimeofday(&end, 0);
			float cgtime = (end.tv_sec - begin.tv_sec) * 1000.0 + (end.tv_usec - begin.tv_usec) / 1000.0;
			printf("\nCG without preconditioning time elapse: %f ms.\n", cgtime);

		    printf("  Test Summary:\n");
		    printf("     Counted total of %d errors\n", nErrors);
		    printf("     qaerr1 = %f \n\n", fabs(qaerr1));
		*/
			// Note this 'exit' will cause abnormal ending of wrapper and fail the call in pb_lsolver.F90!
		    //exit((nErrors == 0 && fabs(qaerr1)<1e-5 ? EXIT_SUCCESS : EXIT_FAILURE));
		}
	#endif // PCG
#endif // CUSPARSE

#ifdef CUSP
	// CUSP library
	#include <cusp/csr_matrix.h>
	#include <cusp/coo_matrix.h>
	#include <cusp/ell_matrix.h>
	#include <cusp/dia_matrix.h>
	#include <cusp/hyb_matrix.h>
	#include <cusp/monitor.h>
	#include <cusp/krylov/cg.h>
	#include <cusp/print.h>
	// For PCG
	#include <cusp/precond/diagonal.h>
	#include <cusp/precond/ainv.h>
	#include <cusp/precond/aggregation/smoothed_aggregation.h>

	/*** CG using CUSP library. Called from pb_lsolver.F90 program ***/
	extern "C" void cuda_cg_wrapper_(float *x, float *b0, float *b1, float *b2, float *b3, float *rhs, int *bwidth, int *xm, int *ym, int *zm, int *maxitn, float *acpt, int *itn, float *residual)
	{
		// CSR matrix parameters
		int xmym = *xm * *ym;
		int N = *xm * *ym * *zm;
		int nz = N + 2 * (N - 1 + N - *xm + N - xmym);
		int *I = NULL, *J = NULL;
		float *val = NULL;
		float *band[4];

		band[0] = b0;
		band[1] = b1;
		band[2] = b2;
		band[3] = b3;

		const float tol = *acpt; // Tolerance
	    int maxiter = *maxitn;

		// Note thrust::device_vector() picks up GPU itself.

		/* Generate CSR matrix A and vector rhs (b) */
		I = (int *)malloc(sizeof(int) * (N + 1));
		J = (int *)malloc(sizeof(int) * nz);
		val = (float *)malloc(sizeof(float) * nz);
		// Caution! DO NOT allocate memory to x anymore! It was already fed by Fortran!

		// Initial approximation of solution
		for (int i = 0; i < N; i++) {
	        x[i] = 0.0;
	    }

		band2csr(I, J, val, N, nz, band, *bwidth, *xm, *ym, *zm);

		// Initialize vectors to device memory first (or cannot assign values to A)
		thrust::device_vector<int> d_I(I, I + N + 1), d_J(J, J + nz);
		thrust::device_vector<float> d_val(val, val + nz);

		// Initialize cusp matrix A on device
		cusp::csr_matrix<int, float, cusp::device_memory> A(N, N, nz);
		A.row_offsets = d_I;
		A.column_indices = d_J;
		A.values = d_val;

		#ifdef CSR
			cusp::csr_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // CSR

		#ifdef COO
			cusp::coo_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // COO

		#ifdef ELL
			cusp::ell_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // ELL

		#ifdef DIA
			cusp::dia_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // DIA

		#ifdef HYB
			cusp::hyb_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // HYB

		// Initialize x & rhs
		thrust::device_vector<float> d_x(x, x + N), d_rhs(rhs, rhs + N);

		// Wrap d_x & d_rhs using array1d
		cusp::array1d<float, cusp::device_memory> cg_x = d_x, cg_b = d_rhs;

		// Iteration control
		cusp::monitor<float> monitor(cg_b, maxiter, tol, 0, false);

		//clock_t start = clock();

		#ifdef CG
			/************* 1. CG **************/
			cusp::krylov::cg(do_A, cg_x, cg_b, monitor);
		#endif // CG

		#ifdef PCG
			/************* 2. PCG **************/
			/**** a) Setup preconditioner ****/
			#ifdef Jacobi
				// 2.1 Diagonal (aka Jacobi) preconditioning CG
				cusp::precond::diagonal<float, cusp::device_memory> M(do_A);
			#endif // Jacobi

			// 2.2.1 AINV preconditioning CG (standard dropping, tolerance .1) -aniv01
			//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> M(A, .1);

			// 2.2.2 AINV preconditioning CG (static dropping, 10 nonzeroes per row) -aniv02
			//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> M(A, 0, 10);

			// 2.2.3 AINV preconditioning CG (novel dropping, Lin strategy, p=2) -aniv03
			//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> M(A, 0, -1, true, 2);
			#ifdef Smooth
				// 2.3 Smoothed aggregation preconditioning CG
				cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> M(do_A);
			#endif // Smooth

			/**** b) Solve ****/
			cusp::krylov::cg(do_A, cg_x, cg_b, monitor, M);
		#endif // PCG

		//clock_t end = clock();
		// Copy back
		cusp::array1d<float, cusp::host_memory> tmp_x = cg_x;
		for (int i = 0; i < N; i++) {
			x[i] = tmp_x[i];
		}

		// Return itn & residual to fortran
			*itn = monitor.iteration_count();
			*residual = monitor.residual_norm();
	/*
		FILE *fp = fopen("time.dat", "w");
		fprintf(fp, "cusp timing: %f\n", (float) (end - start) * 1000 / CLOCKS_PER_SEC);
		fclose(fp);
	*/

	}
#endif // CUSP


/*************** Part 2. PBC (lines 791  ~ end) ******************/

#ifdef CUSPARSE

	/*** CG without preconditioner. Called from pb_lsolver.F90 program ***/
	#ifdef CG
		extern "C" void cuda_cg_wrapper_pbc_(float *x, float *b0, float *b1, float *b2, float *b3, float *b4, float *b5, float *b6, float *rhs, int *bwidth, int *xm, int *ym, int *zm, int *maxitn, float *acpt, int *itn, float *residual)
		{
		    const int maxiter = *maxitn;
			// CSR matrix parameters
			//int xmym = *xm * *ym;
			int N = *xm * *ym * *zm;
			//int nz = N + 2 * (N - 1 + N - *xm + N - xmym);
			int nz = *bwidth * N;
			int *I = NULL, *J = NULL;
			float *val = NULL;
			float *band[7];

			band[0] = b0;
			band[1] = b1;
			band[2] = b2;
			band[3] = b3;
			band[4] = b4;
			band[5] = b5;
			band[6] = b6;

		    int k, *d_col, *d_row; // On device
		    const float tol = *acpt; // Tolerance
		    float r0, r1, alpha, beta; // CG initial parameters
			float mod_b2; // For comparing with residual r1, Ruxi
		    float *d_val, *d_x;
		    float *d_r, *d_p, *d_omega;
			float dot, nalpha;
		    const float constONE = 1.0;
		    const float constZERO = 0.0;

			// Select and use best GPU
			setBestGPU();

			/* Generate CSR matrix A and vector rhs (b) */
			I = (int *)malloc(sizeof(int) * (N + 1));
			J = (int *)malloc(sizeof(int) * nz);
			val = (float *)malloc(sizeof(float) * nz);
			// Caution! DO NOT allocate memory to x anymore! It was already fed by Fortran!

			// Initial approximation of solution
			for (int i = 0; i < N; i++) {
		        x[i] = 0.0;
		    }

		/*	// Check band2csr timing
			struct timeval b2cBegin, b2cEnd;
			gettimeofday(&b2cBegin, 0);
		*/
			band2csr_pbc(I, J, val, N, nz, band, *bwidth, *xm, *ym, *zm);

		/*	gettimeofday(&b2cEnd, 0);
			float b2cTime = (b2cEnd.tv_sec - b2cBegin.tv_sec) * 1000.0 + (b2cEnd.tv_usec - b2cBegin.tv_usec) / 1000.0;
			FILE *b2c_pt = fopen("band2csrTiming", "a");
			fprintf(b2c_pt, "Time spent on band2csr: %f ms\n", b2cTime);
			fclose(b2c_pt);
		*/


			/* Timing begin */
			//struct timeval begin, end;
			//gettimeofday(&begin, 0);

		    /* Create CUBLAS context */
		    cublasHandle_t cublasHandle = 0;
		    cublasStatus_t cublasStatus;
		    cublasStatus = cublasCreate(&cublasHandle);

		    cublasErrorCheck(cublasStatus);

		    /* Create CUSPARSE context */
		    cusparseHandle_t cusparseHandle = 0;
		    cusparseStatus_t cusparseStatus;
		    cusparseStatus = cusparseCreate(&cusparseHandle);

		    cusparseErrorCheck(cusparseStatus);

		    /* Description of the A matrix*/
		    cusparseMatDescr_t descr = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descr);

		    cusparseErrorCheck(cusparseStatus);

		    /* Define the properties of the matrix */
		    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
		    //cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_SYMMETRIC); // Avoid using this, 10x slower - need extra transpose
		    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

		    /* Allocate required memory */
		    cudaErrorCheck(cudaMalloc((void **)&d_col, nz*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_val, nz*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_r, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_p, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_omega, N*sizeof(float)));

		    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_r, rhs, N*sizeof(float), cudaMemcpyHostToDevice);

		    /* CG Algorithm. Reference. Golub and Van Loan, <Matrix Computations> */

			//clock_t start = clock();
		    k = 0;
		    r0 = 0.0;
		    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &mod_b2); // Module b square
			//cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
			r1 = mod_b2;

		    while (r1 > tol * tol * mod_b2 && k <= maxiter) {
		        k++;
		        if (k == 1) {
		            cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
		        }
		        else {
		            beta = r1/r0;
		            cublasSscal(cublasHandle, N, &beta, d_p, 1);
		            cublasSaxpy(cublasHandle, N, &constONE, d_r, 1, d_p, 1) ;
		        }

		        cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &constONE, descr, d_val, d_row, d_col, d_p, &constZERO, d_omega);
		        cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &dot);
		        alpha = r1/dot;
		        cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		        nalpha = -alpha;
		        cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		        r0 = r1;
		        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		    }

		/*	clock_t end = clock();
			FILE *fp = fopen("time-cuda.dat", "w");
			fprintf(fp, "cuda timing: %f\n", (float) (end - start) * 1000 / CLOCKS_PER_SEC);
			fclose(fp);
		*/

		    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

			*itn = k;
			*residual = r1;

		    /* Destroy contexts */
		    cusparseDestroy(cusparseHandle);
		    cublasDestroy(cublasHandle);

		    /* Free device memory */
		    free(I);
		    free(J);
		    free(val);
			// Avoid dupilcate memory deallocation with pb_lsolver.F90
		    //free(x);
		    //free(rhs);
		    cudaFree(d_col);
		    cudaFree(d_row);
		    cudaFree(d_val);
		    cudaFree(d_x);
		    cudaFree(d_r);
		    cudaFree(d_p);
		    cudaFree(d_omega);

		    // clean up all state, flush all profile data
		    cudaDeviceReset();

			// Timing end
		/*	Turn off for efficiency
		 	gettimeofday(&end, 0);
			float cgtime = (end.tv_sec - begin.tv_sec) * 1000.0 + (end.tv_usec - begin.tv_usec) / 1000.0;
			printf("\nCG without preconditioning time elapse: %f ms.\n", cgtime);
		*/
		}
	#endif // CG

	#ifdef PCG
		extern "C" void cuda_cg_wrapper_pbc_(float *x, float *b0, float *b1, float *b2, float *b3, float *b4, float *b5, float *b6, float *rhs, int *bwidth, int *xm, int *ym, int *zm, int *maxitn, float *acpt, int *itn, float *residual)
		{
		    const int maxiter = *maxitn;
			// CSR matrix parameters
			//int xmym = *xm * *ym;
			int N = *xm * *ym * *zm;
			//int nz = N + 2 * (N - 1 + N - *xm + N - xmym);
			int nz = *bwidth * N;
			int *I = NULL, *J = NULL;
			float *val = NULL;
			float *band[7];

			band[0] = b0;
			band[1] = b1;
			band[2] = b2;
			band[3] = b3;
			band[4] = b4;
			band[5] = b5;
			band[6] = b6;

		    int k, *d_col, *d_row; // On device
		    const float tol = *acpt; // Tolerance
		    float r1, alpha, beta; // CG initial parameters
			float mod_b2; // For comparing with residual r1, Ruxi
		    float *d_val, *d_x, *d_valsIncomp;
		    float *d_z1, *d_z2, *d_rm2;
		    float *d_r, *d_p, *d_omega, *d_y;
		    float numerator, denominator, nalpha;
		    const float constONE = 1.0;
		    const float constZERO = 0.0;

			// Select and use best GPU
			setBestGPU();

			/* Generate CSR matrix A and vector rhs (b) */
			I = (int *)malloc(sizeof(int) * (N + 1));
			J = (int *)malloc(sizeof(int) * nz);
			val = (float *)malloc(sizeof(float) * nz);
			// Caution! DO NOT allocate memory to x anymore! It was already fed by Fortran!

			// Initial approximation of solution
			for (int i = 0; i < N; i++) {
		        x[i] = 0.0;
		    }
		/*
			// Check band2csr timing
			struct timeval b2cBegin, b2cEnd;
			gettimeofday(&b2cBegin, 0);
		*/
			band2csr_pbc(I, J, val, N, nz, band, *bwidth, *xm, *ym, *zm);
		/*
			gettimeofday(&b2cEnd, 0);
			float b2cTime = (b2cEnd.tv_sec - b2cBegin.tv_sec) * 1000.0 + (b2cEnd.tv_usec - b2cBegin.tv_usec) / 1000.0;
			printf("Time spent on band2csr: %f ms\n\n", b2cTime);
		*/

			/* Timing begin */
			//struct timeval begin, end;
			//gettimeofday(&begin, 0);

		    /* Create CUBLAS context */
		    cublasHandle_t cublasHandle = 0;
		    cublasStatus_t cublasStatus;
		    cublasStatus = cublasCreate(&cublasHandle);

		    cublasErrorCheck(cublasStatus);

		    /* Create CUSPARSE context */
		    cusparseHandle_t cusparseHandle = 0;
		    cusparseStatus_t cusparseStatus;
		    cusparseStatus = cusparseCreate(&cusparseHandle);

		    cusparseErrorCheck(cusparseStatus);

		    /* Allocate required memory */
		    cudaErrorCheck(cudaMalloc((void **)&d_col, nz*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
		    cudaErrorCheck(cudaMalloc((void **)&d_val, nz*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_r, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_valsIncomp, nz*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_z1, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_z2, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_rm2, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_p, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_omega, N*sizeof(float)));
		    cudaErrorCheck(cudaMalloc((void **)&d_y, N*sizeof(float)));

		    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
		    cudaMemcpy(d_r, rhs, N*sizeof(float), cudaMemcpyHostToDevice);

		#ifdef IC02
		    /*
			 * Preconditioned CG using IC02.
			 * Reference. NVIDIA cuSPARSE documentation, code example under csric02.
			 */

			cusparseMatDescr_t descr_A = 0;
			cusparseMatDescr_t descr_L = 0;
			csric02Info_t info_A  = 0;
			csrsv2Info_t  info_L  = 0;
			csrsv2Info_t  info_Lt = 0;
			int pBufferSize_A;
			int pBufferSize_L;
			int pBufferSize_Lt;
			int pBufferSize;
			void *pBuffer = 0;
			int structural_zero;
			int numerical_zero;
			//const float constONE = 1.;
			const cusparseSolvePolicy_t policy_A  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
			const cusparseSolvePolicy_t policy_L  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
			const cusparseSolvePolicy_t policy_Lt = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
			const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
			const cusparseOperation_t trans_Lt = CUSPARSE_OPERATION_TRANSPOSE;

			// step 1: create a descriptor which contains
			// - matrix A is base-0
			// - matrix L is base-0
			// - matrix L is lower triangular
			// - matrix L has non-unit diagonal
			cusparseCreateMatDescr(&descr_A);
			cusparseSetMatIndexBase(descr_A, CUSPARSE_INDEX_BASE_ZERO);
			cusparseSetMatType(descr_A, CUSPARSE_MATRIX_TYPE_GENERAL);

			cusparseCreateMatDescr(&descr_L);
			cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
			cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
			cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
			cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_NON_UNIT);

			// step 2: create a empty info structure
			// need one info for csric02 and two info's for csrsv2
			cusparseCreateCsric02Info(&info_A);
			cusparseCreateCsrsv2Info(&info_L);
			cusparseCreateCsrsv2Info(&info_Lt);

			// step 3: query how much memory used in csric02 and csrsv2, and allocate the buffer
			cusparseScsric02_bufferSize(cusparseHandle, N, nz, descr_A, d_val, d_row, d_col, info_A, &pBufferSize_A);
			cusparseScsrsv2_bufferSize(cusparseHandle, trans_L, N, nz, descr_L, d_val, d_row, d_col, info_L, &pBufferSize_L);
			cusparseScsrsv2_bufferSize(cusparseHandle, trans_Lt, N, nz, descr_L, d_val, d_row, d_col, info_Lt,&pBufferSize_Lt);

			pBufferSize = std::max(pBufferSize_A, std::max(pBufferSize_L, pBufferSize_Lt));
			// pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
			cudaMalloc((void**)&pBuffer, pBufferSize);

			// Copy A data to IC02 vals as input - Ruxi
			cudaMemcpy(d_valsIncomp, d_val, nz*sizeof(float), cudaMemcpyDeviceToDevice);

			// step 4: perform analysis of incomplete Cholesky on A
			//         perform analysis of triangular solve on L
			//         perform analysis of triangular solve on L'
			// The lower triangular part of A has the same sparsity pattern as L, so
			// we can do analysis of csric02 and csrsv2 simultaneously.

			cusparseScsric02_analysis(cusparseHandle, N, nz, descr_A, d_val, d_row, d_col, info_A, policy_A, pBuffer);
			cusparseStatus = cusparseXcsric02_zeroPivot(cusparseHandle, info_A, &structural_zero);

			if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus) {
			   printf("A(%d,%d) is missing\n", structural_zero, structural_zero);
			}

			cusparseErrorCheck(cusparseStatus);

			cusparseScsrsv2_analysis(cusparseHandle, trans_L, N, nz, descr_L, d_val, d_row, d_col, info_L, policy_L, pBuffer);

			cusparseScsrsv2_analysis(cusparseHandle, trans_Lt, N, nz, descr_L, d_val, d_row, d_col, info_Lt, policy_Lt, pBuffer);

			// step 5: A ~= L * L'
			cusparseScsric02(cusparseHandle, N, nz, descr_A, d_valsIncomp, d_row, d_col, info_A, policy_A, pBuffer);

			cusparseStatus = cusparseXcsric02_zeroPivot(cusparseHandle, info_A, &numerical_zero);
			if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus) {
			   printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
			}

			cusparseErrorCheck(cusparseStatus);

		// =======IC02 iteration Starts=======
		    k = 0;
		    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &mod_b2); // Module b square
		//    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
			r1 = mod_b2;

		    while (r1 > tol * tol * mod_b2 && k <= maxiter) {
				// step 6: solve L * y = r
				cusparseScsrsv2_solve(cusparseHandle, trans_L, N, nz, &constONE, descr_L, d_valsIncomp, d_row, d_col, info_L, d_r, d_y, policy_L, pBuffer);

				// step 7: solve L' * z1 = y
				cusparseScsrsv2_solve(cusparseHandle, trans_Lt, N, nz, &constONE, descr_L, d_valsIncomp, d_row, d_col, info_Lt, d_y, d_z1, policy_Lt, pBuffer);

		        k++;

		        if (k == 1) {
		            cublasScopy(cublasHandle, N, d_z1, 1, d_p, 1);
		        }
		        else {
		            cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		            cublasSdot(cublasHandle, N, d_rm2, 1, d_z2, 1, &denominator);
		            beta = numerator/denominator;
		            cublasSscal(cublasHandle, N, &beta, d_p, 1);
		            cublasSaxpy(cublasHandle, N, &constONE, d_z1, 1, d_p, 1) ;
		        }

		        cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &constONE, descr_A, d_val, d_row, d_col, d_p, &constZERO, d_omega);
		        cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		        cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
		        alpha = numerator / denominator;
		        cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		        cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
		        cublasScopy(cublasHandle, N, d_z1, 1, d_z2, 1);
		        nalpha = -alpha;
		        cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		    }

		// ======= Over =======

		    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

			*itn = k;
			*residual = r1;

		/*
		// Testing 14dd template - Aug 4, 2016
			int num = N;
			float *rx = (float *)malloc(sizeof(float) * num);
			//int *rx = (int *)malloc(sizeof(int) * num);
		    cudaMemcpy(rx, d_z1, num * sizeof(float), cudaMemcpyDeviceToHost);
		    //cudaMemcpy(rx, d_row, num * sizeof(float), cudaMemcpyDeviceToHost);
			FILE *qi = fopen("tmp.dat", "w");
			for (int i = 0; i < num; i++) {
				fprintf(qi, "data: %e\n", rx[i]);
			}
			fclose(qi);
			free(rx);
		// TESTing over
		*/
			// step 8: free resources
			cudaFree(pBuffer);
			cusparseDestroyMatDescr(descr_A);
			cusparseDestroyMatDescr(descr_L);
			cusparseDestroyCsric02Info(info_A);
			cusparseDestroyCsrsv2Info(info_L);
			cusparseDestroyCsrsv2Info(info_Lt);

		#endif // IC02

		#ifdef ILU0
			/*
			 * Preconditioned CG using ILU.
			 * CG Algorithm. Reference. Golub and Van Loan, <Matrix Computations>
			 */

		    /* Description of the A matrix*/
		    cusparseMatDescr_t descr = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descr);

		    cusparseErrorCheck(cusparseStatus);

		    /* Define the properties of the matrix */
		    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
		    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

		    /* create the analysis info object for the A matrix */
		    cusparseSolveAnalysisInfo_t infoA = 0;
		    cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);

		    cusparseErrorCheck(cusparseStatus);

		    /* Perform the analysis for the Non-Transpose case */
		    cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descr, d_val, d_row, d_col, infoA);

		    cusparseErrorCheck(cusparseStatus);

		    /* Copy A data to ILU0 vals as input*/
		    cudaMemcpy(d_valsIncomp, d_val, nz*sizeof(float), cudaMemcpyDeviceToDevice);

		    /* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
		    cusparseStatus = cusparseScsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, d_valsIncomp, d_row, d_col, infoA);

		    cusparseErrorCheck(cusparseStatus);

		    /* Create info objects for the ILU0 preconditioner */
		    cusparseSolveAnalysisInfo_t info_u;
		    cusparseCreateSolveAnalysisInfo(&info_u);

		    cusparseMatDescr_t descrL = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descrL);
			cusparseErrorCheck(cusparseStatus);
		    cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL);
		    cusparseSetMatIndexBase(descrL,CUSPARSE_INDEX_BASE_ZERO);
		    cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
		    cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

		    cusparseMatDescr_t descrU = 0;
		    cusparseStatus = cusparseCreateMatDescr(&descrU);
			cusparseErrorCheck(cusparseStatus);
		    cusparseSetMatType(descrU,CUSPARSE_MATRIX_TYPE_GENERAL);
		    cusparseSetMatIndexBase(descrU,CUSPARSE_INDEX_BASE_ZERO);
		    cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
		    cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
		    cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val, d_row, d_col, info_u);
			cusparseErrorCheck(cusparseStatus);

			// ILU0 Iteration starts
		    k = 0;
		    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &mod_b2); // Module b square
			//cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
			r1 = mod_b2;

		    while (r1 > tol * tol * mod_b2 && k <= maxiter) {
		        // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
		        cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &constONE, descrL, d_valsIncomp, d_row, d_col, infoA, d_r, d_y);

		        cusparseErrorCheck(cusparseStatus);

		        // Back Substitution
		        cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &constONE, descrU, d_valsIncomp, d_row, d_col, info_u, d_y, d_z1);

		        cusparseErrorCheck(cusparseStatus);

		        k++;

		        if (k == 1) {
		            cublasScopy(cublasHandle, N, d_z1, 1, d_p, 1);
		        }
		        else {
		            cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		            cublasSdot(cublasHandle, N, d_rm2, 1, d_z2, 1, &denominator);
		            beta = numerator/denominator;
		            cublasSscal(cublasHandle, N, &beta, d_p, 1);
		            cublasSaxpy(cublasHandle, N, &constONE, d_z1, 1, d_p, 1) ;
		        }

		        cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &constONE, descrU, d_val, d_row, d_col, d_p, &constZERO, d_omega);
		        cublasSdot(cublasHandle, N, d_r, 1, d_z1, 1, &numerator);
		        cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);
		        alpha = numerator / denominator;
		        cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
		        cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
		        cublasScopy(cublasHandle, N, d_z1, 1, d_z2, 1);
		        nalpha = -alpha;
		        cublasSaxpy(cublasHandle, N, &nalpha, d_omega, 1, d_r, 1);
		        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		    }

		    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

			*itn = k;
			*residual = r1;

		    /* Destroy parameters */
		    cusparseDestroySolveAnalysisInfo(infoA);
		    cusparseDestroySolveAnalysisInfo(info_u);

		#endif // ILU0

		    /* Destroy contexts */
		    cusparseDestroy(cusparseHandle);
		    cublasDestroy(cublasHandle);

		    /* Free host memory */
		    free(I);
		    free(J);
		    free(val);
			// Avoid dupilcate memory deallocation with pb_lsolver.F90
		    //free(x);
		    //free(rhs);

			/* Free device memory */
		    cudaFree(d_col);
		    cudaFree(d_row);
		    cudaFree(d_val);
		    cudaFree(d_x);
		    cudaFree(d_y);
		    cudaFree(d_r);
		    cudaFree(d_p);
		    cudaFree(d_omega);
		    cudaFree(d_valsIncomp);
		    cudaFree(d_z1);
		    cudaFree(d_z2);
		    cudaFree(d_rm2);

		    // clean up all state, flush all profile data
		    cudaDeviceReset();

			// Timing end
		/*	Turn off for efficiency
		 	gettimeofday(&end, 0);
			float cgtime = (end.tv_sec - begin.tv_sec) * 1000.0 + (end.tv_usec - begin.tv_usec) / 1000.0;
			printf("\nCG without preconditioning time elapse: %f ms.\n", cgtime);
		*/
			// Note this 'exit' will cause abnormal ending of wrapper and fail the call in pb_lsolver.F90!
		    //exit((nErrors == 0 && fabs(qaerr1)<1e-5 ? EXIT_SUCCESS : EXIT_FAILURE));
		}
	#endif // PCG
#endif // CUSPARSE

#ifdef CUSP
	// CUSP library
	#include <cusp/csr_matrix.h>
	#include <cusp/coo_matrix.h>
	#include <cusp/ell_matrix.h>
	#include <cusp/dia_matrix.h>
	#include <cusp/hyb_matrix.h>
	#include <cusp/monitor.h>
	#include <cusp/krylov/cg.h>
	#include <cusp/print.h>
	// For PCG
	#include <cusp/precond/diagonal.h>
	#include <cusp/precond/ainv.h>
	#include <cusp/precond/aggregation/smoothed_aggregation.h>

	/*** CG using CUSP library. Called from pb_lsolver.F90 program ***/
	extern "C" void cuda_cg_wrapper_pbc_(float *x, float *b0, float *b1, float *b2, float *b3, float *b4, float *b5, float *b6, float *rhs, int *bwidth, int *xm, int *ym, int *zm, int *maxitn, float *acpt, int *itn, float *residual)
	{
		// CSR matrix parameters
		//int xmym = *xm * *ym;
		int N = *xm * *ym * *zm;
		//int nz = N + 2 * (N - 1 + N - *xm + N - xmym);
		int nz = *bwidth * N;
		int *I = NULL, *J = NULL;
		float *val = NULL;
		float *band[7];

		band[0] = b0;
		band[1] = b1;
		band[2] = b2;
		band[3] = b3;
		band[4] = b4;
		band[5] = b5;
		band[6] = b6;

		const float tol = *acpt; // Tolerance
	    int maxiter = *maxitn;

		// Note thrust::device_vector() picks up GPU itself.

		/* Generate CSR matrix A and vector rhs (b) */
		I = (int *)malloc(sizeof(int) * (N + 1));
		J = (int *)malloc(sizeof(int) * nz);
		val = (float *)malloc(sizeof(float) * nz);
		// Caution! DO NOT allocate memory to x anymore! It was already fed by Fortran!
		//x = (float *)malloc(sizeof(float) * N);
	    //rhs = (float *)malloc(sizeof(float) * N);

		// Initial approximation of solution
		for (int i = 0; i < N; i++) {
	        x[i] = 0.0;
	    }

		band2csr_pbc(I, J, val, N, nz, band, *bwidth, *xm, *ym, *zm);

		// Initialize vectors to device memory first (or cannot assign values to A)
		thrust::device_vector<int> d_I(I, I + N + 1), d_J(J, J + nz);
		thrust::device_vector<float> d_val(val, val + nz);

		// Initialize cusp matrix A on device
		cusp::csr_matrix<int, float, cusp::device_memory> A(N, N, nz);
		A.row_offsets = d_I;
		A.column_indices = d_J;
		A.values = d_val;

		#ifdef CSR
			cusp::csr_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // CSR

		#ifdef COO
			cusp::coo_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // COO

		#ifdef ELL
			cusp::ell_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // ELL

		#ifdef DIA
			cusp::dia_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // DIA

		#ifdef HYB
			cusp::hyb_matrix<int, float, cusp::device_memory> do_A(A);
		#endif // HYB

		// Initialize x & rhs
		thrust::device_vector<float> d_x(x, x + N), d_rhs(rhs, rhs + N);

		// Wrap d_x & d_rhs using array1d
		cusp::array1d<float, cusp::device_memory> cg_x = d_x, cg_b = d_rhs;

		// Iteration control
		cusp::monitor<float> monitor(cg_b, maxiter, tol, 0, false);

		//clock_t start = clock();

		#ifdef CG
			/************* 1. CG **************/
			cusp::krylov::cg(do_A, cg_x, cg_b, monitor);
		#endif // CG

		#ifdef PCG
			/************* 2. PCG **************/
			/**** a) Setup preconditioner ****/
			#ifdef Jacobi
				// 2.1 Diagonal (aka Jacobi) preconditioning CG
				cusp::precond::diagonal<float, cusp::device_memory> M(do_A);
			#endif // Jacobi

			// 2.2.1 AINV preconditioning CG (standard dropping, tolerance .1) -aniv01
			//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> M(A, .1);

			// 2.2.2 AINV preconditioning CG (static dropping, 10 nonzeroes per row) -aniv02
			//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> M(A, 0, 10);

			// 2.2.3 AINV preconditioning CG (novel dropping, Lin strategy, p=2) -aniv03
			//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> M(A, 0, -1, true, 2);
			#ifdef Smooth
				// 2.3 Smoothed aggregation preconditioning CG
				cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> M(do_A);
			#endif // Smooth

			/**** b) Solve ****/
			cusp::krylov::cg(do_A, cg_x, cg_b, monitor, M);
		#endif // PCG

		//clock_t end = clock();
		// Copy back
		cusp::array1d<float, cusp::host_memory> tmp_x = cg_x;
		for (int i = 0; i < N; i++) {
			x[i] = tmp_x[i];
		}

		// Return itn & residual to fortran
			*itn = monitor.iteration_count();
			*residual = monitor.residual_norm();
	/*
		FILE *fp = fopen("time.dat", "w");
		fprintf(fp, "cusp timing: %f\n", (float) (end - start) * 1000 / CLOCKS_PER_SEC);
		fclose(fp);
	*/

	}
#endif // CUSP

