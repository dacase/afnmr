/*
 * Wrapped error handling for calling CUDA runtime, cuBLAS and cuSPARSE APIs from host;
 * Maximum performance GPU device pickup. 2017
 * Edit: Added kernel launch error check; updated Arch table with CUDA 9.1. 2018
 *
 * Ruxi Qi, UC Irvine
 */

#ifndef HELPER_CUDA_H
#define HELPER_CUDA_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Avoid lowercased max, confilting with c++ definition
#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif

//******************************
// CUDA runtime API error check
#define cudaErrorCheck(err) { runtimeAPIcheck((err), __FILE__, __LINE__); }
inline void runtimeAPIcheck(cudaError_t code, const char *file, const int line) {
	if (code != cudaSuccess) {
		fprintf(stderr, "CUDA runtime API call failed at %s:%d: %s\n", file, line, cudaGetErrorString(code));
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
}

//******************************
// CUDA kernel launch error check
#define cudaLaunchErrorCheck() { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA launch failed at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
}

//******************************
// cuBLAS API error check
#ifdef CUBLAS_V2_H_ // on when cublas_v2.h is included

#define cublasErrorCheck(err) { cublasSTATUScheck((err), __FILE__, __LINE__); }
static const char *_cublasGetErrorEnum(cublasStatus_t error) {
    switch (error) {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";

        case CUBLAS_STATUS_NOT_SUPPORTED:
            return "CUBLAS_STATUS_NOT_SUPPORTED";

        case CUBLAS_STATUS_LICENSE_ERROR:
            return "CUBLAS_STATUS_LICENSE_ERROR";
    }

    return "<unknown>";
}

inline void cublasSTATUScheck(cublasStatus_t code, const char *file, const int line) {
	if (code != CUBLAS_STATUS_SUCCESS) {
		fprintf(stderr, "cuBLAS API call failed at %s:%d: %s\n", file, line, _cublasGetErrorEnum(code));
		exit(EXIT_FAILURE);
	}
}
#endif // end of CUBLAS_V2_H_

//******************************
// cuSPARSE API error check
#ifdef CUSPARSE_V2_H_ // on when cusparse_v2.h is included

#define cusparseErrorCheck(err) { cusparseSTATUScheck((err), __FILE__, __LINE__); }

static const char *_cusparseGetErrorEnum(cusparseStatus_t error) {
    switch (error) {
        case CUSPARSE_STATUS_SUCCESS:
            return "CUSPARSE_STATUS_SUCCESS";

        case CUSPARSE_STATUS_NOT_INITIALIZED:
            return "CUSPARSE_STATUS_NOT_INITIALIZED";

        case CUSPARSE_STATUS_ALLOC_FAILED:
            return "CUSPARSE_STATUS_ALLOC_FAILED";

        case CUSPARSE_STATUS_INVALID_VALUE:
            return "CUSPARSE_STATUS_INVALID_VALUE";

        case CUSPARSE_STATUS_ARCH_MISMATCH:
            return "CUSPARSE_STATUS_ARCH_MISMATCH";

        case CUSPARSE_STATUS_MAPPING_ERROR:
            return "CUSPARSE_STATUS_MAPPING_ERROR";

        case CUSPARSE_STATUS_EXECUTION_FAILED:
            return "CUSPARSE_STATUS_EXECUTION_FAILED";

        case CUSPARSE_STATUS_INTERNAL_ERROR:
            return "CUSPARSE_STATUS_INTERNAL_ERROR";

        case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";

		// Also put here the return of cusparseXcsric02_zeroPivot
        case CUSPARSE_STATUS_ZERO_PIVOT:
            return "CUSPARSE_STATUS_ZERO_PIVOT";
    }

    return "<unknown>";
}

inline void cusparseSTATUScheck(cusparseStatus_t code, const char *file, const int line) {
	if (code != CUSPARSE_STATUS_SUCCESS) {
		fprintf(stderr, "cuSPARSE API call failed at %s:%d: %s\n", file, line, _cusparseGetErrorEnum(code));
		exit(EXIT_FAILURE);
	}
}
#endif // end of CUSPARSE_V2_H_

//******************************
// Map SM to cores, including Pascal architecture
inline int _MapSMtoCores(int major, int minor) {
    typedef struct {
        int SM;
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] = {
		// Follow CUDA 8.0 & 9.1 official definitions
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
        { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
        { 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
        { 0x53, 128}, // Maxwell Generation (SM 5.3) GM20x class
        { 0x60, 64 }, // Pascal Generation (SM 6.0) GP100 class
        { 0x61, 128}, // Pascal Generation (SM 6.1) GP10x class
        { 0x62, 128}, // Pascal Generation (SM 6.2) GP10x class
		{ 0x70, 64 }, // Volta Generation (SM 7.0) GV100 class
        {   -1, -1 }
    };

    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1) {
		// Convert to hexadecimal
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
            return nGpuArchCoresPerSM[index].Cores;
        }
        index++;
    }

    // Else return previous Cores
    return nGpuArchCoresPerSM[index-1].Cores;
}

//******************************
// Get maximum performance device
inline int getBestGPU() {
    int device;
    int num_devices = 0;
    int best_gpu = 0;
	int sm_per_multiproc = 0;
	int best_SM_arch = 0;
    int prohibited_devices = 0;

    unsigned long long max_compute_perf = 0;
    cudaDeviceProp deviceProp;
    cudaGetDeviceCount(&num_devices);

    cudaErrorCheck(cudaGetDeviceCount(&num_devices));

    if (num_devices == 0) {
        fprintf(stderr, "Error: no CUDA devices found.\n");
        exit(EXIT_FAILURE);
    }
	else {
		// Find best SM device
	    for (device = 0; device < num_devices; device++) {
	        cudaGetDeviceProperties(&deviceProp, device);

	        if (deviceProp.computeMode != cudaComputeModeProhibited) {
	            if (deviceProp.major > 0 && deviceProp.major < 9999) {
	                best_SM_arch = MAX(best_SM_arch, deviceProp.major);
	            }
	        }
	        else {
	            prohibited_devices++;
	        }
	    }

	    if (prohibited_devices == num_devices) {
	    	fprintf(stderr, "Error: all CUDA devices have compute mode prohibited.\n");
	    	exit(EXIT_FAILURE);
	    }

		// Find best GPU
	    for (device = 0; device < num_devices; device++) {
	        cudaGetDeviceProperties(&deviceProp, device);

	        if (deviceProp.computeMode != cudaComputeModeProhibited) {
	            if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
	                sm_per_multiproc = 1;
	            }
	            else {
	                sm_per_multiproc = _MapSMtoCores(deviceProp.major, deviceProp.minor);
	            }

	            unsigned long long compute_perf = (unsigned long long) deviceProp.multiProcessorCount * sm_per_multiproc * deviceProp.clockRate;

	            if (compute_perf  > max_compute_perf) {
	                if (best_SM_arch > 2) {
	                    if (deviceProp.major == best_SM_arch) {
	                        max_compute_perf = compute_perf;
	                        best_gpu = device;
	                    }
	                }
	                else {
	                    max_compute_perf = compute_perf;
	                    best_gpu = device;
	                }
	            }
	        }
	    }
	}
    return best_gpu;
}

//******************************
// Set best GPU for host calls
inline void setBestGPU() {
	cudaDeviceProp deviceProp;
	int devID= 0;
	devID = getBestGPU();

	cudaErrorCheck(cudaSetDevice(devID));
	cudaErrorCheck(cudaGetDeviceProperties(&deviceProp, devID));

	int version = (deviceProp.major * 0x10 + deviceProp.minor);
	// Future version should feed in the parameter of Amber min_SM_requre
	if (version < 0x20) {
		printf("PBSA requires CUDA compute capability >= 2.0\n");
		cudaDeviceReset();
		exit(EXIT_SUCCESS);
	}
}

#endif
