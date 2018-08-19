#ifndef CUDA_MG_WRAPPER_H_
#define CUDA_MG_WRAPPER_H_

#ifdef __cplusplus
extern "C" void init_param_c_(int *nx, int *ny, int *nz, int *p_maxitn, int *p_bcopt, float *p_accept, float *p_pbkappa, float *p_epsout, float *p_h, float *p_wsor);
extern "C" void allocate_array_cuda_(int *solvopt);
extern "C" void deallocate_array_cuda_();
extern "C" void init_array_cuda_(int *solvopt, float *epsx, float *epsy, float *epsz, float *p_bv, float *p_iv, float *p_xs);
__global__ void init_vector_kernel(float *vec, int m);
__global__ void copy_vector_kernel(float *vec, float *vec_f, int m);
__global__ void inv_vector_kernel(float *vec, float *inv, int m);
__global__ void feedepsintoam_kernel(int lxm, int lym, int lzm, float *am1, float *am2, float *am3, float *eps1, float *eps2, float *eps3);
__host__ void restrict_eps_map(float *epsx, float *epsy, float *epsz, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr);
__host__ __device__ float r_map_exp_x(float *epsx, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float r_map_exp_y(float *epsy, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float r_map_exp_z(float *epsz, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float hmav(float a, float b);
__global__ void restrict_eps_map_kernel(float *epsx, float *epsy, float *epsz, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr);
__host__ void set_am_ad(float *epsx, float *epsy, float *epsz, float *iv, float *lam1, float *lam2, float *lam3, float *lad, float *lbz, int xn, int yn, int zn, float lfactor, float epsout);
//__global__ void set_am_ad_kernel(float *epsx, float *epsy, float *epsz, float *iv, float *lam1, float *lam2, float *lam3, float *lad, float *lbz, int xn,int yn, int zn, float lfactor, float epsout);
__global__ void set_am_ad_kernel_head(float *epsx, float *epsy, float *epsz, float *lam1, float *lam2, float *lam3, int xnynzn);
__global__ void set_am_ad_kernel_body(float *lam1, float *lam2, float *lam3, float *lad, float *lbz, float *iv, int xn, int yn, int zn, float lfactor, float epsout);
__global__ void set_am_ad_kernel_tail(float *lam1, float *lam2, float *lam3, int xn, int yn, int zn);
__host__ void restrict_v(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr);
__global__ void restrict_v_kernel(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr);
__host__ __device__ int f_id(int i, int j, int k, int nx, int ny);
__host__ void restrict_cuda_(int level);
__host__ void interpolate(int level);
//__global__ void interpolate_kernel(float *v, int xn, int yn, int zn, float *vi, int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float *lbz, float epsout);
__global__ void interpolate_kernel_head(int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float epsout);
__global__ void interpolate_kernel_body(float *v, int xn, int yn, int zn, float *vi, int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float *lbz);

__host__ void ipl_chain_h(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__host__ void ipl_chain2_h(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__device__ void ipl_chain_d(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__device__ void ipl_chain2_d(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);

__host__ __device__ float ipl_comp1(float v, int l, float *lbz, float *am_1, int xnyn, int xnynzn, int shift_1);
__host__ __device__ float ipl_comp2(float v, int l, float *lbz, float *am_1, float *am_2, int xnyn, int xnynzn, int shift_1, int shift_2);
__host__ __device__ float ipl_comp3(float v, int l, float *lbz, float *am_1, float *am_2, float *am_3, int xnyn, int xnynzn, int shift_1, int shift_2, int shift_3);
__host__ void interpolate_cuda(int level);
__host__ void relax(int level, int ncyc);
__host__ void relax_cuda(int level, int ncyc);
__global__ void solver_black_kernel(float *phi, float *epsi, float *epsj, float *epsk, float *repsc, float *rho, float wsor, int xm, int ym, int zm);
__global__ void solver_red_kernel(float *phi, float *epsi, float *epsj, float *epsk, float *repsc, float *rho, float wsor, int xm, int ym, int zm);
__global__ void residue_kernel(float *phi, float *epsi, float *epsj, float *epsk, float *epsc, float *rho, int xm, int ym, int zm, float *res);
__host__ void VCycle(int level);
extern "C" void pb_mg_cuda_(float *phi_f, float *xs_f);
//Defined in F code
//extern "C" int  get_itn_(); 
//extern "C" float get_inorm();
//extern "C" void get_norm_();

#else

void init_param_c_(int *nx, int *ny, int *nz, int *p_maxitn, int *p_bcopt, float *p_accept, float *p_pbkappa, float *p_epsout, float *p_h, float *p_wsor);
void allocate_array_cuda_(int *solvopt);
void deallocate_array_cuda_();
void init_array_cuda_(int *solvopt, float *epsx, float *epsy, float *epsz, float *p_bv, float *p_iv, float *p_xs);
__global__ void init_vector_kernel(float *vec, int m);
__global__ void copy_vector_kernel(float *vec, float *vec_f, int m);
__global__ void inv_vector_kernel(float *vec, float *inv, int m);
__global__ void feedepsintoam_kernel(int lxm, int lym, int lzm, float *am1, float *am2, float *am3, float *eps1, float *eps2, float *eps3);
__host__ void restrict_eps_map(float *epsx, float *epsy, float *epsz, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr);
__host__ __device__ float r_map_exp_x(float *epsx, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float r_map_exp_y(float *epsy, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float r_map_exp_z(float *epsz, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float hmav(float a, float b);
__global__ void restrict_eps_map_kernel(float *epsx, float *epsy, float *epsz, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr);
__host__ void set_am_ad(float *epsx, float *epsy, float *epsz, float *iv, float *lam1, float *lam2, float *lam3, float *lad, float *lbz, int xn, int yn, int zn, float lfactor, float epsout);
//__global__ void set_am_ad_kernel(float *epsx, float *epsy, float *epsz, float *iv, float *lam1, float *lam2, float *lam3, float *lad, float *lbz, int xn,int yn, int zn, float lfactor, float epsout);
__global__ void set_am_ad_kernel_head(float *epsx, float *epsy, float *epsz, float *lam1, float *lam2, float *lam3, int xnynzn);
__global__ void set_am_ad_kernel_body(float *lam1, float *lam2, float *lam3, float *lad, float *lbz, float *iv, int xn, int yn, int zn, float lfactor, float epsout);
__global__ void set_am_ad_kernel_tail(float *lam1, float *lam2, float *lam3, int xn, int yn, int zn);

__host__ void restrict_v(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr);
__global__ void restrict_v_kernel(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr);
__host__ __device__ int f_id(int i, int j, int k, int nx, int ny);
__host__ void restrict_cuda_(int level);
__host__ void interpolate(int level);
//__global__ void interpolate_kernel(float *v, int xn, int yn, int zn, float *vi, int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float *lbz, float epsout);
__global__ void interpolate_kernel_head(int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float epsout);
__global__ void interpolate_kernel_body(float *v, int xn, int yn, int zn, float *vi, int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float *lbz);

__host__ void ipl_chain_h(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__host__ void ipl_chain2_h(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__device__ void ipl_chain_d(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__device__ void ipl_chain2_d(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);

__host__ __device__ float ipl_comp1(float v, int l, float *lbz, float *am_1, int xnyn, int xnynzn, int shift_1);
__host__ __device__ float ipl_comp2(float v, int l, float *lbz, float *am_1, float *am_2, int xnyn, int xnynzn, int shift_1, int shift_2);
__host__ __device__ float ipl_comp3(float v, int l, float *lbz, float *am_1, float *am_2, float *am_3, int xnyn, int xnynzn, int shift_1, int shift_2, int shift_3);
__host__ void interpolate_cuda(int level);
__host__ void relax(int level, int ncyc);
__host__ void relax_cuda(int level, int ncyc);
__global__ void solver_black_kernel(float *phi, float *epsi, float *epsj, float *epsk, float *repsc, float *rho, float wsor, int xm, int ym, int zm);
__global__ void solver_red_kernel(float *phi, float *epsi, float *epsj, float *epsk, float *repsc, float *rho, float wsor, int xm, int ym, int zm);
__global__ void residue_kernel(float *phi, float *epsi, float *epsj, float *epsk, float *epsc, float *rho, int xm, int ym, int zm, float *res);
__host__ void VCycle(int level);
void pb_mg_cuda_(float *phi_f, float *xs_f);
void pb_mg_cuda_(float *phi_f, float *xs_f);
//float get_itn_(int *itn);
//float get_inorm();
//float get_norm_(float *norm);


#endif //__cplusplus
#endif //CUDA_MG_WRAPPER_H_
