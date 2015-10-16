#ifndef NDFT2_2D_CUH
#define NDFT2_2D_CUH

extern "C" void NDFT2_2D_GPU(cublasHandle_t, const double * __restrict__, const double * __restrict__, const double * __restrict__, 
							 const double * __restrict__, const double2_ * __restrict__, double2_ * __restrict__, const int, const int, 
							 const int, const int);

#endif
