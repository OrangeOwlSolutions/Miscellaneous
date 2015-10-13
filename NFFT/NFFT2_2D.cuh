#ifndef NFFT2_2D_CUH
#define NFFT2_2D_CUH

extern	cufftHandle	NFFT2_2D_GPUplan;
void	Calculate_cuFFT_plan_Z2Z_NFFT2_2D(const int, const int);
void	Calculate_cuFFT_plan_C2C_NFFT2_2D(const int, const int);
void	Destroy_cuFFT_plan_NFFT2_2D();
void	NFFT2_2D_GPU(float2 * __restrict__, const float2 * __restrict__, const float * __restrict__, const float * __restrict__, 
					const int, const int, const int);
void	NFFT2_2D_GPU(double2 * __restrict__, const double2 * __restrict__, const double * __restrict__, const double * __restrict__, 
	                 const int, const int, const int);

#endif
