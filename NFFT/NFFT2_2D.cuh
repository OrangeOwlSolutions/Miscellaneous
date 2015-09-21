#ifndef NFFT2_2D_CUH
#define NFFT2_2D_CUH

extern "C" void NFFT2_2D_GPU(double2 * __restrict__, const double2 * __restrict__, const double * __restrict__, const double * __restrict__, 
	                         const int, const int, const int);

#endif
