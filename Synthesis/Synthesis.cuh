#ifndef __RA_SYNTHESIS_CUH__
#define __RA_SYNTHESIS_CUH__

template <class T>
thrust::pair<thrust::pair<T *, T*>, T*> defineSpectralQuantities(const T, const T, const T, const T, const T, int *, int *);

void calculateApertureField(const float  * __restrict__, const float  * __restrict__, const float  * __restrict__,
											 const float  * __restrict__, float2_  * __restrict__, const float , const float , const float ,
											 const float , const int);

void calculateApertureField(const double * __restrict__, const double * __restrict__, const double * __restrict__,
											 const double * __restrict__, double2_ * __restrict__, const double, const double, const double, 
											 const double, const int);

float2_ * raFarFieldCalculation(const float * __restrict__,    const float * __restrict__, 
							    const float * __restrict__, const float * __restrict__, const float * __restrict__,   
								const float * __restrict__, const float * __restrict__,
								const float * __restrict__,
								const int, const int, const int,
								const cublasHandle_t, 
								const float, const float, const float,
								const float, const float, const float, 
								const float, const float,
								const int, const int,
								const int,  const int);

double2_ * raFarFieldCalculation(const double * __restrict__,    const double * __restrict__, 
								 const double * __restrict__, const double * __restrict__, const double * __restrict__,   
								 const double * __restrict__, const double * __restrict__,
								 const double * __restrict__,
								 const int, const int, const int,
								 const cublasHandle_t, 
								 const double, const double, const double,
								 const double, const double, const double, 
								 const double, const double,
								 const int, const int,
								 const int,  const int);

float  raFunctionalCalculation(const float2_  * __restrict__, const float  * __restrict__, const float  * __restrict__, const int, const int);
double raFunctionalCalculation(const double2_ * __restrict__, const double * __restrict__, const double * __restrict__, const int, const int);

#endif
