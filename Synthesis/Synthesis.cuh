#ifndef __RA_SYNTHESIS_CUH__
#define __RA_SYNTHESIS_CUH__

template <class T>
thrust::pair<thrust::pair<T *, T*>, T*> defineSpectralQuantities(const T, const T, const T, const T, const T, int *, int *);

template <class T>
thrust::pair<thrust::pair<T *, T*>, T*> h_defineSpectralQuantities(const T, const T, const T, const T, const T, int *, int *);

void calculateApertureField(const float  * __restrict__, const float  * __restrict__, const float  * __restrict__,
											 const float  * __restrict__, float2_  * __restrict__, const float , const float , const float ,
											 const float , const int);

void calculateApertureField(const double * __restrict__, const double * __restrict__, const double * __restrict__,
											 const double * __restrict__, double2_ * __restrict__, const double, const double, const double, 
											 const double, const int);

float minimumDistance(const float * __restrict__, const float * __restrict__, const int, const int);

//float2_ * raFarFieldCalculation(float *);
//thrust::pair<float2_ *, float> raFarFieldCalculation(float *);
thrust::tuple<float2_ *, float, float> raFarFieldCalculation(float *);
float2_ * raFarFieldCalculationSaving(float *);

double2_ * raFarFieldCalculation(double *);

float  raFunctionalCalculation(const float2_  * __restrict__);
double raFunctionalCalculation(const double2_ * __restrict__);

//float raCostFunctional(float *);
//thrust::pair<float, float> raCostFunctional(float *);
thrust::tuple<float, float, float> raCostFunctional(float *);

#endif
