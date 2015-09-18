#ifndef __RA_SYNTHESIS_CUH__
#define __RA_SYNTHESIS_CUH__

void calculateApertureField(const float  * __restrict__, const float  * __restrict__, const float  * __restrict__,
											 const float  * __restrict__, float2_  * __restrict__, const float , const float , const float ,
											 const float , const int);

void calculateApertureField(const double * __restrict__, const double * __restrict__, const double * __restrict__,
											 const double * __restrict__, double2_ * __restrict__, const double, const double, const double, 
											 const double, const int);

#endif
