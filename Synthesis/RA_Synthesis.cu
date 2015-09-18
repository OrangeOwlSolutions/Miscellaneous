#include "Utilities.cuh"
#include "BBComplex.h"
#include "RA_Synthesis.cuh"


#define BLOCKSIZE_APERTUREFIELD	256

/*****************************************************************************************/
/* FUNCTIONS TO CALCULATE THE APERTURE FIELD OVER THE REFLECTARRAY SURFACE - DOUBLE CASE */
/*****************************************************************************************/
__global__ void calculateApertureFieldKernel(const double * __restrict__ d_X, const double * __restrict__ d_Y, const double * __restrict__ d_Z,
											 const double * __restrict__ d_phases, double2_ * __restrict__ d_excitations, const double alpha, 
											 const double beta, const double mfact, const double scale_factor, const int N) {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (tid < N) {
		
		const double2_		im_unit(0., 1.);				// Immaginary unit

		double x_third		= -d_X[tid];
		double y_third		= d_Y[tid] * cos(alpha) + d_Z[tid] * sin(alpha);
		double z_third		= d_Y[tid] * sin(alpha) - d_Z[tid] * cos(alpha);
	
		double r			= sqrt(x_third * x_third + y_third * y_third + z_third * z_third);
		double theta		= asin(sqrt(x_third * x_third + y_third * y_third) / r);

		d_excitations[tid]	= pow(cos(theta), mfact) * exp(-im_unit * beta * r) * exp(im_unit * d_phases[tid]) / (scale_factor * r);	 
	
	}

}

void calculateApertureField(const double * __restrict__ d_X, const double * __restrict__ d_Y, const double * __restrict__ d_Z,
											 const double * __restrict__ d_phases, double2_ * __restrict__ d_excitations, const double alpha,
											 const double beta, const double mfact, const double scale_factor, const int N) {

	calculateApertureFieldKernel<<<iDivUp(N, BLOCKSIZE_APERTUREFIELD), BLOCKSIZE_APERTUREFIELD>>>(d_X, d_Y, d_Z, d_phases, d_excitations, alpha, beta, mfact, scale_factor, N);

}

/****************************************************************************************/
/* FUNCTIONS TO CALCULATE THE APERTURE FIELD OVER THE REFLECTARRAY SURFACE - FLOAT CASE */
/****************************************************************************************/
__global__ void calculateApertureFieldKernel(const float  * __restrict__ d_X, const float  * __restrict__ d_Y, const float  * __restrict__ d_Z,
											 const float  * __restrict__ d_phases, float2_ * __restrict__ d_excitations, const float alpha, 
											 const float beta, const float mfact, const float scale_factor, const int N) {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (tid < N) {
		
		const float2_		im_unit(0.f, 1.f);				// Immaginary unit

		float x_third		= -d_X[tid];
		float y_third		= d_Y[tid] * cos(alpha) + d_Z[tid] * sin(alpha);
		float z_third		= d_Y[tid] * sin(alpha) - d_Z[tid] * cos(alpha);
	
		float r				= sqrt(x_third * x_third + y_third * y_third + z_third * z_third);
		float theta			= asin(sqrt(x_third * x_third + y_third * y_third) / r);

		d_excitations[tid]	= pow(cos(theta), mfact) * exp(-im_unit * beta * r) * exp(im_unit * d_phases[tid]) / (scale_factor * r);	 
	
	}

}

void calculateApertureField(const float  * __restrict__ d_X, const float  * __restrict__ d_Y, const float  * __restrict__ d_Z,
											 const float  * __restrict__ d_phases, float2_  * __restrict__ d_excitations, const float alpha,
											 const float beta, const float mfact, const float  scale_factor, const int N) {

	calculateApertureFieldKernel<<<iDivUp(N, BLOCKSIZE_APERTUREFIELD), BLOCKSIZE_APERTUREFIELD>>>(d_X, d_Y, d_Z, d_phases, d_excitations, alpha, beta, mfact, scale_factor, N);

}
