#include "Utilities.cuh"
#include "BBComplex.h"
#include "Matlab_like.cuh"
#include "Synthesis.cuh"
#include "NDFT2_2D.cuh"
#include "NFFT2_2D.cuh"
#include "InputOutput.cuh"

#include <thrust\device_vector.h>
#include <thrust\reduce.h>
#include <thrust\transform_reduce.h>
#include <thrust\pair.h>

#define DEBUG

#define pi 3.141592653589793238463

/********************************/
/* DEFINING SPECTRAL QUANTITIES */
/********************************/
#define BLOCKSIZE_FILTER	256

template <class T>
__global__ void filterConstructionKernel(T * __restrict__ d_Filter, const T * __restrict__ d_U, const T * __restrict__ d_V, const T beta, const int N) {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N) d_Filter[tid] = static_cast<T>((d_U[tid] * d_U[tid] + d_V[tid] * d_V[tid]) <= beta * beta);

}

template <class T>
thrust::pair<thrust::pair<T *, T*>, T*> defineSpectralQuantities(const T uMax, const T vMax, const T aPrime, const T bPrime, const T beta, int *Nu, int *Nv) {

	Nu[0] = (int)(uMax / (pi / aPrime)) + 1;
	Nv[0] = (int)(vMax / (pi / bPrime)) + 1;

	T *d_u_discrete = colon(-static_cast<T>(Nu[0]), static_cast<T>(1), static_cast<T>(Nu[0] - 1));
	T *d_v_discrete = colon(-static_cast<T>(Nv[0]), static_cast<T>(1), static_cast<T>(Nu[0] - 1));

	thrust::pair<T *, T *> d_UV_discrete = meshgrid(d_u_discrete, 2 * Nu[0], d_v_discrete, 2 * Nv[0]);
	T *d_U_discrete = d_UV_discrete.first;
	T *d_V_discrete = d_UV_discrete.second;

	T *d_u;	gpuErrchk(cudaMalloc(&d_u, (2 * Nu[0]) * sizeof(T)));
	T *d_v;	gpuErrchk(cudaMalloc(&d_v, (2 * Nv[0]) * sizeof(T)));

	gpuErrchk(cudaMemcpy(d_u, d_u_discrete, (2 * Nu[0]) * sizeof(T), cudaMemcpyDeviceToDevice));
	gpuErrchk(cudaMemcpy(d_v, d_v_discrete, (2 * Nv[0]) * sizeof(T), cudaMemcpyDeviceToDevice));
	
	vectorMulConstant(d_u, static_cast<T>(pi / aPrime), (2 * Nu[0]));
	vectorMulConstant(d_v, static_cast<T>(pi / bPrime), (2 * Nv[0]));

	thrust::pair<T *, T *> d_UV = meshgrid(d_u, 2 * Nu[0], d_v, 2 * Nv[0]);
	T *d_U = d_UV.first;
	T *d_V = d_UV.second;

	T * d_Filter; gpuErrchk(cudaMalloc(&d_Filter, (2 * Nu[0]) * (2 * Nv[0]) * sizeof(T)));
	filterConstructionKernel<<<iDivUp((2 * Nu[0]) * (2 * Nv[0]), BLOCKSIZE_FILTER), BLOCKSIZE_FILTER>>>(d_Filter, d_U, d_V, beta, (2 * Nu[0]) * (2 * Nv[0]));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	gpuErrchk(cudaFree(d_u_discrete));
	gpuErrchk(cudaFree(d_v_discrete));
	gpuErrchk(cudaFree(d_u));
	gpuErrchk(cudaFree(d_v));
	gpuErrchk(cudaFree(d_U));
	gpuErrchk(cudaFree(d_V));

	return thrust::make_pair(thrust::make_pair(d_U_discrete, d_V_discrete), d_Filter);

}

template thrust::pair<thrust::pair<float  *, float  *>, float  *> defineSpectralQuantities(const float  uMax, const float  vMax, const float  aPrime, const float  bPrime, const float  beta, int *, int *);
template thrust::pair<thrust::pair<double *, double *>, double *> defineSpectralQuantities(const double uMax, const double vMax, const double aPrime, const double bPrime, const double beta, int *, int *);

/*****************************************************************************************/
/* FUNCTIONS TO CALCULATE THE APERTURE FIELD OVER THE REFLECTARRAY SURFACE - DOUBLE CASE */
/*****************************************************************************************/
#define BLOCKSIZE_APERTUREFIELD	256

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

/**************************************/
/* REFLECTARRAY FAR-FIELD CALCULATION */
/**************************************/
#define BLOCKSIZE_FAR_FIELD_FILTERING	256
#define NFFT

__global__ void farFieldFiltering(float2_ * __restrict__ d_far_field, const float * __restrict__ d_Filter, const int N) {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N) d_far_field[tid] = d_Filter[tid] * d_far_field[tid];

}

__global__ void farFieldFiltering(double2_ * __restrict__ d_far_field, const double * __restrict__ d_Filter, const int N) {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N) d_far_field[tid] = d_Filter[tid] * d_far_field[tid];

}

float2_ * raFarFieldCalculation(const float * __restrict__ d_Coeff_Zernike,    const float * __restrict__ d_ZERNIKE, 
			 					 const float * __restrict__ d_Coeff_Lagrange_x, const float * __restrict__ d_Coeff_Lagrange_y, const float * __restrict__ d_LEG,   
								 const float * __restrict__ d_U_discrete, const float * __restrict__ d_V_discrete,
								 const float * __restrict__ d_Filter,
								 const int Num_unknowns_phases, const int Num_unknowns_x, const int Num_unknowns_y,
								 const cublasHandle_t handle, 
								 const float feedCenterX, const float feedCenterY, const float feedCenterZ,
								 const float alfa, const float beta, const float mfact, 
								 const float a_prime, const float b_prime,
								 const int M_x, const int M_y,
								 const int Nu,  const int Nv) {
	
	// --- Calculating aperture phase distribution
	float *d_phases;			gpuErrchk(cudaMalloc(&d_phases, M_x * M_y * sizeof(float)));
	linearCombination(d_Coeff_Zernike, d_ZERNIKE, d_phases, Num_unknowns_phases, M_x * M_y, handle);

	//float *h_phases = (float *)malloc(M_x * M_y * sizeof(float));
	//gpuErrchk(cudaMemcpy(h_phases, d_phases, M_x * M_y * sizeof(float), cudaMemcpyDeviceToHost));
	//for (int k=0; k<M_x * M_y; k++) printf("%i %f\n", k, h_phases[k]);

	// --- Calculating the non-uniform element locations
	float *d_X;	gpuErrchk(cudaMalloc(&d_X, M_x * M_y * sizeof(float)));	linearCombination(d_Coeff_Lagrange_x, d_LEG, d_X, Num_unknowns_x * Num_unknowns_y, M_x * M_y, handle);
	float *d_Y;	gpuErrchk(cudaMalloc(&d_Y, M_x * M_y * sizeof(float)));	linearCombination(d_Coeff_Lagrange_y, d_LEG, d_Y, Num_unknowns_x * Num_unknowns_y, M_x * M_y, handle);
	float *d_Z;	gpuErrchk(cudaMalloc(&d_Z, M_x * M_y * sizeof(float)));	gpuErrchk(cudaMemset(d_Z, 0, M_x * M_y * sizeof(float)));
	vectorAddConstant(d_X, feedCenterX, M_x * M_y);
	vectorAddConstant(d_Y, feedCenterY, M_x * M_y);
	vectorAddConstant(d_Z, feedCenterZ, M_x * M_y);

	//float *h_X = (float *)malloc(M_x * M_y * sizeof(float));
	//gpuErrchk(cudaMemcpy(h_X, d_X, M_x * M_y * sizeof(float), cudaMemcpyDeviceToHost));
	//for (int k=0; k<M_x * M_y; k++) printf("%i %f\n", k, h_X[k]);

	// --- Calculating the aperture field
	float2_ *d_excitations;		gpuErrchk(cudaMalloc(&d_excitations, M_x * M_y * sizeof(float2_)));
	calculateApertureField(d_X, d_Y, d_Z, d_phases, d_excitations, alfa, beta, mfact, 1., M_x * M_y);

	// --- Calculating the far field
	float *d_X_scaled;		gpuErrchk(cudaMalloc(&d_X_scaled, M_x * M_y * sizeof(float)));		gpuErrchk(cudaMemcpy(d_X_scaled, d_X, M_x * M_y * sizeof(float), cudaMemcpyDeviceToDevice));
	float *d_Y_scaled;		gpuErrchk(cudaMalloc(&d_Y_scaled, M_x * M_y * sizeof(float)));		gpuErrchk(cudaMemcpy(d_Y_scaled, d_Y, M_x * M_y * sizeof(float), cudaMemcpyDeviceToDevice));

	vectorMulConstant<float>(d_X_scaled, static_cast<float>(2 * Nu * (pi / a_prime)) / (static_cast<float>(2) * pi), M_x * M_y);
	vectorMulConstant<float>(d_Y_scaled, static_cast<float>(2 * Nv * (pi / b_prime)) / (static_cast<float>(2) * pi), M_x * M_y);

	float2_ *d_far_field;			gpuErrchk(cudaMalloc(&d_far_field, (2 * Nu) * (2 * Nv) * sizeof(float2_)));
	Calculate_cuFFT_plan_C2C_NFFT2_2D(2 * Nu, 2 * Nv);
#ifdef NFFT
	NFFT2_2D_GPU(reinterpret_cast<float2 *>(d_far_field), reinterpret_cast<float2 *>(d_excitations), d_X_scaled, d_Y_scaled, 2 * Nv, 2 * Nu, M_x * M_y);
#else
	NDFT2_2D_GPU(handle, d_X_scaled, d_Y_scaled, d_U_discrete, d_V_discrete, d_excitations, d_far_field, 2 * Nu, 2 * Nv, M_x * M_y, (2 * Nu) * (2 * Nv));
#endif

	farFieldFiltering<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_FAR_FIELD_FILTERING), BLOCKSIZE_FAR_FIELD_FILTERING>>>(d_far_field, d_Filter, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	return d_far_field;

}

double2_ * raFarFieldCalculation(const double * __restrict__ d_Coeff_Zernike,    const double * __restrict__ d_ZERNIKE, 
			 					 const double * __restrict__ d_Coeff_Lagrange_x, const double * __restrict__ d_Coeff_Lagrange_y, const double * __restrict__ d_LEG,   
								 const double * __restrict__ d_U_discrete, const double * __restrict__ d_V_discrete,
								 const double * __restrict__ d_Filter,
								 const int Num_unknowns_phases, const int Num_unknowns_x, const int Num_unknowns_y,
								 const cublasHandle_t handle, 
								 const double feedCenterX, const double feedCenterY, const double feedCenterZ,
								 const double alfa, const double beta, const double mfact, 
								 const double a_prime, const double b_prime,
								 const int M_x, const int M_y,
								 const int Nu,  const int Nv) {
	
	// --- Calculating aperture phase distribution
	double *d_phases;			gpuErrchk(cudaMalloc(&d_phases, M_x * M_y * sizeof(double)));
	linearCombination(d_Coeff_Zernike, d_ZERNIKE, d_phases, Num_unknowns_phases, M_x * M_y, handle);

	// --- Calculating the non-uniform element locations
	double *d_X;	gpuErrchk(cudaMalloc(&d_X, M_x * M_y * sizeof(double)));	linearCombination(d_Coeff_Lagrange_x, d_LEG, d_X, Num_unknowns_x * Num_unknowns_y, M_x * M_y, handle);
	double *d_Y;	gpuErrchk(cudaMalloc(&d_Y, M_x * M_y * sizeof(double)));	linearCombination(d_Coeff_Lagrange_y, d_LEG, d_Y, Num_unknowns_x * Num_unknowns_y, M_x * M_y, handle);
	double *d_Z;	gpuErrchk(cudaMalloc(&d_Z, M_x * M_y * sizeof(double)));	gpuErrchk(cudaMemset(d_Z, 0, M_x * M_y * sizeof(double)));
	vectorAddConstant(d_X, feedCenterX, M_x * M_y);
	vectorAddConstant(d_Y, feedCenterY, M_x * M_y);
	vectorAddConstant(d_Z, feedCenterZ, M_x * M_y);

	// --- Calculating the aperture field
	double2_ *d_excitations;		gpuErrchk(cudaMalloc(&d_excitations, M_x * M_y * sizeof(double2_)));
	calculateApertureField(d_X, d_Y, d_Z, d_phases, d_excitations, alfa, beta, mfact, 1., M_x * M_y);

	// --- Calculating the far field
	double *d_X_scaled;		gpuErrchk(cudaMalloc(&d_X_scaled, M_x * M_y * sizeof(double)));		gpuErrchk(cudaMemcpy(d_X_scaled, d_X, M_x * M_y * sizeof(double), cudaMemcpyDeviceToDevice));
	double *d_Y_scaled;		gpuErrchk(cudaMalloc(&d_Y_scaled, M_x * M_y * sizeof(double)));		gpuErrchk(cudaMemcpy(d_Y_scaled, d_Y, M_x * M_y * sizeof(double), cudaMemcpyDeviceToDevice));

	vectorMulConstant(d_X_scaled, static_cast<double>(2 * Nu) * (pi / a_prime) / (static_cast<double>(2) * pi), M_x * M_y);
	vectorMulConstant(d_Y_scaled, static_cast<double>(2 * Nv) * (pi / b_prime) / (static_cast<double>(2) * pi), M_x * M_y);

	double2_ *d_far_field;			gpuErrchk(cudaMalloc(&d_far_field, (2 * Nu) * (2 * Nv) * sizeof(double2_)));
	Calculate_cuFFT_plan_Z2Z_NFFT2_2D(2 * Nu, 2 * Nv);
#ifdef NFFT
	NFFT2_2D_GPU(reinterpret_cast<double2 *>(d_far_field), reinterpret_cast<double2 *>(d_excitations), d_X_scaled, d_Y_scaled, 2 * Nv, 2 * Nu, M_x * M_y);
#else
	NDFT2_2D_GPU(handle, d_X_scaled, d_Y_scaled, d_U_discrete, d_V_discrete, d_excitations, d_far_field, 2 * Nu, 2 * Nv, M_x * M_y, (2 * Nu) * (2 * Nv));
#endif

	farFieldFiltering<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_FAR_FIELD_FILTERING), BLOCKSIZE_FAR_FIELD_FILTERING>>>(d_far_field, d_Filter, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	return d_far_field;

}

/***************************************/
/* REFLECTARRAY FUNCTIONAL CALCULATION */
/***************************************/
#define BLOCKSIZE_ABSKERNEL		256

__global__ void absKernel(float * __restrict__ d_far_field_abs, const float2_ * __restrict__ d_far_field, const int N)  {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N) d_far_field_abs[tid] = sqrt(d_far_field[tid].c.x * d_far_field[tid].c.x + d_far_field[tid].c.y * d_far_field[tid].c.y);

}

__global__ void absKernel(double * __restrict__ d_far_field_abs, const double2_ * __restrict__ d_far_field, const int N)  {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N) d_far_field_abs[tid] = sqrt(d_far_field[tid].c.x * d_far_field[tid].c.x + d_far_field[tid].c.y * d_far_field[tid].c.y);

}

#define BLOCKSIZE_PROJECTION	256

template <class T>
__global__ void evaluateProjection(T * __restrict__ d_far_field_abs_projected, const T * __restrict__ d_far_field_abs, 
	                               const T * __restrict__ d_Internal_Coverage, const T * __restrict__ d_External_Coverage, const int N)  {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N) {

		if (d_far_field_abs[tid] * d_far_field_abs[tid] > d_External_Coverage[tid]) d_far_field_abs_projected[tid] = sqrt(d_External_Coverage[tid]);
		else 
			if (d_far_field_abs[tid] * d_far_field_abs[tid] < d_Internal_Coverage[tid]) d_far_field_abs_projected[tid] = sqrt(d_Internal_Coverage[tid]);
			else d_far_field_abs_projected[tid] = d_far_field_abs[tid];

	}
