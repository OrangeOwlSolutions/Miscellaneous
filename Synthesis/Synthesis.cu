// --- WARNING: float AND double VERSIONS ARE NOT ALIGNED. float VERSIONS ARE THE NEWEST!

#include <float.h>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/pair.h>
#include <thrust/tuple.h>

#include "Utilities.cuh"
#include "BBComplex.h"
#include "Matlab_like.cuh"
#include "Synthesis.cuh"
#include "NDFT2_2D.cuh"
#include "NFFT2_2D.cuh"
#include "InputOutput.cuh"

#define DEBUG
#define DEBUG_SAVING
//#define DEBUG_PRINTING

#define pi 3.141592653589793238463

extern const int	M_x;							// Number of reflectarray elements along the x-axis
extern const int	M_y;							// Number of reflectarray elements along the y-axis
extern const int	Num_unknowns_x;					// Number of unknowns for the element positions along the x-axis
extern const int	Num_unknowns_y;					// Number of unknowns for the element positions along the y-axis
extern const int	Num_unknowns_phases;			// Number of unknowns for the phase representation
			 int	Nu;
			 int    Nv;
extern const float	feed_center_x_f;								
extern const float	feed_center_y_f;
extern const float	feed_center_z_f;
extern const float  alfa_f;							// Feed illumination angle
extern const float  beta_f;							// Wavenumber
extern const float  mfact_f;                        // Feed pattern: cos^mfact(theta)
extern const float  a_prime_f;
extern const float  b_prime_f;
extern const double	feed_center_x;								
extern const double	feed_center_y;
extern const double	feed_center_z;
extern const double beta;							// Wavenumber
extern const double alfa;							// Feed illumination angle
extern const double mfact;                          // Feed pattern: cos^mfact(theta)
extern const double a_prime;
extern const double b_prime;
cublasHandle_t cublasHandleSynthesis; 

extern float  *d_U_discrete_f;
extern float  *d_V_discrete_f;
extern float  *d_Filter_f;
extern float  *d_LEG_f;
extern float  *d_ZERNIKE_f;
extern float  *d_Internal_Coverage_f;
extern float  *d_External_Coverage_f;
extern double *d_U_discrete;
extern double *d_V_discrete;
extern double *d_Filter;
extern double *d_LEG;
extern double *d_ZERNIKE;
extern double *d_Internal_Coverage;
extern double *d_External_Coverage;

/**************************************/
/* DEFINING SPECTRAL QUANTITIES - GPU */
/**************************************/
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

/**************************************/
/* DEFINING SPECTRAL QUANTITIES - CPU */
/**************************************/
template <class T>
thrust::pair<thrust::pair<T *, T*>, T*> h_defineSpectralQuantities(const T uMax, const T vMax, const T aPrime, const T bPrime, const T beta, int *Nu, int *Nv) {

	Nu[0] = (int)(uMax / (pi / aPrime)) + 1;
	Nv[0] = (int)(vMax / (pi / bPrime)) + 1;

	T *h_u_discrete = h_colon(-static_cast<T>(Nu[0]), static_cast<T>(1), static_cast<T>(Nu[0] - 1));
	T *h_v_discrete = h_colon(-static_cast<T>(Nv[0]), static_cast<T>(1), static_cast<T>(Nu[0] - 1));

	thrust::pair<T *, T *> h_UV_discrete = h_meshgrid(h_u_discrete, 2 * Nu[0], h_v_discrete, 2 * Nv[0]);
	T *h_U_discrete = h_UV_discrete.first;
	T *h_V_discrete = h_UV_discrete.second;

	T *h_u = (T *)malloc((2 * Nu[0]) * sizeof(T));
	T *h_v = (T *)malloc((2 * Nv[0]) * sizeof(T));

	memcpy(h_u, h_u_discrete, (2 * Nu[0]) * sizeof(T));
	memcpy(h_v, h_v_discrete, (2 * Nv[0]) * sizeof(T));
	
	h_vectorMulConstant(h_u, static_cast<T>(pi / aPrime), (2 * Nu[0]));
	h_vectorMulConstant(h_v, static_cast<T>(pi / bPrime), (2 * Nv[0]));

	thrust::pair<T *, T *> h_UV = h_meshgrid(h_u, 2 * Nu[0], h_v, 2 * Nv[0]);
	T *h_U = h_UV.first;
	T *h_V = h_UV.second;

	T * h_Filter = (T *)malloc((2 * Nu[0]) * (2 * Nv[0]) * sizeof(T));
	for (int i = 0; i < (2 * Nu[0]) * (2 * Nv[0]); i++) h_Filter[i] = static_cast<T>((h_U[i] * h_U[i] + h_V[i] * h_V[i]) <= beta * beta);

	free(h_u_discrete);
	free(h_v_discrete);
	free(h_u);
	free(h_v);
	free(h_U);
	free(h_V);

	return thrust::make_pair(thrust::make_pair(h_U_discrete, h_V_discrete), h_Filter);

}

template thrust::pair<thrust::pair<float  *, float  *>, float  *> h_defineSpectralQuantities(const float  uMax, const float  vMax, const float  aPrime, const float  bPrime, const float  beta, int *, int *);
template thrust::pair<thrust::pair<double *, double *>, double *> h_defineSpectralQuantities(const double uMax, const double vMax, const double aPrime, const double bPrime, const double beta, int *, int *);

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

/****************************************************************************************/
/* KERNEL FUNCTION TO COMPUTE THE RECIPROCAL DISTANCES FOR MINIMUM DISTANCE CALCULATION */
/****************************************************************************************/
#define BLOCKSIZE_MIN_RECIPROCAL_DISTANCES 256

__global__ void computeMinReciprocalDistances(const float * __restrict__ d_X, const float * __restrict__ d_Y, float * __restrict__ d_distances,
                                              const int p, const int M_x, const int M_y) {

    const int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < M_x * M_y) {
		if (i == p)		d_distances[i] = FLT_MAX;
		else			d_distances[i] = sqrt((d_X[i] - d_X[p]) * (d_X[i] - d_X[p]) + (d_Y[i] - d_Y[p]) * (d_Y[i] - d_Y[p]));
	}
}

/*********************************************************/
/* RECIPROCAL DISTANCES FOR MINIMUM DISTANCE CALCULATION */
/*********************************************************/
float minimumDistance(const float * __restrict__ d_X, const float * __restrict__ d_Y, const int M_x, const int M_y) {

    thrust::device_vector<float> d_distances(M_x * M_x);

	float min_distance = FLT_MAX;
	for (int p = 0; p < M_x * M_y; p++) {
		computeMinReciprocalDistances<<<iDivUp(M_x * M_y, BLOCKSIZE_MIN_RECIPROCAL_DISTANCES), BLOCKSIZE_MIN_RECIPROCAL_DISTANCES>>>(d_X, d_Y, thrust::raw_pointer_cast(d_distances.data()), p, M_x, M_y);
#ifdef DEBUG
		gpuErrchk(cudaPeekAtLastError());
		gpuErrchk(cudaDeviceSynchronize());
#endif
		float temp_min_distance = thrust::reduce(d_distances.begin(), d_distances.end(), FLT_MAX, thrust::minimum<float>());
		if (temp_min_distance < min_distance) min_distance = temp_min_distance;
	}
	
	return min_distance;
}

/****************************************************************************************/
/* KERNEL FUNCTION TO COMPUTE THE RECIPROCAL DISTANCES FOR MAXIMUM DISTANCE CALCULATION */
/****************************************************************************************/
#define BLOCKSIZE_MAX_RECIPROCAL_DISTANCES_X 16
#define BLOCKSIZE_MAX_RECIPROCAL_DISTANCES_Y 16

__global__ void computeMaxReciprocalDistances(const float * __restrict__ d_X, const float * __restrict__ d_Y, float * __restrict__ d_distances,
                                              const int M_x, const int M_y) {

    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    const int j = threadIdx.y + blockIdx.y * blockDim.y;

    float distance = 0.f;
	float temp_distance;

	if ((i > 0) && (i < M_x - 1) && (j > 0) && (j < M_y - 1)) {
        // --- Left distance
		distance      = sqrt((d_X[j * M_x + i] - d_X[j * M_x + i - 1]) * (d_X[j * M_x + i] - d_X[j * M_x + i - 1]) + 
			                 (d_Y[j * M_x + i] - d_Y[j * M_x + i - 1]) * (d_Y[j * M_x + i] - d_Y[j * M_x + i - 1]));
        // --- Right distance
        temp_distance = sqrt((d_X[j * M_x + i] - d_X[j * M_x + i + 1]) * (d_X[j * M_x + i] - d_X[j * M_x + i + 1]) + 
			                 (d_Y[j * M_x + i] - d_Y[j * M_x + i + 1]) * (d_Y[j * M_x + i] - d_Y[j * M_x + i + 1]));
		if (temp_distance > distance) distance = temp_distance;
        
        // --- Up distance
        temp_distance = sqrt((d_X[j * M_x + i] - d_X[(j + 1) * M_x + i]) * (d_X[j * M_x + i] - d_X[(j + 1) * M_x + i]) + 
			                 (d_Y[j * M_x + i] - d_Y[(j + 1) * M_x + i]) * (d_Y[j * M_x + i] - d_Y[(j + 1) * M_x + i]));
		if (temp_distance > distance) distance = temp_distance;
		
        // --- Down distance
        temp_distance = sqrt((d_X[j * M_x + i] - d_X[(j - 1) * M_x + i]) * (d_X[j * M_x + i] - d_X[(j - 1) * M_x + i]) + 
			                 (d_Y[j * M_x + i] - d_Y[(j - 1) * M_x + i]) * (d_Y[j * M_x + i] - d_Y[(j - 1) * M_x + i]));
		if (temp_distance > distance) distance = temp_distance;

		d_distances[j * M_x + i] = distance;
	}
}

/*********************************************************/
/* RECIPROCAL DISTANCES FOR MINIMUM DISTANCE CALCULATION */
/*********************************************************/
float maximumDistance(const float * __restrict__ d_X, const float * __restrict__ d_Y, const int M_x, const int M_y) {

    thrust::device_vector<float> d_distances(M_x * M_x, 0.f);

	dim3 GridDim(iDivUp(M_x, BLOCKSIZE_MAX_RECIPROCAL_DISTANCES_X), iDivUp(M_y, BLOCKSIZE_MAX_RECIPROCAL_DISTANCES_Y));
	dim3 BlockDim(BLOCKSIZE_MAX_RECIPROCAL_DISTANCES_X, BLOCKSIZE_MAX_RECIPROCAL_DISTANCES_Y);

	computeMaxReciprocalDistances<<<GridDim, BlockDim>>>(d_X, d_Y, thrust::raw_pointer_cast(d_distances.data()), M_x, M_y);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	float max_distance = thrust::reduce(d_distances.begin(), d_distances.end(), 0.f, thrust::maximum<float>());
	
	return max_distance;
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

//thrust::pair<float2_ *, float> raFarFieldCalculation(float *d_in) {
thrust::tuple<float2_ *, float, float> raFarFieldCalculation(float *d_in) {
	
	float *d_Coeff_Zernike = d_in;
	float *d_Coeff_Lagrange_x = &d_in[Num_unknowns_phases];
	float *d_Coeff_Lagrange_y = &d_in[Num_unknowns_phases + Num_unknowns_x * Num_unknowns_y];
									 
	// --- Calculating aperture phase distribution
	float *d_phases;			gpuErrchk(cudaMalloc(&d_phases, M_x * M_y * sizeof(float)));
	linearCombination(d_Coeff_Zernike, d_ZERNIKE_f, d_phases, Num_unknowns_phases, M_x * M_y, cublasHandleSynthesis);

	// --- Calculating the non-uniform element locations
	float *d_X;	gpuErrchk(cudaMalloc(&d_X, M_x * M_y * sizeof(float)));	linearCombination(d_Coeff_Lagrange_x, d_LEG_f, d_X, Num_unknowns_x * Num_unknowns_y, M_x * M_y, cublasHandleSynthesis);
	float *d_Y;	gpuErrchk(cudaMalloc(&d_Y, M_x * M_y * sizeof(float)));	linearCombination(d_Coeff_Lagrange_y, d_LEG_f, d_Y, Num_unknowns_x * Num_unknowns_y, M_x * M_y, cublasHandleSynthesis);
	float *d_Z;	gpuErrchk(cudaMalloc(&d_Z, M_x * M_y * sizeof(float)));	gpuErrchk(cudaMemset(d_Z, 0, M_x * M_y * sizeof(float)));
	vectorAddConstant(d_X, feed_center_x_f, M_x * M_y);
	vectorAddConstant(d_Y, feed_center_y_f, M_x * M_y);
	vectorAddConstant(d_Z, feed_center_z_f, M_x * M_y);

	float min_distance = minimumDistance(d_X, d_Y, M_x, M_y);
	float max_distance = maximumDistance(d_X, d_Y, M_x, M_y);

	// --- Calculating the aperture field
	float2_ *d_excitations;		gpuErrchk(cudaMalloc(&d_excitations, M_x * M_y * sizeof(float2_)));
	calculateApertureField(d_X, d_Y, d_Z, d_phases, d_excitations, alfa_f, beta_f, mfact_f, 1., M_x * M_y);

	// --- Calculating the far field
	float *d_X_scaled;		gpuErrchk(cudaMalloc(&d_X_scaled, M_x * M_y * sizeof(float)));		gpuErrchk(cudaMemcpy(d_X_scaled, d_X, M_x * M_y * sizeof(float), cudaMemcpyDeviceToDevice));
	float *d_Y_scaled;		gpuErrchk(cudaMalloc(&d_Y_scaled, M_x * M_y * sizeof(float)));		gpuErrchk(cudaMemcpy(d_Y_scaled, d_Y, M_x * M_y * sizeof(float), cudaMemcpyDeviceToDevice));

	vectorMulConstant<float>(d_X_scaled, static_cast<float>(2 * Nu * (pi / a_prime_f)) / (static_cast<float>(2) * pi), M_x * M_y);
	vectorMulConstant<float>(d_Y_scaled, static_cast<float>(2 * Nv * (pi / b_prime_f)) / (static_cast<float>(2) * pi), M_x * M_y);

	float2_ *d_far_field;			gpuErrchk(cudaMalloc(&d_far_field, (2 * Nu) * (2 * Nv) * sizeof(float2_)));
	//Calculate_cuFFT_plan_C2C_NFFT2_2D(2 * Nu, 2 * Nv);
#ifdef NFFT
	NFFT2_2D_GPU(reinterpret_cast<float2 *>(d_far_field), reinterpret_cast<float2 *>(d_excitations), d_X_scaled, d_Y_scaled, 2 * Nv, 2 * Nu, M_x * M_y);
#else
	NDFT2_2D_GPU(cublasHandleSynthesis, d_X_scaled, d_Y_scaled, d_U_discrete_f, d_V_discrete_f, d_excitations, d_far_field, 2 * Nu, 2 * Nv, M_x * M_y, (2 * Nu) * (2 * Nv));
#endif

	farFieldFiltering<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_FAR_FIELD_FILTERING), BLOCKSIZE_FAR_FIELD_FILTERING>>>(d_far_field, d_Filter_f, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	gpuErrchk(cudaFree(d_phases));
	gpuErrchk(cudaFree(d_excitations));
	gpuErrchk(cudaFree(d_X)); gpuErrchk(cudaFree(d_Y)); gpuErrchk(cudaFree(d_Z));
	gpuErrchk(cudaFree(d_X_scaled)); gpuErrchk(cudaFree(d_Y_scaled)); 
	
	//return thrust::make_pair(d_far_field, min_distance);
	return thrust::make_tuple(d_far_field, min_distance, max_distance);
}

float2_ * raFarFieldCalculationSaving(float *d_in) {
	
	float *d_Coeff_Zernike = d_in;
	float *d_Coeff_Lagrange_x = &d_in[Num_unknowns_phases];
	float *d_Coeff_Lagrange_y = &d_in[Num_unknowns_phases + Num_unknowns_x * Num_unknowns_y];
									 
	// --- Calculating aperture phase distribution
	float *d_phases;			gpuErrchk(cudaMalloc(&d_phases, M_x * M_y * sizeof(float)));
	linearCombination(d_Coeff_Zernike, d_ZERNIKE_f, d_phases, Num_unknowns_phases, M_x * M_y, cublasHandleSynthesis);

	saveGPUrealtxt(d_phases,	"/home/angelo/cuda-workspace/ParticleSwarmSynthesis/Release/PHASES_CUDA.txt", M_x * M_y);

	// --- Calculating the non-uniform element locations
	float *d_X;	gpuErrchk(cudaMalloc(&d_X, M_x * M_y * sizeof(float)));	linearCombination(d_Coeff_Lagrange_x, d_LEG_f, d_X, Num_unknowns_x * Num_unknowns_y, M_x * M_y, cublasHandleSynthesis);
	float *d_Y;	gpuErrchk(cudaMalloc(&d_Y, M_x * M_y * sizeof(float)));	linearCombination(d_Coeff_Lagrange_y, d_LEG_f, d_Y, Num_unknowns_x * Num_unknowns_y, M_x * M_y, cublasHandleSynthesis);
	float *d_Z;	gpuErrchk(cudaMalloc(&d_Z, M_x * M_y * sizeof(float)));	gpuErrchk(cudaMemset(d_Z, 0, M_x * M_y * sizeof(float)));
	vectorAddConstant(d_X, feed_center_x_f, M_x * M_y);
	vectorAddConstant(d_Y, feed_center_y_f, M_x * M_y);
	vectorAddConstant(d_Z, feed_center_z_f, M_x * M_y);

	float min_distance = minimumDistance(d_X, d_Y, M_x, M_y);
	printf("Minimum interelement distance = %f\n", min_distance);

	float max_distance = maximumDistance(d_X, d_Y, M_x, M_y);
	printf("Maximum interelement distance = %f\n", max_distance);

	saveGPUrealtxt(d_X,	"/home/angelo/cuda-workspace/ParticleSwarmSynthesis/Release/X_CUDA.txt", M_x * M_y);
	saveGPUrealtxt(d_Y,	"/home/angelo/cuda-workspace/ParticleSwarmSynthesis/Release/Y_CUDA.txt", M_x * M_y);
	saveGPUrealtxt(d_Z,	"/home/angelo/cuda-workspace/ParticleSwarmSynthesis/Release/Z_CUDA.txt", M_x * M_y);

	// --- Calculating the aperture field
	float2_ *d_excitations;		gpuErrchk(cudaMalloc(&d_excitations, M_x * M_y * sizeof(float2_)));
	calculateApertureField(d_X, d_Y, d_Z, d_phases, d_excitations, alfa_f, beta_f, mfact_f, 1., M_x * M_y);

	saveGPUcomplextxt(d_excitations,	"/home/angelo/cuda-workspace/ParticleSwarmSynthesis/Release/Excitations_CUDA.txt", M_x * M_y);

	// --- Calculating the far field
	float *d_X_scaled;		gpuErrchk(cudaMalloc(&d_X_scaled, M_x * M_y * sizeof(float)));		gpuErrchk(cudaMemcpy(d_X_scaled, d_X, M_x * M_y * sizeof(float), cudaMemcpyDeviceToDevice));
	float *d_Y_scaled;		gpuErrchk(cudaMalloc(&d_Y_scaled, M_x * M_y * sizeof(float)));		gpuErrchk(cudaMemcpy(d_Y_scaled, d_Y, M_x * M_y * sizeof(float), cudaMemcpyDeviceToDevice));

	vectorMulConstant<float>(d_X_scaled, static_cast<float>(2 * Nu * (pi / a_prime_f)) / (static_cast<float>(2) * pi), M_x * M_y);
	vectorMulConstant<float>(d_Y_scaled, static_cast<float>(2 * Nv * (pi / b_prime_f)) / (static_cast<float>(2) * pi), M_x * M_y);

	float2_ *d_far_field;			gpuErrchk(cudaMalloc(&d_far_field, (2 * Nu) * (2 * Nv) * sizeof(float2_)));
	//Calculate_cuFFT_plan_C2C_NFFT2_2D(2 * Nu, 2 * Nv);
#ifdef NFFT
	NFFT2_2D_GPU(reinterpret_cast<float2 *>(d_far_field), reinterpret_cast<float2 *>(d_excitations), d_X_scaled, d_Y_scaled, 2 * Nv, 2 * Nu, M_x * M_y);
#else
	NDFT2_2D_GPU(cublasHandleSynthesis, d_X_scaled, d_Y_scaled, d_U_discrete_f, d_V_discrete_f, d_excitations, d_far_field, 2 * Nu, 2 * Nv, M_x * M_y, (2 * Nu) * (2 * Nv));
#endif

	farFieldFiltering<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_FAR_FIELD_FILTERING), BLOCKSIZE_FAR_FIELD_FILTERING>>>(d_far_field, d_Filter_f, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	return d_far_field;

}

double2_ * raFarFieldCalculation(double *d_in) {
	
	double *d_Coeff_Zernike = d_in;
	double *d_Coeff_Lagrange_x = &d_in[Num_unknowns_phases];
	double *d_Coeff_Lagrange_y = &d_in[Num_unknowns_phases + Num_unknowns_x * Num_unknowns_y];

	// --- Calculating aperture phase distribution
	double *d_phases;			gpuErrchk(cudaMalloc(&d_phases, M_x * M_y * sizeof(double)));
	linearCombination(d_Coeff_Zernike, d_ZERNIKE, d_phases, Num_unknowns_phases, M_x * M_y, cublasHandleSynthesis);

	// --- Calculating the non-uniform element locations
	double *d_X;	gpuErrchk(cudaMalloc(&d_X, M_x * M_y * sizeof(double)));	linearCombination(d_Coeff_Lagrange_x, d_LEG, d_X, Num_unknowns_x * Num_unknowns_y, M_x * M_y, cublasHandleSynthesis);
	double *d_Y;	gpuErrchk(cudaMalloc(&d_Y, M_x * M_y * sizeof(double)));	linearCombination(d_Coeff_Lagrange_y, d_LEG, d_Y, Num_unknowns_x * Num_unknowns_y, M_x * M_y, cublasHandleSynthesis);
	double *d_Z;	gpuErrchk(cudaMalloc(&d_Z, M_x * M_y * sizeof(double)));	gpuErrchk(cudaMemset(d_Z, 0, M_x * M_y * sizeof(double)));
	vectorAddConstant(d_X, feed_center_x, M_x * M_y);
	vectorAddConstant(d_Y, feed_center_y, M_x * M_y);
	vectorAddConstant(d_Z, feed_center_z, M_x * M_y);

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
	NDFT2_2D_GPU(cublasHandleSynthesis, d_X_scaled, d_Y_scaled, d_U_discrete, d_V_discrete, d_excitations, d_far_field, 2 * Nu, 2 * Nv, M_x * M_y, (2 * Nu) * (2 * Nv));
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

}

#define BLOCKSIZE_SQUARED_DIFFERENCE	256

template <class T>
__global__ void evaluateSquaredDifference(const T * __restrict__ d_far_field_abs_projected, T * __restrict__ d_far_field_abs, const int N)  {

	const int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N) {
		
		d_far_field_abs[tid] = (d_far_field_abs[tid] * d_far_field_abs[tid] - d_far_field_abs_projected[tid] * d_far_field_abs_projected[tid]) * 
							   (d_far_field_abs[tid] * d_far_field_abs[tid] - d_far_field_abs_projected[tid] * d_far_field_abs_projected[tid]);

	}

}

template <class T> 
struct AbsFourth { __host__ __device__ T operator()(const T& x) const { return x * x * x * x; } };

float raFunctionalCalculation(const float2_ * __restrict__ d_far_field) {

	float *d_far_field_abs;	gpuErrchk(cudaMalloc(&d_far_field_abs, (2 * Nu) * (2 * Nv) * sizeof(float)));
	absKernel<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_ABSKERNEL), BLOCKSIZE_ABSKERNEL>>>(d_far_field_abs, d_far_field, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	
	thrust::device_ptr<float> d_far_field_abs_device_pointer = thrust::device_pointer_cast(d_far_field_abs);
	float scale_factor = thrust::reduce(d_far_field_abs_device_pointer, d_far_field_abs_device_pointer + (2 * Nu) * (2 * Nv), static_cast<float>(0), thrust::maximum<float>());
		
	vectorMulConstant(d_far_field_abs, static_cast<float>(1) / scale_factor, (2 * Nu) * (2 * Nv));

	float *d_far_field_abs_projected;	gpuErrchk(cudaMalloc(&d_far_field_abs_projected, (2 * Nu) * (2 * Nv) * sizeof(float)));
	evaluateProjection<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_PROJECTION), BLOCKSIZE_PROJECTION>>>(d_far_field_abs_projected, d_far_field_abs, 
		                                             d_Internal_Coverage_f, d_External_Coverage_f, (2 * Nu) * (2 * Nv));
	
	evaluateSquaredDifference<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_SQUARED_DIFFERENCE), BLOCKSIZE_SQUARED_DIFFERENCE>>>(d_far_field_abs_projected, d_far_field_abs, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	thrust::device_ptr<float> d_far_field_abs_projected_device_pointer = thrust::device_pointer_cast(d_far_field_abs_projected);
	float Numerator	= thrust::reduce(d_far_field_abs_device_pointer, d_far_field_abs_device_pointer + (2 * Nu) * (2 * Nv));
	float Denominator	= thrust::transform_reduce(d_far_field_abs_projected_device_pointer, d_far_field_abs_projected_device_pointer + (2 * Nu) * (2 * Nv),
                                  AbsFourth<float>(), static_cast<float>(0), thrust::plus<float>());
	
#ifdef DEBUG_PRINTING
	printf("Numerator = %f; Denominator = %f; Functional = %f\n", Numerator, Denominator, Numerator / Denominator);
#endif
	
	gpuErrchk(cudaFree(d_far_field_abs));
	gpuErrchk(cudaFree(d_far_field_abs_projected));
	
	return Numerator / Denominator;

}

double raFunctionalCalculation(const double2_ * __restrict__ d_far_field) {

	double *d_far_field_abs;	gpuErrchk(cudaMalloc(&d_far_field_abs, (2 * Nu) * (2 * Nv) * sizeof(double)));
	absKernel<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_ABSKERNEL), BLOCKSIZE_ABSKERNEL>>>(d_far_field_abs, d_far_field, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	
	thrust::device_ptr<double> d_far_field_abs_device_pointer = thrust::device_pointer_cast(d_far_field_abs);
	double scale_factor = thrust::reduce(d_far_field_abs_device_pointer, d_far_field_abs_device_pointer + (2 * Nu) * (2 * Nv), static_cast<double>(0), thrust::maximum<double>());
		
	vectorMulConstant(d_far_field_abs, static_cast<double>(1) / scale_factor, (2 * Nu) * (2 * Nv));

	double *d_far_field_abs_projected;	gpuErrchk(cudaMalloc(&d_far_field_abs_projected, (2 * Nu) * (2 * Nv) * sizeof(double)));
	evaluateProjection<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_PROJECTION), BLOCKSIZE_PROJECTION>>>(d_far_field_abs_projected, d_far_field_abs, 
		                                             d_Internal_Coverage, d_External_Coverage, (2 * Nu) * (2 * Nv));
	
	evaluateSquaredDifference<<<iDivUp((2 * Nu) * (2 * Nv), BLOCKSIZE_SQUARED_DIFFERENCE), BLOCKSIZE_SQUARED_DIFFERENCE>>>(d_far_field_abs_projected, d_far_field_abs, (2 * Nu) * (2 * Nv));
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	thrust::device_ptr<double> d_far_field_abs_projected_device_pointer = thrust::device_pointer_cast(d_far_field_abs_projected);
	double Numerator	= thrust::reduce(d_far_field_abs_device_pointer, d_far_field_abs_device_pointer + (2 * Nu) * (2 * Nv));
	double Denominator	= thrust::transform_reduce(d_far_field_abs_projected_device_pointer, d_far_field_abs_projected_device_pointer + (2 * Nu) * (2 * Nv),
                                  AbsFourth<double>(), static_cast<double>(0), thrust::plus<double>());
	
#ifdef DEBUG_PRINTING
	printf("Numerator = %f; Denominator = %f; Functional = %f\n", Numerator, Denominator, Numerator / Denominator);
#endif

	return Numerator / Denominator;

}

/********************************/
/* REFLECTARRAY COST FUNCTIONAL */
/********************************/
//float raCostFunctional(float *d_in) {
//
//	float2_ *d_far_field = raFarFieldCalculation(d_in);
//
//	return raFunctionalCalculation(d_far_field);
//
//}

//thrust::pair<float, float> raCostFunctional(float *d_in) {
thrust::tuple<float, float, float> raCostFunctional(float *d_in) {

	//thrust::pair<float2_ *, float> d_raFFCalc = raFarFieldCalculation(d_in);
	//	
	//float2_ *d_far_field	= d_raFFCalc.first;
	//float    min_distance	= d_raFFCalc.second;

	thrust::tuple<float2_ *, float, float> d_raFFCalc = raFarFieldCalculation(d_in);
		
	//float2_ *d_far_field	= d_raFFCalc.first;
	//float    min_distance	= d_raFFCalc.second;
	float2_ *d_far_field	= d_raFFCalc.get<0>();
	float    min_distance	= d_raFFCalc.get<1>();
	float    max_distance	= d_raFFCalc.get<2>();

	float Functional = raFunctionalCalculation(d_far_field);

	gpuErrchk(cudaFree(d_far_field));
	
	//return thrust::make_pair(Functional, min_distance);
	return thrust::make_tuple(Functional, min_distance, max_distance);

}
