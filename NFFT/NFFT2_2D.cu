#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

#include "Utilities.cuh"
#include "Bessel.cuh"
#include "cuFFT_auxiliary.cuh"
#include "NFFT2_2D.cuh"
#include "InputOutput.cuh"

cufftHandle	NFFT2_2D_GPUplan;

#define BLOCKSIZE_INTERPOLATION	 256
#define BLOCK_SIZE_x			  16
#define BLOCK_SIZE_y			  16
#define BLOCKSIZE_BESSEL		  32

#define DEBUG

#define cc 2
#define K 6

#define IDX2R(i,j,N) (((i)*(N))+(j))

#define pi_double	3.141592653589793238463

__constant__ double alpha=(2.-1./cc)*pi_double-0.01;

__constant__ int constant1_GPU;
__constant__ int constant2_GPU;
__constant__ int constant3_GPU;
__constant__ int constant4_GPU;
__constant__ int constant5_GPU;
__constant__ int constant6_GPU;
__constant__ int constant7_GPU;
__constant__ int constant8_GPU;
__constant__ int constant9_GPU;
__constant__ int constant10_GPU;
__constant__ float constant11_GPU_f;
__constant__ float constant12_GPU_f;
__constant__ double constant11_GPU;
__constant__ double constant12_GPU;
__constant__ int constant13_GPU;
__constant__ int constant14_GPU;

/**************************/
/* cuFFT PLAN CALCULATION */
/**************************/
void Calculate_cuFFT_plan_C2C_NFFT2_2D(const int N1, const int N2) { cufftSafeCall(cufftPlan2d(&NFFT2_2D_GPUplan, cc*N1, cc*N2, CUFFT_C2C)); }
void Calculate_cuFFT_plan_Z2Z_NFFT2_2D(const int N1, const int N2) { cufftSafeCall(cufftPlan2d(&NFFT2_2D_GPUplan, cc*N1, cc*N2, CUFFT_Z2Z)); }

/*************************/
/* FFTW PLAN CALCULATION */
/*************************/
void Destroy_cuFFT_plan_NFFT2_2D() { cufftSafeCall(cufftDestroy(NFFT2_2D_GPUplan)); }

/***********************************************/
/* MODIFIED BESSEL FUNCTION CALCULATION KERNEL */
/***********************************************/
template<class T>
__global__ void Kernel_Bessel(T * __restrict__ Bessel_vector, const int N)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	if (i<N) {
        T xi = (static_cast<T>(2*pi_double)*(i-(N/2)))/(cc*N);
		Bessel_vector[i] = static_cast<T>(1)/(bessi0(static_cast<T>(K)*sqrt(static_cast<T>(alpha*alpha)-xi*xi)));
	}
}

/**************************/
/* DECIMATION AND SCALING */
/**************************/
__global__ void Decimation_and_Scaling(const float2* __restrict__ data, float2* __restrict__ result, const float* __restrict__ Bessel_vector_x, const float* __restrict__ Bessel_vector_y, const int N1, const int N2)
{
    int i = threadIdx.y + blockDim.y * blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;

    if((((i >= constant1_GPU)  && (i < constant2_GPU)) && ((j >= constant3_GPU)  && (j < constant4_GPU))))
	{
		float a = Bessel_vector_x[i-constant1_GPU]*Bessel_vector_y[j-constant3_GPU];

		result[IDX2R(i-constant1_GPU,j-constant3_GPU,N2)].x=data[IDX2R(i,j,cc*N2)].x*a;
		result[IDX2R(i-constant1_GPU,j-constant3_GPU,N2)].y=data[IDX2R(i,j,cc*N2)].y*a;
	}
}

__global__ void Decimation_and_Scaling(const double2* __restrict__ data, double2* __restrict__ result, const double* __restrict__ Bessel_vector_x, const double* __restrict__ Bessel_vector_y, const int N1, const int N2)
{
    int i = threadIdx.y + blockDim.y * blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;

    if((((i >= constant1_GPU)  && (i < constant2_GPU)) && ((j >= constant3_GPU)  && (j < constant4_GPU))))
	{
		double a = Bessel_vector_x[i-constant1_GPU]*Bessel_vector_y[j-constant3_GPU];

		result[IDX2R(i-constant1_GPU,j-constant3_GPU,N2)].x=data[IDX2R(i,j,cc*N2)].x*a;
		result[IDX2R(i-constant1_GPU,j-constant3_GPU,N2)].y=data[IDX2R(i,j,cc*N2)].y*a;
	}
}

/*****************************************************************************************/
/* KERNEL FUNCTION TO CALCULATE SERIES TERMS FOR INTERPOLATION USING DYNAMIC PARALLELISM */
/*****************************************************************************************/
__global__ void series_terms(float2 temp_data, float2* __restrict__ result, const float r_cc_points1, const float cc_diff1, const float r_cc_points2, const float cc_diff2, const int N1, const int N2)
{
	int m = threadIdx.x;
	int n = threadIdx.y;

    float tempd, phi_cap;

	float P = K*K-(cc_diff1-(m-K))*(cc_diff1-(m-K));

	if(P<0.) {tempd=rsqrt(-P); phi_cap = (static_cast<float>(1./pi_double))*((sin(alpha/tempd))*tempd);  }
	else if(P>0.f) {tempd=rsqrt(P); phi_cap = static_cast<float>(1./pi_double)*((sinh(alpha/tempd))*tempd); }
	else phi_cap = static_cast<float>(alpha/pi_double);

	P = K*K-(cc_diff2-(n-K))*(cc_diff2-(n-K));

	if(P<0.) {tempd=rsqrt(-P); phi_cap = phi_cap*static_cast<float>(1./pi_double)*((sin(alpha/tempd))*tempd);  }
	else if(P>0.f) {tempd=rsqrt(P); phi_cap = phi_cap*static_cast<float>(1./pi_double)*((sinh(alpha/tempd))*tempd); }
	else phi_cap = static_cast<float>(phi_cap*alpha/pi_double);

	int PP1 = modulo((r_cc_points1+(m-K)+N1*cc/2),(cc*N1));
	int PP2 = modulo((r_cc_points2+(n-K)+N2*cc/2),(cc*N2));

	atomicAdd(&result[IDX2R(PP1,PP2,cc*N2)].x,temp_data.x*phi_cap);
	atomicAdd(&result[IDX2R(PP1,PP2,cc*N2)].y,temp_data.y*phi_cap);

}

__global__ void series_terms(double2 temp_data, double2* __restrict__ result, const double r_cc_points1, const double cc_diff1, const double r_cc_points2, const double cc_diff2, const int N1, const int N2)
{
	int m = threadIdx.x;
	int n = threadIdx.y;

    double tempd, phi_cap;

	double P = K*K-(cc_diff1-(m-K))*(cc_diff1-(m-K));

	if(P<0.) {tempd=rsqrt(-P); phi_cap = (1./pi_double)*((sin(alpha/tempd))*tempd);  }
	else if(P>0.) {tempd=rsqrt(P); phi_cap = (1./pi_double)*((sinh(alpha/tempd))*tempd); }
	else phi_cap = alpha/pi_double;

	P = K*K-(cc_diff2-(n-K))*(cc_diff2-(n-K));

	if(P<0.) {tempd=rsqrt(-P); phi_cap = phi_cap*(1./pi_double)*((sin(alpha/tempd))*tempd);  }
	else if(P>0.) {tempd=rsqrt(P); phi_cap = phi_cap*(1./pi_double)*((sinh(alpha/tempd))*tempd); }
	else phi_cap = phi_cap*alpha/pi_double;

	int PP1 = modulo((r_cc_points1+(m-K)+N1*cc/2),(cc*N1));
	int PP2 = modulo((r_cc_points2+(n-K)+N2*cc/2),(cc*N2));

	atomicAdd(&result[IDX2R(PP1,PP2,cc*N2)].x,temp_data.x*phi_cap);
	atomicAdd(&result[IDX2R(PP1,PP2,cc*N2)].y,temp_data.y*phi_cap);

}

/************************/
/* INTERPOLATION 2D NED */
/************************/
// --- Code using dynamic parallelism
__global__ void Interpolation_NFFT2_2D_GPUKernel(const float2* __restrict__ data, float2* __restrict__ result, const float* __restrict__ x, const float* __restrict__ y, const int N1, const int N2, int M)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	float cc_points1=cc*x[i];
	float r_cc_points1=rint(cc_points1);				// It is the mu in Fourmont's paper
	const float cc_diff1 = cc_points1-r_cc_points1;

	float cc_points2=cc*y[i];
	float r_cc_points2=rint(cc_points2);				// It is the mu in Fourmont's paper
	const float cc_diff2 = cc_points2-r_cc_points2;

#	if __CUDA_ARCH__ >= 350
	float2 temp_data = data[i];

	dim3 dimBlock(13,13); dim3 dimGrid(1,1);

	if(i<M) series_terms<<<dimGrid,dimBlock>>>(temp_data,result,r_cc_points1,cc_diff1,r_cc_points2,cc_diff2,N1,N2);

#   else
	int PP1, PP2;
	float P, tempd;

	float phi_cap1, phi_cap2;

	if(i<M) {

		for(int m=0; m<constant9_GPU; m++) {

			P = constant10_GPU-(cc_diff1-(m-K))*(cc_diff1-(m-K));

			PP1 = modulo((r_cc_points1+(m-K)+N1*cc/2),constant13_GPU);

			if(P<0.) {tempd=rsqrt(-P); phi_cap1 = constant11_GPU_f*((sin(alpha/tempd))*tempd);  }
			else if(P>0.f) {tempd=rsqrt(P); phi_cap1 = constant11_GPU_f*((sinh(alpha/tempd))*tempd); }
			else phi_cap1 = constant12_GPU_f;

			for(int n=0; n<constant9_GPU; n++) {

				P = constant10_GPU-(cc_diff2-(n-K))*(cc_diff2-(n-K));

				PP2 = modulo((r_cc_points2+(n-K)+N2*cc/2),constant14_GPU);

				if(P<0.f) {tempd=rsqrt(-P); phi_cap2 = phi_cap1*constant11_GPU_f*((sin(alpha/tempd))*tempd);  }
				else if(P>0.f) {tempd=rsqrt(P); phi_cap2 = phi_cap1*constant11_GPU_f*((sinh(alpha/tempd))*tempd); }
				else phi_cap2 = phi_cap1*constant12_GPU_f;

				atomicAdd(&result[IDX2R(PP1,PP2,constant14_GPU)].x,data[i].x*phi_cap2);
				atomicAdd(&result[IDX2R(PP1,PP2,constant14_GPU)].y,data[i].y*phi_cap2);

			}
		}
	}
#   endif

}

// --- Code using dynamic parallelism
__global__ void Interpolation_NFFT2_2D_GPUKernel(const double2* __restrict__ data, double2* __restrict__ result, const double* __restrict__ x, const double* __restrict__ y, const int N1, const int N2, int M)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;

	double cc_points1=cc*x[i];
	double r_cc_points1=rint(cc_points1);				// It is the mu in Fourmont's paper
	const double cc_diff1 = cc_points1-r_cc_points1;

	double cc_points2=cc*y[i];
	double r_cc_points2=rint(cc_points2);				// It is the mu in Fourmont's paper
	const double cc_diff2 = cc_points2-r_cc_points2;

#	if __CUDA_ARCH__ >= 350
	double2 temp_data = data[i];

	dim3 dimBlock(13,13); dim3 dimGrid(1,1);

	if(i<M) series_terms<<<dimGrid,dimBlock>>>(temp_data,result,r_cc_points1,cc_diff1,r_cc_points2,cc_diff2,N1,N2);

#   else
	int PP1, PP2;
	double P, tempd;

	double phi_cap1, phi_cap2;

	if(i<M) {

		for(int m=0; m<constant9_GPU; m++) {

			P = constant10_GPU-(cc_diff1-(m-K))*(cc_diff1-(m-K));

			PP1 = modulo((r_cc_points1+(m-K)+N1*cc/2),constant13_GPU);

			if(P<0.) {tempd=rsqrt(-P); phi_cap1 = constant11_GPU*((sin(alpha/tempd))*tempd);  }
			else if(P>0.) {tempd=rsqrt(P); phi_cap1 = constant11_GPU*((sinh(alpha/tempd))*tempd); }
			else phi_cap1 = constant12_GPU;

			for(int n=0; n<constant9_GPU; n++) {

				P = constant10_GPU-(cc_diff2-(n-K))*(cc_diff2-(n-K));

				PP2 = modulo((r_cc_points2+(n-K)+N2*cc/2),constant14_GPU);

				if(P<0.) {tempd=rsqrt(-P); phi_cap2 = phi_cap1*constant11_GPU*((sin(alpha/tempd))*tempd);  }
				else if(P>0.) {tempd=rsqrt(P); phi_cap2 = phi_cap1*constant11_GPU*((sinh(alpha/tempd))*tempd); }
				else phi_cap2 = phi_cap1*constant12_GPU;

				atomicAdd(&result[IDX2R(PP1,PP2,constant14_GPU)].x,data[i].x*phi_cap2);
				atomicAdd(&result[IDX2R(PP1,PP2,constant14_GPU)].y,data[i].y*phi_cap2);

			}
		}
	}
#   endif

}

// --- Code not using dynamic parallelism
//__global__ void Interpolation_NFFT2_2D_GPUKernel(const double2* __restrict__ data, double2* __restrict__ result, const double* __restrict__ x, const double* __restrict__ y, const int N1, const int N2, int M)
//{
//	int i = threadIdx.x + blockDim.x * blockIdx.x;
//
//	double cc_points1=cc*x[i];
//	double r_cc_points1=rint(cc_points1);
//	const double cc_diff1 = cc_points1-r_cc_points1;
//
//	double cc_points2=cc*y[i];
//	double r_cc_points2=rint(cc_points2);
//	const double cc_diff2 = cc_points2-r_cc_points2;
//
//	int PP1, PP2;
//	double P, tempd;
//
//	double phi_cap1, phi_cap2;
//
//	if(i<M) {
//
//		for(int m=0; m<constant9_GPU; m++) {
//
//			P = constant10_GPU-(cc_diff1-(m-K))*(cc_diff1-(m-K));
//
//			PP1 = modulo((r_cc_points1+(m-K)+N1*cc/2),constant13_GPU);
//
//			if(P<0.) {tempd=rsqrt(-P); phi_cap1 = constant11_GPU*((sin(alpha/tempd))*tempd);  }
//			else if(P>0.) {tempd=rsqrt(P); phi_cap1 = constant11_GPU*((sinh(alpha/tempd))*tempd); }
//			else phi_cap1 = constant12_GPU;
//
//			for(int n=0; n<constant9_GPU; n++) {
//
//				P = constant10_GPU-(cc_diff2-(n-K))*(cc_diff2-(n-K));
//
//				PP2 = modulo((r_cc_points2+(n-K)+N2*cc/2),constant14_GPU);
//
//				if(P<0.) {tempd=rsqrt(-P); phi_cap2 = phi_cap1*constant11_GPU*((sin(alpha/tempd))*tempd);  }
//				else if(P>0.) {tempd=rsqrt(P); phi_cap2 = phi_cap1*constant11_GPU*((sinh(alpha/tempd))*tempd); }
//				else phi_cap2 = phi_cap1*constant12_GPU;
//
//				atomicAdd(&result[IDX2R(PP1,PP2,constant14_GPU)].x,data[i].x*phi_cap2);
//				atomicAdd(&result[IDX2R(PP1,PP2,constant14_GPU)].y,data[i].y*phi_cap2);
//
//			}
//		}
//	}
//}

/***************************/
/* NUFFT NED 2D EVALUATION */
/***************************/
void NFFT2_2D_GPU(float2 * __restrict__ result, const float2 * __restrict__ data, const float * __restrict__ x, const float * __restrict__ y, const int N1, const int N2, const int M)
{

	float alfa_CPU=static_cast<float>((2.-1./cc)*pi_double-0.01);
	int constant1_CPU = (cc-1)*N1/2;								gpuErrchk(cudaMemcpyToSymbol(constant1_GPU,  &constant1_CPU,  sizeof(int)));
	int constant2_CPU = (cc+1)*N1/2;								gpuErrchk(cudaMemcpyToSymbol(constant2_GPU,  &constant2_CPU,  sizeof(int)));
	int constant3_CPU = (cc-1)*N2/2;								gpuErrchk(cudaMemcpyToSymbol(constant3_GPU,  &constant3_CPU,  sizeof(int)));
	int constant4_CPU = (cc+1)*N2/2;								gpuErrchk(cudaMemcpyToSymbol(constant4_GPU,  &constant4_CPU,  sizeof(int)));
	int constant5_CPU = (cc-1)*N1/2-N1/2;							gpuErrchk(cudaMemcpyToSymbol(constant5_GPU,  &constant5_CPU,  sizeof(int)));
	int constant6_CPU = (cc-1)*N2/2-N2/2;							gpuErrchk(cudaMemcpyToSymbol(constant6_GPU,  &constant6_CPU,  sizeof(int)));
	int constant7_CPU = 2.*pi_double/(cc*N1);						gpuErrchk(cudaMemcpyToSymbol(constant7_GPU,  &constant7_CPU,  sizeof(int)));
	int constant8_CPU = 2.*pi_double/(cc*N2);						gpuErrchk(cudaMemcpyToSymbol(constant8_GPU,  &constant8_CPU,  sizeof(int)));
	int constant9_CPU = 2*K+1;										gpuErrchk(cudaMemcpyToSymbol(constant9_GPU,  &constant9_CPU,  sizeof(int)));
	int constant10_CPU = K*K;										gpuErrchk(cudaMemcpyToSymbol(constant10_GPU, &constant10_CPU, sizeof(int)));
	float constant11_CPU = static_cast<float>(1./pi_double);		gpuErrchk(cudaMemcpyToSymbol(constant11_GPU_f, &constant11_CPU, sizeof(float)));
	float constant12_CPU = static_cast<float>(alfa_CPU/pi_double);	gpuErrchk(cudaMemcpyToSymbol(constant12_GPU_f, &constant12_CPU, sizeof(float)));
	int constant13_CPU = cc*N1;										gpuErrchk(cudaMemcpyToSymbol(constant13_GPU, &constant13_CPU, sizeof(int)));
	int constant14_CPU = cc*N2;										gpuErrchk(cudaMemcpyToSymbol(constant14_GPU, &constant14_CPU, sizeof(int)));

	/* CALCULATION OF BESSEL FUNCTIONS */
	float* Bessel_vector_x; gpuErrchk(cudaMalloc((void **)&Bessel_vector_x,sizeof(float)*N1));
	float* Bessel_vector_y; gpuErrchk(cudaMalloc((void **)&Bessel_vector_y,sizeof(float)*N2));

	dim3 dimBlock01(BLOCKSIZE_BESSEL,1); dim3 dimGrid01(iDivUp(N1,BLOCKSIZE_BESSEL));
	Kernel_Bessel<<<dimGrid01,dimBlock01>>>(Bessel_vector_x, N1);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	dim3 dimBlock02(BLOCKSIZE_BESSEL,1); dim3 dimGrid02(iDivUp(N2,BLOCKSIZE_BESSEL));
	Kernel_Bessel<<<dimGrid02,dimBlock02>>>(Bessel_vector_y, N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* ALLOCATIONS AND INITIALIZATIONS */
	cufftComplex *temp_result; gpuErrchk(cudaMalloc((void **)&temp_result,sizeof(cufftComplex)*cc*N1*cc*N2));
	gpuErrchk(cudaMemset(temp_result,0,sizeof(cufftComplex)*cc*N1*cc*N2));

	/* INTERPOLATION */
	dim3 dimBlock1(BLOCKSIZE_INTERPOLATION,1); dim3 dimGrid1(iDivUp(M,BLOCKSIZE_INTERPOLATION));
#	if __CUDA_ARCH__ >= 350
	gpuErrchk(cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, dimGrid1.x*BLOCKSIZE_INTERPOLATION));
#   endif
	Interpolation_NFFT2_2D_GPUKernel<<<dimGrid1,dimBlock1>>>(data,temp_result,x,y,N1,N2,M);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* FFTSHIFT 2D */
    dim3 dimBlock2(BLOCK_SIZE_x,BLOCK_SIZE_y);
	dim3 dimGrid2((cc*N1)/BLOCK_SIZE_x + ((cc*N1)%BLOCK_SIZE_x == 0 ? 0:1),(cc*N2)/BLOCK_SIZE_y + ((cc*N2)%BLOCK_SIZE_y == 0 ? 0:1));
	fftshift_2D<<<dimGrid2,dimBlock2>>>(temp_result,cc*N1,cc*N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* FFT */
    cufftSafeCall(cufftExecC2C(NFFT2_2D_GPUplan, temp_result, temp_result, CUFFT_FORWARD));

	/* FFTSHIFT 2D */
	fftshift_2D<<<dimGrid2,dimBlock2>>>(temp_result,cc*N1,cc*N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* DECIMATION AND SCALING */
    dim3 dimBlock3(BLOCK_SIZE_x,BLOCK_SIZE_y);
	dim3 dimGrid3((cc*N2)/BLOCK_SIZE_x + ((cc*N2)%BLOCK_SIZE_x == 0 ? 0:1),(cc*N1)/BLOCK_SIZE_y + ((cc*N1)%BLOCK_SIZE_y == 0 ? 0:1));
	Decimation_and_Scaling<<<dimGrid3,dimBlock3>>>(temp_result,result,Bessel_vector_x,Bessel_vector_y,N1,N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	gpuErrchk(cudaFree(Bessel_vector_x));
	gpuErrchk(cudaFree(Bessel_vector_y));
	gpuErrchk(cudaFree(temp_result));
}

void NFFT2_2D_GPU(double2 * __restrict__ result, const double2 * __restrict__ data, const double * __restrict__ x, const double * __restrict__ y, const int N1, const int N2, const int M)
{

	double alfa_CPU=(2.-1./cc)*pi_double-0.01;
	int constant1_CPU = (cc-1)*N1/2;			gpuErrchk(cudaMemcpyToSymbol(constant1_GPU,  &constant1_CPU,  sizeof(int)));
	int constant2_CPU = (cc+1)*N1/2;			gpuErrchk(cudaMemcpyToSymbol(constant2_GPU,  &constant2_CPU,  sizeof(int)));
	int constant3_CPU = (cc-1)*N2/2;			gpuErrchk(cudaMemcpyToSymbol(constant3_GPU,  &constant3_CPU,  sizeof(int)));
	int constant4_CPU = (cc+1)*N2/2;			gpuErrchk(cudaMemcpyToSymbol(constant4_GPU,  &constant4_CPU,  sizeof(int)));
	int constant5_CPU = (cc-1)*N1/2-N1/2;		gpuErrchk(cudaMemcpyToSymbol(constant5_GPU,  &constant5_CPU,  sizeof(int)));
	int constant6_CPU = (cc-1)*N2/2-N2/2;		gpuErrchk(cudaMemcpyToSymbol(constant6_GPU,  &constant6_CPU,  sizeof(int)));
	int constant7_CPU = 2.*pi_double/(cc*N1);	gpuErrchk(cudaMemcpyToSymbol(constant7_GPU,  &constant7_CPU,  sizeof(int)));
	int constant8_CPU = 2.*pi_double/(cc*N2);	gpuErrchk(cudaMemcpyToSymbol(constant8_GPU,  &constant8_CPU,  sizeof(int)));
	int constant9_CPU = 2*K+1;					gpuErrchk(cudaMemcpyToSymbol(constant9_GPU,  &constant9_CPU,  sizeof(int)));
	int constant10_CPU = K*K;					gpuErrchk(cudaMemcpyToSymbol(constant10_GPU, &constant10_CPU, sizeof(int)));
	double constant11_CPU = 1./pi_double;		gpuErrchk(cudaMemcpyToSymbol(constant11_GPU, &constant11_CPU, sizeof(double)));
	double constant12_CPU = alfa_CPU/pi_double;	gpuErrchk(cudaMemcpyToSymbol(constant12_GPU, &constant12_CPU, sizeof(double)));
	int constant13_CPU = cc*N1;					gpuErrchk(cudaMemcpyToSymbol(constant13_GPU, &constant13_CPU, sizeof(int)));
	int constant14_CPU = cc*N2;					gpuErrchk(cudaMemcpyToSymbol(constant14_GPU, &constant14_CPU, sizeof(int)));

	/* CALCULATION OF BESSEL FUNCTIONS */
	double* Bessel_vector_x; gpuErrchk(cudaMalloc((void **)&Bessel_vector_x,sizeof(double)*N1));
	double* Bessel_vector_y; gpuErrchk(cudaMalloc((void **)&Bessel_vector_y,sizeof(double)*N2));

	dim3 dimBlock01(BLOCKSIZE_BESSEL,1); dim3 dimGrid01(iDivUp(N1,BLOCKSIZE_BESSEL));
	Kernel_Bessel<<<dimGrid01,dimBlock01>>>(Bessel_vector_x, N1);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	dim3 dimBlock02(BLOCKSIZE_BESSEL,1); dim3 dimGrid02(iDivUp(N2,BLOCKSIZE_BESSEL));
	Kernel_Bessel<<<dimGrid02,dimBlock02>>>(Bessel_vector_y, N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* ALLOCATIONS AND INITIALIZATIONS */
	cufftDoubleComplex *temp_result; gpuErrchk(cudaMalloc((void **)&temp_result,sizeof(cufftDoubleComplex)*cc*N1*cc*N2));
	gpuErrchk(cudaMemset(temp_result,0,sizeof(cufftDoubleComplex)*cc*N1*cc*N2));

	/* INTERPOLATION */
	dim3 dimBlock1(BLOCKSIZE_INTERPOLATION,1); dim3 dimGrid1(iDivUp(M,BLOCKSIZE_INTERPOLATION));
#	if __CUDA_ARCH__ >= 350
	gpuErrchk(cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, dimGrid1.x*BLOCKSIZE_INTERPOLATION));
#   endif
	Interpolation_NFFT2_2D_GPUKernel<<<dimGrid1,dimBlock1>>>(data,temp_result,x,y,N1,N2,M);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* FFTSHIFT 2D */
    dim3 dimBlock2(BLOCK_SIZE_x,BLOCK_SIZE_y);
	dim3 dimGrid2((cc*N1)/BLOCK_SIZE_x + ((cc*N1)%BLOCK_SIZE_x == 0 ? 0:1),(cc*N2)/BLOCK_SIZE_y + ((cc*N2)%BLOCK_SIZE_y == 0 ? 0:1));
	fftshift_2D<<<dimGrid2,dimBlock2>>>(temp_result,cc*N1,cc*N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* FFT */
    cufftSafeCall(cufftExecZ2Z(NFFT2_2D_GPUplan, temp_result, temp_result, CUFFT_FORWARD));

	/* FFTSHIFT 2D */
	fftshift_2D<<<dimGrid2,dimBlock2>>>(temp_result,cc*N1,cc*N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	/* DECIMATION AND SCALING */
    dim3 dimBlock3(BLOCK_SIZE_x,BLOCK_SIZE_y);
	dim3 dimGrid3((cc*N2)/BLOCK_SIZE_x + ((cc*N2)%BLOCK_SIZE_x == 0 ? 0:1),(cc*N1)/BLOCK_SIZE_y + ((cc*N1)%BLOCK_SIZE_y == 0 ? 0:1));
	Decimation_and_Scaling<<<dimGrid3,dimBlock3>>>(temp_result,result,Bessel_vector_x,Bessel_vector_y,N1,N2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	gpuErrchk(cudaFree(Bessel_vector_x));
	gpuErrchk(cudaFree(Bessel_vector_y));
	gpuErrchk(cudaFree(temp_result));
}
