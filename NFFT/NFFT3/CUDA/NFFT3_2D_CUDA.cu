#include<cstdlib>
#include<stdio.h>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<time.h>
#include<algorithm>
#include<assert.h>
#include<fstream>
#include<string>
#include<cfloat>

#include<cufft.h>
#include<device_launch_parameters.h>

#include <thrust\device_vector.h>
#include <thrust\reduce.h>
#include <thrust\device_ptr.h>

#include "NUFFT3_2D_CUDA.cuh"
#include "TimingGPU.cuh"
#include "Utilities.cuh"

#define PI 3.141592653589793238463

#define BLOCKSIZE                   512
#define BLOCKSIZE_SPATIAL           32
#define BLOCKSIZE_SPECTRAL          32
#define BLOCKSIZE_X_FFTSHIFT        32
#define BLOCKSIZE_Y_FFTSHIFT        32
#define BLOCKSIZEX_DP				32
#define BLOCKSIZEY_DP				32

//#define DYNAMIC_PARALLELISM

/***************************/
/* CIRCHSHIFT AND FFTSHIFT */
/***************************/
__global__ void circAndFFTShiftsKernel(double2 * __restrict__ out, const double2 * __restrict__ in, const int xdim, const int ydim, const int xshift, const int yshift)
{

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if ((i < xdim) && (j < ydim)) {

		int ii = (i + xshift) % xdim;
		int jj = (j + yshift) % ydim;

		out[jj * xdim + ii].x = (1. - 2 * ((i + j) & 1)) * in[j * xdim + i].x;
		out[jj * xdim + ii].y = (1. - 2 * ((i + j) & 1)) * in[j * xdim + i].y;

	}

}

/*****************************************************/
/* SPATIAL CONVOLUTION KERNEL WITH ATOMIC OPERATIONS */
/*****************************************************/
__global__ void summationDynamicParallelism(double2 * __restrict__ d_f_tau, const double fact, 
											const double real_part_temp, const double imag_part_temp,  
											const double x_temp, const double y_temp, const double t1, 
											const double Dx, const double Dy, const double b, const int mm_offset,
											const int nn_offset, const int msp, const int Mrx, const int Mry) {

	const int tidx = threadIdx.x + blockIdx.x * blockDim.x;
	const int tidy = threadIdx.y + blockIdx.y * blockDim.y;

	if ((tidx >= (2 * msp + 1)) || (tidy >= (2 * msp + 1))) return;
	
	const int mm = tidx - msp;
	const int nn = tidy - msp;

	int j = (mm + mm_offset);

	double d_n_sj = j - (int)Mry / 2;

	int i = nn + nn_offset;

	int tid = j * Mrx + i;

	double d_n_si = i - (int)Mrx / 2;

	double exp_temp = fact * exp(b * (Dy * Dy * d_n_sj * d_n_sj + Dx * Dx * d_n_si * d_n_si) - t1 * (((x_temp - d_n_si) * (x_temp - d_n_si) + (y_temp - d_n_sj) * (y_temp - d_n_sj))));

	atomicAdd(&d_f_tau[tid].x, real_part_temp * exp_temp);
	atomicAdd(&d_f_tau[tid].y, imag_part_temp * exp_temp);

}

__global__ void spatialConvolutionKernelAtomic(const double2 * __restrict__ d_f, const double * __restrict__ d_x, const double * __restrict__ d_y, double2 * __restrict__ d_f_tau, const double t1, const int len_in,
	const int msp, const double xb, const double yb, const double sb, const double tb, const double Dx, const double Dy, const double b, const int Mrx, const int Mry) {

	int kk = threadIdx.x + blockIdx.x * blockDim.x;

	if (kk < len_in) {

		double fact = 1. / (4.* PI *b);
		fact = fact * fact;

		double arg = sb * d_x[kk] + tb * d_y[kk];

		double real_part_temp = d_f[kk].x * cos(arg) + d_f[kk].y * sin(arg);
		double imag_part_temp = d_f[kk].y * cos(arg) - d_f[kk].x * sin(arg);

		double x_temp = (d_x[kk] - xb) / Dx;
		double y_temp = (d_y[kk] - yb) / Dy;

		int nn_offset = floor(Mrx / 2 + x_temp);
		int mm_offset = floor(Mry / 2 + y_temp);

#ifdef DYNAMIC_PARALLELISM
		dim3 dimBlock(BLOCKSIZEX_DP, BLOCKSIZEY_DP); dim3 dimGrid(iDivUp(2 * msp + 1, BLOCKSIZEX_DP), iDivUp(2 * msp + 1, BLOCKSIZEY_DP));

		summationDynamicParallelism <<< dimGrid, dimBlock >>> (d_f_tau, fact, real_part_temp, imag_part_temp,
			x_temp, y_temp, t1, Dx, Dy, b, mm_offset, nn_offset,
			msp, Mrx, Mry);
#else
		for (int mm = -msp; mm <= msp; mm++) {

			int j = (mm + mm_offset);

			double d_n_sj = j - (int)Mry / 2;

			for (int nn = -msp; nn < msp; nn++) {

				int i = nn + nn_offset;

				int tid = j * Mrx + i;

				double d_n_si = i - (int)Mrx / 2;

				double exp_temp = fact * exp(b * (Dy * Dy * d_n_sj * d_n_sj + Dx * Dx * d_n_si * d_n_si) - t1 * (((x_temp - d_n_si) * (x_temp - d_n_si) + (y_temp - d_n_sj) * (y_temp - d_n_sj))));

				atomicAdd(&d_f_tau[tid].x, real_part_temp * exp_temp);
				atomicAdd(&d_f_tau[tid].y, imag_part_temp * exp_temp);
			}
		}
#endif

	}
}

/*****************************************************/
/* SPECTRAL CONVOLUTION KERNEL AND POST-COMPENSATION */
/*****************************************************/
__global__ void spectralConvolutionKernelCompensation(double2 * __restrict__ d_F_stm, const double2 * __restrict__ d_F_st, const double * __restrict d_s, const double * __restrict d_t, const double sb, const double tb,
													  const int msp, const double t1, const int len_out, const double Ds, const double Dt, const double Dx, const double Dy, const double xb, const double yb, const double b,
													  const int Mrx, const int Mry) {

	int kk = threadIdx.x + blockIdx.x * blockDim.x;

	if (kk < len_out) {

		int d_nsn = (int)floor(Mrx / 2 + (d_s[kk] - sb) / Ds);
		double d_difs = (Mrx / 2 + (d_s[kk] - sb) / Ds) - d_nsn;
		int d_ntn = (int)floor(Mry / 2 + (d_t[kk] - tb) / Dt);
		double d_dift = (Mry / 2 + (d_t[kk] - tb) / Dt) - d_ntn;

		double temp1 = Dy * (d_t[kk] - tb); temp1 = temp1 * temp1;
		double temp2 = Dx * (d_s[kk] - sb); temp2 = temp2 * temp2;
		double temp = exp(b * (temp1 + temp2));

		double2 d_temp4; d_temp4.x = 0.; d_temp4.y = 0.;

		for (int nn = d_nsn - msp; nn <= d_nsn + msp; nn++){
			for (int mm = d_ntn - msp; mm <= d_ntn + msp; mm++){

				int i = mm - d_ntn + msp;
				int j = nn - d_nsn + msp;
				double d_n_pi = i - msp;
				double d_n_pj = j - msp;

				double d_temp3 = exp(-t1 * ((d_dift - d_n_pi) * (d_dift - d_n_pi) + (d_difs - d_n_pj)*(d_difs - d_n_pj)));

				d_temp4.x = d_temp4.x + d_F_st[mm * Mrx + nn].x * d_temp3;
				d_temp4.y = d_temp4.y + d_F_st[mm * Mrx + nn].y * d_temp3;

			}
		}

		double arg = (d_s[kk] - sb) * xb + (d_t[kk] - tb) * yb;
		double cosarg = cos(arg);
		double sinarg = sin(arg);

		d_F_stm[kk].x = (d_temp4.x * cosarg + d_temp4.y * sinarg) * temp;
		d_F_stm[kk].y = (d_temp4.y * cosarg - d_temp4.x * sinarg) * temp;

	}
}

/******************/
/* SCALING KERNEL */
/******************/
__global__ void scalingKernel(double * __restrict__ d_inout1, double * __restrict__ d_inout2, const double a, const double b, const int N) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid >= N) return;

	d_inout1[tid] = a * d_inout1[tid];
	d_inout2[tid] = b * d_inout2[tid];

}

/**********************/
/* ROUTINE NUFFT-3 2D */
/**********************/
void  NFFT2_2D_GPU(double * __restrict__ d_x, double * __restrict__ d_y, double * __restrict__ d_s, double * __restrict__ d_t, double2 *d_f, double2 * __restrict__ d_F_stm, const double eps, const int len_in, const int len_out) {

	//TimingGPU timerGPU;
	//std::ofstream timingFile;
	//timingFile.open("timingFile.txt");

	/************************************/
	/* SETTING THE ALGORITHM PARAMETERS */
	/************************************/
	thrust::device_ptr<double> dev_ptr_x = thrust::device_pointer_cast(d_x);
	thrust::device_ptr<double> dev_ptr_y = thrust::device_pointer_cast(d_y);
	thrust::device_ptr<double> dev_ptr_s = thrust::device_pointer_cast(d_s);
	thrust::device_ptr<double> dev_ptr_t = thrust::device_pointer_cast(d_t);

	double Max_x = thrust::reduce(dev_ptr_x, dev_ptr_x + len_in, -FLT_MAX, thrust::maximum<double>());
	double min_x = thrust::reduce(dev_ptr_x, dev_ptr_x + len_in,  FLT_MAX, thrust::minimum<double>());
	double xb = (Max_x + min_x) / 2;
	double X1 = fabs(Max_x - xb);

	double Max_y = thrust::reduce(dev_ptr_y, dev_ptr_y + len_in, -FLT_MAX, thrust::maximum<double>());
	double min_y = thrust::reduce(dev_ptr_y, dev_ptr_y + len_in,  FLT_MAX, thrust::minimum<double>());
	double yb = (Max_y + min_y) / 2;
	double Y1 = fabs(Max_y - yb);

	double Max_s = thrust::reduce(dev_ptr_s, dev_ptr_s + len_out, -FLT_MAX, thrust::maximum<double>());
	double min_s = thrust::reduce(dev_ptr_s, dev_ptr_s + len_out,  FLT_MAX, thrust::minimum<double>());
	double sb = (Max_s + min_s) / 2;
	double S = fabs(Max_s - sb);

	double Max_t = thrust::reduce(dev_ptr_t, dev_ptr_t + len_out, -FLT_MAX, thrust::maximum<double>());
	double min_t = thrust::reduce(dev_ptr_t, dev_ptr_t + len_out,  FLT_MAX, thrust::minimum<double>());
	double tb = (Max_t + min_t) / 2;
	double T = fabs(Max_t - tb);

	double R = 2.2;
	double msp = 2 + 2 * (R * R) * (-log(eps / 76)) / (PI * (R * R - 2)); // --- 2*msp are the numbers of non negligible samples of the Gaussian bell function

	int Mrx = 2 * ceil((X1 * S / PI) * R * R + R * msp);				  // --- Number of points who satisfies condition Ca and Cb
	int Mry = 2 * ceil((Y1 * T / PI) * R * R + R * msp);

	double Dx = 2 * PI / (Mrx); double Ds = 1;
	double Dy = 2 * PI / (Mry); double Dt = 1;

	double kk_x = X1 / (PI / 2);        double kk_y = Y1 / (PI / 2);

	scalingKernel << <iDivUp(len_in, BLOCKSIZE), BLOCKSIZE >> >(d_x, d_y, 1. / kk_x, 1. / kk_y, len_in);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	scalingKernel << <iDivUp(len_out, BLOCKSIZE), BLOCKSIZE >> >(d_s, d_t, kk_x, kk_y, len_out);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	Max_x = Max_x / kk_x;
	min_x = min_x / kk_x;
	xb = (Max_x + min_x) / 2;

	Max_s = Max_s * kk_x;
	min_s = min_s * kk_x;
	sb = (Max_s + min_s) / 2;

	Max_y = Max_y / kk_y;
	min_y = min_y / kk_y;
	yb = (Max_y + min_y) / 2;

	Max_t = Max_t * kk_y;
	min_t = min_t * kk_y;
	tb = (Max_t + min_t) / 2;

	// --- Precompute Convolution Constants
	double b = (msp - 2) / (4 * PI); double t1 = 1 / (4 * b); msp = floor(msp - 2);

	/*******************************************************************/
	/* STEP 1: PRE-SCALING AND GAUSSIAN GRIDDING IN THE SPATIAL DOMAIN */
	/*******************************************************************/

	double2 *d_f_tau;   gpuErrchk(cudaMalloc((void**)&d_f_tau, Mrx * Mry * sizeof(double2)));

	//timerGPU.StartCounter();
	
#ifdef DYNAMIC_PARALLELISM
	gpuErrchk(cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, BLOCKSIZE_SPATIAL * (2 * msp + 1) * (2 * msp + 1)));
#endif

	gpuErrchk(cudaMemset(d_f_tau, 0, Mrx * Mry * sizeof(double2)));
	spatialConvolutionKernelAtomic << <iDivUp(len_in, BLOCKSIZE_SPATIAL), BLOCKSIZE_SPATIAL >> >(d_f, d_x, d_y, d_f_tau, t1, len_in, msp, xb, yb, sb, tb, Dx, Dy, b, Mrx, Mry);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	//timingFile << "spatialConvolutionKernelAtomic " << timerGPU.GetCounter() << "\n";

	/***********************/
	/* STEP 2: STANDARD FFT*/
	/***********************/
	//timerGPU.StartCounter();
	cufftHandle plan;
	cufftPlan2d(&plan, Mry, Mrx, CUFFT_Z2Z);
	//cufftExecZ2Z(plan, d_f_tau, d_F_stx, CUFFT_FORWARD);
	cufftExecZ2Z(plan, d_f_tau, d_f_tau, CUFFT_FORWARD);
	cufftDestroy(plan);
	//timingFile << "FFT " << timerGPU.GetCounter() << "\n";

	//timerGPU.StartCounter();
	double2 *d_F_stx2; gpuErrchk(cudaMalloc((void**)&d_F_stx2, Mrx * Mry * sizeof(cuDoubleComplex)));

	dim3 GridSize_FFTSHIFT(iDivUp(Mrx, BLOCKSIZE_X_FFTSHIFT), iDivUp(Mry, BLOCKSIZE_Y_FFTSHIFT));
	dim3 BlockSize_FFTSHIFT(BLOCKSIZE_X_FFTSHIFT, BLOCKSIZE_Y_FFTSHIFT);
	circAndFFTShiftsKernel << <GridSize_FFTSHIFT, BlockSize_FFTSHIFT >> >(d_F_stx2, d_f_tau, Mrx, Mry, Mrx / 2, Mry / 2);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	//timingFile << "FFT shift " << timerGPU.GetCounter() << "\n";

	/*************************************************************************/
	/*STEP 3: GAUSSIAN GRIDDING IN THE SPECTRAL DOMAIN AND POST-COMPENSATION */
	/*************************************************************************/
	//timerGPU.StartCounter();
	spectralConvolutionKernelCompensation << <iDivUp(len_out, BLOCKSIZE_SPECTRAL), BLOCKSIZE_SPECTRAL >> >(d_F_stm, d_F_stx2, d_s, d_t, sb, tb, msp, t1, len_out, Ds, Dt, Dx, Dy, xb, yb, b, Mrx, Mry);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif
	//timingFile << "spectralConvolutionKernelCompensation " << timerGPU.GetCounter() << "\n";

	/***********************************/
	/* FREEING THE ALLOCATED VARIABLES */
	/***********************************/
	gpuErrchk(cudaFree(d_f_tau));
	gpuErrchk(cudaFree(d_F_stx2));

}
