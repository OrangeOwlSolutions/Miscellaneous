#include "stdafx.h"
#pragma warning(disable:4996)

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<complex>
#include<time.h>
#include<algorithm>
#include<cfloat>
#include<vector>
#include<fstream>

#include<omp.h>

#include "mkl_dfti.h"

#include "NUFFT3_2D_C.h"
#include "TimingCPU.h"

#define PI 3.141592653589793238463

#include<fftw3.h>
#include<fstream>
#include<string>

#include<thrust/host_vector.h>
#include<thrust/reduce.h>
#include<thrust/device_ptr.h>
#include <thrust/iterator/constant_iterator.h>

/********************************/
/* SPATIAL CONVOLUTION FUNCTION */
/********************************/
void spatialConvolutionFunction(const fftw_complex * __restrict f, const double * __restrict x, const double * __restrict y, fftw_complex * __restrict f_tau, const double t1, const int len_in, const int msp,
	                            const double xb, const double yb, const double sb, const double tb, const double Dx, const double Dy, const double b, const int Mrx, const int Mry) {

	omp_set_num_threads(omp_get_max_threads());
	#pragma omp parallel for
	for (int kk = 0; kk < len_in; kk++) {

		double fact = 1. / (4.* PI *b);
		fact = fact * fact;

		double arg = sb * x[kk] + tb * y[kk];

		double real_part_temp = f[kk][0] * cos(arg) + f[kk][1] * sin(arg);
		double imag_part_temp = f[kk][1] * cos(arg) - f[kk][0] * sin(arg);

		double x_temp = (x[kk] - xb) / Dx;
		double y_temp = (y[kk] - yb) / Dy;

		int nn_offset = floor(Mrx / 2 + x_temp);
		int mm_offset = floor(Mry / 2 + y_temp);

		double difx = (Mrx / 2 + (x[kk] - xb) / Dx) - nn_offset;
		double dify = (Mry / 2 + (y[kk] - yb) / Dy) - mm_offset;

		for (int mm = -msp; mm <= msp; mm++) {

			int j = (mm + mm_offset);

			double d_n_sj = j - (int)Mry / 2;

			for (int nn = -msp; nn <= msp; nn++) {

				int i = nn + nn_offset;

				double d_n_si = i - (int)Mrx / 2;

				double exp_temp = fact * exp(b * (Dy * Dy * d_n_sj * d_n_sj + Dx * Dx * d_n_si * d_n_si) - t1 * ((x_temp - d_n_si) * (x_temp - d_n_si) + (y_temp - d_n_sj) * (y_temp - d_n_sj)));

				int index = j * Mrx + i;

				#pragma omp atomic
				f_tau[index][0] += real_part_temp * exp_temp;
				#pragma omp atomic
				f_tau[index][1] += imag_part_temp * exp_temp;
			}
		}
	}
}

//void spatialConvolutionFunction(const fftw_complex * __restrict f, const double * __restrict x, const double * __restrict y, fftw_complex * __restrict f_tau, const double t1, const int len_in, const int msp,
//	const double xb, const double yb, const double sb, const double tb, const double Dx, const double Dy, const double b, const int Mrx, const int Mry) {
//
//	for (int j = 0; j < Mry; j++)
//		for (int i = 0; i < Mrx; i++) {
//
//		int tid = j * Mrx + i;
//
//		double fact = 1. / (4.* PI *b);
//		fact = fact * fact;
//
//		double2 temp_sum;
//
//		temp_sum.x = 0.;
//		temp_sum.y = 0.;
//
//		double2 temp_sum_partial;
//
//		double d_n_si = i - (int)Mrx / 2;
//		double d_n_sj = j - (int)Mry / 2;
//
//		for (int kk = 0; kk < len_in; kk++) {
//
//			double x_temp = x[kk];
//			double y_temp = y[kk];
//
//			if ((abs(x_temp - d_n_si) <= msp) && (abs(y_temp - d_n_sj) <= msp)) {
//
//				double exp_temp = exp(b * (Dy * Dy * d_n_sj * d_n_sj + Dx * Dx * d_n_si * d_n_si) - t1 * (((x_temp - d_n_si) * (x_temp - d_n_si) + (y_temp - d_n_sj) * (y_temp - d_n_sj))));
//
//				temp_sum.x = temp_sum.x + f[kk][0] * exp_temp;
//				temp_sum.y = temp_sum.y + f[kk][1] * exp_temp;
//
//			}
//
//		}
//		
//		f_tau[tid][0] = temp_sum.x;
//		f_tau[tid][1] = temp_sum.y;
//
//	}
//
//}


/***************************/
/* CIRCHSHIFT AND FFTSHIFT */
/***************************/
void circAndFFTShiftsFunction(fftw_complex * __restrict out, const fftw_complex * __restrict in, const int xdim, const int ydim, const int xshift, const int yshift) {

	omp_set_nested(1);
	int j;
	#pragma omp parallel for private(j) num_threads(4)
	for (int i = 0; i < xdim; i++) {
		#pragma omp parallel for num_threads(4)
		for (j = 0; j < ydim; j++) {

			int ii = (i + xshift) % xdim;
			int jj = (j + yshift) % ydim;

			out[jj * xdim + ii][0] = (1. - 2 * ((i + j) & 1)) * in[j * xdim + i][0];
			out[jj * xdim + ii][1] = (1. - 2 * ((i + j) & 1)) * in[j * xdim + i][1];
		}
	}
}

/*******************************************************/
/* SPECTRAL CONVOLUTION FUNCTION AND POST-COMPENSATION */
/*******************************************************/
void spectralConvolutionFunctionCompensation(fftw_complex * __restrict F_stm, const fftw_complex * __restrict F_st,	const double * __restrict s, const double * __restrict t, const double sb, const double tb,	const int msp,
											 const double t1, const int len_out, const double Ds, const double Dt, const double Dx, const double Dy, const double xb, const double yb, const double b, const int Mrx, const int Mry) {

	omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for
	for (int kk = 0; kk < len_out; kk++) {

		int nsn = (int)floor(Mrx / 2 + (s[kk] - sb) / Ds);
		double difs = (Mrx / 2 + (s[kk] - sb) / Ds) - nsn;
		int ntn = (int)floor(Mry / 2 + (t[kk] - tb) / Dt);
		double dift = (Mry / 2 + (t[kk] - tb) / Dt) - ntn;

		fftw_complex temp4; temp4[0] = 0.; temp4[1] = 0.;

		double temp1 = Dy * (t[kk] - tb); temp1 = temp1 * temp1;
		double temp2 = Dx * (s[kk] - sb); temp2 = temp2 * temp2;
		double temp = exp(b * (temp1 + temp2));

		for (int nn = nsn - msp; nn <= nsn + msp; nn++){
			for (int mm = ntn - msp; mm <= ntn + msp; mm++){

				int i = mm - ntn + msp;
				int j = nn - nsn + msp;
				double n_pi = i - msp;
				double n_pj = j - msp;

				double temp3 = exp(-t1 * ((dift - n_pi) * (dift - n_pi) + (difs - n_pj) * (difs - n_pj)));

				temp4[0] = temp4[0] + F_st[mm * Mrx + nn][0] * temp3;
				temp4[1] = temp4[1] + F_st[mm * Mrx + nn][1] * temp3;
			}
		}

		double arg = (s[kk] - sb) * xb + (t[kk] - tb) * yb;
		
		F_stm[kk][0] = (temp4[0] * cos(arg) + temp4[1] * sin(arg)) * temp;
		F_stm[kk][1] = (temp4[1] * cos(arg) - temp4[0] * sin(arg)) * temp;
	}
}

/**********************/
/* ROUTINE NUFFT-3 2D */
/**********************/
void NUFFT_2D_CPU(const double * __restrict x, const double * __restrict y,	const double * __restrict s, const double * __restrict t, const fftw_complex * __restrict f, fftw_complex * __restrict F,
				  const double eps, const int len_in, const int len_out) {

	//TimingCPU timerCPU;
	//std::ofstream timingFile;
	//timingFile.open("timingFile.txt");

	/************************************/
	/* SETTING THE ALGORITHM PARAMETERS */
	/************************************/
	thrust::host_vector<double> h_x2(x, x + len_in);
	thrust::host_vector<double> h_y2(y, y + len_in);
	thrust::host_vector<double> h_s2(s, s + len_out);
	thrust::host_vector<double> h_t2(t, t + len_out);

	double Max_x = thrust::reduce(h_x2.begin(), h_x2.end(), -FLT_MAX, thrust::maximum<double>());
	double min_x = thrust::reduce(h_x2.begin(), h_x2.end(),  FLT_MAX, thrust::minimum<double>());
	double xb = (Max_x + min_x) / 2;
	double X1 = fabs(Max_x - xb);

	double Max_y = thrust::reduce(h_y2.begin(), h_y2.end(), -FLT_MAX, thrust::maximum<double>());
	double min_y = thrust::reduce(h_y2.begin(), h_y2.end(),  FLT_MAX, thrust::minimum<double>());
	double yb = (Max_y + min_y) / 2;
	double Y1 = fabs(Max_y - yb);

	double Max_s = thrust::reduce(h_s2.begin(), h_s2.end(), -FLT_MAX, thrust::maximum<double>());
	double min_s = thrust::reduce(h_s2.begin(), h_s2.end(),  FLT_MAX, thrust::minimum<double>());
	double sb = (Max_s + min_s) / 2;
	double S = fabs(Max_s - sb);

	double Max_t = thrust::reduce(h_t2.begin(), h_t2.end(), -FLT_MAX, thrust::maximum<double>());
	double min_t = thrust::reduce(h_t2.begin(), h_t2.end(),  FLT_MAX, thrust::minimum<double>());
	double tb = (Max_t + min_t) / 2;
	double T = fabs(Max_t - tb);

	double R = 2.2;
	double msp = 2 + 2 * (R * R) * (-log(eps / 76)) / (PI * (R * R - 2)); // --- 2*msp are the numbers of non negligible samples of the Gaussian bell function

	int Mrx = 2 * ceil((X1 * S / PI) * R * R + R * msp);				  // --- Number of points who satisfies condition Ca and Cb
	int Mry = 2 * ceil((Y1 * T / PI) * R * R + R * msp);

	double Dx = 2 * PI / (Mrx); double Ds = 1;
	double Dy = 2 * PI / (Mry); double Dt = 1;

	double kk_x = X1 / (PI / 2);		double kk_y = Y1 / (PI / 2);

	thrust::transform(h_x2.begin(), h_x2.end(), thrust::make_constant_iterator(kk_x), h_x2.begin(), thrust::divides<double>());
	thrust::transform(h_y2.begin(), h_y2.end(), thrust::make_constant_iterator(kk_y), h_y2.begin(), thrust::divides<double>());
	thrust::transform(h_s2.begin(), h_s2.end(), thrust::make_constant_iterator(kk_x), h_s2.begin(), thrust::multiplies<double>());
	thrust::transform(h_t2.begin(), h_t2.end(), thrust::make_constant_iterator(kk_y), h_t2.begin(), thrust::multiplies<double>());

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

	fftw_complex *f_tau = (fftw_complex *)malloc(Mrx * Mry * sizeof(fftw_complex));
	memset(f_tau, 0, Mrx * Mry * sizeof(fftw_complex));
	//timerCPU.StartCounter();
	spatialConvolutionFunction(f, thrust::raw_pointer_cast(h_x2.data()), thrust::raw_pointer_cast(h_y2.data()), f_tau, t1, len_in, msp, xb, yb, sb, tb, Dx, Dy, b, Mrx, Mry);
	//timingFile << "timing spatialConvolutionFunction " << timerCPU.GetCounter() << "\n";

	/***********************/
	/* STEP 2: STANDARD FFT*/
	/***********************/
	//timerCPU.StartCounter();
	fftw_complex *F_stx = (fftw_complex*)fftw_malloc(Mrx * Mry * sizeof(fftw_complex));
	//fftw_plan p = fftw_plan_dft_2d(Mry, Mrx, f_tau, f_tau, FFTW_FORWARD, FFTW_ESTIMATE);
	//fftw_execute(p);
	//fftw_destroy_plan(p);	
	DFTI_DESCRIPTOR_HANDLE pMKL;
	MKL_LONG len[2]; len[0] = Mry; len[1] = Mrx;
	DftiCreateDescriptor(&pMKL, DFTI_DOUBLE, DFTI_COMPLEX, 2, len);
	DftiCommitDescriptor(pMKL);
	DftiComputeForward(pMKL, f_tau);
	DftiFreeDescriptor(&pMKL);
	//timingFile << "timing FFT " << timerCPU.GetCounter() << "\n";
	
	//timerCPU.StartCounter();
	fftw_complex *F_stx2 = (fftw_complex *)malloc(Mrx * Mry * sizeof(fftw_complex));
	circAndFFTShiftsFunction(F_stx2, f_tau, Mrx, Mry, Mrx / 2, Mry / 2);
	//timingFile << "timing post FFT " << timerCPU.GetCounter() << "\n";

	/*************************************************************************/
	/*STEP 3: GAUSSIAN GRIDDING IN THE SPECTRAL DOMAIN AND POST-COMPENSATION */
	/*************************************************************************/
	//timerCPU.StartCounter();
	spectralConvolutionFunctionCompensation(F, F_stx2, thrust::raw_pointer_cast(h_s2.data()), thrust::raw_pointer_cast(h_t2.data()), sb, tb, msp, t1, len_out, Ds, Dt, Dx, Dy, xb, yb, b, Mrx, Mry);
	//timingFile << "spectralConvolutionFunctionCompensation " << timerCPU.GetCounter() << "\n";

	/***********************************/
	/* FREEING THE ALLOCATED VARIABLES */
	/***********************************/
	free(f_tau);
	free(F_stx2);
}
