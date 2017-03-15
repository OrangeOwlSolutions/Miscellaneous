#include <iostream>
#include <complex>
#include <string.h>
#include <fstream>

#include "cuda.h"
#include <cuda_runtime.h>

#include "mex.h" 

#include "NUFFT3_2D_CUDA.cuh"

#include "Utilities.cuh"
#include "TimingGPU.cuh"

/****************/
/* MEX FUNCTION */
/****************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	TimingGPU timerGPU;
	std::ofstream timingFile;
	timingFile.open("timingFileGPU.txt");
	
	/**********/
	/* INPUTS */
	/**********/
	double *x = (double *)mxGetData(prhs[0]);
	double *y = (double *)mxGetData(prhs[1]);
	double *s = (double *)mxGetData(prhs[2]);
	double *t = (double *)mxGetData(prhs[3]);

	const int len_in	= mxGetN(prhs[0]);
	const int len_out	= mxGetN(prhs[2]);

	double* fr = (double *)mxGetPr(prhs[4]);
	double* fi = (double *)mxGetPi(prhs[4]);

	double2 *f = (double2 *)malloc(len_in * sizeof(double2));
	for (int i = 0; i<len_in; i++) {
		f[i].x = fr[i];
		if (fi != NULL) f[i].y = fi[i];
		else            f[i].y = 0.;
	}

	double* eps = (double*)mxGetData(prhs[5]);

	/***********/
	/* OUTPUTS */
	/***********/
	plhs[0] = mxCreateDoubleMatrix(1, len_out, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, len_out, mxREAL);

	/***************************************/
	/* MOVING THE DATA FROM HOST TO DEVICE */
	/***************************************/
	double *d_x;        gpuErrchk(cudaMalloc((void**)&d_x, len_in*sizeof(double)));
	double *d_y;        gpuErrchk(cudaMalloc((void**)&d_y, len_in*sizeof(double)));
	gpuErrchk(cudaMemcpy(d_x, x, len_in*sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_y, y, len_in*sizeof(double), cudaMemcpyHostToDevice));

	double *d_s;        gpuErrchk(cudaMalloc((void**)&d_s, len_out*sizeof(double)));
	double *d_t;        gpuErrchk(cudaMalloc((void**)&d_t, len_out*sizeof(double)));
	gpuErrchk(cudaMemcpy(d_s, s, len_out*sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_t, t, len_out*sizeof(double), cudaMemcpyHostToDevice));

	double2* d_f;       gpuErrchk(cudaMalloc((void**)&d_f, len_in*sizeof(double2)));
	gpuErrchk(cudaMemcpy(d_f, f, len_in*sizeof(double2), cudaMemcpyHostToDevice));

	/****************************************/
	/* ALLOCATING THE SPACE FOR THE RESULTS */
	/****************************************/
	double *Fr = (double *)mxGetPr(plhs[0]);
	double *Fi = (double *)mxGetPr(plhs[1]);

	double2 *F = (double2 *)malloc(len_out * sizeof(double2));

	double2* d_F_stm; gpuErrchk(cudaMalloc((void**)&d_F_stm, len_out * sizeof(double2)));

	/***************/
	/* COMPUTATION */
	/***************/
	timerGPU.StartCounter();
	NFFT2_2D_GPU(d_x, d_y, d_s, d_t, d_f, d_F_stm, eps[0], len_in, len_out);
	timingFile << "Timing " << len_in << " " << len_out << " " << timerGPU.GetCounter() << "\n";

	/******************************************/
	/* MOVING THE RESULTS FROM DEVICE TO HOST */
	/******************************************/
	gpuErrchk(cudaMemcpy(F, d_F_stm, len_out * sizeof(double2), cudaMemcpyDeviceToHost));

	for (int i = 0; i < len_out; i++) {
		Fr[i] = F[i].x;
		Fi[i] = F[i].y;
	}

	timingFile.close();

	/************************************/
	/* FREEING HOST AND DEVICE MEMORIES */
	/************************************/
	free(f);
	gpuErrchk(cudaFree(d_x));
	gpuErrchk(cudaFree(d_y));
	gpuErrchk(cudaFree(d_s));
	gpuErrchk(cudaFree(d_t));
	gpuErrchk(cudaFree(d_f));
	gpuErrchk(cudaFree(d_F_stm));
}


