#include "stdafx.h" 

#include <complex>
#include <string.h>
#include <fstream>

#include "mex.h" 

#include "fftw3.h"

#include "NUFFT3_2D_C.h"

#include "TimingCPU.h"

/****************/
/* MEX FUNCTION */
/****************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	TimingCPU timerCPU;
	std::ofstream timingFile;
	timingFile.open("timingFileCPU.txt");

	/**********/
	/* INPUTS */
	/**********/
	double *x = (double *)mxGetData(prhs[0]);
	double *y = (double *)mxGetData(prhs[1]);
	double *s = (double *)mxGetData(prhs[2]);
	double *t = (double *)mxGetData(prhs[3]);

	const int len_in  = mxGetN(prhs[0]);
	const int len_out = mxGetN(prhs[2]);

	double* fr = (double *)mxGetPr(prhs[4]);
	double* fi = (double *)mxGetPi(prhs[4]);

	fftw_complex *f = (fftw_complex *)malloc(len_in * sizeof(fftw_complex));
	for (int i = 0; i<len_in; i++)
	{
		f[i][0] = fr[i];
		f[i][1] = fi[i];
	}

	double* eps = (double*)mxGetData(prhs[5]);

	/***********/
	/* OUTPUTS */
	/***********/
	plhs[0] = mxCreateDoubleMatrix(1, len_out, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, len_out, mxREAL);

	/****************************************/
	/* ALLOCATING THE SPACE FOR THE RESULTS */
	/****************************************/
	double *Fr = (double *)mxGetPr(plhs[0]);
	double *Fi = (double *)mxGetPr(plhs[1]);

	fftw_complex *F = (fftw_complex *)malloc(len_out * sizeof(fftw_complex));

	/***************/
	/* COMPUTATION */
	/***************/
	timerCPU.StartCounter();
	NUFFT_2D_CPU(x, y, s, t, f, F, eps[0], len_in, len_out);
	timingFile << "Timing " << len_in << " " << len_out << " " << timerCPU.GetCounter() << "\n";

	/***************************/
	/* REARRANGING THE RESULTS */
	/***************************/
	for (int i = 0; i < len_out; i++)
	{
		Fr[i] = F[i][0];
		Fi[i] = F[i][1];
	}

	timingFile.close();
	
	/******************/
	/* FREEING MEMORY */
	/******************/
	free(f);
}


