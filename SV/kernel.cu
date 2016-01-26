#define maxNrows				20
#define maxNcols				20
#define DEBUG
//#define DEBUG_SAVE
#define TEMPLATE

#include <iostream>
#include <time.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "Utilities.cuh"
#include "InputOutput.cuh"
#include "TimingGPU.cuh"

#define real_bool	1				// --- 1 is for double, 0 for float
#if (real_bool == 1)
	#define real double
	#define real_MAX DBL_MAX;
#else
	#define real float
	#define real_MAX FLT_MAX;
#endif

#define numGPUs		1

#include "sv.h"

/********/
/* MAIN */
/********/
int main() {

    const unsigned int Nrows			= 6;
	const unsigned int Ncols			= 6;
	const unsigned int numMatrices		= 5000;

	const real		   tol				= 0.0000001;

    const int blockSizeBidiagonalize	= 128;
    const int blockSizeExtractDiagonals	= 128;
    const int blockSizeRearrange		= 128;
    const int blockSizeSturm			= 32;

    srand(time(NULL));

#ifdef DEBUG
	TimingGPU timerGPU;

	double	rearrangeTime				= 0.,
			bidiagonalizationTime		= 0.,
			tridiagAndBisectionTime		= 0.,
			totalTime					= 0.,
			hostToDeviceTime			= 0.,
			totalTimeWithTransfers		= 0.,
			deviceToHostTime			= 0.,
			totalTransfersTime			= 0.,
			totalTimeTest				= 0.;
#endif

    real *inputMatrices; cudaMallocHost(&inputMatrices, Nrows * Ncols * numMatrices * numGPUs * sizeof(real));

#ifdef DEBUG
	real *singularValuesHost; cudaMallocHost(&singularValuesHost, Ncols * numMatrices * numGPUs * sizeof(real));
#endif

    svdPlan<real> plan[numGPUs];
    for (int k = 0; k < numGPUs; k++) createPlan(plan[k], Nrows, Ncols, numMatrices, k);

	const unsigned int numExecutions = 1;
    for (unsigned k = 0; k < numExecutions; k++) {

        // --- Generate random matrices
	    srand(k);
		int N = 5;
	    for (int h = 0; h < Nrows * Ncols * numMatrices * numGPUs; h++)
	            inputMatrices[h]= 2000 * ((double)rand() / (double)(RAND_MAX - 0.2) * (1. / (double)N)) + 100 * ((double)rand() / (double)(RAND_MAX - 0.2) *
				                (1. / (double)N)) + 24;

#ifdef DEBUG_SAVE
		saveCPUrealtxt(inputMatrices, "/home/angelo/cuda-workspace/SVD/Release/inputMatrices.txt", Nrows * Ncols * numMatrices * numGPUs);
#endif

		// --- Compute batched SVD
#ifdef DEBUG
		timerGPU.StartCounter();
#endif
#ifdef TEMPLATE
		my_svd<real, numMatrices, Nrows, Ncols, blockSizeBidiagonalize, blockSizeExtractDiagonals, blockSizeRearrange, blockSizeSturm>(plan, inputMatrices, rearrangeTime, bidiagonalizationTime, tridiagAndBisectionTime, hostToDeviceTime, tol);
#else
		my_svd<real, blockSizeSturm>(plan, inputMatrices, rearrangeTime, bidiagonalizationTime, tridiagAndBisectionTime, hostToDeviceTime, numMatrices, Nrows, Ncols, tol, blockSizeBidiagonalize, blockSizeExtractDiagonals, blockSizeRearrange);
#endif
#ifdef DEBUG
    	totalTimeTest		+= timerGPU.GetCounter();
#endif

#ifdef DEBUG
		timerGPU.StartCounter();
		for (int k = 0; k < numGPUs; k++) {
			gpuErrchk(cudaSetDevice(k));
			gpuErrchk(cudaMemcpyAsync(singularValuesHost + k * Ncols * numMatrices, plan[k].singularValues, Ncols * numMatrices * sizeof(real), cudaMemcpyDeviceToHost));
		}

		deviceToHostTime	+= timerGPU.GetCounter();
#endif
#ifdef DEBUG_SAVE
		saveCPUrealtxt(singularValuesHost, "/home/angelo/cuda-workspace/SVD/Release/singularValues.txt", Ncols * numMatrices * numGPUs);
#endif
	}

    rearrangeTime			/= (numExecutions * numMatrices);
	bidiagonalizationTime	/= (numExecutions * numMatrices);
	tridiagAndBisectionTime	/= (numExecutions * numMatrices);
	hostToDeviceTime		/= (numExecutions * numMatrices);
	deviceToHostTime		/= (numExecutions * numMatrices);
	totalTimeTest			/= (numExecutions * numMatrices);

    totalTimeWithTransfers	= bidiagonalizationTime + tridiagAndBisectionTime + hostToDeviceTime + deviceToHostTime;
    totalTime				= bidiagonalizationTime + tridiagAndBisectionTime;
    totalTransfersTime		= hostToDeviceTime   + deviceToHostTime;

	std::cout << std::scientific << "Nrows \t\t\t\t: " << Nrows << "\n";
	std::cout << std::scientific << "Ncols \t\t\t\t: " << Ncols << "\n";
	std::cout << std::scientific << "numMatrices \t\t\t: " << numMatrices << "\n";
	std::cout << std::scientific << "Rearrange time \t\t\t: " << rearrangeTime << "\n";
	std::cout << std::scientific << "Bidiagonalization time \t\t: " << bidiagonalizationTime << "\n";
	std::cout << std::scientific << "Tridiag and Bisection time \t: " << tridiagAndBisectionTime << "\n";
	std::cout << std::scientific << "Total transfers time \t\t: " << totalTransfersTime << "\n";
	std::cout << std::scientific << "Total time \t\t\t: " << totalTime << "\n";
	std::cout << std::scientific << "Total time with transfers\t: " << totalTimeWithTransfers << "\n";
	std::cout << std::scientific << "Total time test\t\t\t: " << totalTimeTest << "\n";

    for (int k = 0; k < numGPUs; k++) destroyPlan(plan[k], k);

	gpuErrchk(cudaDeviceReset());

    std::cout << "finish\n";

    return 0;
}
