#ifndef __SV_H__
#define __SV_H__

#include "step1.cuh"
#include "steps2_and_3.cuh"
#include <time.h>
#include <iomanip>

#include "TimingGPU.cuh"

/**********************/
/* SVD PLAN STRUCTURE */
/**********************/
template<class T>
struct svdPlan {
	T *dev_mat_opt;
    T *dev_mat;
	T *dev_diag;		// --- Diagonal of bidiagonal
	T *dev_supdiag;		// --- Off-diagonal of bidiagonal
	T *alpha;			// --- Diagonal of tridiagonal
	T *beta;			// --- Off-diagonal of tridiagonal
	T *singularValues;	// --- Singular values
};

/*********************/
/* SVD PLAN CREATION */
/*********************/
template<class T>
void createPlan(svdPlan<T>& plan, unsigned int Nrows, unsigned int Ncols, unsigned int numMatrices, unsigned int gpuID) {

    // --- Device allocations
    gpuErrchk(cudaSetDevice(gpuID));
	gpuErrchk(cudaMalloc(&(plan.dev_mat_opt),		Nrows	*	Ncols		*	numMatrices * sizeof(T)));
	gpuErrchk(cudaMalloc(&(plan.dev_mat),			Nrows   *   Ncols		*   numMatrices * sizeof(T)));
 	gpuErrchk(cudaMalloc(&(plan.dev_diag),			            Ncols		*   numMatrices * sizeof(T)));
    gpuErrchk(cudaMalloc(&(plan.dev_supdiag),		            (Ncols - 1) *	numMatrices * sizeof(T)));
    gpuErrchk(cudaMalloc(&(plan.alpha),							Ncols       *   numMatrices * sizeof(T)));
	gpuErrchk(cudaMalloc(&(plan.beta),				            (Ncols - 1) *   numMatrices * sizeof(T)));
	gpuErrchk(cudaMalloc(&(plan.singularValues),	            Ncols 		*   numMatrices * sizeof(T)));

}

/************************/
/* SVD PLAN DESTRUCTION */
/************************/
template<class real_type>
void destroyPlan(svdPlan<real_type>& plan, unsigned int gpuID) {

    // --- Device deallocations
    gpuErrchk(cudaSetDevice(gpuID));
    gpuErrchk(cudaFree(plan.dev_mat_opt));
	gpuErrchk(cudaFree(plan.dev_mat));
	gpuErrchk(cudaFree(plan.dev_diag));
    gpuErrchk(cudaFree(plan.dev_supdiag));
    gpuErrchk(cudaFree(plan.beta));
	gpuErrchk(cudaFree(plan.alpha));
	gpuErrchk(cudaFree(plan.singularValues));
}

/*********************************/
/* HOST-SIDE COMPUTATION ROUTINE */
/*********************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols, const int blockSizeBidiagonalize,
         const int blockSizeExtractDiagonals, const int blockSizeRearrange, const int blockSizeSturm>
void my_svd(svdPlan<T> *plan, T *inputMatrices, double &rearrangeTime, double &bidiagonalizationTime,
	        double &tridiagAndBisectionTime, double &hostToDeviceTime, const T tol) {
#else
template<class T, const int blockSizeSturm>
void my_svd(svdPlan<T> *plan, T *inputMatrices, double &rearrangeTime, double &bidiagonalizationTime,
	        double &tridiagAndBisectionTime, double &hostToDeviceTime, const unsigned int numMatrices, const unsigned int Nrows,
			const unsigned int Ncols, const T tol, const int blockSizeBidiagonalize, const int blockSizeExtractDiagonals, const int blockSizeRearrange) {
#endif

#ifdef DEBUG
	TimingGPU timer;
#endif

	// --- Host -> Device memory transfers
#ifdef DEBUG
	timer.StartCounter();
#endif
	for (int k = 0; k < numGPUs; k++) {
		gpuErrchk(cudaSetDevice(k));
		gpuErrchk(cudaMemcpyAsync(plan[k].dev_mat, inputMatrices + k * Nrows * Ncols * numMatrices, Nrows * Ncols * numMatrices * sizeof(T), cudaMemcpyHostToDevice));
	}
#ifdef DEBUG
	hostToDeviceTime	+= timer.GetCounter();
#endif

	// --- Reorganize the input matrix
#ifdef DEBUG
	timer.StartCounter();
#endif
#ifdef TEMPLATE
	for (int k = 0; k < numGPUs; k++) {
		gpuErrchk(cudaSetDevice(k));
		rearrange<T, numMatrices, Nrows, Ncols, blockSizeRearrange>(plan[k].dev_mat, plan[k].dev_mat_opt);
	}
#else
	for (int k = 0; k < numGPUs; k++) {
		gpuErrchk(cudaSetDevice(k));
		rearrange(plan[k].dev_mat, plan[k].dev_mat_opt, numMatrices, Nrows, Ncols, blockSizeRearrange);
	}
#endif
#ifdef DEBUG
	rearrangeTime	+= timer.GetCounter();
#endif

	// --- Compute bidiagonalization
#ifdef DEBUG
	timer.StartCounter();
#endif
#ifdef TEMPLATE
	for (int k = 0; k < numGPUs; k++) {
		gpuErrchk(cudaSetDevice(k));
		bidiagonalize<T, numMatrices, Nrows, Ncols, blockSizeBidiagonalize, blockSizeExtractDiagonals>(plan[k].dev_mat_opt, plan[k].dev_diag, plan[k].dev_supdiag);
	}
#else
	for (int k = 0; k < numGPUs; k++) {
		gpuErrchk(cudaSetDevice(k));
		bidiagonalize(plan[k].dev_mat_opt, plan[k].dev_diag, plan[k].dev_supdiag, numMatrices, Nrows, Ncols, blockSizeBidiagonalize, blockSizeExtractDiagonals);
	}
#endif
#ifdef DEBUG
	bidiagonalizationTime += timer.GetCounter();
#endif

	// --- Compute tridiagonal matrix and find eigenvalues
#ifdef DEBUG
	timer.StartCounter();
#endif
#ifdef TEMPLATE
	for (int k = 0; k < numGPUs; k++) {
		gpuErrchk(cudaSetDevice(k));
		tridiagonalizeAndFindEigenvalues<T, numMatrices, Ncols, blockSizeSturm>(plan[k].dev_diag, plan[k].dev_supdiag, plan[k].alpha, plan[k].beta, plan[k].singularValues, tol);
	}
#else
	for (int k = 0; k < numGPUs; k++) {
		gpuErrchk(cudaSetDevice(k));
		tridiagonalizeAndFindEigenvalues<T, blockSizeSturm>(plan[k].dev_diag, plan[k].dev_supdiag, plan[k].alpha, plan[k].beta, plan[k].singularValues, numMatrices, Ncols, tol);
	}
#endif
#ifdef DEBUG
	tridiagAndBisectionTime	+= timer.GetCounter();
#endif

}

#endif
