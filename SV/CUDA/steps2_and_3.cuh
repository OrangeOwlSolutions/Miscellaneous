#ifndef __STEPS2_AND_3_H__
#define __STEPS2_AND_3_H__

#include <cfloat>
#include <cstdio>

//#define DEBUG_STURM_BISECTION

/*****************************/
/* CUSTOMIZED fabs FUNCTIONS */
/*****************************/
__host__ __device__ float  fabsT(float  x) { return fabsf(x); }
__host__ __device__ double fabsT(double x) { return fabs (x); }

/*****************************/
/* CUSTOMIZED fmax FUNCTIONS */
/*****************************/
__host__ __device__ float  fmaxT(float  x, float  y) { return fmaxf(x, y); }
__host__ __device__ double fmaxT(double x, double y) { return fmax (x, y); }

/**********************************************/
/* COMPUTE TRIDIAGONAL MATRIX DEVICE FUNCTION */
/**********************************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Ncols>
__host__ __device__ __forceinline__ void computeTridiagonalMatrix(const T *d, const T *e, T * __restrict__ alpha, T * __restrict__ beta) {
#else
template<class T>
__host__ __device__ void computeTridiagonalMatrix(const T *d, const T *e, T * __restrict__ alpha, T * __restrict__ beta,
	                                              const unsigned int numMatrices, const unsigned int Ncols){
#endif

    alpha[0]	= d[0] * d[0];
	beta[0]		= d[0] * e[0];

#pragma unroll
	for (unsigned int i = 1; i < Ncols - 1; i++) {
		alpha[i * numMatrices] = e[(i - 1) * numMatrices] * e[(i - 1) * numMatrices] + d[i * numMatrices] * d[i * numMatrices];
        beta [i * numMatrices] = d[i * numMatrices] * e[i * numMatrices];
    }

	if (Ncols > 1) alpha[(Ncols - 1) * numMatrices] = e[(Ncols - 2) * numMatrices] * e[(Ncols - 2) * numMatrices] + d[(Ncols - 1) * numMatrices] * d[(Ncols - 1) * numMatrices];

}

/*****************/
/* COMPUTE PIVOT */
/*****************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Ncols>
__host__ __device__ __forceinline__ void computePivot(const T * __restrict__ b, T * __restrict__ pivot) {
#else
template<class T>
__host__ __device__ __forceinline__ void computePivot(const T * __restrict__ b, T * __restrict__ pivot, const unsigned int numMatrices, const unsigned int Ncols) {
#endif

    const unsigned int length = Ncols - 1;
    T currentMax = b[0];
    for (unsigned int i = 1; i < length; i++) {
        T temp = b[i * numMatrices];
        if (temp > currentMax) currentMax = temp;
    }

    *pivot = (static_cast<T>(2) * fmaxT(static_cast<T>(1), currentMax)) / real_MAX;
}

//#ifdef TEMPLATE
//template<const unsigned int numMatrices, const unsigned int Ncols>
//__host__ __device__ __forceinline__ void computePivot(const double * __restrict__ b, double * __restrict__ pivot) {
//#else
//__host__ __device__ __forceinline__ void computePivot(const double * __restrict__ b, double * __restrict__ pivot, const unsigned int numMatrices, const unsigned int Ncols) {
//#endif
//
//    const unsigned int length = Ncols - 1;
//    double currentMax = b[0];
//    for (unsigned int i = 1; i < length; i++) {
//        double temp = b[i * numMatrices];
//        if (temp > currentMax) currentMax = temp;
//    }
//
//    *pivot = (static_cast<double>(2) * fmax(static_cast<double>(1), currentMax)) / DBL_MAX;
//}
//
//__host__ __device__ __forceinline__ void computePivot(const float * __restrict__ b, float * __restrict__ pivot, const unsigned int Ncols, const unsigned int numMatrices) {
//
//    const unsigned int length = Ncols - 1;
//    float currentMax = b[0];
//    for (unsigned int i = 1; i < length; i++) {
//        float temp = b[i * numMatrices];
//        if (temp > currentMax) currentMax = temp;
//    }
//
//    *pivot = (static_cast<float>(2) * fmaxf(static_cast<float>(1), currentMax)) / FLT_MAX;
//}

/***********************************/
/* COMPUTE INITIAL SEARCH INTERVAL */
/***********************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Ncols>
__host__ __device__ __forceinline__ void computeInitialInterval(const T * __restrict__ d, const T * __restrict__ b, T &alpha,
																T &beta) {
#else
template<class T>
__host__ __device__ __forceinline__ void computeInitialInterval(const T * __restrict__ d, const T * __restrict__ b, T &alpha,
																T &beta, const unsigned int numMatrices, const unsigned int Ncols) {
#endif

    T temp1 = static_cast<T>(0);
    T temp2 = static_cast<T>(0);

    alpha	= d[0] - fabsT(b[0]);
    beta	= fabsT(d[0]) + fabsT(b[0]);

	temp1 = d[(Ncols - 1) * numMatrices] - fabsT(b[(Ncols - 2) * numMatrices]);
    temp2 = fabsT(d[(Ncols - 1) * numMatrices]) + fabsT(b[(Ncols - 2) * numMatrices]);

	if (temp1 < alpha) alpha = temp1;
    if (temp2 > beta)  beta	 = temp2;

	for (unsigned int i = 1; i < Ncols - 1; i++) {
        temp1 = d[i * numMatrices] - fabsT(b[(i - 1) * numMatrices]) - fabsT(b[i * numMatrices]);
        temp2 = d[i * numMatrices] + fabsT(b[(i - 1) * numMatrices]) + fabsT(b[i * numMatrices]);
        if (temp1 < alpha) alpha = temp1;
        if (temp2 > beta)  beta  = temp2;
    }

}

//#ifdef TEMPLATE
//template<const unsigned int numMatrices, const unsigned int Ncols>
//__host__ __device__ __forceinline__ void computeInitialInterval(const double * __restrict__ d, const double * __restrict__ b, double &alpha,
//																double &beta) {
//#else
//__host__ __device__ __forceinline__ void computeInitialInterval(const double * __restrict__ d, const double * __restrict__ b, double &alpha,
//																double &beta, const unsigned int numMatrices, const unsigned int Ncols) {
//#endif
//
//    double temp1 = static_cast<double>(0);
//    double temp2 = static_cast<double>(0);
//
//    alpha	= d[0] - fabs(b[0]);
//    beta	= fabs(d[0]) + fabs(b[0]);
//
//	temp1 = d[(Ncols - 1) * numMatrices] - fabs(b[(Ncols - 2) * numMatrices]);
//    temp2 = fabs(d[(Ncols - 1) * numMatrices]) + fabs(b[(Ncols - 2) * numMatrices]);
//
//	if (temp1 < alpha) alpha = temp1;
//    if (temp2 > beta)  beta	 = temp2;
//
//	for (unsigned int i = 1; i < Ncols - 1; i++) {
//        temp1 = d[i * numMatrices] - fabs(b[(i - 1) * numMatrices]) - fabs(b[i * numMatrices]);
//        temp2 = d[i * numMatrices] + fabs(b[(i - 1) * numMatrices]) + fabs(b[i * numMatrices]);
//        if (temp1 < alpha) alpha = temp1;
//        if (temp2 > beta)  beta  = temp2;
//    }
//
//}
//
//#ifdef TEMPLATE
//template<const unsigned int numMatrices, const unsigned int Ncols>
//__host__ __device__ __forceinline__ void computeInitialInterval(const float * __restrict__ d, const float * __restrict__ b, float &alpha,
//																float &beta) {
//#else
//__host__ __device__ __forceinline__ void computeInitialInterval(const float * __restrict__ d, const float * __restrict__ b, float &alpha,
//																float &beta, const unsigned int numMatrices, const unsigned int Ncols) {
//#endif
//
//    float temp1 = static_cast<float>(0);
//    float temp2 = static_cast<float>(0);
//
//    alpha	= d[0] - fabsf(b[0]);
//    beta	= fabsf(d[0]) + fabsf(b[0]);
//
//	temp1 = d[(Ncols - 1) * numMatrices] - fabsf(b[(Ncols - 2) * numMatrices]);
//    temp2 = fabsf(d[(Ncols - 1) * numMatrices]) + fabsf(b[(Ncols - 2) * numMatrices]);
//
//	if (temp1 < alpha) alpha = temp1;
//    if (temp2 > beta)  beta	 = temp2;
//
//	for (unsigned int i = 1; i < Ncols - 1; i++) {
//        temp1 = d[i * numMatrices] - fabsf(b[(i - 1) * numMatrices]) - fabsf(b[i * numMatrices]);
//        temp2 = d[i * numMatrices] + fabsf(b[(i - 1) * numMatrices]) + fabsf(b[i * numMatrices]);
//        if (temp1 < alpha) alpha = temp1;
//        if (temp2 > beta)  beta  = temp2;
//    }
//
//}

/*************/
/* CHECKSIGN */
/*************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Ncols>
__host__ __device__ void chcksign(const T * __restrict__ d, const T * __restrict__ b, const T pivot, const T x, unsigned int &numChanges) {
#else
template<class T>

__host__ __device__ void chcksign(const T * __restrict__ d, const T * __restrict__ b, const T pivot, const T x, unsigned int &numChanges,
	                              const unsigned int numMatrices, const unsigned int Ncols) {
#endif

    numChanges = 0;

    T q = d[0] - x;
    if (fabsT(q) <= pivot) q = -pivot;
    if (q < static_cast<T>(0)) numChanges = numChanges + 1;

    for (unsigned int i = 2; i <= Ncols; i++) {
        // --- This order of operation preserve monotonicity
		q = (d[(i - 1) * numMatrices] - ((b[(i - 2) * numMatrices] * b[(i - 2) * numMatrices]) / q)) - x;
        if (fabsT(q) <= pivot) q = -pivot;
        if (q < static_cast<T>(0)) numChanges = numChanges + 1;
    }
}

/************************************/
/* STURM METHOD BASED ON BISECTIONS */
/************************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Ncols, const int blockSizeSturm>
__global__ void sturmBisection(const T * __restrict__ d, const T * __restrict__ b, T * __restrict__ dd, T * __restrict__ bb,
                               T * __restrict__ singularValues, const T tol){
#else
template<class T, const int blockSizeSturm>
__global__ void sturmBisection(const T * __restrict__ d, const T * __restrict__ b, T * __restrict__ dd, T * __restrict__ bb,
                               T * __restrict__ singularValues, const unsigned int numMatrices, const unsigned int Ncols, const T tol) {
#endif

    int tidx = threadIdx.x + blockDim.x * blockIdx.x;
    int tidy = threadIdx.y;

    T alpha, beta, pivot;

    __shared__ T alpha_shared [blockSizeSturm];
    __shared__ T beta_shared  [blockSizeSturm];
    __shared__ T pivmin_shared[blockSizeSturm];

    if (tidx < numMatrices && tidy == 0) {
#ifdef TEMPLATE
        computeTridiagonalMatrix<T, numMatrices, Ncols>(d + tidx, b + tidx, dd + tidx, bb + tidx);
#else
        computeTridiagonalMatrix(d + tidx, b + tidx, dd + tidx, bb + tidx, numMatrices, Ncols);
#endif
#ifdef TEMPLATE
	    computePivot<T, numMatrices, Ncols>(bb + tidx, &pivmin_shared[threadIdx.x]);
#else
	    computePivot<T>(bb + tidx, &pivmin_shared[threadIdx.x], numMatrices, Ncols);
#endif
#ifdef TEMPLATE
		computeInitialInterval<T, numMatrices, Ncols>(dd + tidx, bb + tidx, alpha_shared[threadIdx.x], beta_shared[threadIdx.x]);
#else
		computeInitialInterval(dd + tidx, bb + tidx, alpha_shared[threadIdx.x], beta_shared[threadIdx.x], numMatrices, Ncols);
#endif
    }
    __syncthreads();

    if (tidx < numMatrices && tidy < Ncols) {

        unsigned int eigIndex = tidy % Ncols;
        //unsigned int eigIndex = tidy & (Ncols - 1);

        unsigned int	numChanges = 0;
        T				c		   = 0;

        alpha	= alpha_shared [threadIdx.x];
        beta	= beta_shared  [threadIdx.x];
        pivot	= pivmin_shared[threadIdx.x];

        T dist	= fabsT(beta - alpha); // --- Initial amplitude of the search interval
        T s		= fabsT(beta) + fabsT(alpha);

        while (dist > tol * s) { // --- Quarteroni
        //while(dist > max(tol, DBL_EPSILON * max(abs(alpha), abs(beta)))) { // --- Demmel - Dhillon

            c = (alpha + beta) / static_cast<T>(2);

#ifdef TEMPLATE
			chcksign<T, numMatrices, Ncols>(dd + tidx, bb + tidx, pivot, c, numChanges);
#else
			chcksign(dd + tidx, bb + tidx, pivot, c, numChanges, numMatrices, Ncols);
#endif

            if (numChanges > (Ncols - (eigIndex + 1)))	beta  = c;
            else										alpha = c;

            dist	= fabsT(beta - alpha);
            s		= fabsT(beta) + fabsT(alpha);
        }

        singularValues[eigIndex + tidx * Ncols] = sqrt(fabsT(c));
        //printf("%f\n", singularValues[eigIndex + tidx * Ncols]);
    }
}

/*********************************************************************/
/* TRIDIAGONALIZE AND FIND THE EIGENVALUES OF THE TRIDIAGONAL MATRIX */
/*********************************************************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Ncols, const int blockSizeSturm>
void tridiagonalizeAndFindEigenvalues(const T * __restrict__ d, const T * __restrict__ b, T * __restrict__ dd, T * __restrict__ bb,
	                                  T * __restrict__ singularValues, const T tol) {
#else
template<class T, const int blockSizeSturm>
void tridiagonalizeAndFindEigenvalues(const T * __restrict__ d, const T * __restrict__ b, T * __restrict__ dd, T * __restrict__ bb,
	                                  T * __restrict__ singularValues, const unsigned int numMatrices, const unsigned int Ncols, const T tol) {
#endif

	const unsigned int NUM_BLOCKS_STURM = iDivUp(numMatrices, blockSizeSturm);
    dim3 BLOCKGRID(blockSizeSturm, Ncols);

#ifdef TEMPLATE
	gpuErrchk(cudaFuncSetCacheConfig(sturmBisection<T, numMatrices, Ncols, blockSizeSturm>,	cudaFuncCachePreferL1));
#else
	gpuErrchk(cudaFuncSetCacheConfig(sturmBisection<T, blockSizeSturm>,						cudaFuncCachePreferL1));
#endif

#ifdef TEMPLATE
	sturmBisection<T, numMatrices, Ncols, blockSizeSturm><<<NUM_BLOCKS_STURM, BLOCKGRID>>>(d, b, dd, bb, singularValues, tol);
#else
	sturmBisection<T, blockSizeSturm><<<NUM_BLOCKS_STURM, BLOCKGRID>>>(d, b, dd, bb, singularValues, numMatrices, Ncols, tol);
#endif
#ifdef DEBUG_STURM_BISECTION
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
#endif

}

#endif
