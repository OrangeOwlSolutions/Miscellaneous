#ifndef __STEP1_H__
#define __STEP1_H__

//#define DEBUG_REARRANGE
//#define DEBUG_BIDIAGONALIZE
//#define DEBUG_EXTRACT_DIAGONALS

/*************************************/
/* EXTRACT DIAGONALS KERNEL FUNCTION */
/*************************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols>
__global__ void extractDiagonals(const T * __restrict__ inputMatrices, T * __restrict__ d, T * __restrict__ e) {
#else
template<class T>
__global__ void extractDiagonals(const T * __restrict__ inputMatrices, T * __restrict__ d, T * __restrict__ e,
	                             const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols) {
#endif

	int tid = threadIdx.x + blockDim.x * blockIdx.x;

	if (tid < numMatrices) {

		//for (unsigned int i = 0; i < Ncols; i++)     d[tid + i * numMatrices] = inputMatrices[tid + i * numMatrices * (Ncols + 1)];
		//for (unsigned int i = 0; i < Ncols - 1; i++) e[tid + i * numMatrices] = inputMatrices[tid + numMatrices + i * numMatrices * (Ncols + 1)];
#pragma unroll
		for (unsigned int i = 0; i < Ncols - 1; i++) {
			e[tid + i * numMatrices] = inputMatrices[tid + numMatrices + i * numMatrices * (Ncols + 1)];
			d[tid + i * numMatrices] = inputMatrices[tid + i * numMatrices * (Ncols + 1)];
		}
		d[tid + (Ncols - 1) * numMatrices] = inputMatrices[tid + (Ncols - 1) * numMatrices * (Ncols + 1)];

	}
}

/************************************/
/* COMPUTE LEFT HOUSEHOLDER VECTORS */
/************************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols>
__device__ __forceinline__ void computeLeftHouseholderVectors(const T * __restrict__ inputMatrices, T * __restrict__ outputLeftHourseholderVectors, T &betaq,
											  const unsigned int tid, const unsigned int iterCounter) {
#else
template<class T>
__device__ __forceinline__ void computeLeftHouseholderVectors(const T * __restrict__ inputMatrices, T * __restrict__ outputLeftHourseholderVectors, T &betaq,
											  const unsigned int tid, const unsigned int iterCounter, const unsigned int numMatrices,
											  const unsigned int Nrows, const unsigned int Ncols) {
#endif

	T x0, norm2 = static_cast<T>(0), mu, v0;
#pragma unroll
	for (unsigned int i = 0; i < iterCounter ; i++) outputLeftHourseholderVectors[i] = static_cast<T>(0);

#pragma unroll
	for (unsigned int i = iterCounter + 1; i < Nrows; i++) {
		const T buffer = inputMatrices[tid + iterCounter * numMatrices + i * numMatrices * Ncols];
		outputLeftHourseholderVectors[i] = buffer;
        // --- Evaluate the squared norm of input vector x
        norm2 = norm2 + buffer * buffer;
	}

	x0 = inputMatrices[tid + (iterCounter * numMatrices) * (1 + Ncols)];

	if (norm2 == static_cast<T>(0)) betaq = static_cast<T>(0);
    else {
		mu = sqrt(x0 * x0 + norm2);

		if (x0 <= 0)	v0 = x0 - mu;
        else			v0 = -norm2 / (x0 + mu);

        T temp = static_cast<T>(1) / v0;
        betaq = static_cast<T>(2) * ((v0 * v0) / (norm2 + (v0 * v0)));
        outputLeftHourseholderVectors[iterCounter] = static_cast<T>(1);
#pragma unroll
		for(unsigned int s = iterCounter + 1; s < Nrows ; s++) outputLeftHourseholderVectors[s] *= temp;
	}
}

/*************************************/
/* COMPUTE RIGHT HOUSEHOLDER VECTORS */
/*************************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols>
__device__ __forceinline__ void computeRightHouseholderVectors(const T * __restrict__ inputVectors, T * __restrict__ outputRightHouseholderVectors, T &betap, const unsigned int tid,
											   const unsigned int iterCounter) {
#else
template<class T>
__device__ __forceinline__ void computeRightHouseholderVectors(const T * __restrict__ inputVectors, T * __restrict__ outputRightHouseholderVectors, T &betap, const unsigned int tid,
											   const unsigned int iterCounter, const unsigned int numMatrices, const unsigned int Nrows,
											   const unsigned int Ncols) {
#endif

	T x0, norm2 = static_cast<T>(0), mu, v0;
#pragma unroll
	for (unsigned int i = 0; i < iterCounter + 1 ; i++) outputRightHouseholderVectors[i] = static_cast<T>(0);

#pragma unroll
	for (unsigned int i = iterCounter + 2; i < Ncols ; i++) {
		const T buffer = inputVectors[i];
    	outputRightHouseholderVectors[i] = buffer;
    	// --- Evaluate the squared norm of input vector x
    	norm2 = norm2 + buffer * buffer;
	}

    x0 = inputVectors[iterCounter + 1];

    if(norm2 == 0.) betap = static_cast<T>(0);
    else {
		mu = sqrt(x0 * x0 + norm2);

		if (x0 <= 0)	v0 = x0 - mu;
    	else			v0 = -norm2 / (x0 + mu);

        T temp = static_cast<T>(1) / v0;
    	betap = static_cast<T>(2) * ((v0 * v0) / (norm2 + v0 * v0));
    	outputRightHouseholderVectors[iterCounter + 1] = static_cast<T>(1);
#pragma unroll
		for (unsigned int s = iterCounter + 2; s < Ncols ; s++) outputRightHouseholderVectors[s] *= temp;
	}
}

/*************************************/
/* BIDIAGONALIZATION KERNEL FUNCTION */
/*************************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols>
__global__ void bidiagonalizeKernel(T *inputMatrices, const unsigned int maxIter) {
#else
template<class T>
__global__ void bidiagonalizeKernel(T *inputMatrices, const unsigned int maxIter, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols) {
#endif

    int tid = threadIdx.x + blockDim.x * blockIdx.x;

#ifdef TEMPLATE
	T leftHouseholderVectors[Nrows];
    T rightHouseholderVectors[Ncols];
    T inputVectorsForRightHouseholderReflections[Ncols];
    T wi[Nrows];
    T xi[Ncols];
    T zi[Ncols];
#else
	T leftHouseholderVectors[maxNrows];
    T rightHouseholderVectors[maxNcols];
    T inputVectorsForRightHouseholderReflections[maxNcols];
    T wi[maxNrows];
    T xi[maxNcols];
    T zi[maxNcols];
#endif
    T betau = static_cast<T>(0), betav = static_cast<T>(0);

#pragma unroll
	for (unsigned int iterCounter = 0; iterCounter < maxIter; iterCounter++) {

		if (iterCounter < Ncols - 2) {

			if (tid < numMatrices) {

				// --- Calculate left Householder vectors
#ifdef TEMPLATE
				computeLeftHouseholderVectors<T, numMatrices, Nrows, Ncols>(inputMatrices, leftHouseholderVectors, betau, tid, iterCounter);
#else
				computeLeftHouseholderVectors(inputMatrices, leftHouseholderVectors, betau, tid, iterCounter, numMatrices, Nrows, Ncols);
#endif
				// --- Calculate xi and input vectors for right Householder reflections
				{
#pragma unroll
					for (unsigned int j = 0; j < Ncols; j++) {
						T sum = static_cast<T>(0);
#pragma unroll
						for (unsigned int i = 0; i < Nrows; i++) sum += -betau * leftHouseholderVectors[i] * inputMatrices[tid + i * numMatrices * Ncols + j * numMatrices];
						inputVectorsForRightHouseholderReflections[j] = sum + inputMatrices[tid + j*numMatrices + iterCounter*numMatrices*Ncols];
						xi[j] = -sum;
					}
				}

				// --- Calculate right Householder vectors
				//computeRightHouseholderVectors(inputVectorsForRightHouseholderReflections, rightHouseholderVectors, betav, tid, iterCounter, numMatrices, Nrows, Ncols);
#ifdef TEMPLATE
				computeRightHouseholderVectors<T, numMatrices, Nrows, Ncols>(inputVectorsForRightHouseholderReflections, rightHouseholderVectors, betav, tid, iterCounter);
#else
				computeRightHouseholderVectors(inputVectorsForRightHouseholderReflections, rightHouseholderVectors, betav, tid, iterCounter, numMatrices, Nrows, Ncols);
#endif
				// --- Calculate wi
				{
#pragma unroll
					for (unsigned int i = 0; i < Nrows; i++) {
						T sum = static_cast<T>(0);
#pragma unroll
						for (unsigned int j = 0; j < Ncols; j++) sum += betav * rightHouseholderVectors[j] * inputMatrices[tid + i * numMatrices * Ncols + j * numMatrices];
						wi[i] = sum;
					}
				}

				// --- Calculate zi
				{
					//T temp	= static_cast<T>(0);
					T sum	= static_cast<T>(0);
#pragma unroll
					for(unsigned int i = 0; i < Ncols; i++) sum += xi[i] * rightHouseholderVectors[i];
					//temp = sum;
					//for (unsigned int i = 0; i < Ncols; i++) zi[i] = xi[i] - betav * temp * rightHouseholderVectors[i];
#pragma unroll
					for (unsigned int i = 0; i < Ncols; i++) zi[i] = xi[i] - betav * sum * rightHouseholderVectors[i];
				}

				// --- Matrix update
#pragma unroll
				for (unsigned int r = 0; r < Nrows; r++)
#pragma unroll
					for(unsigned int c = 0; c < Ncols; c++)
						inputMatrices[tid + r * numMatrices * Ncols + c * numMatrices] = -leftHouseholderVectors[r] * zi[c]
						                                                                 -wi[r] * rightHouseholderVectors[c] +
																						  inputMatrices[tid + r * numMatrices * Ncols + c * numMatrices];
			}

		}

		if (iterCounter >= (Ncols - 2)) {

			if (tid < numMatrices) {

				// --- Calculate the Householder vectors
				//computeLeftHouseholderVectors(inputMatrices, leftHouseholderVectors, betau, tid, numMatrices, Nrows, Ncols, iterCounter);
#ifdef TEMPLATE
				computeLeftHouseholderVectors<T, numMatrices, Nrows, Ncols>(inputMatrices, leftHouseholderVectors, betau, tid, iterCounter);
#else
				computeLeftHouseholderVectors(inputMatrices, leftHouseholderVectors, betau, tid, iterCounter, numMatrices, Nrows, Ncols);
#endif
				// --- Calculate the xi
#pragma unroll
				for (unsigned int j = 0; j < Ncols; j++) {
					T sum = static_cast<T>(0);
#pragma unroll
					for (unsigned int i = 0; i < Nrows; i++) sum += leftHouseholderVectors[i] * inputMatrices[tid + i * numMatrices * Ncols + j * numMatrices];
					xi[j] = sum;
				}

				// --- Matrix update
#pragma unroll
				for (unsigned int r = 0; r < Nrows; r++)
#pragma unroll
					for (unsigned int c = 0; c < Ncols; c++)
						inputMatrices[tid + r * numMatrices * Ncols + c * numMatrices] = -betau * xi[c] * leftHouseholderVectors[r] +
							                                                              inputMatrices[tid + r * numMatrices * Ncols + c * numMatrices];
			}
		}
	}
}

/*****************/
/* BIDIAGONALIZE */
/*****************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols, const int blockSizeBidiagonalize,
         const int blockSizeExtractDiagonals>
void bidiagonalize(T *inputMatrices, T *d, T *e) {
#else
template<class T>
void bidiagonalize(T *inputMatrices, T *d, T *e, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols,
				   const int blockSizeBidiagonalize, const int blockSizeExtractDiagonals) {
#endif

    const unsigned int maxIter = (Nrows == Ncols ? Ncols - 1 : Ncols);

#ifdef TEMPLATE
	gpuErrchk(cudaFuncSetCacheConfig(bidiagonalizeKernel<T, numMatrices, Nrows, Ncols>,	cudaFuncCachePreferL1));
	gpuErrchk(cudaFuncSetCacheConfig(extractDiagonals<T, numMatrices, Nrows, Ncols>,	cudaFuncCachePreferL1));
#else
	gpuErrchk(cudaFuncSetCacheConfig(bidiagonalizeKernel<T>,	cudaFuncCachePreferL1));
	gpuErrchk(cudaFuncSetCacheConfig(extractDiagonals<T>,		cudaFuncCachePreferL1));
#endif

	// --- ??? Questo ciclo for si può portare nel kernel?
	//for (unsigned int iterCounter = 0; iterCounter < maxIter ; iterCounter++) {
 //
 //       //bidiagonalizeKernel<<<iDivUp(numMatrices, BLOCKSIZE), BLOCKSIZE>>>(inputMatrices, numMatrices, Nrows, Ncols, iterCounter);
 //       bidiagonalizeKernel<T, numMatrices, Nrows, Ncols><<<iDivUp(numMatrices, BLOCKSIZE), BLOCKSIZE>>>(inputMatrices, iterCounter);
	//	#ifdef DEBUG
	//	    gpuErrchk(cudaPeekAtLastError());
	//	    gpuErrchk(cudaDeviceSynchronize());
 //       #endif
	//}

#ifdef TEMPLATE
	bidiagonalizeKernel<T, numMatrices, Nrows, Ncols><<<iDivUp(numMatrices, blockSizeBidiagonalize), blockSizeBidiagonalize>>>(inputMatrices, maxIter);
#else
	bidiagonalizeKernel<T><<<iDivUp(numMatrices, blockSizeBidiagonalize), blockSizeBidiagonalize>>>(inputMatrices, maxIter, numMatrices, Nrows, Ncols);
#endif
#ifdef DEBUG_BIDIAGONALIZE
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

#ifdef TEMPLATE
	extractDiagonals<T, numMatrices, Nrows, Ncols><<<iDivUp(numMatrices, blockSizeExtractDiagonals), blockSizeExtractDiagonals>>>(inputMatrices, d, e);
#else
	extractDiagonals<<<iDivUp(numMatrices, blockSizeExtractDiagonals), blockSizeExtractDiagonals>>>(inputMatrices, d, e, numMatrices, Nrows, Ncols);
#endif

#ifdef DEBUG_EXTRACT_DIAGONALS
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
#endif
}

/*******************************/
/* MATRIX REARRANGEMENT KERNEL */
/*******************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols>
__global__ void rearrangeKernel(const T * __restrict__ inputMatrices, T * __restrict__ outputMatrices) {
#else
template<class T>
__global__ void rearrangeKernel(const T * __restrict__ inputMatrices, T * __restrict__ outputMatrices, const unsigned int numMatrices,
	                            const unsigned int Nrows, const unsigned int Ncols) {
#endif

    int tid = threadIdx.x + blockDim.x * blockIdx.x;

    if (tid < numMatrices){

#pragma unroll
		for (unsigned int i = 0; i < Nrows; i++)
#pragma unroll
			for (unsigned int j = 0; j < Ncols; j++)
				outputMatrices[tid + j * numMatrices + i * numMatrices * Ncols] = inputMatrices[tid * Nrows * Ncols + j * Nrows + i];

    }
}

/************************/
/* MATRIX REARRANGEMENT */
/************************/
#ifdef TEMPLATE
template<class T, const unsigned int numMatrices, const unsigned int Nrows, const unsigned int Ncols, const int blockSizeRearrange>
void rearrange(const T * __restrict__ inputMatrices, T * __restrict__ outputMatrices) {
#else
template<class T>
void rearrange(const T * __restrict__ inputMatrices, T * __restrict__ outputMatrices, const unsigned int numMatrices, const unsigned int Nrows,
	           const unsigned int Ncols, const int blockSizeRearrange) {
#endif

#ifdef TEMPLATE
	rearrangeKernel<T, numMatrices, Nrows, Ncols><<<iDivUp(numMatrices, blockSizeRearrange), blockSizeRearrange>>>(inputMatrices, outputMatrices);
#else
	rearrangeKernel<T><<<iDivUp(numMatrices, blockSizeRearrange), blockSizeRearrange>>>(inputMatrices, outputMatrices, numMatrices, Nrows, Ncols);
#endif
#ifdef DEBUG_REARRANGE
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
#endif
}

#endif
