#ifndef __SVD_CUH__
#define __SVD_CUH__

#include <time.h>

//#include "my_utils.h"
//#include "bidiag.cuh"
//#include "sturm.cuh"

#include "Utilities.cuh"

//#define DEBUG_SVD_D

#define NUM_THREADS_PER_BLOCK			128  
#define NUM_THREADS_PER_BLOCK_STURM		32

/**********************/
/* SVD PLAN STRUCTURE */
/**********************/
template<class T>
struct svd_plan{
	T* dev_mat_opt;
    T* dev_mat;
    T* dev_mat_copy; 		 
	T* dev_diag;		// --- Diagonal of bidiagonal		 
	T* dev_supdiag;		// --- Off-diagonal of bidiagonal     
	T* dev_diag_redux;   
	T *alpha;			// --- Diagonal of tridiagonal
	T *beta;			// --- Off-diagonal of tridiagonal
};

/*********************/
/* SVD PLAN CREATION */
/*********************/
//template<class T> void create_plan(svd_plan<T> &, unsigned int, unsigned int, unsigned int);
template<class T>
void create_plan(svd_plan<T> &plan, unsigned int Nrows, unsigned int Ncols, unsigned int batch_size){

     // --- Device allocations
     gpuErrchk(cudaMalloc(&(plan.dev_mat_opt),		Nrows * Ncols       * batch_size * sizeof(T)));
	 gpuErrchk(cudaMalloc(&(plan.dev_mat),			Nrows * Ncols       * batch_size * sizeof(T)));
	 gpuErrchk(cudaMalloc(&(plan.dev_mat_copy),		Nrows * Ncols       * batch_size * sizeof(T)));
 	 gpuErrchk(cudaMalloc(&(plan.dev_diag),		            Ncols       * batch_size * sizeof(T)));
     gpuErrchk(cudaMalloc(&(plan.dev_supdiag),		        (Ncols - 1) * batch_size * sizeof(T)));
     gpuErrchk(cudaMalloc(&(plan.dev_diag_redux),		    (Ncols / 2) * batch_size * sizeof(T)));
     gpuErrchk(cudaMalloc(&(plan.alpha),				    Ncols       * batch_size * sizeof(T)));
	 gpuErrchk(cudaMalloc(&(plan.beta),				        (Ncols-1)	* batch_size * sizeof(T)));
	 #ifdef DEBUG_SVD_D
         gpuErrchk(cudaMalloc(&(plan.dev_mat_opt),		Nrows * Ncols       * batch_size * sizeof(T)));
	     gpuErrchk(cudaMalloc(&(plan.dev_mat),			Nrows * Ncols       * batch_size * sizeof(T)));
     	 gpuErrchk(cudaMalloc(&(plan.dev_diag),			        Ncols       * batch_size * sizeof(T)));
         gpuErrchk(cudaMalloc(&(plan.dev_supdiag),		        (Ncols - 1) * batch_size * sizeof(T)));
         gpuErrchk(cudaMalloc(&(plan.dev_diag_redux),	        (Ncols / 2) * batch_size * sizeof(T)));
         gpuErrchk(cudaMalloc(&(plan.alpha),				    Ncols	    * batch_size * sizeof(T)));
	     gpuErrchk(cudaMalloc(&(plan.beta),					    (Ncols - 1) * batch_size * sizeof(T)));
	#endif     	 

}
/************************/
/* SVD PLAN DESTRUCTION */
/************************/
template<class real_type>
void destroy_plan(svd_plan<real_type>& plan){
    
     gpuErrchk(cudaFree(plan.dev_mat_opt));
	 gpuErrchk(cudaFree(plan.dev_mat));
	 gpuErrchk(cudaFree(plan.dev_mat_copy));
 	 gpuErrchk(cudaFree(plan.dev_diag));
     gpuErrchk(cudaFree(plan.dev_supdiag));
     gpuErrchk(cudaFree(plan.dev_diag_redux));
     gpuErrchk(cudaFree(plan.beta));
	 gpuErrchk(cudaFree(plan.alpha));	
}

/*******************************/
/* DEVICE SIDE SVD CALCULATION */
/*******************************/
void deviceSideSVD(thrust::device_vector<double> &, thrust::device_vector<double> &, svd_plan<double> &, unsigned int,  unsigned int, unsigned int);

/*********************************************/
/* COMPUTE BIDIAGONALIZATION KERNEL FUNCTION */
/*********************************************/
template<class real_type>
__global__ void transform_matrix(real_type *in_mat, real_type *out_mat, unsigned int batchSize, unsigned int Nrows, unsigned int Ncols){

    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    
    if (tid < batchSize)
        for (unsigned int num = 0; num < Nrows; num++)
            for (unsigned int j = 0; j < Ncols; j++) out_mat[tid + j*batchSize + num*batchSize*Ncols] = in_mat[tid*Nrows*Ncols + j*Nrows + num];
}

/************************/
/* MATRIX REARRANGEMENT */
/************************/
template<class real_type>
void transform(real_type *in_mat, real_type *out_mat, unsigned int Nrows, unsigned int Ncols, unsigned int batchSize) {
    
    unsigned int NUM_BLOCKS = iDivUp(batchSize, NUM_THREADS_PER_BLOCK);
    transform_matrix<<<NUM_BLOCKS, NUM_THREADS_PER_BLOCK>>>(in_mat, out_mat, batchSize, Nrows, Ncols);
    
}

/*********************************************/
/* COMPUTE BIDIAGONALIZATION KERNEL FUNCTION */
/*********************************************/
//template<class real_type>
//__global__ void build_bidiag(real_type *in_mat, unsigned int batchSize, unsigned int Nrows, unsigned int Ncols, unsigned int step_index) {
//							   
//    int tid_x = threadIdx.x + blockDim.x * blockIdx.x;
//   
//	real_type x0, norm2 = 0., mu, v0;	// --- Local parameters
//	real_type out[20];					// --- Householder vector
//    real_type y[20];
//    real_type beta; 
//
//    if(tid_x < batchSize) {
//    
//        for (unsigned int i = 0; i < step_index; i++) out[i] = 0.f;
//
//        for (unsigned int i = step_index + 1; i < Nrows; i++) {
//            const real_type buffer = in_mat[tid_x + step_index * batchSize + i * batchSize * Ncols]; 
//	        out[i] = buffer;
//	        // --- Evaluate the norm^2 of input vector x
//	        norm2 = norm2 + buffer * buffer;
//	    }
//	        
//        x0 = in_mat[tid_x + (step_index * batchSize) * (1 + Ncols)];
//
//        if(norm2 == 0.) beta = 0.;
//        else {
//	        mu = sqrt((x0 * x0) + norm2);
//
//	        if (x0 <= 0) v0 = x0 - mu;
//	        else v0 = -norm2 / (x0 + mu);
//
//	        real_type temp = 1. / v0;
//	        beta = 2. * ((v0 * v0) / (norm2 + (v0 * v0)));
//	        out[step_index] = 1.;
//	        for (unsigned int s = step_index + 1; s < Nrows; s++) out[s] *= temp;
//        }
//        
//		{
//			real_type sum = 0.f;
//			for (unsigned int j = 0; j < Ncols; j++) {
//				sum = 0.f; 
//				for (unsigned int i = 0; i < Nrows; i++)  sum += out[i] * in_mat[tid_x + i * batchSize * Ncols + j * batchSize];
//				y[j] = sum;
//			}
//		}   
//
//        // --- Matrix update
//        for (unsigned int r = 0; r < Nrows; r++)
//            for (unsigned int c = 0; c < Ncols; c++) in_mat[tid_x + r * batchSize * Ncols + c * batchSize] = 
//				 -beta * y[c] * out[r] + in_mat[tid_x + r * batchSize * Ncols + c * batchSize];
//      
//	}							   
//	
//	__syncthreads(); // --- Waiting for the matrix update 
//    
//    if (step_index < Ncols - 2) {
//    
//	    // --- Build the right Householder vector
//        if (tid_x < batchSize) {
//
//		    // --- Initializing norm2
//		    norm2 = 0.;
//                    
//            for (unsigned int i = 0; i < step_index + 1; i++) out[i] = 0.f;
//			
//		    for (unsigned int i = step_index + 2; i < Ncols; i++) {
//		        const real_type buffer = in_mat[tid_x + i * batchSize + step_index * batchSize * Ncols];
//			    out[i] = buffer;
//			    // --- Evaluate the norm^2 of input vector x
//			    norm2 = norm2 + buffer * buffer;
//		    }
//		
//		    x0 = in_mat[tid_x + (step_index + 1) * batchSize + step_index * batchSize * Ncols];
//
//		    if (norm2 == 0.) beta = 0.;
//		    else
//		    {
//			    mu = sqrt((x0 * x0) + norm2);
//
//			    if (x0 <= 0) v0 = x0 - mu;
//			    else v0 = -norm2 / (x0 + mu);
//
//			    real_type temp = 1. / v0;
//			    beta = 2. * ((v0 * v0) / (norm2 + (v0 * v0)));
//			    out[step_index + 1] = 1.;
//			    for(unsigned int s = step_index + 2; s < Ncols; s++) out[s] *= temp;
//		    }
//
//			// -- Evaluate the product A * (right Householder vector)
//			{
//				real_type sum = 0.f;
//				for (unsigned int i = 0; i < Nrows; i++) {
//					sum = 0.f; 
//					for (unsigned int j = 0; j < Ncols; j++) sum += out[j] * in_mat[tid_x + i * batchSize * Ncols + j * batchSize];
//					y[i] = sum;
//				}
//			}
//
//			{
//				// --- Matrix update
//				for (unsigned int r = 0; r < Nrows; r++)
//					for(unsigned int c = 0; c < Ncols; c++) in_mat[tid_x + r*batchSize*Ncols + c*batchSize] = - beta*y[r]*out[c] + in_mat[tid_x + r*batchSize*Ncols + c*batchSize];
//			}   
//	     
//		}
//	}
//}

/*********************************/
/* BUILD LEFT HOUSEHOLDER VECTOR */
/*********************************/
template<class real_type>
__device__ void build_bidiag_left(const real_type *in_mat, real_type *out, real_type &betaq, unsigned int tid_x, unsigned int batchSize,
                                unsigned int Nrows, unsigned int Ncols, unsigned int step_index) {

	real_type x0, norm2 = 0.f, mu, v0;
    for (unsigned int i = 0; i < step_index; i++) out[i] = 0.f;

	for (unsigned int i = step_index + 1; i < Nrows ; i++) {
		const real_type buffer = in_mat[tid_x + step_index * batchSize + i * batchSize * Ncols]; 
		out[i] = buffer;
		// --- Evaluate the norm^2 of input vector x
		norm2 = norm2 + buffer * buffer;
	}
                
	x0 = in_mat[tid_x + (step_index * batchSize) * (1 + Ncols)];

	if (norm2 == 0.) betaq = 0.;
	else {
		mu = sqrt( (x0 * x0) + norm2 );
		
		if (x0 <= 0) v0 = x0 - mu;
		else v0 = -norm2/(x0 + mu);

        real_type temp = 1. / v0;
        betaq = 2. * ((v0 * v0) / (norm2 + (v0 * v0)));
        out[step_index] = 1.;
        for(unsigned int s = step_index + 1; s < Nrows ; s++) out[s] *= temp;
	}
}

/**********************************/
/* BUILD RIGHT HOUSEHOLDER VECTOR */
/**********************************/
template<class real_type>
__device__ void build_bidiag_right(const real_type *pseudo_v, real_type *v, real_type &betap, unsigned int tid_x, unsigned int batchSize,
								   unsigned int Nrows, unsigned int Ncols, unsigned int step_index){

	real_type x0, norm2 = 0.f, mu, v0;
    for (unsigned int i = 0; i < step_index + 1; i++) v[i] = 0.f;      
	            
	for (unsigned int i = step_index + 2; i < Ncols ; i++) {
		const real_type buffer = pseudo_v[i];
    	v[i] = buffer;
    	// --- Evaluate the norm^2 of input vector x
    	norm2 = norm2 + buffer * buffer;
	}	
             
    x0 = pseudo_v[step_index+1];
         
    if (norm2 == 0.) betap = 0.;
	else {
		mu = sqrt((x0 * x0) + norm2);

		if (x0 <= 0) v0 = x0 - mu;
    	else v0 = -norm2 / (x0 + mu);
         
		real_type temp = 1. / v0;
		betap = 2. * ((v0 * v0) / (norm2 + (v0 * v0)));
		v[step_index + 1] = 1.;
		for (unsigned int s = step_index + 2; s < Ncols; s++) v[s] *= temp;
	}
}

/*****************************************************************/
/* COMPUTE BIDIAGONALIZATION KERNEL FUNCTION - OPTIMIZED VERSION */
/*****************************************************************/
template<class real_type>
__global__ void build_bidiag_opt(real_type *in_mat, unsigned int batchSize, unsigned int Nrows, unsigned int Ncols, unsigned int step_index) {
							   
    int tid_x = threadIdx.x + blockDim.x * blockIdx.x;
	
	real_type out[20];		// --- Householder vector
    real_type pseudo_v[20];
    real_type v[20];
    real_type x1[20];
    real_type w[20];
    real_type z[20];
    real_type betaq = 0.f, betap = 0.f; 
    
    if (step_index < Ncols - 2) {

       if (tid_x < batchSize) {
        
            // --- Build the left Householder vector
            build_bidiag_left(in_mat, out, betaq, tid_x, batchSize, Nrows, Ncols, step_index);
           
			// --- Build x1 and pseudo_v
			{
				real_type sum = 0.f;
				for (unsigned int j = 0; j < Ncols; j++) {
					sum = 0.f; 
					for (unsigned int i = 0; i < Nrows; i++) sum += -betaq*out[i]*in_mat[tid_x + i*batchSize*Ncols + j*batchSize];
					pseudo_v[j] = sum + in_mat[tid_x + j * batchSize + step_index * batchSize * Ncols];
					x1[j] = -sum;
				}
			}

			// --- Build the right Householder vector  
			build_bidiag_right(pseudo_v, v, betap, tid_x, batchSize, Nrows, Ncols, step_index);
       
			// --- Build w 
			{
				real_type sum = 0.f;
				for (unsigned int i = 0; i < Nrows; i++) {
					sum = 0.f; 
					for (unsigned int j = 0; j < Ncols; j++) sum += betap*v[j]*in_mat[tid_x + i*batchSize*Ncols + j*batchSize];
					w[i] = sum;
				}
			} 
         
			{
				real_type temp = 0.f;
				real_type sum = 0.f;
				for (unsigned int i = 0; i < Ncols; i++) sum += x1[i]*v[i];
				temp = sum;   
				for (unsigned int i = 0; i < Ncols; i++) z[i] = x1[i] - (betap*temp*v[i]);
			}

        
			// --- Matrix update
			for (unsigned int r = 0; r < Nrows; r++) 
				for (unsigned int c = 0; c < Ncols; c++) in_mat[tid_x + r*batchSize*Ncols + c*batchSize] = - out[r]*z[c] - w[r]*v[c] + in_mat[tid_x + r*batchSize*Ncols + c*batchSize];
       }							   
	     
	}
    
    if (step_index >= (Ncols - 2)) {
    
        if (tid_x < batchSize) {
            
			build_bidiag_left(in_mat, out, betaq, tid_x, batchSize, Nrows, Ncols, step_index);
            
            real_type sum = 0.f;
            for (unsigned int j = 0; j < Ncols; j++) {
                sum = 0.f; 
                for (unsigned int i = 0; i < Nrows; i++) sum += out[i]*in_mat[tid_x + i*batchSize*Ncols + j*batchSize];
                x1[j] = sum;
            }   
            
            // --- Matrix update
            for (unsigned int r = 0; r < Nrows; r++)
				for(unsigned int c = 0; c < Ncols; c++) 
					in_mat[tid_x + r*batchSize*Ncols + c*batchSize] = - betaq*x1[c]*out[r] + in_mat[tid_x + r*batchSize*Ncols + c*batchSize];
            
        }
    }
}
 
/**************************************/
/* SELECT DIAGONAL AND UPPER DIAGONAL */ 
/**************************************/
template<class real_type>
__global__ void select_bidiag(real_type *mat, real_type *d, real_type *e, unsigned int Nrows, unsigned int Ncols, unsigned int batchSize) {

	int tid = threadIdx.x + blockDim.x * blockIdx.x;

	if (tid < batchSize){
		
		for(unsigned int i = 0; i < Ncols; i++) d[tid + i * batchSize] = mat[tid + i * batchSize * (Ncols + 1)];
		for(unsigned int j = 0; j < Ncols - 1; j++) e[tid + j * batchSize] = mat[tid + batchSize + j * batchSize * (Ncols + 1)];
	}
}

/*****************************/
/* COMPUTE BIDIAGONALIZATION */
/*****************************/
template<class real_type>
void bidiag(real_type *mat, real_type *d, real_type *e, unsigned int Nrows, unsigned int Ncols, unsigned int batchSize) {

    unsigned int NUM_BLOCKS = iDivUp(batchSize, NUM_THREADS_PER_BLOCK);
 
    unsigned int max_iter = (Nrows == Ncols ? Ncols - 1 : Ncols);
	unsigned int step_index = 0;

	for(step_index = 0; step_index < max_iter; step_index++) { 
        
        build_bidiag_opt<<<NUM_BLOCKS,NUM_THREADS_PER_BLOCK>>>(mat, batchSize, Nrows, Ncols, step_index);
        #ifdef DEBUG_SVD_D
		    gpuErrchk( cudaPeekAtLastError() );
		    gpuErrchk( cudaDeviceSynchronize() );
        #endif
	}

	select_bidiag<<<NUM_BLOCKS,NUM_THREADS_PER_BLOCK>>>(mat, d, e, Nrows, Ncols, batchSize);
	#ifdef DEBUG
	    gpuErrchk( cudaPeekAtLastError() );
	    gpuErrchk( cudaDeviceSynchronize() );
	#endif    
}

/*************************************************/
/* BUILD TRIDIAGONAL MATRIX FROM DIAGONAL MATRIX */
/*************************************************/
template<class real_type>
__host__ __device__ void build_tri(const real_type *d, const real_type *e, real_type *alpha, real_type *beta, unsigned int length, 
								  unsigned int batchSize){
    
    // --- Build alpha
    for (unsigned int i = 0; i < length; i++) {
        if (i == 0) alpha[i * batchSize] = d[i * batchSize] * d[i * batchSize];
        else alpha[i*batchSize] = (e[(i - 1) * batchSize] * e[(i - 1) * batchSize]) + d[i * batchSize] * d[i * batchSize];
    }
    
    // --- Build beta
    for (unsigned int i = 0; i < length - 1; i++) beta[i*batchSize] = d[i*batchSize]*e[i*batchSize];        
}

/********************************************/
/* BUILD TRIDIAGONAL MATRIX FROM BIDIAGONAL */
/********************************************/
template<class real_type>
__host__ __device__ void build_pivmin(const real_type *b, real_type *pivmin, unsigned int cols, unsigned int batchSize) {
    
    const unsigned int length = cols - 1;
    real_type max_temp = b[0];
    for (unsigned int i = 1; i < length; i++) {
        real_type temp = b[i * batchSize];
        if (temp > max_temp) max_temp = temp;
    } 
    
    *pivmin = (2 * fmax(1, max_temp)) / DBL_MAX;       
}

/*******************************************************************/
/* CHCKSIGN COMPUTES THE NUMBER OF SIGN CHANGES IN STURM SEQUENCES */
/*******************************************************************/
template<class real_type>
__host__ __device__ void chcksign_pivmin(const real_type *d, const real_type *b, const real_type pivmin, real_type x, unsigned int &nch, 
										 unsigned int length, unsigned int batchSize){
    
    unsigned int n = length;
    nch = 0;
    
    real_type q = d[0] - x;
    if(fabs(q) <= pivmin) q = -pivmin;
    if(q < 0.0f) nch = nch + 1;
    
    for (unsigned int i = 2; i <= n; i++) {
        // -- This order of operation preserve monotonicity 
		q = (d[(i-1) * batchSize]  - ((b[(i - 2) * batchSize] * b[(i - 2) * batchSize]) / q) ) - x; 
        if (fabs(q) <= pivmin) q = -pivmin;
        if (q < 0.0f) nch = nch + 1;        
    }
}

/*******************************************************************************************/
/* USE GERSHGORING DISK TO COMPUTE THE INTERVAL WITHIN WHICH THE EIGENVALUES ARE CONTAINED */
/*******************************************************************************************/
template<class real_type>
__host__ __device__ void bound(const real_type *d, const real_type *b, real_type &alpha, real_type &beta, unsigned int length, unsigned int batchSize) {
    
    real_type temp = 0;
    unsigned int n = length;
    
    alpha = d[0] - fabs(b[0]);
    temp = d[(n - 1) * batchSize] - fabs(b[(n - 2) * batchSize]);
    if (temp < alpha) alpha = temp;
    for (unsigned int i = 1; i < n - 1; i++) {
        temp = d[i * batchSize] - fabs(b[(i - 1) * batchSize]) - fabs(b[i * batchSize]);
        if (temp < alpha) alpha = temp;
    }  
    beta = fabs(d[0]) + fabs(b[0]);
    temp = fabs(d[(n - 1) * batchSize]) + fabs(b[(n - 2) * batchSize]);
    if (temp > beta) beta = temp;
    for (unsigned int i = 1; i < n - 1; i++) {
        temp = d[i * batchSize] + fabs(b[(i - 1) * batchSize]) + fabs(b[i * batchSize]);
        if (temp > beta) beta = temp;
    }          
}

/******************************************/
/* GIVENS METHOD BASED ON STURM SEQUENCES */
/******************************************/
template<class real_type>
__global__ void givsturm_parallel_pivmin(real_type *d, real_type *b, real_type *dd, real_type *bb, real_type *eig, const unsigned int Ncols, 
										 const unsigned int batchSize) {

    int tid_x = threadIdx.x + blockDim.x * blockIdx.x;
    int tid_y = threadIdx.y;
    
    real_type tol = 0.0000000000001;
    
    __shared__ real_type alpha_shared[NUM_THREADS_PER_BLOCK_STURM];
    __shared__ real_type beta_shared[NUM_THREADS_PER_BLOCK_STURM];
    __shared__ real_type pivmin_shared[NUM_THREADS_PER_BLOCK_STURM];
    
    const unsigned int n = Ncols;
    
    if (tid_x < batchSize && tid_y == 0) {
       build_tri(d + tid_x, b + tid_x, dd + tid_x, bb + tid_x, n, batchSize);
       build_pivmin(bb + tid_x, &pivmin_shared[threadIdx.x], Ncols, batchSize);
       bound(dd + tid_x, bb + tid_x, alpha_shared[threadIdx.x], beta_shared[threadIdx.x], Ncols, batchSize);
    }
    
    __syncthreads();
    
    if(tid_x < batchSize && tid_y < Ncols) {

        unsigned int eig_index = tid_y % n;
           
        real_type c = 0;
        unsigned int nch = 0;
                     
        real_type alpha, beta, pivmin;
        
        alpha	= alpha_shared[threadIdx.x];
        beta	= beta_shared[threadIdx.x];
        pivmin	= pivmin_shared[threadIdx.x];
        
        real_type dist = fabs(beta - alpha);		// --- Initial amplitude of the search interval along x, different from thread to thread
        real_type s = fabs(beta) + fabs(alpha);

		while(dist > tol * s) { // --- Quarternoni's criterion

			c = (alpha + beta)/2;

			chcksign_pivmin<real_type>(dd + tid_x, bb + tid_x, pivmin, c, nch, n, batchSize); 
                        
            if (nch > (n - ((eig_index) + 1))) beta = c;
            else alpha = c;
            
            dist = fabs(beta - alpha);
            s = fabs(beta) + fabs(alpha);  
                      
        }  
        
        eig[eig_index + tid_x*n] = sqrt(fabs(c));
    }
}

/********************/
/* GIVSTURM WRAPPER */
/********************/
template<class real_type>
void givsturm_par(real_type *d, real_type *b, real_type *dd, real_type *bb, real_type *eig, unsigned int cols, const unsigned int batchSize) {
                                          
    const unsigned int NUM_THREADS_Y = cols;
    const unsigned int NUM_BLOCKS_STURM = iDivUp(batchSize, NUM_THREADS_PER_BLOCK_STURM);
    dim3 BLOCKGRID(NUM_THREADS_PER_BLOCK_STURM, NUM_THREADS_Y);
                
    givsturm_parallel_pivmin<<<NUM_BLOCKS_STURM, BLOCKGRID>>>(d, b, dd, bb, eig, cols, batchSize);
                                       
}

/*****************************/
/* SINGULAR VALUES SELECTION */
/*****************************/
template<class real_type> 
__host__ __device__ void select_diag_test(real_type *input, real_type *output, unsigned int Ncols) {

	for (int i = 0; i < Ncols / 2; i++) output[i] = input[2 * i];
}

template<class real_type>
__global__ void select_sing(real_type *input, real_type *output, unsigned int Ncols, unsigned int batchSize) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < batchSize) select_diag_test(input + Ncols * tid, output + (Ncols / 2) * tid, Ncols);
}

/************************************/
/* DEVICE SIDE SVD COMPUTATION CORE */
/************************************/
template<class real_type>
void my_svd(svd_plan<real_type>& plan, real_type *mat, real_type *sing, unsigned int Nrows, unsigned int Ncols, unsigned int batchSize) {

    // -- Reorganize the input matrix
	transform(mat, plan.dev_mat_opt, Nrows, Ncols, batchSize);
    #ifdef DEBUG_SVD_D
        gpuErrchk(cudaPeekAtLastError());
        gpuErrchk(cudaDeviceSynchronize());
    #endif	
    gpuErrchk(cudaDeviceSynchronize());
    
    // -- Compute bidiagonalization
	bidiag<real_type>(plan.dev_mat_opt, plan.dev_diag, plan.dev_supdiag, Nrows, Ncols, batchSize);
    #ifdef DEBUG_SVD_D
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
    #endif	
	gpuErrchk(cudaDeviceSynchronize());
		
    givsturm_par<real_type>(plan.dev_diag, plan.dev_supdiag, plan.alpha, plan.beta, sing, Ncols, batchSize);
    #ifdef DEBUG_SVD_D
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
    #endif 
    gpuErrchk(cudaDeviceSynchronize());   

    unsigned int NUM_BLOCKS = iDivUp(batchSize, NUM_THREADS_PER_BLOCK);
	// --- Takes every two singular values
	select_sing<<<NUM_BLOCKS, NUM_THREADS_PER_BLOCK>>>(sing, plan.dev_diag_redux, Ncols, batchSize); 
	gpuErrchk(cudaMemcpy(sing, plan.dev_diag_redux, (Ncols / 2) * batchSize * sizeof(real_type), cudaMemcpyDeviceToDevice));
}

void deviceSideSVD(thrust::device_vector<double> &A, thrust::device_vector<double> &S, svd_plan<double> &plan, unsigned int Nrows, 
				   unsigned int Ncols, unsigned int batchSize) {

	my_svd(plan, thrust::raw_pointer_cast(A.data()), thrust::raw_pointer_cast(S.data()), Nrows, Ncols, batchSize);
}

//void deviceSideSVD(double *A, double *S, svd_plan<double> &plan, unsigned int Nrows, 
//				   unsigned int Ncols, unsigned int batchSize) {
//
//	my_svd(plan, A, S, Nrows, Ncols, batchSize);
//}

#endif
