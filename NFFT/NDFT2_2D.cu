#include "BBComplex.h"
#include "Utilities.cuh"
#include "InputOutput.cuh"

#define BLOCKSIZE_NUDFT2_2D_X	16
#define BLOCKSIZE_NUDFT2_2D_Y	16

#define DEBUG

#define pi 3.141592653589793238463

/*************************/
/* KERNEL MATRIX FILLING */
/*************************/
__global__ void Kernel_Matrix_Filling(const double * __restrict__ d_X, const double * __restrict__ d_Y, const double * __restrict__ d_u, 
									  const double * __restrict__ d_v, double2_ * __restrict__ d_Kernel_Matrix, const int Nu, const int Nv, 
									  const int M, const int N)
{
	const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int tidy = blockIdx.y * blockDim.y + threadIdx.y;
    
	// --- Evaluates the matrix filling index
	const int tid = tidy * M + tidx;
    
 	const double2_		im_unit(0., 1.);				// Immaginary unit
	
	if (tidx < M && tidy < N) {
		d_Kernel_Matrix[tid] = exp(-2. * pi * im_unit * ((d_u[tidy] * d_X[tidx]) / static_cast<double>(Nu) 
			                                           + (d_v[tidy] * d_Y[tidx]) / static_cast<double>(Nv))); 
	}

}

/************/
/* NDFT2 2D */
/************/
extern "C" void NDFT2_2D_GPU(cublasHandle_t handle, const double * __restrict__ d_X, const double * __restrict__ d_Y, const double * __restrict__ d_u, 
							 const double * __restrict__ d_v, double2_ * __restrict__ d_in, double2_ * __restrict__ d_out, 
							 const int Nu, const int Nv, const int M, const int N) {

	// --- N:		length of d_u and d_v
	// --- M:		length of d_X and d_Y

	double2_ *d_Kernel_Matrix;		gpuErrchk(cudaMalloc(&d_Kernel_Matrix, M * N * sizeof(double2_)));
								 
	// --- Filling the kernel matrix   
	dim3 dimBlock(BLOCKSIZE_NUDFT2_2D_X, BLOCKSIZE_NUDFT2_2D_Y);
	dim3 dimGrid(iDivUp(M, BLOCKSIZE_NUDFT2_2D_X), iDivUp(N, BLOCKSIZE_NUDFT2_2D_Y));
   
	Kernel_Matrix_Filling <<<dimGrid, dimBlock>>>(d_X, d_Y, d_u, d_v, d_Kernel_Matrix, Nu, Nv, M, N);
#ifdef DEBUG
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());
#endif

	saveGPUcomplextxt(d_Kernel_Matrix,	"C:\\Users\\angelo\\Documents\\CEM\\ParticleSwarm\\ParticleSwarmSynthesis\\ParticleSwarmSynthesisMatlab\\Far_Field_Kernel.txt", M * N);

	// --- Matrix multiplication
	double2 alpha;	alpha.x = 1.; alpha.y = 0.;
	double2 beta;	beta.x  = 0.; beta.y  = 0.;
   
	//cublasSafeCall(cublasZgemv('n', 'n',dimu,1,dimx,alpha,d_Kernel_Matrix,dimu,d_eccitazioni,dimx,beta,d_Campo_scalato,dimu));	

	double2 *d_Kernel_Matrix_pointer = reinterpret_cast<double2 *>(d_Kernel_Matrix);
	double2 *d_in_pointer			 = reinterpret_cast<double2 *>(d_in);

	cublasSafeCall(cublasZgemv(handle, CUBLAS_OP_T, M, N, &alpha, d_Kernel_Matrix_pointer, M, d_in_pointer, 1, &beta, reinterpret_cast<double2 *>(d_out), 1));

	// --- Freeing device memory
	gpuErrchk(cudaFree(d_Kernel_Matrix));

}
