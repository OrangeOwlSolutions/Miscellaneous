#include "svd.h"

#include <Eigen/Dense>
#include <iostream>

/*******************************/
/* SVD CALCULATION USING EIGEN */
/*******************************/
void eigen_svd(const double *A, double *S, unsigned int Nrows, unsigned int Ncols, unsigned int batchSize) {
    
    Eigen::MatrixXf Hgj(Nrows, Ncols);
    
    for (unsigned int j = 0; j < batchSize; j++) {
    
        for (unsigned int r = 0; r < Nrows; r++) 
            for (unsigned int c = 0; c < Ncols; c++) Hgj(r, c) = A[j * (Nrows * Ncols) + c * Nrows + r];

		//printf("\n\nMatrix nr %i\n", j);
		//for (size_t r = 0; r < Nrows; r++) {
		//	for (size_t c = 0; c < Ncols; c++) {
		//		std::cout << Hgj(r, c) << " ";  
		//	}
		//	std::cout << std::endl;
		//}
        Eigen::VectorXf singularv = Hgj.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).singularValues(); 

        for (unsigned int i = 0; i < Ncols / 2; i++) S[j * Ncols / 2 +i] = singularv(2 * i);
    }
}

/*****************************/
/* HOST SIDE SVD CALCULATION */
/*****************************/
extern "C" void hostSideSVD(const thrust::host_vector<double> &A, thrust::host_vector<double> &S, unsigned int Nrows, unsigned int Ncols, unsigned int batchSize) {           
	
	eigen_svd(thrust::raw_pointer_cast(A.data()), thrust::raw_pointer_cast(S.data()), Nrows, Ncols, batchSize);
}
