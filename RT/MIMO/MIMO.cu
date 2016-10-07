#include <math.h>  

#include <cuda.h>
#include <cuda_runtime.h>

#include <thrust\iterator\discard_iterator.h>

#include "MIMO.cuh"
#include "Utilities.cuh"

/************************************************/
/* EVALUATE THE FROBENIUS NORM OF A REAL VECTOR */
/************************************************/
double frobeniusNorm(const thrust::host_vector<double>::iterator start, int size) {

    double tempSum = 0.;
    
    for (int index = 0; index < size; index++) tempSum += fabs(start[index]) * fabs(start[index]);

    return sqrt(tempSum);        
}

/*************************************************/
/* NORMALIZE A REAL VECTOR BY ITS FROBENIUS NORM */
/*************************************************/
void normalizeVector(thrust::host_vector<double> &matrix, int num_TX, int num_RX) {

    double frobNorm		= 0.;
    double scale_factor = 0.;
    
    frobNorm		= frobeniusNorm(matrix.begin(), num_TX * num_RX);
    scale_factor	= (frobNorm / sqrt((double)(num_TX * num_RX)));

	if (frobNorm != 0) for (int k = 0; k < num_TX; k++) for (int j = 0; j < num_RX; j++) matrix[j + k * num_TX] = (matrix[j + k * num_TX] / scale_factor);
 
}

/***************************/
/* INTEGER DIVISION STRUCT */
/***************************/
struct idivision{
    
	unsigned int num_TX;
    
    idivision(unsigned int NTX_) : num_TX(NTX_) {}
    
    __host__ __device__ unsigned int operator()(unsigned int i) const { return i / num_TX; }
};

/***********************************/
/* COMPUTE CHANNEL CAPACITY STRUCT */
/***********************************/
struct compute_capacity : public thrust::unary_function<double, double>
{
    double SNR_linear_over_NTX;

    compute_capacity(double SNR_linear_over_NTX_) : SNR_linear_over_NTX(SNR_linear_over_NTX_) {}
        
	__host__ __device__ double operator()(double x) const { return log2(1 + (SNR_linear_over_NTX) * (x * x) ); }
};

/************************************************/
/* COMPUTE CHANNEL CAPACITY - HOST SIDE VERSION */
/************************************************/
void computeCapacityHost(const thrust::host_vector<double> &S, thrust::host_vector<double> &C_h, double SNR_linear, unsigned int num_TX) {
    
	thrust::host_vector<unsigned int> keys(S.size());
    
    thrust::transform(thrust::counting_iterator<unsigned int>(0), thrust::counting_iterator<unsigned int>(0) + S.size(),
                      keys.begin(), idivision(num_TX));    
    
	thrust::reduce_by_key(keys.begin(), keys.end(), thrust::make_transform_iterator(S.begin(), compute_capacity(SNR_linear / num_TX)),                         
                         thrust::make_discard_iterator(), C_h.begin());
}

/**************************************************/
/* COMPUTE CHANNEL CAPACITY - DEVICE SIDE VERSION */
/**************************************************/
void computeCapacityDevice(const thrust::device_vector<double> &S, thrust::device_vector<double> &C_d, double SNR_linear, unsigned int num_TX) {
    
	thrust::device_vector<unsigned int> keys(S.size());
	
    thrust::transform(thrust::counting_iterator<unsigned int>(0), thrust::counting_iterator<unsigned int>(0) + S.size(),
                      keys.begin(), idivision(num_TX));    
    
	thrust::reduce_by_key(keys.begin(), keys.end(), thrust::make_transform_iterator(S.begin(), compute_capacity(SNR_linear / num_TX)),                         
                          thrust::make_discard_iterator(), C_d.begin());
}

/*************************************************************/
/* COMPUTE MAXIMUM CHANNEL CAPACITY PER GROUP - HOST VERSION */
/*************************************************************/
void findMaxChannelCapacityPerGroupHost(const thrust::host_vector<double> &C_h, thrust::host_vector<double> &Cg_h, unsigned int Nconf) {
    
	thrust::host_vector<unsigned int> keys(C_h.size());
    thrust::transform(thrust::counting_iterator<unsigned int>(0), thrust::counting_iterator<unsigned int>(0) + C_h.size(),
                      keys.begin(), idivision(Nconf));  
    
    thrust::equal_to<unsigned int> binary_pred;
    thrust::maximum<double> binary_op;
    
    thrust::reduce_by_key(keys.begin(), keys.end(), C_h.begin(), thrust::make_discard_iterator(), Cg_h.begin(), binary_pred, binary_op);
                     
}

/***************************************************************/
/* COMPUTE MAXIMUM CHANNEL CAPACITY PER GROUP - DEVICE VERSION */
/***************************************************************/
void findMaxChannelCapacityPerGroupDevice(const thrust::device_vector<double> &C_d, thrust::device_vector<double> &Cg_d, unsigned int Nconf) {
    
	thrust::device_vector<unsigned int> keys(C_d.size());
    thrust::transform(thrust::counting_iterator<unsigned int>(0), thrust::counting_iterator<unsigned int>(0) + C_d.size(),
                      keys.begin(), idivision(Nconf));  
    
    thrust::equal_to<unsigned int> binary_pred;
    thrust::maximum<double> binary_op;
    
    thrust::reduce_by_key(keys.begin(), keys.end(), C_d.begin(), thrust::make_discard_iterator(), Cg_d.begin(), binary_pred, binary_op);
                          
}
