#ifndef __MIMO_CUH__
#define __MIMO_CUH__

#include <thrust\host_vector.h>
#include <thrust\device_vector.h>

#include "cfloat.h"

void normalizeVector(thrust::host_vector<double> &, int, int);

void computeCapacityHost(const thrust::host_vector<double> &, thrust::host_vector<double> &, double, unsigned int);
void computeCapacityDevice(const thrust::device_vector<double> &, thrust::device_vector<double> &, double, unsigned int);

void findMaxChannelCapacityPerGroupHost  (const thrust::host_vector  <double> &, thrust::host_vector  <double> &, unsigned int);
void findMaxChannelCapacityPerGroupDevice(const thrust::device_vector<double> &, thrust::device_vector<double> &, unsigned int);

#endif
