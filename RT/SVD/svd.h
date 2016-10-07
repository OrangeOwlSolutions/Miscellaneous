#ifndef __SVD_H__
#define __SVD_H__

#include <thrust\host_vector.h>

extern "C" void hostSideSVD(const thrust::host_vector<double> &, thrust::host_vector<double> &, unsigned int, unsigned int, unsigned int);

#endif
