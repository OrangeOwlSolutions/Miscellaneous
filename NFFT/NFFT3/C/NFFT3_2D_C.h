#ifndef __NUFFT_3T_2DC_H__
#define __NUFFT_3T_2DC_H__
#include <complex>
#include "fftw3.h"

void  NUFFT_2D_CPU(const double * __restrict,			// x
	const double * __restrict,			// y
	const double * __restrict,			// s
	const double * __restrict,			// t 
	const fftw_complex * __restrict,		// f
	fftw_complex * __restrict,			// F
	const double,						// eps
	const int,							// len_in
	const int);							// len_out


#endif
