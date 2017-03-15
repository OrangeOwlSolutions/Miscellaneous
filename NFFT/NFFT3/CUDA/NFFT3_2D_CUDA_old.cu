//#include <thrust\adjacent_difference.h>
//#include <thrust\execution_policy.h>
//#include <thrust\sort.h>
//#include <thrust\gather.h>
//#include <thrust\iterator\constant_iterator.h>
//#include <thrust\binary_search.h>
//
////__global__ void CircshiftKernel(double2 * __restrict  out, //double2  __restrict * out,
////	const double2 * __restrict in,
////	const int xdim, const int ydim,
////	const int xshift, const int yshift)
////{
////
////	int i = threadIdx.x + blockIdx.x * blockDim.x;
////	int j = threadIdx.y + blockIdx.y * blockDim.y;
////
////	if ((i < xdim) && (j < ydim)) {
////
////		int ii = (i + xshift) % xdim;
////		int jj = (j + yshift) % ydim;
////
////		//out[jj * xdim + ii].x = in[j * xdim + i].x;
////		//out[jj * xdim + ii].y = in[j * xdim + i].y;
////		out[jj * xdim + ii].x = -(1. - 2 * ((i + j) & 1)) * in[j * xdim + i].x;
////		out[jj * xdim + ii].y = -(1. - 2 * ((i + j) & 1)) * in[j * xdim + i].y;
////
////	}
////
////}
////
////__global__ void FftshiftKernel(double2 * __restrict data,
////	const int N1,
////	const int N2) {
////
////	int i = threadIdx.x + blockIdx.x * blockDim.x;
////	int j = threadIdx.y + blockIdx.y * blockDim.y;
////
////	if ((i < N1) && (j < N2)) {
////
////		data[j*N1 + i].x *= -(1. - 2 * ((i + j) & 1));
////		data[j*N1 + i].y *= -(1. - 2 * ((i + j) & 1));
////	}
////
////}
//
///************************/
///* MORTON CODE ENCODING */
///************************/
//// "Insert" a 0 bit after each of the 16 low bits of x
//// --- ^ xor
//// --- & and
//__host__ __device__ uint32_t Part1By1(uint32_t x){
//	x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
//	x = (x ^ (x << 8)) & 0x00ff00ff;  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
//	x = (x ^ (x << 4)) & 0x0f0f0f0f;  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
//	x = (x ^ (x << 2)) & 0x33333333;  // x = --fe --dc --ba --98 --76 --54 --32 --10
//	x = (x ^ (x << 1)) & 0x55555555;  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
//	return x;
//}
//
//__host__ __device__ uint32_t encode_morton2d(uint32_t x, uint32_t y)
//{
//	return (Part1By1(y) << 1) + Part1By1(x);
//}
//
///********************************************************/
///* SPATIAL CONVOLUTION KERNEL WITHOUT ATOMIC OPERATIONS */
///********************************************************/
////__host__ __device__ double2 addContributions(const double * __restrict__ d_x, const double * __restrict__ d_y, const double2 * __restrict d_f, double fact, double2 temp_sum, const double xb, const double yb, const double sb,
////	const double tb, const double Dx, const double Dy, const int msp, const double d_n_si, const double d_n_sj, const double t1, const double b, const int num_elements) {
////
////	for (int kk = 0; kk < num_elements; kk++) {
////
////		double x_temp = (d_x[kk] - xb) / Dx;
////		double y_temp = (d_y[kk] - yb) / Dy;
////
////		if ((abs(x_temp - d_n_si) <= msp) && (abs(y_temp - d_n_sj) <= msp)) {
////
////			double arg = sb * d_x[kk] + tb * d_y[kk];
////
////			double real_part_temp = d_f[kk].x * cos(arg) + d_f[kk].y * sin(arg);
////			double imag_part_temp = d_f[kk].y * cos(arg) - d_f[kk].x * sin(arg);
////
////			double exp_temp = fact * exp(b * (Dy * Dy * d_n_sj * d_n_sj + Dx * Dx * d_n_si * d_n_si) - t1 * (((x_temp - d_n_si) * (x_temp - d_n_si) + (y_temp - d_n_sj) * (y_temp - d_n_sj))));
////
////			temp_sum.x = temp_sum.x + real_part_temp * exp_temp;
////			temp_sum.y = temp_sum.y + imag_part_temp * exp_temp;
////
////		}
////
////	}
////
////	return temp_sum;
////
////}
//
//struct complex_sum {
//
//	__device__ double2 operator()(const double2 &a, const double2 &b) {
//
//		double2 result;
//		result.x = a.x + b.x;
//		result.y = a.y + b.y;
//
//		return result;
//	}
//};
//
//struct transf
//{
//	int msp;
//	double t1, d_n_si, d_n_sj;
//
//	__device__ transf(double t1_, double d_n_si_, double d_n_sj_, int msp_) : t1(t1_),
//		d_n_si(d_n_si_), d_n_sj(d_n_sj_), msp(msp_) { }
//
//	__device__ double2 operator()(thrust::tuple<double, double, double2> t)
//	{
//		double x_temp = thrust::get<0>(t);
//		double y_temp = thrust::get<1>(t);
//
//		double exp_temp = exp(-t1 * (((x_temp - d_n_si) * (x_temp - d_n_si) + (y_temp - d_n_sj) * (y_temp - d_n_sj))));
//
//		double2 temp_sum;
//		temp_sum.x = thrust::get<2>(t).x;
//		temp_sum.y = thrust::get<2>(t).y;
//
//		double test = ((abs(x_temp - d_n_si) <= msp) && (abs(y_temp - d_n_sj) <= msp));
//
//		temp_sum.x = test * temp_sum.x * exp_temp;
//		temp_sum.y = test * temp_sum.y * exp_temp;
//
//		return temp_sum;
//
//	}
//
//};
//
//__device__ double2 addContributions(const double * __restrict__ d_x, const double * __restrict__ d_y, const double2 * __restrict d_f, double fact, double2 temp_sum, const double xb, const double yb, const double sb,
//	const double tb, const double Dx, const double Dy, const int msp, const double d_n_si, const double d_n_sj, const double t1, const double b, const int num_elements) {
//
//	//auto begin = thrust::make_zip_iterator(thrust::make_tuple(thrust::device_pointer_cast(d_x), thrust::device_pointer_cast(d_y), thrust::device_pointer_cast(d_f)));
//	//auto end   = thrust::make_zip_iterator(thrust::make_tuple(thrust::device_pointer_cast(d_x), thrust::device_pointer_cast(d_y), thrust::device_pointer_cast(d_f))) + num_elements;
//
//	//double2 initial_result;
//	//initial_result.x = 0.;
//	//initial_result.y = 0.;
//
//	//double2 result = thrust::transform_reduce(thrust::device, begin, end, transf(t1, d_n_si, d_n_sj, msp), temp_sum, complex_sum());
//	//
//	//return result;
//
//	for (int kk = 0; kk < num_elements; kk++) {
//
//		double x_temp = d_x[kk];
//		double y_temp = d_y[kk];
//
//		if ((abs(x_temp - d_n_si) <= msp) && (abs(y_temp - d_n_sj) <= msp)) {
//
//			double exp_temp = exp(-t1 * (((x_temp - d_n_si) * (x_temp - d_n_si) + (y_temp - d_n_sj) * (y_temp - d_n_sj))));
//
//			temp_sum.x = temp_sum.x + d_f[kk].x * exp_temp;
//			temp_sum.y = temp_sum.y + d_f[kk].y * exp_temp;
//
//		}
//
//	}
//
//	return temp_sum;
//
//}
//
//__global__ void spatialConvolutionKernelNoAtomic(double2 * __restrict__ d_f_tau, const double2 * __restrict__ d_f, const double * __restrict__ d_x, const double * __restrict__ d_y, const int * __restrict__ d_cumsum,
//	const double sb, const double tb, const int len_in,
//	const int Mrx, const int Mry, const double b, const double t1, const double xb, const double yb, const double Dx, const double Dy, const double offsetx, const double offsety, const int num_bins, const int msp) {
//
//	int i = threadIdx.x + blockIdx.x * blockDim.x;
//	int j = threadIdx.y + blockIdx.y * blockDim.y;
//
//	int tid = j * Mrx + i;
//
//	if ((i < Mrx) && (j < Mry)) {
//
//		double fact = 1. / (4.* PI *b);
//		fact = fact * fact;
//
//		double2 temp_sum;
//
//		temp_sum.x = 0.;
//		temp_sum.y = 0.;
//
//		double d_n_si = i - (int)Mrx / 2;
//		double d_n_sj = j - (int)Mry / 2;
//
//		uint32_t indicesx = floor((d_n_si + offsetx) / (2 * msp + 1));
//		uint32_t indicesy = floor((d_n_sj + offsety) / (2 * msp + 1));
//
//		uint32_t d_code = encode_morton2d(indicesx, indicesy);				// --- Center
//		uint32_t d_code_W = encode_morton2d(indicesx - 1, indicesy);				// --- West
//		uint32_t d_code_E = encode_morton2d(indicesx + 1, indicesy);				// --- East
//		uint32_t d_code_S = encode_morton2d(indicesx, indicesy - 1);			// --- South
//		uint32_t d_code_N = encode_morton2d(indicesx, indicesy + 1);			// --- North
//		uint32_t d_code_NE = encode_morton2d(indicesx + 1, indicesy + 1);			// --- North-East
//		uint32_t d_code_NW = encode_morton2d(indicesx - 1, indicesy + 1);			// --- North-West
//		uint32_t d_code_SE = encode_morton2d(indicesx + 1, indicesy - 1);			// --- South-East
//		uint32_t d_code_SW = encode_morton2d(indicesx - 1, indicesy - 1);			// --- South-West
//
//		if ((d_code < num_bins) && ((d_cumsum[d_code] - d_cumsum[d_code - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code - 1], d_y + d_cumsum[d_code - 1], d_f + d_cumsum[d_code - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code] - d_cumsum[d_code - 1]);
//		}
//		if ((d_code_W < num_bins) && ((d_cumsum[d_code_W] - d_cumsum[d_code_W - 1]) > 0) && (d_code_W > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_W - 1], d_y + d_cumsum[d_code_W - 1], d_f + d_cumsum[d_code_W - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_W] - d_cumsum[d_code_W - 1]);
//		}
//		if ((d_code_E < num_bins) && ((d_cumsum[d_code_E] - d_cumsum[d_code_E - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_E - 1], d_y + d_cumsum[d_code_E - 1], d_f + d_cumsum[d_code_E - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_E] - d_cumsum[d_code_E - 1]);
//		}
//		if ((d_code_S < num_bins) && ((d_cumsum[d_code_S] - d_cumsum[d_code_S - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_S - 1], d_y + d_cumsum[d_code_S - 1], d_f + d_cumsum[d_code_S - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_S] - d_cumsum[d_code_S - 1]);
//		}
//		if ((d_code_N < num_bins) && ((d_cumsum[d_code_N] - d_cumsum[d_code_N - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_N - 1], d_y + d_cumsum[d_code_N - 1], d_f + d_cumsum[d_code_N - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_N] - d_cumsum[d_code_N - 1]);
//		}
//		if ((d_code_NE < num_bins) && ((d_cumsum[d_code_NE] - d_cumsum[d_code_NE - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_NE - 1], d_y + d_cumsum[d_code_NE - 1], d_f + d_cumsum[d_code_NE - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_NE] - d_cumsum[d_code_NE - 1]);
//		}
//		if ((d_code_NW < num_bins) && ((d_cumsum[d_code_NW] - d_cumsum[d_code_NW - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_NW - 1], d_y + d_cumsum[d_code_NW - 1], d_f + d_cumsum[d_code_NW - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_NW] - d_cumsum[d_code_NW - 1]);
//		}
//		if ((d_code_SE < num_bins) && ((d_cumsum[d_code_SE] - d_cumsum[d_code_SE - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_SE - 1], d_y + d_cumsum[d_code_SE - 1], d_f + d_cumsum[d_code_SE - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_SE] - d_cumsum[d_code_SE - 1]);
//		}
//		if ((d_code_SW < num_bins) && ((d_cumsum[d_code_SW] - d_cumsum[d_code_SW - 1]) > 0)) {
//			temp_sum = addContributions(d_x + d_cumsum[d_code_SW - 1], d_y + d_cumsum[d_code_SW - 1], d_f + d_cumsum[d_code_SW - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_SW] - d_cumsum[d_code_SW - 1]);
//		}
//
//		d_f_tau[tid].x = fact * exp(b * (Dy * Dy * d_n_sj * d_n_sj + Dx * Dx * d_n_si * d_n_si)) * temp_sum.x;
//		d_f_tau[tid].y = fact * exp(b * (Dy * Dy * d_n_sj * d_n_sj + Dx * Dx * d_n_si * d_n_si)) * temp_sum.y;
//
//	}
//}
//
//
////void spatialConvolutionKernelNoAtomicCPU(double2 * __restrict__ d_f_tau, const double2 * __restrict d_f, const double * __restrict__ d_x, const double * __restrict__ d_y, const int * __restrict__ d_cumsum,
////	const double sb, const double tb, const int len_in,
////	const int Mrx, const int Mry, const double b, const double t1, const double xb, const double yb, const double Dx, const double Dy, const double offsetx, const double offsety, const int num_bins, const int msp) {
////
////	for (int j = 0; j < Mry; j++)
////		for (int i = 0; i < Mrx; i++) {
////
////		int tid = j * Mrx + i;
////
////		double fact = 1. / (4.* PI *b);
////		fact = fact * fact;
////
////		double2 temp_sum;
////
////		temp_sum.x = 0.;
////		temp_sum.y = 0.;
////
////		double2 temp_sum_partial;
////
////		double d_n_si = i - (int)Mrx / 2;
////		double d_n_sj = j - (int)Mry / 2;
////
////		uint32_t indicesx = floor((d_n_si + offsetx) / (2 * msp + 1));
////		uint32_t indicesy = floor((d_n_sj + offsety) / (2 * msp + 1));
////
////		uint32_t d_code = encode_morton2d(indicesx, indicesy);						// --- Center
////		uint32_t d_code_W = encode_morton2d(indicesx - 1, indicesy);				// --- West
////		uint32_t d_code_E = encode_morton2d(indicesx + 1, indicesy);				// --- East
////		uint32_t d_code_S = encode_morton2d(indicesx, indicesy - 1);				// --- South
////		uint32_t d_code_N = encode_morton2d(indicesx, indicesy + 1);				// --- North
////		uint32_t d_code_NE = encode_morton2d(indicesx + 1, indicesy + 1);			// --- North-East
////		uint32_t d_code_NW = encode_morton2d(indicesx - 1, indicesy + 1);			// --- North-West
////		uint32_t d_code_SE = encode_morton2d(indicesx + 1, indicesy - 1);			// --- South-East
////		uint32_t d_code_SW = encode_morton2d(indicesx - 1, indicesy - 1);			// --- South-West
////
////		if ((d_code < num_bins) && ((d_cumsum[d_code] - d_cumsum[d_code - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code - 1], d_y + d_cumsum[d_code - 1], d_f + d_cumsum[d_code - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code] - d_cumsum[d_code - 1]);
////		}
////		if ((d_code_W < num_bins) && ((d_cumsum[d_code_W] - d_cumsum[d_code_W - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_W - 1], d_y + d_cumsum[d_code_W - 1], d_f + d_cumsum[d_code_W - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_W] - d_cumsum[d_code_W - 1]);
////		}
////		if ((d_code_E < num_bins) && ((d_cumsum[d_code_E] - d_cumsum[d_code_E - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_E - 1], d_y + d_cumsum[d_code_E - 1], d_f + d_cumsum[d_code_E - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_E] - d_cumsum[d_code_E - 1]);
////		}
////		if ((d_code_S < num_bins) && ((d_cumsum[d_code_S] - d_cumsum[d_code_S - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_S - 1], d_y + d_cumsum[d_code_S - 1], d_f + d_cumsum[d_code_S - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_S] - d_cumsum[d_code_S - 1]);
////		}
////		if ((d_code_N < num_bins) && ((d_cumsum[d_code_N] - d_cumsum[d_code_N - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_N - 1], d_y + d_cumsum[d_code_N - 1], d_f + d_cumsum[d_code_N - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_N] - d_cumsum[d_code_N - 1]);
////		}
////		if ((d_code_NE < num_bins) && ((d_cumsum[d_code_NE] - d_cumsum[d_code_NE - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_NE - 1], d_y + d_cumsum[d_code_NE - 1], d_f + d_cumsum[d_code_NE - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_NE] - d_cumsum[d_code_NE - 1]);
////		}
////		if ((d_code_NW < num_bins) && ((d_cumsum[d_code_NW] - d_cumsum[d_code_NW - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_NW - 1], d_y + d_cumsum[d_code_NW - 1], d_f + d_cumsum[d_code_NW - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_NW] - d_cumsum[d_code_NW - 1]);
////		}
////		if ((d_code_SE < num_bins) && ((d_cumsum[d_code_SE] - d_cumsum[d_code_SE - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_SE - 1], d_y + d_cumsum[d_code_SE - 1], d_f + d_cumsum[d_code_SE - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_SE] - d_cumsum[d_code_SE - 1]);
////		}
////		if ((d_code_SW < num_bins) && ((d_cumsum[d_code_SW] - d_cumsum[d_code_SW - 1]) > 0)) {
////			temp_sum = addContributions(d_x + d_cumsum[d_code_SW - 1], d_y + d_cumsum[d_code_SW - 1], d_f + d_cumsum[d_code_SW - 1], fact, temp_sum, xb, yb, sb, tb, Dx, Dy, msp, d_n_si, d_n_sj, t1, b, d_cumsum[d_code_SW] - d_cumsum[d_code_SW - 1]);
////		}
////
////		d_f_tau[tid].x = temp_sum.x;
////		d_f_tau[tid].y = temp_sum.y;
////
////		}
////
////}
//
///***********************/
///* PARTITIONING KERNEL */
///***********************/
////__global__ void partitionKernel(const double * __restrict__ d_x, const double * __restrict__ d_y, uint32_t * __restrict__ d_code, const double offsetx, const double offsety, const double xb, const double yb,
////	const double Dx, const double Dy, const int msp, const int N) {
////
////	const int tid = threadIdx.x + blockIdx.x * blockDim.x;
////
////	if (tid > N) return;
////
////	double d_x_temp = (d_x[tid] - xb) / Dx + offsetx;
////	double d_y_temp = (d_y[tid] - yb) / Dy + offsety;
////
////	uint32_t indicesx = floor(d_x_temp / (2 * msp + 1));
////	uint32_t indicesy = floor(d_y_temp / (2 * msp + 1));
////
////	d_code[tid] = encode_morton2d(indicesx, indicesy);
////
////}
//
//__global__ void partitionShiftKernel(double * __restrict__ d_x, double * __restrict__ d_y, double2 * __restrict__ d_f, uint32_t * __restrict__ d_code, const double offsetx, const double offsety,
//	const double xb, const double yb, const double sb, const double tb, const double Dx, const double Dy, const int msp, const int N) {
//
//	const int tid = threadIdx.x + blockIdx.x * blockDim.x;
//
//	if (tid > N) return;
//
//	double d_x_temp1 = (d_x[tid] - xb) / Dx;
//	double d_y_temp1 = (d_y[tid] - yb) / Dy;
//
//	double d_x_temp = d_x_temp1 + offsetx;
//	double d_y_temp = d_y_temp1 + offsety;
//
//	uint32_t indicesx = floor(d_x_temp / (2 * msp + 1));
//	uint32_t indicesy = floor(d_y_temp / (2 * msp + 1));
//
//	d_code[tid] = encode_morton2d(indicesx, indicesy);
//
//	double arg = sb * d_x[tid] + tb * d_y[tid];
//
//	double real_part_temp = d_f[tid].x * cos(arg) + d_f[tid].y * sin(arg);
//	double imag_part_temp = d_f[tid].y * cos(arg) - d_f[tid].x * sin(arg);
//
//	d_f[tid].x = real_part_temp;
//	d_f[tid].y = imag_part_temp;
//
//	d_x[tid] = d_x_temp1;
//	d_y[tid] = d_y_temp1;
//}
//
////gpuErrchk(cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 32768));
//
////thrust::device_vector<double> d_x_sorted(len_in);
////thrust::device_vector<double> d_y_sorted(len_in);
////thrust::device_vector<double2> d_f_sorted(len_in);
//
////thrust::device_vector<uint32_t> d_code(len_in);
//
////	double offsetx = std::max(0.5 * Mrx, (Max_x - xb) / Dx) + 2. * (double)(2 * msp + 1);
////	double offsety = std::max(0.5 * Mry, (Max_y - yb) / Dy) + 2. * (double)(2 * msp + 1);
////	
////	timerGPU.StartCounter();
////	//partitionKernel << <iDivUp(len_in, BLOCKSIZE), BLOCKSIZE >> >(d_x, d_y, thrust::raw_pointer_cast(d_code.data()), offsetx, offsety, xb, yb, Dx, Dy, msp, len_in);
////	partitionShiftKernel << <iDivUp(len_in, BLOCKSIZE), BLOCKSIZE >> >(d_x, d_y, d_f, thrust::raw_pointer_cast(d_code.data()), offsetx, offsety, xb, yb, sb, tb, Dx, Dy, msp, len_in);
////#ifdef DEBUG
////	gpuErrchk(cudaPeekAtLastError());
////	gpuErrchk(cudaDeviceSynchronize());
////#endif
////	//timingFile << "Partition kernel " << timerGPU.GetCounter() << "\n";
////
////	// --- Initialize indices vector to [0, 1, 2, ...]
////	//timerGPU.StartCounter();
////	thrust::counting_iterator<int> iter(0);
////	thrust::device_vector<int> indices(len_in);
////	thrust::copy(iter, iter + indices.size(), indices.begin());
////	
////		// --- First, sort the keys and indices by the keys
////		thrust::sort_by_key(d_code.begin(), d_code.end(), indices.begin());
////	
////		// --- Now reorder the ID arrays using the sorted indices
////		thrust::gather(indices.begin(), indices.end(), thrust::device_pointer_cast(d_x), d_x_sorted.begin());
////		thrust::gather(indices.begin(), indices.end(), thrust::device_pointer_cast(d_y), d_y_sorted.begin());
////		thrust::gather(indices.begin(), indices.end(), thrust::device_pointer_cast(d_f), d_f_sorted.begin());
////	
////		thrust::device_vector<int> d_cumsum;
////	
////		// --- The number of d_cumsum bins is equal to the maximum value plus one
////		int num_bins = d_code.back() + 1;
////	
////		// --- Resize d_cumsum storage
////		d_cumsum.resize(num_bins);
////	
////		// --- Find the end of each bin of values - Cumulative d_cumsum
////		thrust::counting_iterator<int> search_begin(0);
////		thrust::upper_bound(d_code.begin(), d_code.end(), search_begin, search_begin + num_bins, d_cumsum.begin());
////	//timingFile << "Parte thrust " << timerGPU.GetCounter() << "\n";
////
////	//timerGPU.StartCounter();
////	dim3 GridSize_SpatialScaling(iDivUp(Mrx, BLOCKSIZE_X_SPATIAL_SCALING), iDivUp(Mry, BLOCKSIZE_Y_SPATIAL_SCALING));
////	dim3 BlockSize_SpatialScaling(BLOCKSIZE_X_SPATIAL_SCALING, BLOCKSIZE_Y_SPATIAL_SCALING);
////	spatialConvolutionKernelNoAtomic << <GridSize_SpatialScaling, BlockSize_SpatialScaling >> >(d_f_tau, thrust::raw_pointer_cast(d_f_sorted.data()), thrust::raw_pointer_cast(d_x_sorted.data()), thrust::raw_pointer_cast(d_y_sorted.data()), thrust::raw_pointer_cast(d_cumsum.data()), sb, tb, len_in, Mrx, Mry, b, t1, xb, yb, Dx, Dy, offsetx, offsety, num_bins, msp);
////#ifdef DEBUG
////	gpuErrchk(cudaPeekAtLastError());
////	gpuErrchk(cudaDeviceSynchronize());
////#endif
////	//timingFile << "Parte kernel " << timerGPU.GetCounter() << "\n";
//
////timingFile << timerGPU.GetCounter();
//
