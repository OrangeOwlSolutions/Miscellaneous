#ifndef __DEVICE_SIMPLE_TRACER_H__
#define __DEVICE_SIMPLE_TRACER_H__

#include "ray_list.h"
#include "geointer.h"

void device_simple_trace(device_ray_list::R0_R1_iterator begin, device_ray_list::R0_R1_iterator end, 
	                     thrust::device_vector<int4>::const_iterator faces_cbegin, thrust::device_vector<int4>::const_iterator faces_cend,
						 thrust::device_vector<float4>::const_iterator verts_cbegin, thrust::device_vector<geointer>::iterator geointer_begin);

#endif
