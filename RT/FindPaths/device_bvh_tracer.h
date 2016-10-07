#ifndef __DEVICE_BVH_TRACER_H__
#define __DEVICE_BVH_TRACER_H__

#include "ray_list.h"
#include "geointer.h"
#include "aabb.h"

void device_bvh_trace(device_ray_list::R0_R1_iterator begin, device_ray_list::R0_R1_iterator end, thrust::device_vector<int4>::const_iterator faces_begin, 
	                  thrust::device_vector<float4>::const_iterator verts_begin, thrust::device_vector<int>::const_iterator rmap_begin,
					  thrust::device_vector<int2>::const_iterator N0_begin, thrust::device_vector<float4>::const_iterator N1_begin,
					  thrust::device_vector<float4>::const_iterator N2_begin, thrust::device_vector<float4>::const_iterator N3_begin,    
					  thrust::device_vector<geointer>::iterator geointer_begin, const aabb &scn_aabb);

#endif
