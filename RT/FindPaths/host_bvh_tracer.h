#ifndef __HOST_BVH_TRACER_H__
#define __HOST_BVH_TRACER_H__

#include "ray_list.h"
#include "geointer.h"
#include "aabb.h"

struct  host_compute_geo_inters_with_bvh {

    // --- Iterators for the bvh    
    thrust::host_vector<int4>::const_iterator     faces_begin;
    thrust::host_vector<float4>::const_iterator   verts_begin;
    
    thrust::host_vector<int>::const_iterator      rmap_begin;

	thrust::host_vector<int2>::const_iterator     N0_begin;
    thrust::host_vector<float4>::const_iterator   N1_begin;
    thrust::host_vector<float4>::const_iterator   N2_begin;
    thrust::host_vector<float4>::const_iterator   N3_begin;
    
    aabb										  scn_aabb;
    
    
    host_compute_geo_inters_with_bvh(thrust::host_vector<int4>::const_iterator faces_begin_, thrust::host_vector<float4>::const_iterator verts_begin_,
									 thrust::host_vector<int>::const_iterator rmap_begin_, thrust::host_vector<int2>::const_iterator N0_begin_,
									 thrust::host_vector<float4>::const_iterator N1_begin_, thrust::host_vector<float4>::const_iterator N2_begin_,
									 thrust::host_vector<float4>::const_iterator N3_begin_, aabb scn_aabb_);
    
    // --- Functions to extract near or far distance 
    inline float near_fun(const float2 &near_far);    
    inline float far_fun (const float2 &near_far);
    
    geointer operator()(const vec3 &orig, const vec3 &dir);
    
    geointer operator()(const R0_R1_tuple &r);

};

void host_bvh_trace(host_ray_list::R0_R1_iterator begin, host_ray_list::R0_R1_iterator end, thrust::host_vector<int4>::const_iterator faces_begin,    
					thrust::host_vector<float4>::const_iterator verts_begin, thrust::host_vector<int>::const_iterator rmap_begin, 
					thrust::host_vector<int2>::const_iterator N0_begin, thrust::host_vector<float4>::const_iterator N1_begin,
					thrust::host_vector<float4>::const_iterator N2_begin, thrust::host_vector<float4>::const_iterator N3_begin,
					thrust::host_vector<geointer>::iterator geointer_begin, const aabb &scn_aabb);

#endif

