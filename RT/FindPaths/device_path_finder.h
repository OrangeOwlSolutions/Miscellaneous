//#ifndef __DEVICE_PATH_FINDER_H__
//#define __DEVICE_PATH_FINDER_H__
//
//#include <vector_types.h>			// --- CUDA definition of float4
//#include <vector_functions.h>		// --- make_float4
//#include <thrust/host_vector.h>
//
//#include "path_finder.h"
//#include "aabb.h"
//#include "RXTX.h"
//#include "ray_path.h"
//
//pathfinder_p device_create_path_finder(
//    const size_t                        device_idx,    
//    const thrust::host_vector<float4>   &verts,
//    const thrust::host_vector<int4>     &faces,
//    const thrust::host_vector<float4>   &normals_and_k1,
//    const thrust::host_vector<float4>   &x1_and_k2,    
//    const thrust::host_vector<int>      &remapping_table,
//    const thrust::host_vector<int2>     &N0,
//    const thrust::host_vector<float4>   &N1,
//    const thrust::host_vector<float4>   &N2,
//    const thrust::host_vector<float4>   &N3,
//    const aabb                          &bb,
//    const size_t                        max_cnt = 20
//);
//
//void device_find_paths(pathfinder_p);
//
//void device_reset(pathfinder_p, std::vector<TX> *txs, std::vector<RX> *rxs);
//
//void device_append_new_paths(pathfinder_p, std::vector<path> &global_path_list);
//
//void device_destroy_path_finder(pathfinder_p);
//
//#endif 

#ifndef __DEVICE_PATH_FINDER_H__
#define __DEVICE_PATH_FINDER_H__

#include <vector_types.h>			// --- CUDA definition of float4
#include <vector_functions.h>		// --- make_float4

#include <curand_kernel.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "path_finder.h"
#include "aabb.h"
#include "RXTX.h"
#include "ray_path.h"
#include "ray_list.h"
#include "geointer.h"

/*********************/
/* PATHFINDER STRUCT */
/*********************/
struct device_pathfinder: public pathfinder {

	thrust::device_vector<float4>		verts;
	thrust::device_vector<int4>			faces;
	thrust::device_vector<float4>		normals_and_k1;
	thrust::device_vector<float4>		x1_and_k2; 
	thrust::device_vector<int>			remapping_table;
	thrust::device_vector<int2>			N0;
	thrust::device_vector<float4>		N1;
	thrust::device_vector<float4>		N2;
	thrust::device_vector<float4>		N3;
	aabb								bb;
	
	std::vector<TX>						*txs;				// --- Transmitter locations. Set in the device_reset function.
	std::vector<RX>						*rxs;				// --- Receiver locations. Set in the device_reset function.

	size_t								device_idx;			// --- Device identification

	std::vector<path>					path_list;

	device_ray_list						global_ray_list;	
	device_ray_list						starting_ray_list;	// --- Starting positions of the rays
	device_ray_list						saved_ray_list;		// --- Saved ray list
  
	thrust::device_vector<geointer>		ginters;			
	thrust::device_vector<float>		dinters;			// --- Distance of each ray from the origin
  
	thrust::device_vector<curandState>	rnd_states;

	thrust::device_vector<TX>			dtxs;				// --- Detectors (used by the generator of new rays) ???
  
	size_t								max_cnt;
  
	// --- Pathfinder constructor
	device_pathfinder(const size_t device_idx_, const thrust::host_vector<float4> &verts_, const thrust::host_vector<int4> &faces_,
					  const thrust::host_vector<float4> &normals_and_k1_, const thrust::host_vector<float4> &x1_and_k2_, 
					  const thrust::host_vector<int> &remapping_table_, const thrust::host_vector<int2> &N0_, const thrust::host_vector<float4> &N1_,
					  const thrust::host_vector<float4> &N2_, const thrust::host_vector<float4> &N3_, const aabb &bb_, const size_t max_cnt_);
  
};

pathfinder_p device_create_path_finder(
    const size_t                        device_idx,    
    const thrust::host_vector<float4>   &verts,
    const thrust::host_vector<int4>     &faces,
    const thrust::host_vector<float4>   &normals_and_k1,
    const thrust::host_vector<float4>   &x1_and_k2,    
    const thrust::host_vector<int>      &remapping_table,
    const thrust::host_vector<int2>     &N0,
    const thrust::host_vector<float4>   &N1,
    const thrust::host_vector<float4>   &N2,
    const thrust::host_vector<float4>   &N3,
    const aabb                          &bb,
    const size_t                        max_cnt = 20
);

void device_find_paths(pathfinder_p);

void device_reset(pathfinder_p, std::vector<TX> *txs, std::vector<RX> *rxs);

void device_append_new_paths(pathfinder_p, std::vector<path> &global_path_list);

void device_destroy_path_finder(pathfinder_p);

#endif 
