#ifndef __HOST_PATH_FINDER_H__
#define __HOST_PATH_FINDER_H__

#include <vector_types.h> //cuda definition of float4
#include <vector_functions.h> //make_float4

#include <thrust/host_vector.h>

#include "path_finder.h"
#include "aabb.h"
#include "RXTX.h"
#include "ray_path.h"
#include "edge.h"

pathfinder_p host_create_path_finder(    
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
    const std::vector<edge>             &edge_list
);

void host_find_paths(pathfinder_p);

void host_append_new_paths(pathfinder_p, std::vector<path> &global_path_list);

void host_reset(pathfinder_p, std::vector<TX> *txs, std::vector<RX> *rxs);

void host_destroy_path_finder(pathfinder_p);

#endif
