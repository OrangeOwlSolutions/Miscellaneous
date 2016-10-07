#ifndef __BVH_NODE_H__
#define __BVH_NODE_H__

#ifdef __CUDACC__
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else 
#define HOST_DEVICE_INLINE  inline
#endif

/***************************/
/* CHECK IF A NODE IS LEAF */
/***************************/
// --- By convention, a leaf node has a negative address in its parent node
HOST_DEVICE_INLINE bool is_leaf(int node_idx) { return node_idx < 0; }

// --- Returns index of first child node
HOST_DEVICE_INLINE int get_c0_idx(const int2 &N0) { return N0.x; }

// --- Returns index of first child node
HOST_DEVICE_INLINE int get_c1_idx(const int2 &N0) { return N0.y; }

// --- Finds the coordinates of the "minimum point" of the first child AABB
HOST_DEVICE_INLINE vec3 get_c0_min_vec_from_N1_N3(const float4 &N1, const float4 &N3) { return make_vec3(N1.x, N1.z, N3.x); }

// --- Finds the coordinates of the "maximum point" of the first child AABB
HOST_DEVICE_INLINE vec3 get_c0_max_vec_from_N1_N3(const float4 &N1, const float4 &N3) { return make_vec3(N1.y, N1.w, N3.y); }

// --- Finds the coordinates of the "minimum point" of the second child AABB
HOST_DEVICE_INLINE vec3 get_c1_min_vec_from_N2_N3(const float4 &N2, const float4 &N3) { return make_vec3(N2.x, N2.z, N3.z); }

// --- Finds the coordinates of the "maximum point" of the second child AABB
HOST_DEVICE_INLINE vec3 get_c1_max_vec_from_N2_N3(const float4 &N2, const float4 &N3) { return make_vec3(N2.y, N2.w, N3.w); }

#endif
