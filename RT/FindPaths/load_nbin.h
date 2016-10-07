//#ifndef __LOAD_NBIN_H__
//#define __LOAD_NBIN_H__
//
//#include <iostream>
//#include <string>
//
//#include <vector_types.h>
//#include <vector_functions.h> //make_float4
//
//#include <thrust\host_vector.h>
//
//#include "packed_data_types.h"
//#include "aabb.h"
//
///***************/
///* MESH STRUCT */
///***************/
//struct nvmesh{
//   
//    thrust::host_vector<float4>		verts;				// --- List of the mesh vertices. Each vertex is a flat4 - (x, y, z, 1)
//    thrust::host_vector<int4>		faces;				// --- List of the mesh faces. Each face (triangle) is identified by the indices of the three
//														//     vertices and by a fourth integer, namely is like (a, b, c, idx). idx identifies the sub-mesh
//														//     the triangle belongs to and is used by the code-based approach to separate the wavefronts.
//    thrust::host_vector<int>		remapping_table;	// --- Maps the triangles in the BVH leaves to the original triangle ordering.
//    thrust::host_vector<float4>		normals_and_k1;		// --- The three components of the normals to the vertices and the first principal curvature k1 
//														//     at the same vertices.
//    thrust::host_vector<float4>		x1_and_k2;			// --- Principal direction at the vertices and second principal curvature k2 at the same vertices.
//    
//    // --- BVH nodes. The following fields contain information on the BVH nodes (either inner or leaves). The vectors N0, N1, N2 and N3 have the same
//	//     lengths equal to the number of BVH nodes. Each node is represented by an int2 and by 3 float4s. The contained information depend on the kind
//	//     of node, mamely, if it is internal or leaf. In the case of an internal node, assume c0 and c1 are the two child nodes. The following node format
//	//     is used. N1, N2 and N3 contain the coordinates of the two extremals of the child nodes. N0 reports the indices of the two child nodes and tells 
//	//     if the two child nodes are internal or leaves. In particular, by convention, if the node is internal, its index is stored unchanged, while if it
//	//     is a leaf the its bits are reversed. In this way, a positive node index denotes internal node, while a negative node index denotes a leaf node.
//	//     For a leaf node, only the triangle list contained in the node is stored (its bounding box is stored in the parent node). To this end, only the
//	//     starting and post-end indices of the remapping table are stored in N0. In other words, the portion of the remapping table containing the list
//	//     of the triangles belonging to the node is stored. In this case, N1, N2 and N3 are not used.
//
//    //  NVIDIA NODE FORMAT
//    //  N0:int2  [c0.inner or ~c0.leaf , c1.inner or ~c1.leaf]
//    //  N1:real4 [c0.lo.x, c0.hi.x, c0.lo.y, c0.hi.y]
//    //  N2:real4 [c1.lo.x, c1.hi.x, c1.lo.y, c1.hi.y]
//    //  N3:real4 [c0.lo.z, c0.hi.z, c1.lo.z, c1.hi.z]
//
//    thrust::host_vector<int2>		N0;
//    thrust::host_vector<float4>		N1;
//    thrust::host_vector<float4>		N2;
//    thrust::host_vector<float4>		N3;    
//    
//    aabb							bb;					// --- Bounding box of the whole scene
//    
//    float min_len;
//	
//	vec3 get_normal(const int ID, float u, float v) const; 
//};
//
//
//bool load_model_nbin(const std::string filename, nvmesh &m);
//
//std::ostream &operator << (std::ostream &os, const nvmesh &m);
//
//#endif
#ifndef __LOAD_NBIN_H__
#define __LOAD_NBIN_H__

#include <iostream>
#include <string>

#include <vector_types.h>
#include <vector_functions.h> //make_float4

#include <thrust\host_vector.h>

#include "packed_data_types.h"
#include "aabb.h"

/***************/
/* MESH STRUCT */
/***************/
struct nvmesh{
   
    thrust::host_vector<float4>		verts;				// --- List of the mesh vertices. Each vertex is a flat4 - (x, y, z, 1)
    thrust::host_vector<int4>		faces;				// --- List of the mesh faces. Each face (triangle) is identified by the indices of the three
														//     vertices and by a fourth integer, namely is like (a, b, c, idx). idx identifies the sub-mesh
														//     the triangle belongs to and is used by the code-based approach to separate the wavefronts.
    thrust::host_vector<int>		remapping_table;	// --- Maps the triangles in the BVH leaves to the original triangle ordering.
    thrust::host_vector<float4>		normals_and_k1;		// --- The three components of the normals to the vertices and the first principal curvature k1 
														//     at the same vertices.
    thrust::host_vector<float4>		x1_and_k2;			// --- Principal direction at the vertices and second principal curvature k2 at the same vertices.
    
    // --- BVH nodes. The following fields contain information on the BVH nodes (either inner or leaves). The vectors N0, N1, N2 and N3 have the same
	//     lengths equal to the number of BVH nodes. Each node is represented by an int2 and by 3 float4s. The contained information depend on the kind
	//     of node, mamely, if it is internal or leaf. In the case of an internal node, assume c0 and c1 are the two child nodes. The following node format
	//     is used. N1, N2 and N3 contain the coordinates of the two extremals of the child nodes. N0 reports the indices of the two child nodes and tells 
	//     if the two child nodes are internal or leaves. In particular, by convention, if the node is internal, its index is stored unchanged, while if it
	//     is a leaf the its bits are reversed. In this way, a positive node index denotes internal node, while a negative node index denotes a leaf node.
	//     For a leaf node, only the triangle list contained in the node is stored (its bounding box is stored in the parent node). To this end, only the
	//     starting and post-end indices of the remapping table are stored in N0. In other words, the portion of the remapping table containing the list
	//     of the triangles belonging to the node is stored. In this case, N1, N2 and N3 are not used.

    //  NVIDIA NODE FORMAT
    //  N0:int2  [c0.inner or ~c0.leaf , c1.inner or ~c1.leaf]
    //  N1:real4 [c0.lo.x, c0.hi.x, c0.lo.y, c0.hi.y]
    //  N2:real4 [c1.lo.x, c1.hi.x, c1.lo.y, c1.hi.y]
    //  N3:real4 [c0.lo.z, c0.hi.z, c1.lo.z, c1.hi.z]

    thrust::host_vector<int2>		N0;
    thrust::host_vector<float4>		N1;
    thrust::host_vector<float4>		N2;
    thrust::host_vector<float4>		N3;    
    
    aabb							bb;					// --- Bounding box of the whole scene
    
    float min_len;
	
	vec3 get_normal(const int ID, float u, float v) const; 

    void get_normal_and_curvature(const int ID, float u, float  v, vec3 &normal_v, vec3 &U1_v, vec3 &U2_v, float &a1, float &a2) const;
};


bool load_model_nbin(const std::string filename, nvmesh &m);

std::ostream &operator << (std::ostream &os, const nvmesh &m);

#endif
