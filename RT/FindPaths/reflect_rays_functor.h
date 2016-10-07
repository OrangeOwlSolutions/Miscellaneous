#ifndef __REFLECT_RAYS_FUNCTOR__
#define __REFLECT_RAYS_FUNCTOR__

#include "vec3.h"
#include "geointer.h"
#include "geometry.h"
#include "crc24.h"

#ifdef __CUDACC__
    #define HOST_DEVICE_INLINE  __host__ __device__ inline
#else 
    #define HOST_DEVICE_INLINE  inline
#endif

/************************************************/
/* THRUST FUNCTOR TO COMPUTE THE REFLECTED RAYS */
/************************************************/
template<typename float4_iterator, typename int4_iterator > 
struct compute_reflected_ray {        
    
	float4_iterator  verts_begin;
    int4_iterator    faces_begin;  
    float4_iterator  normals_and_k1_begin;
    
    // --- Constructor
	compute_reflected_ray(float4_iterator verts_begin_, int4_iterator faces_begin_, float4_iterator normals_and_k1_begin_) : 
						  verts_begin(verts_begin_), faces_begin(faces_begin_), normals_and_k1_begin(normals_and_k1_begin_) {}
   
    // --- From the normals at the triangle vertices determines the normal at the intersection centers by interpolation
    HOST_DEVICE_INLINE vec3 get_normal(const geointer &gi) {
        
		const int ID = gi.ID;
        
        // --- Get the vertex indices of the triangle
        const int4 face   = faces_begin[ID];
        const uint v1_idx = face.x;
        const uint v2_idx = face.y;
        const uint v3_idx = face.z;
        
        // --- Get the normals
        const float4 data_l2_v1 = normals_and_k1_begin[v1_idx];
        const float4 data_l2_v2 = normals_and_k1_begin[v2_idx];
        const float4 data_l2_v3 = normals_and_k1_begin[v3_idx];
        
        // --- Normals at each vertex
        const vec3 normal_v1  = make_vec3(data_l2_v1.x,data_l2_v1.y,data_l2_v1.z);
        const vec3 normal_v2  = make_vec3(data_l2_v2.x,data_l2_v2.y,data_l2_v2.z);
        const vec3 normal_v3  = make_vec3(data_l2_v3.x,data_l2_v3.y,data_l2_v3.z);                
        
        // --- Barycentric coordinate of the intersection point
        const float u = gi.u;
        const float v = gi.v;
        
        // --- Perform interpolation
        const vec3 normal_v = normalize(normal_v1 + u * (normal_v2 - normal_v1) + v * (normal_v3 - normal_v1));
        
        return normal_v;
                                                       
    }
    
	// --- return the reflected ray    
	HOST_DEVICE_INLINE ray_tuple operator ()(const ray_tuple &r, const geointer &gi) { 

        // --- Extract raw data from the tuple
        const float4 R0 = thrust::get<0>(r);	// --- R0: [orig.x, orig.y, orig.z, length]
        const float4 R1 = thrust::get<1>(r);	// --- R1: [dir_v.x, dir_v.y, dir_v.z, distance to the receiver]
        const uint4  R2 = thrust::get<2>(r);	// --- R2: [c0, c1, c2, c3]
        const int4 face = faces_begin[gi.ID];
    
        const float  t  = gi.t;           
        
        const vec3  origin_i = get_orig_from_R0(R0);	// --- Ray origin
        const vec3  dir_v    = get_dir_v_from_R1(R1);   // --- Ray direction     
        
        float  len			 = get_len_from_R0(R0);		// --- Ray length (overall)
        float  dist_from_det = get_dist_from_R1(R1);	// --- Ray distance from detector
  
        // --- The origin of the reflected ray is the point of intersection, namely Qr
        vec3 origin_r		 = origin_i + dir_v * t;
        
        // --- Updates the overall ray distance: distance from the origin of the incident ray to the Qr
        len  += dist(origin_i, origin_r);
        
        // --- Unit normal at Qr
        vec3 normal_v = get_normal(gi);
               
        // --- Direction of the reflected ray
        const vec3 dir_r_v = normalize(dir_v - 2*dot(normal_v, dir_v) * normal_v);
        
        // --- Detach the origin of the reflected ray from the surface 
        // --- Divergenge factor
        const float s_detach = 0.001f; // --- A small value
        origin_r = origin_r + dir_r_v * s_detach;
        
        // --- Costructs the reflected ray
        const float4 R0r		= make_R0_from_orig_and_len(origin_r, len);            // --- R0: [orig.x, orig.y, orig.z, length]
        const float4 R1r		= make_R1_from_dir_v_and_dist(dir_r_v, dist_from_det); // --- [dir_v.x, dir_v.y, dir_v.z, distance to the receiver]
        const uint new_crc24	= crc_24(R2.x, face.w); 
        const uint4 R2r			= {new_crc24, R2.y, R2.z + 1, R2.w};  
       
        return thrust::make_tuple(R0r,R1r,R2r);
        
    }
};

/**********************/
/* MAKE REFLECTED RAY */
/**********************/
template< typename float4_iterator, typename int4_iterator>        
inline compute_reflected_ray<float4_iterator, int4_iterator> make_compute_reflected_ray(float4_iterator verts_begin, int4_iterator faces_begin,  
																						float4_iterator normals_and_k1_begin) {
    
	return compute_reflected_ray<float4_iterator,int4_iterator>(verts_begin, faces_begin, normals_and_k1_begin);
}

#undef HOST_DEVICE_INLINE

#endif
