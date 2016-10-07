#ifndef __AABB_H__
#define __AABB_H__

#include "vec3.h"

#ifdef __CUDACC__
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else 
#include <cmath>
#define HOST_DEVICE_INLINE  inline
#endif

/***************/
/* AABB STRUCT */
/***************/
struct aabb {
    
	bool empty;
    vec3 min_vec;
    vec3 max_vec;
    
	// --- Constructs an empty AABB: the two vertices coincide
    HOST_DEVICE_INLINE
    aabb() : empty(true), min_vec(make_vec3(0,0,0)), max_vec(make_vec3(0,0,0)) {};

    // --- Constructs a non-empty AABB
	HOST_DEVICE_INLINE
    aabb(const vec3 &min_vec_, const vec3 &max_vec_) : empty(false), min_vec(min_vec_), max_vec(max_vec_) {}
    
    // --- Finds the envelope of two AABBs
	HOST_DEVICE_INLINE
    void grow(const aabb &bb) { grow(bb.min_vec); grow(bb.max_vec); }
    
    HOST_DEVICE_INLINE
    void  grow(const vec3 &v) {
        if(!empty) {
            min_vec = cwisemin(min_vec, v);
            max_vec = cwisemax(max_vec, v);
        } else {
            min_vec = v;
            max_vec = v;
            empty = false;
        }     
    }
    
    // --- Magnifies the AABB
	HOST_DEVICE_INLINE
    void  scale(float factor) { *this = this->scaled(factor); }

    HOST_DEVICE_INLINE
    aabb  scaled(float factor) const {
        vec3 center		= get_center();
        vec3 half_size	= get_half_size();
        return aabb(center - half_size * factor, center + half_size * factor);
    }

	// --- Gets the center of the AABB
    HOST_DEVICE_INLINE
    vec3  get_center() const { return (max_vec + min_vec) * 0.5f; }

	// --- Gets the half-size of the AABB
    HOST_DEVICE_INLINE
    vec3 get_half_size() const { return (max_vec - min_vec) * 0.5f; }

	// --- Returns if the AABB is empty
    HOST_DEVICE_INLINE
    bool is_empty() const { return empty; }

	// --- Returns if the AABB is planar
    HOST_DEVICE_INLINE
    bool is_planar() const {
        vec3 s = max_vec - min_vec;
        return  almost_equal(s.x, 0.0f) || almost_equal(s.y, 0.0f) || almost_equal(s.z, 0.0f) ;
    }

	// --- Sets the extremals of the AABB
    HOST_DEVICE_INLINE
    void set(const vec3 &min_vec_, const vec3 &max_vec_) {
		min_vec = min_vec_; 
        max_vec = max_vec_;    
    }

	// --- Gets the overall surface area of the AABB
    HOST_DEVICE_INLINE
    float get_surface_area() const {
        vec3 side = max_vec - min_vec;
        return 2.0f * (side.x * side.y + side.x * side.z + side.y * side.z);
    }

	// --- Gets the overall volume of the AABB
    HOST_DEVICE_INLINE
    float get_vol() const {
        vec3 side = max_vec - min_vec;
        return side.x * side.y * side.z;
    } 
};

// --- Prints the AABB
std::ostream &operator<<(std::ostream &os, const aabb & bb);

#endif
