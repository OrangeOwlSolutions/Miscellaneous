#ifndef __RAY_LIST_H__
#define __RAY_LIST_H__


#include "vec3.h"

#include <vector_types.h> //cuda definition of float4
#include <vector_functions.h> //make_float4
#include "types_utils.h"

#ifdef __CUDACC__
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else 
#define HOST_DEVICE_INLINE  inline

#endif

#ifdef _WIN32 || _WIN64
typedef unsigned int uint;
#endif

// Description
//  R0: [orig.x, orig.y, orig.z, length]
//  R1: [dir_v.x,dir_v.y,dir_v.z, distance to the receiver]
//  R2: [c0,c1,c2,c3]  [code, rx, numdiff | numref, tx]

typedef thrust::tuple<float4, float4, uint4> ray_tuple;
typedef thrust::tuple<float4, float4>        R0_R1_tuple;

/*******************/
/* RAY LIST STRUCT */
/*******************/
// --- The ray list class is shared by host and device and is interfaced by a generic "container". The container specifies into a thrust::host_vector for
//     the host and into a thrust::device_vector for the device. The class is specified in the host_ray_list and device_ray_list below.
template <template <typename, typename> class container>
class ray_list {

private:
    
    typedef typename container<uint4,  typename type_utils::get_containter_allocator<container, uint4> ::allocator>::iterator u4_iterator;

	typedef typename container<float4, typename type_utils::get_containter_allocator<container, float4>::allocator>::iterator f4_iterator;
    
    typedef thrust::tuple<f4_iterator, f4_iterator, u4_iterator> iterator_tuple;
    
public:
  
	// --- R0: [orig.x, orig.y, orig.z, length]
	//         (orig.x, orig.y, orig.z) is the ray origin.
	//         length is the overall ray length which is updated at each reflection (this is done in reflect_rays_functor.h, operator()). knowing length
	//         lets knowing the ray length once it hits the receiver.
	container<float4, typename type_utils::get_containter_allocator<container, float4>::allocator> R0;

    // --- R1: [dir_v.x, dir_v.y, dir_v.z, distance to the receiver]
	//         (dir_v.x, dir_v.y, dir_v.z) is the ray unit vector.
	//         "distance to the receiver" is si used when the ray intersects the receiver sphere and denotes the distance between ray and receiver, meant
	//         as the distance between a point (the receiver center) and a line (the ray).
	container<float4, typename type_utils::get_containter_allocator<container, float4>::allocator> R1;

	// --- R2: [c0, c1, c2, c3]  [code, rx, numdiff | numref, tx]
	//         R2 identifies the different paths and enables their unique identification to separate duplicates. In particular:
	//         c0:      unsigned 32 bit which is updated with the cyclic redundancy check algorithm. The update occurs based on the previous value of c0
	//					and by the facet ID of the reflecting the ray. The facet ID is the fourth value in the faces field of nvmesh. This approach enables
	//                  to separate wavefronts following different paths. The update is operated by the crc_24 function.
	//         c1:	    unsigned 32 bit which contains the receiver's code on which the ray hits. The receivers are assigned a progressive number, 
	//                  as long as they are introduced in the scene. 
	//		   c2:      32 bits: the first 16 (less significant) contain the number of reflections, while the last 16 bits (most significant) contain the
	//                  number of diffractions undergone by the ray during its path.
	//         c3:		unsigned 32 bit which contains the transmitter's code. The transmitters are assigned a progressive number, as long as they are
	//                  introduced in the scene.
	container<uint4,  typename type_utils::get_containter_allocator<container, uint4> ::allocator> R2;
	    
    typedef thrust::zip_iterator<iterator_tuple> ray_list_iterator;

    typedef thrust::zip_iterator<thrust::tuple<f4_iterator, f4_iterator>> R0_R1_iterator;
    
public:    
    
	// --- Default constructor
	explicit  ray_list() :  R0(0), R1(0), R2(0) {}
    
	// --- Constructor: creates a list of N rays
    explicit  ray_list(size_t N) : R0(N), R1(N), R2(N) {}
    
    // --- Copy-constructor
    template <template <typename ,typename> class other_container>
    ray_list(const ray_list<other_container>& v) { R0 = v.R0; R1 = v.R1; R2 = v.R2; }
    
    // --- Operator =
	template <template <typename ,typename > class other_container>
    ray_list& operator = (const ray_list<other_container>& v) { R0 = v.R0; R1 = v.R1; R2 = v.R2; return *this; }

    // --- Returns the number of rays
	size_t size() const { return R0.size(); }

    // --- Resize the rays array
	void resize(size_t N) { R0.resize(N); R1.resize(N); R2.resize(N); }

    // --- Reserve space for the rays array
	void reserve(size_t N) { R0.reserve(N); R1.reserve(N); R2.reserve(N); }

    // --- Ray list iterator begin
    ray_list_iterator begin() { return thrust::make_zip_iterator(thrust::make_tuple(R0.begin(), R1.begin(), R2.begin())); }

    // --- Ray list end
	ray_list_iterator end() { return begin() + size() ; }
    
    // --- (R0, R1) iterator begin
	R0_R1_iterator begin_R0_R1() { return thrust::make_zip_iterator(thrust::make_tuple(R0.begin(), R1.begin())); }

    // --- (R0, R1) iterator end
	R0_R1_iterator end_R0_R1()   { return begin_R0_R1() + size(); }

};

// --- Specifies the ray list class into a host ray list and a device ray list.
typedef ray_list<thrust::host_vector>	host_ray_list;
typedef ray_list<thrust::device_vector> device_ray_list;

// --- Return ray origin
HOST_DEVICE_INLINE vec3 get_orig_from_R0(const float4 &R0) { return make_vec3(R0.x, R0.y, R0.z); }

// --- Return ray length
HOST_DEVICE_INLINE float get_len_from_R0(const float4 &R0) { return R0.w; }

// --- Return ray direction
HOST_DEVICE_INLINE vec3 get_dir_v_from_R1(const float4 &R1) { return make_vec3(R1.x, R1.y, R1.z); }

// --- Return ray distance from detector
HOST_DEVICE_INLINE float get_dist_from_R1(const float4 &R1) { return R1.w; }

// --- Make R0
HOST_DEVICE_INLINE float4 make_R0_from_orig_and_len(const vec3& orig, const float len) { return make_float4(orig.x, orig.y, orig.z, len); }

// --- Make R1
HOST_DEVICE_INLINE float4 make_R1_from_dir_v_and_dist(const vec3& dir_v, const float dist) { return make_float4(dir_v.x, dir_v.y, dir_v.z, dist); }

HOST_DEVICE_INLINE uint4 make_R2(const uint code, const uint rx, const uint num_diff, const uint num_ref, const uint tx) {
    return make_uint4(code,rx,num_diff << 16 | num_ref,tx);
}

#undef HOST_DEVICE_INLINE

#endif
