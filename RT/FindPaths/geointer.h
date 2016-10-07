#ifndef __GEO_INTERSECTION_H__
#define __GEO_INTERSECTION_H__

#ifdef __CUDACC__
#define ALING_FOUR  __builtin_align__(16)
//#define ALING_FOUR
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else 
#ifdef __linux__
#define ALIGN_FOUR __attribute__ ((aligned (16)))
#define HOST_DEVICE_INLINE  inline
#elif _WIN32 || _WIN64
#define ALING_FOUR __declspec(align(16))
#define ALING_FOUR 
#define HOST_DEVICE_INLINE  inline
#endif
#endif

/*******************/
/* GEOINTER STRUCT */
/*******************/
struct ALING_FOUR geointer {
    int		ID;					// --- ID of the intersection triangle
    float	t;					// --- t is the path length (actually, the distance) from the ray origin to the intersection point. 
								//     If t < 0, no intersection occurs.
    float	u;					// --- (u, v) baricentric coordinates of the intersection point
    float	v;
};

// --- Make geointer
HOST_DEVICE_INLINE geointer make_geointer(int ID, float t, float u, float v) { geointer g = {ID, t, u, v}; return g; }

/********************************/
/* IS VALID INTERSECTION STRUCT */
/********************************/
// ??? Penso non sia mai usata
struct is_valid_intersection{

	HOST_DEVICE_INLINE bool operator()(const geointer & g) const { return g.t >= 0; }
    
};

#undef ALING_FOUR
#undef HOST_DEVICE_INLINE

#endif
