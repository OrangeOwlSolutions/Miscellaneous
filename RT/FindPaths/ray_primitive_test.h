#ifndef __RAY_PRIMITIVE_TEST_H__
#define __RAY_PRIMITIVE_TEST_H__

#include <algorithm>

#include "vec3.h"

#ifdef __CUDACC__
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else 
#define HOST_DEVICE_INLINE  inline
#endif

/****************************************************************************/
/* FINDS A "MEDIAN" DISTANCE BETWEEN THE RAY ORIGIN AND THE DETECTOR CENTER */
/****************************************************************************/
HOST_DEVICE_INLINE float ray_sphere_median(const vec3 &orig, const vec3 &dir, const vec3 &center, const float &radius) {

	const vec3 a = orig - center;

	const float a_dot_dir   = dot(a, dir);
	const float dir_norm_sq = dot(dir, dir);
	const float delta = a_dot_dir * a_dot_dir - dir_norm_sq * (dot(a, a) - radius * radius);    
	
	if (delta <0) return -INF_R;
	
	float t_med = -a_dot_dir / dir_norm_sq;
	
    return t_med;
}


//
HOST_DEVICE_INLINE
float ray_plane_nearest(const vec3 &orig,
                        const vec3 &dir,
                        const vec3 &Q,
                        const vec3 &n
                       ){
    
    float den = dot(n,dir);
    if(den != 0){
        float t = dot(n,Q-orig)/den;
        return t;
    }
    
    return -1;
}

/******************************************/
/* MöLLER–TRUMBORE INTERSECTION ALGORITHM */
/******************************************/
HOST_DEVICE_INLINE float3 ray_triangle_nearest_with_uv(const vec3 &orig, const vec3 &dir, const vec3 &v1, const vec3 &v2, const vec3 &v3) {
    
	const float3 miss = make_float3(-1.f, 0.f, 0.f);
    
    vec3 e1, e2;  // --- Edge1, Edge2
    vec3 P, Q, T;
    float det, inv_det, u, v;
    float t;
 
    // --- Find vectors for two edges sharing V1
    e1 =  v2 - v1;
    e2 =  v3 - v1;

	// --- Begin calculating determinant - also used to calculate u parameter
    P = cross(dir, e2);
    
	// --- If determinant is close to zero, ray lies in the plane of the triangle
    det = dot(e1, P);
    
	// --- NOT CULLING
    if(det > -EPS_R && det < EPS_R) return miss;

    inv_det = 1.f / det;
 
    // --- Calculate distance from V1 to ray origin
    T = orig - v1;
	//printf("T.x %f T.y %f T.z %f\n", T.x, T.y, T.z);
 
    // --- Calculate u parameter and test bound
    u = dot(T, P) * inv_det;
	//printf("u %f\n", u);
    
	// --- The intersection lies outside of the triangle
    if (u < 0.f || u > 1.f) return miss;
 
	// --- Prepare to test v parameter
    Q = cross(T, e1);
	//printf("Q.x %f Q.y %f Q.z %f\n", Q.x, Q.y, Q.z);
 
    // --- Calculate V parameter and test bound
    v = dot(dir, Q) * inv_det;
    
	// --- The intersection lies outside of the triangle
    if(v < 0.f || u + v  > 1.f) return miss;
 
    t = dot(e2, Q) * inv_det;
	//printf("t %f\n", t);
 
    // --- Ray intersection  
	if(t < 0.f) return miss;
 
    return make_float3(t, u, v);
  
}

/**************************************************************/
/* FINDS NEAR AND FAR INTERSECTIONS BETWEEN A RAY AND AN AABB */
/**************************************************************/
HOST_DEVICE_INLINE float2 ray_aabb_near_far(const vec3 &orig, const vec3 &dir, const vec3 &min_vec, const vec3 &max_vec) {

	// --- Improved Smith's method, see "An Efficient and Robust Ray–Box Intersection Algorithm"
    
	//printf("Here ray_aabb_near_far\n");
	
	float tnear	= -INF_R;
    float tfar  =  INF_R;

	//printf("tnear %f tfar %f\n", tnear, tfar);
	//printf("orig.x %f orig.y %f orig.z %f\n", orig.x, orig.y, orig.z);
	//printf("dir.x %f dir.y %f dir.z %f\n", dir.x, dir.y, dir.z);
	//printf("min_vec.x %f min_vec.y %f min_vec.z %f\n", min_vec.x, min_vec.y, min_vec.z);
	//printf("max_vec.x %f max_vec.y %f max_vec.z %f\n", max_vec.x, max_vec.y, max_vec.z);

	{	// --- Check the x axis
        float tx1, tx2;
        float divx = 1 / dir.x;
        if(divx >= 0.0) {
            tx1 = (min_vec.x - orig.x) * divx;
            tx2 = (max_vec.x - orig.x) * divx;
        } else {
            tx1 = (max_vec.x - orig.x) * divx;
            tx2 = (min_vec.x - orig.x) * divx;
        }
        // --- Note that the sign direction condition garantees that tx1 < tx2
        //tnear = max(tnear, min(tx1, tx2));
        //tfar  = min(tfar,  max(tx1, tx2));    
#ifdef __linux__
        tnear = fmaxf(tnear, tx1);
        tfar  = fminf(tfar,  tx2);
#elif _WIN32 || _WIN64
	#ifdef __CUDACC__
		tnear = max(tnear, tx1);
        tfar  = min(tfar,  tx2);
	#else
		tnear = std::max(tnear, tx1);
        tfar  = std::min(tfar,  tx2);
	#endif
#endif
		//printf("tnear %f tfar %f\n", tnear, tfar);
        // --- DEBUG_MSG("tnear ",tnear," tfar ", tfar,"\n");        
    }
    
	{	// --- Check the y axis
        float ty1, ty2;
        float divy = 1 / dir.y;
        if(divy >= 0.0){
            ty1 = (min_vec.y - orig.y) * divy;
            ty2 = (max_vec.y - orig.y) * divy;
        } else {
            ty1 = (max_vec.y - orig.y) * divy;
            ty2 = (min_vec.y - orig.y) * divy;
        }   
        // --- Note that the sign direction condition garantees that ty1 < ty2
        //tnear = max(tnear, min(ty1, ty2));
        //tfar  = min(tfar,  max(ty1, ty2));
#ifdef __linux__
        tnear = fmaxf(tnear, ty1);
        tfar  = fminf(tfar,  ty2);
#elif _WIN32 || _WIN64
	#ifdef __CUDACC__
        tnear = max(tnear, ty1);
        tfar  = min(tfar,  ty2);
	#else
        tnear = std::max(tnear, ty1);
        tfar  = std::min(tfar,  ty2);
	#endif
#endif
        // --- DEBUG_MSG("tnear ",tnear," tfar ", tfar,"\n");        
    }
    
    {	// --- Check the z axis
        float tz1, tz2;
        float divz = 1 / dir.z;        
        if(divz >= 0.0){
            tz1 = (min_vec.z - orig.z) * divz;
            tz2 = (max_vec.z - orig.z) * divz;
        } else {
            tz1 = (max_vec.z - orig.z) * divz;
            tz2 = (min_vec.z - orig.z) * divz;
        }   
        // --- Note that the sign direction condition garantees that tz1 < tz2       
		//tnear = max(tnear, min(tz1, tz2));
		//tfar  = min(tfar,  max(tz1, tz2));
#ifdef __linux__
        tnear = fmaxf(tnear, tz1);
        tfar  = fminf(tfar,  tz2);
#elif _WIN32 || _WIN64
	#ifdef __CUDACC__
        tnear = max(tnear, tz1);
        tfar  = min(tfar,  tz2);
	#else
        tnear = std::max(tnear, tz1);
        tfar  = std::min(tfar,  tz2);
	#endif
#endif        
        // --- DEBUG_MSG("tnear ",tnear," tfar ", tfar,"\n");
    }
                
    return make_float2(tnear,tfar);                    
}


#undef HOST_DEVICE_INLINE

#endif
