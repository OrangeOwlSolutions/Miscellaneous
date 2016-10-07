#ifndef __VEC3_H__
#define __VEC3_H__

#include <iostream>

#include "math_and_physical_constants.h"

#include <cuda.h>
#include <cuda_runtime.h>

#define maximum(a, b) (((a) > (b)) ? (a) : (b))
#define minimum(a, b) (((a) < (b)) ? (a) : (b))

#ifdef __CUDACC__
#define HOST_DEVICE_INLINE __host__ __device__ inline
typedef float3 vec3;
#else
#include <cmath>
#define HOST_DEVICE_INLINE  inline
//typedef float3 vec3;
struct vec3 {
	float x,y,z;
};
//class vec3{
//public:
//	float x,y,z;
//};
#endif


HOST_DEVICE_INLINE
vec3 make_vec3(float x,float y, float z){
    vec3 v = {x,y,z};
    return v;
}

HOST_DEVICE_INLINE
vec3 operator+(const vec3 &a,const vec3 &b){
    return make_vec3(a.x+b.x,a.y+b.y,a.z+b.z);
}

HOST_DEVICE_INLINE
vec3 operator-(const vec3 &a,const vec3 &b){
    return make_vec3(a.x-b.x,a.y-b.y,a.z-b.z);
}

HOST_DEVICE_INLINE
vec3 operator+(const vec3 &a){
    return a;
}

HOST_DEVICE_INLINE
vec3 operator-(const vec3 &a){
    return make_vec3(-a.x,-a.y,-a.z);
}


HOST_DEVICE_INLINE
vec3 operator*(const vec3 &a,const vec3 &b){
    return make_vec3(a.x*b.x,a.y*b.y,a.z*b.z);
}

HOST_DEVICE_INLINE
vec3 operator*(const float a,const vec3 &b){
    return make_vec3(a*b.x,a*b.y,a*b.z);
}

HOST_DEVICE_INLINE
vec3 operator*(const vec3 &a,const float b){
    return make_vec3(a.x*b,a.y*b,a.z*b);
}

HOST_DEVICE_INLINE
vec3 operator/(const vec3 &a,const vec3 &b){
    return make_vec3(a.x/b.x,a.y/b.y,a.z/b.z);
}

// --- Scalar product between two vec3's
HOST_DEVICE_INLINE float dot(const vec3 &a, const vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// --- Vector product
HOST_DEVICE_INLINE vec3 cross(const vec3 &a, const vec3 &b) { return make_vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }

// --- Normalize a vec3
HOST_DEVICE_INLINE
vec3 normalize(const vec3 &a){
    #ifdef __CUDACC__
    return a*rsqrtf(dot(a,a));
    #else
    return a*(1.0f/sqrtf(dot(a,a)));
    #endif
}

HOST_DEVICE_INLINE
float norm(const vec3 &a){
    return sqrtf(dot(a,a));
}

// --- Distance between two vec3 vectors
HOST_DEVICE_INLINE float dist(const vec3 &a, const vec3 &b) { return norm(b - a); }


HOST_DEVICE_INLINE
bool almost_equal(float x, float y){
#ifdef __linux__
	float max_x_y_one = fmaxf(fmaxf(1.0f,fabsf(x)),fabsf(y)) ;
#elif _WIN32 || _WIN64
	//float max_x_y_one = std::max(std::max(1.0f,fabsf(x)),fabsf(y)) ;
	float max_x_y_one = maximum(maximum(1.0f, fabsf(x)), fabsf(y)) ;
#endif
    return fabsf(x - y) <= max_x_y_one*EPS_R ;
}

HOST_DEVICE_INLINE
bool almost_equal(const vec3 &a, const vec3 &b){
    return almost_equal(a.x,b.x) &&
           almost_equal(a.y,b.y) &&
           almost_equal(a.z,b.z);
}

HOST_DEVICE_INLINE
float fminf(const vec3 &a){
#ifdef __linux__
	return fminf(fminf(a.x,a.y),a.z);
#elif _WIN32 || _WIN64
    //return std::min(std::min(a.x,a.y),a.z);
    return minimum(minimum(a.x, a.y), a.z);
#endif
}

HOST_DEVICE_INLINE
float fmaxf(const vec3 &a){
#ifdef __linux__
    return fmaxf(fmaxf(a.x,a.y),a.z);
#elif _WIN32 || _WIN64
	//return std::max(std::max(a.x,a.y),a.z);
	return maximum(maximum(a.x, a.y), a.z);
#endif
}

/********************************************************************************/
/* TAKES TWO POINTS AND RETURNS THE POINT HAVING THE MINIMUM OF EACH COORDINATE */
/********************************************************************************/
HOST_DEVICE_INLINE
vec3 cwisemin(const vec3 &a,const vec3 &b){
#ifdef __linux__
    return make_vec3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
#elif _WIN32 || _WIN64
    //return make_vec3(std::min(a.x,b.x), std::min(a.y,b.y), std::min(a.z,b.z));
    return make_vec3(minimum(a.x,b.x), minimum(a.y,b.y), minimum(a.z,b.z));
#endif
}

/********************************************************************************/
/* TAKES TWO POINTS AND RETURNS THE POINT HAVING THE MAXIMUM OF EACH COORDINATE */
/********************************************************************************/
HOST_DEVICE_INLINE
vec3 cwisemax(const vec3 &a,const vec3 &b){
#ifdef __linux__
    return make_vec3(fmaxf(a.x,b.x),fmaxf(a.y,b.y),fmaxf(a.z,b.z));
#elif _WIN32 || _WIN64
    //return make_vec3(std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z));
    return make_vec3(maximum(a.x,b.x), maximum(a.y,b.y), maximum(a.z,b.z));
#endif
}



std::ostream &operator << (std::ostream &os ,const vec3 & v);

#undef HOST_DEVICE_INLINE

#endif
