#ifndef __RAY_TRACING_COMPLEX_VEC3_H__
#define __RAY_TRACING_COMPLEX_VEC3_H__


#include <iostream>
#include "cfloat.h"
#include "vec3.h"

#ifdef __CUDACC__
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else 
#define HOST_DEVICE_INLINE  inline
#endif


// this structure isn't 128-bit aligned
struct cvec3{
    cfloat x;
    cfloat y;
    cfloat z;
};



HOST_DEVICE_INLINE
cvec3 make_cvec3(const cfloat &x, const cfloat & y, const cfloat &z){
    cvec3 cv;
    cv.x = x;
    cv.y = y;
    cv.z = z;
    return cv;
}


HOST_DEVICE_INLINE
cvec3 operator+(const cvec3 &a,const cvec3 &b){
    return make_cvec3(a.x+b.x,a.y+b.y,a.z+b.z);
}

HOST_DEVICE_INLINE
cvec3 operator-(const cvec3 &a,const cvec3 &b){
    return make_cvec3(a.x-b.x,a.y-b.y,a.z-b.z);
}


HOST_DEVICE_INLINE
cvec3 operator+(const cvec3 &a){
    return a;
}

HOST_DEVICE_INLINE
cvec3 operator-(const cvec3 &a){
    return make_cvec3(-a.x,-a.y,-a.z); 
}


HOST_DEVICE_INLINE
cvec3 operator*(const cvec3 &a,const cvec3 &b){
    return make_cvec3(a.x*b.x,a.y*b.y,a.z*b.z);
}

HOST_DEVICE_INLINE
cvec3 operator*(const cfloat &a,const cvec3 &b){
    return make_cvec3(a*b.x,a*b.y,a*b.z);
}

HOST_DEVICE_INLINE
cvec3 operator*(const cvec3 &a,const cfloat &b){
    return make_cvec3(a.x*b,a.y*b,a.z*b);
}

HOST_DEVICE_INLINE
cvec3 operator*(const float a,const cvec3 &b){
    return make_cvec3(a*b.x,a*b.y,a*b.z);
}

HOST_DEVICE_INLINE
cvec3 operator*(const cvec3 &a,const float b){
    return make_cvec3(a.x*b,a.y*b,a.z*b);
}

HOST_DEVICE_INLINE
cfloat edot(const cvec3 &a,const cvec3 &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

HOST_DEVICE_INLINE
cfloat edot(const cvec3 &a,const vec3 &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

HOST_DEVICE_INLINE
cfloat edot(const vec3 &a,const cvec3 &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}


HOST_DEVICE_INLINE
cvec3 cross(const vec3 &a,const cvec3 &b){
    return make_cvec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}


HOST_DEVICE_INLINE
cvec3 cross(const cvec3 &a,const vec3 &b){
    return make_cvec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}


HOST_DEVICE_INLINE
cvec3 conjugate(const cvec3 &a){
    return make_cvec3(conjugate(a.x),conjugate(a.y),conjugate(a.z));
}

HOST_DEVICE_INLINE
cfloat dot(const cvec3 &a,const cvec3 &b){
    return edot(a,conjugate(b));
}

HOST_DEVICE_INLINE
float  norm(const cvec3 &a){
    return sqrtf( REAL(a.x)*REAL(a.x) + 
                 REAL(a.y)*REAL(a.y) +
                 REAL(a.z)*REAL(a.z) +
                 IMAG(a.x)*IMAG(a.x) + 
                 IMAG(a.y)*IMAG(a.y) +
                 IMAG(a.z)*IMAG(a.z));
}


HOST_DEVICE_INLINE
cvec3 operator*(const vec3 &a,const cfloat &b){
    return make_cvec3(a.x*b,a.y*b,a.z*b);
}

HOST_DEVICE_INLINE
cvec3 operator*(const cfloat &a,const vec3 &b){
    return make_cvec3(a*b.x,a*b.y,a*b.z);
}

HOST_DEVICE_INLINE
cvec3 operator+(const vec3 &a,const cvec3 &b){
    return make_cvec3(make_cfloat(a.x,0)+b.x,
                      make_cfloat(a.y,0)+b.y,
                      make_cfloat(a.z,0)+b.z);
}

HOST_DEVICE_INLINE
cvec3 operator+(const cvec3 &a,const vec3 &b){
    return make_cvec3(make_cfloat(b.x,0)+a.x,
                      make_cfloat(b.y,0)+a.y,
                      make_cfloat(b.z,0)+a.z);

}




HOST_DEVICE_INLINE
cvec3 operator/(const cvec3 &a,const float b){
    return make_cvec3(a.x/b,a.y/b,a.z/b);
}

HOST_DEVICE_INLINE
cvec3 operator/(const cvec3 &a,const cfloat &b){
    return make_cvec3(a.x/b,a.y/b,a.z/b);
}

std::ostream &operator << (std::ostream &os ,const cvec3 & v);


#define cvec3_zero  (make_cvec3(make_cfloat(0,0),make_cfloat(0,0),make_cfloat(0,0)))

#undef HOST_DEVICE_INLINE

#endif
