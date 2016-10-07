//#ifndef __RAY_TRACING_COMPLEX_MATH_H__
//#define __RAY_TRACING_COMPLEX_MATH_H__
//
//// --- Complex float library
//
////#ifdef __linux__ || _WIN32 || _WIN64
//#ifndef __CUDACC__
//#include <iostream>
//#include <cmath>
//#endif
//
//#include "math_and_physical_constants.h"
//
//#define maximum(a, b) (((a) > (b)) ? (a) : (b))
//#define minimum(a, b) (((a) < (b)) ? (a) : (b))
//
//#ifdef __CUDACC__
//#define ALING_TWO_REAL __builtin_align__(8)
//#define HOST_DEVICE_INLINE __host__ __device__ inline
//#else 
//#ifdef __linux__
//#define ALING_TWO_REAL __attribute__ ((aligned (8)))
//#elif _WIN32 || _WIN64
//#define ALING_TWO_REAL __declspec(align(8))
//#endif
//#define HOST_DEVICE_INLINE  inline
//#endif
//
//
//struct ALING_TWO_REAL cfloat {
//    float x;
//    float y;
//};
//
//
//#ifdef __linux__
//HOST_DEVICE_INLINE
//#endif
//float &REAL(cfloat &);
////float &REAL(cfloat &a){
////    return a.x;
////}
//
//
//HOST_DEVICE_INLINE
//float &IMAG(cfloat &a){
//    return a.y;
//}
//
//// --- Real part of a complex number
//HOST_DEVICE_INLINE const float &REAL(const cfloat &a) { return a.x; }
//
//HOST_DEVICE_INLINE
//const float &IMAG(const cfloat &a){
//    return a.y;
//}
//
//
//HOST_DEVICE_INLINE
//cfloat make_cfloat(float r,float i){
//    cfloat c = {r,i};
//    return c;    
//}
//
//HOST_DEVICE_INLINE
//cfloat operator+(const cfloat &a,const cfloat &b){
//    return make_cfloat(REAL(a)+REAL(b),IMAG(a)+IMAG(b));
//}
//
//
//HOST_DEVICE_INLINE
//cfloat operator-(const cfloat &a,const cfloat &b){
//    return make_cfloat(REAL(a)-REAL(b),IMAG(a)-IMAG(b));
//}
//
//
//HOST_DEVICE_INLINE
//cfloat operator-(const cfloat &a){
//    return make_cfloat(-REAL(a),-IMAG(a));
//}
//
//HOST_DEVICE_INLINE
//cfloat operator+(const cfloat &a){
//    return a;
//}
//
//HOST_DEVICE_INLINE
//cfloat operator*(const cfloat &a,const cfloat &b){
//    return make_cfloat(REAL(a)*REAL(b)-IMAG(a)*IMAG(b),
//                      REAL(a)*IMAG(b)+IMAG(a)*REAL(b));
//}
//
//HOST_DEVICE_INLINE
//cfloat operator*(const float &a,const cfloat &b){
//    return make_cfloat(a*REAL(b),a*IMAG(b));
//}
//
//HOST_DEVICE_INLINE
//cfloat operator*(const cfloat &a,const float &b){
//    return make_cfloat(b*REAL(a),b*IMAG(a));
//}
//
//HOST_DEVICE_INLINE
//cfloat operator/(const cfloat &a,const float &b){
//    return make_cfloat(REAL(a)/b,IMAG(a)/b);
//}
//
//
//HOST_DEVICE_INLINE
//cfloat operator/(const float &a,const cfloat &b){
//    const float den = hypotf(REAL(b),IMAG(b)); 
//    const float sq_den = den*den; 
//    const float m = a/sq_den;
//
//    return make_cfloat(REAL(b)*m,-IMAG(b)*m);
//}
//
//HOST_DEVICE_INLINE
//cfloat operator/(const cfloat &a,const cfloat &b){
//    cfloat c;
//    const float den = hypot(REAL(b),IMAG(b)); 
//    const float sq_den = den*den; 
//    REAL(c) = (REAL(a)*REAL(b)+IMAG(a)*IMAG(b))/sq_den;
//    IMAG(c) = (IMAG(a)*REAL(b)-REAL(a)*IMAG(b))/sq_den; 
//    return c;
//}
//
//
//HOST_DEVICE_INLINE
//cfloat conjugate(const cfloat &a){
//    return make_cfloat(REAL(a),-IMAG(a));
//}
//
//HOST_DEVICE_INLINE
//cfloat expf(const cfloat &a){
//    return expf(REAL(a))*make_cfloat(cosf(IMAG(a)),sinf(IMAG(a)));
//}
//
//HOST_DEVICE_INLINE
//cfloat sqrtf(const cfloat &z){
//    const float rr = REAL(z);
//    const float ii = IMAG(z);
//    const float zero = 0.0f;
//    const float half = 0.5f;
//    
//    float x = fabsf(rr); 
//    float y = fabsf(ii);
//
//    if ((x == zero)&(y == zero)) 
//        return make_cfloat(zero,zero);
//
//    float temp = hypotf(x,y)+x;
//    x = sqrtf(half * temp);
//    y = sqrtf(half * y * y / temp);
//    if (ii >= zero){
//        if (rr >= 0.f) 
//            return make_cfloat(x, y); 
//        else 
//            return make_cfloat(y, x); 
//    }else{
//        if (rr >= 0.f) 
//            return make_cfloat(x, -y); 
//        else 
//            return make_cfloat(y, -x); 
//    }
//}
//
//HOST_DEVICE_INLINE
//float fabsf(const cfloat &a){
//
//    return hypotf(REAL(a),IMAG(a));
//}
//
//
//HOST_DEVICE_INLINE
//bool almost_equal(const cfloat &x, const cfloat &y,const float eps){
//#ifdef __linux__
//	float max_x_y_one = fmaxf(fmaxf(1.0,fabsf(x)),fabsf(y)) ;
//#elif _WIN32 || _WIN64
//	//float max_x_y_one = std::max(std::max(1.0f,fabsf(x)),fabsf(y)) ;
//	float max_x_y_one = maximum(maximum(1.0f, fabsf(x)), fabsf(y)) ;
//#endif
//    return fabsf(x - y) <= max_x_y_one*eps ;
//}
//
//#ifdef __linux__ || _WIN32 || _WIN64
//std::ostream & operator<<(std::ostream &os, const cfloat &cf);
//#endif
//
//#undef ALING_TWO_REAL
//#undef HOST_DEVICE_INLINE
//
//#endif
//
//
#ifndef __RAY_TRACING_COMPLEX_MATH_H__
#define __RAY_TRACING_COMPLEX_MATH_H__

// --- Complex float library

//#ifdef __linux__ || _WIN32 || _WIN64
#ifndef __CUDACC__
#include <iostream>
#include <cmath>
#endif

#include "math_and_physical_constants.h"

#define maximum(a, b) (((a) > (b)) ? (a) : (b))
#define minimum(a, b) (((a) < (b)) ? (a) : (b))

#ifdef __CUDACC__
#define ALING_TWO_REAL __builtin_align__(8)
#define HOST_DEVICE_INLINE __host__ __device__ inline
#else 
#ifdef __linux__
#define ALING_TWO_REAL __attribute__ ((aligned (8)))
#elif _WIN32 || _WIN64
#define ALING_TWO_REAL __declspec(align(8))
#endif
#define HOST_DEVICE_INLINE  inline
#endif


struct ALING_TWO_REAL cfloat {
    float x;
    float y;
};


#ifdef __linux__
HOST_DEVICE_INLINE
#endif
float &REAL(cfloat &);
//float &REAL(cfloat &a){
//    return a.x;
//}


HOST_DEVICE_INLINE
float &IMAG(cfloat &a){
    return a.y;
}

// --- Real part of a complex number
HOST_DEVICE_INLINE const float &REAL(const cfloat &a) { return a.x; }

HOST_DEVICE_INLINE
const float &IMAG(const cfloat &a){
    return a.y;
}


HOST_DEVICE_INLINE
cfloat make_cfloat(float r,float i){
    cfloat c = {r,i};
    return c;    
}

HOST_DEVICE_INLINE
cfloat operator+(const cfloat &a,const cfloat &b){
    return make_cfloat(REAL(a)+REAL(b),IMAG(a)+IMAG(b));
}

HOST_DEVICE_INLINE
cfloat operator+(const float &a,const cfloat &b){
    return make_cfloat(a+REAL(b),IMAG(b));
}

HOST_DEVICE_INLINE
cfloat operator+(const cfloat &a,const float &b){
    return make_cfloat(REAL(a)+b,IMAG(a));
}

HOST_DEVICE_INLINE
cfloat operator-(const cfloat &a,const cfloat &b){
    return make_cfloat(REAL(a)-REAL(b),IMAG(a)-IMAG(b));
}

HOST_DEVICE_INLINE
cfloat operator-(const float &a,const cfloat &b){
    return make_cfloat(a-REAL(b),IMAG(b));
}


HOST_DEVICE_INLINE
cfloat operator-(const cfloat &a,const float &b){
    return make_cfloat(REAL(a)-b,IMAG(a));
}

HOST_DEVICE_INLINE
cfloat operator-(const cfloat &a){
    return make_cfloat(-REAL(a),-IMAG(a));
}

HOST_DEVICE_INLINE
cfloat operator+(const cfloat &a){
    return a;
}

HOST_DEVICE_INLINE
cfloat operator*(const cfloat &a,const cfloat &b){
    return make_cfloat(REAL(a)*REAL(b)-IMAG(a)*IMAG(b),
                      REAL(a)*IMAG(b)+IMAG(a)*REAL(b));
}

HOST_DEVICE_INLINE
cfloat operator*(const float &a,const cfloat &b){
    return make_cfloat(a*REAL(b),a*IMAG(b));
}

HOST_DEVICE_INLINE
cfloat operator*(const cfloat &a,const float &b){
    return make_cfloat(b*REAL(a),b*IMAG(a));
}

HOST_DEVICE_INLINE
cfloat operator/(const cfloat &a,const float &b){
    return make_cfloat(REAL(a)/b,IMAG(a)/b);
}


HOST_DEVICE_INLINE
cfloat operator/(const float &a,const cfloat &b){
    const float den = hypotf(REAL(b),IMAG(b)); 
    const float sq_den = den*den; 
    const float m = a/sq_den;

    return make_cfloat(REAL(b)*m,-IMAG(b)*m);
}

HOST_DEVICE_INLINE
cfloat operator/(const cfloat &a,const cfloat &b){
    cfloat c;
    const float den = hypot(REAL(b),IMAG(b)); 
    const float sq_den = den*den; 
    REAL(c) = (REAL(a)*REAL(b)+IMAG(a)*IMAG(b))/sq_den;
    IMAG(c) = (IMAG(a)*REAL(b)-REAL(a)*IMAG(b))/sq_den; 
    return c;
}


HOST_DEVICE_INLINE
cfloat conjugate(const cfloat &a){
    return make_cfloat(REAL(a),-IMAG(a));
}

HOST_DEVICE_INLINE
cfloat expf(const cfloat &a){
    return expf(REAL(a))*make_cfloat(cosf(IMAG(a)),sinf(IMAG(a)));
}

HOST_DEVICE_INLINE
cfloat sqrtf(const cfloat &z){
    const float rr = REAL(z);
    const float ii = IMAG(z);
    const float zero = 0.0f;
    const float half = 0.5f;
    
    float x = fabsf(rr); 
    float y = fabsf(ii);

    if ((x == zero)&(y == zero)) 
        return make_cfloat(zero,zero);

    float temp = hypotf(x,y)+x;
    x = sqrtf(half * temp);
    y = sqrtf(half * y * y / temp);
    if (ii >= zero){
        if (rr >= 0.f) 
            return make_cfloat(x, y); 
        else 
            return make_cfloat(y, x); 
    }else{
        if (rr >= 0.f) 
            return make_cfloat(x, -y); 
        else 
            return make_cfloat(y, -x); 
    }
}

HOST_DEVICE_INLINE
float fabsf(const cfloat &a){

    return hypotf(REAL(a),IMAG(a));
}


HOST_DEVICE_INLINE
bool almost_equal(const cfloat &x, const cfloat &y,const float eps){
#ifdef __linux__
	float max_x_y_one = fmaxf(fmaxf(1.0,fabsf(x)),fabsf(y)) ;
#elif _WIN32 || _WIN64
	//float max_x_y_one = std::max(std::max(1.0f,fabsf(x)),fabsf(y)) ;
	float max_x_y_one = maximum(maximum(1.0f, fabsf(x)), fabsf(y)) ;
#endif
    return fabsf(x - y) <= max_x_y_one*eps ;
}

#ifdef __linux__ || _WIN32 || _WIN64
std::ostream & operator<<(std::ostream &os, const cfloat &cf);
#endif

#undef ALING_TWO_REAL
#undef HOST_DEVICE_INLINE

#endif


