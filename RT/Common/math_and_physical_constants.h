#ifndef __MATH_AND_PHYSICAL_CONSTANTS_H__
#define __MATH_AND_PHYSICAL_CONSTANTS_H__

/**************************/
/* MATHEMATICAL CONSTANTS */
/**************************/
#ifdef __CUDACC__
    #include <math_constants.h> // CUDA CONSTANTS DEFINES

    #define PI_R         CUDART_PI_F
    #define SQRT_TWO_R   CUDART_SQRT_TWO_F
    #define SQRT_HALH_R  CUDART_SQRT_HALF_F
    //#define INF_R        (1.0f/0.0f) 
    #define INF_R        0x7ff0000000000000		// --- Credit to https://devtalk.nvidia.com/default/topic/456829/how-to-assign-infinity-to-variables-in-cuda-code-/
    //#define INF_R        fCUDART_INF_F
    #define EPS_R        1.19209290E-07F
    #define NAN_R        CUDART_NAN_F

#else
    #include <limits>

    #define PI_R         3.14159265358979323846f
    #define SQRT_TWO_R   1.414213562f
    #define SQRT_HALH_R  0.707106781f
    #define INF_R		 std::numeric_limits<float>::infinity()
    #define EPS_R		 1.19209290E-07F
    #define NAN_R		 std::numeric_limits<float>::signaling_NaN()

#endif

/**********************/
/* PHYSICAL CONSTANTS */
/**********************/
// --- Free space dielectric constant  [F/m]
#define EPSILON_0 8.854187817E-12

// --- Free space magnetic permiability [H/m]
#define MU_0      (4E-7 * PI_R)

// --- Free space characteristic impedence [Ohm]
#define ETA_0     376.7303134749689

#endif
