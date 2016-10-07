//#ifndef __RAY_PATH_H__
//#define __RAY_PATH_H__
//
//#include <vector>
//#include <iostream>
//
//#include <vector_types.h>		// --- CUDA definition of float4
//#include <vector_functions.h>	// --- make_float4
//
//#include "vec3.h"
//
//#define EMISSION	0			// --- The point is the starting point of the path, namely, the trasmitter location.
//#define REFLECTION	1			// --- The point is an arrival point of the path, namely, the receiver center.
//#define DETECTION	2			// --- At that point, a reflection occurs, and then the point is on the surface.
//#define DIFFRACTION 3			// --- At that point, a diffraction occurs, and then the point is on the surface edge.
//
///****************************/
///* STRUCT PATH INTERSECTION */
///****************************/
//// --- Contains the different points composing a path. Each point is identified by the type (according to the same convention for inter_point) and the 
////     coordinates in pos. The other fields are used by the electromagnetic tool.
//struct path_inters {
//    
//    int type;
//    vec3 pos;
//    
//    int ID;			// --- Denotes the intersected triangle index in the case of reflection and the edge index in the case of diffraction
//    float u, v;		// --- Honogeneous coordinates of the intersection point in the case of reflection.
//
//};
//
///***************/
///* STRUCT PATH */
///***************/
//// --- Such structure is partially filled by host_path_finder and device_path_finder which update only the fields code, rx, tx, orig and dir_v.
////     Subsequently, the last CPU tracing will find the rays and fill len_dist and intersections. This holds with one exception. For direct rays and
////     rays with single reflection (computed with image theory) which are created on the CPU, the data structure is completely filled. In this case,
////     the rays are not retraced.
//struct path {
//
//    uint4 code;				// --- Code which uniquely identifies each path.
//    int tx;					// --- Transmitter index (contained also in code).
//    int rx;					// --- Receiver index (contained also in code).
//    
//    vec3 orig;				// --- Ray origin (corresponds to the transmitter position).
//    vec3 dir_v;				// --- Ray direction at the transmitter position.
//    
//    float len;				// --- Path length
//    float dist;				// --- Distance of the last path ray from the receiver center.
//     
//    bool good;				// --- True if the path hits the receiver.
//    
//    std::vector<path_inters> intersections;
//							// --- List of all the points (path_inters) of the path comprised between the origin (tx) and the end (the closest point to the
//							//     receiver center.
//
//};
//
//struct compare_ray_path { bool operator() (const path &p1,const path &p2); };
//
//struct is_ray_path_equal { bool operator() (const path &p1, const path &p2); };
//
//struct is_ray_path_not_good { bool operator() (const path &p); };
//
//std::ostream &operator<<(std::ostream &os, const path &p);
//
//#endif 
#ifndef __RAY_PATH_H__
#define __RAY_PATH_H__

#include <vector>
#include <iostream>

#include <vector_types.h>		// --- CUDA definition of float4
#include <vector_functions.h>	// --- make_float4

#include "vec3.h"

#define EMISSION	0			// --- The point is the starting point of the path, namely, the trasmitter location.
#define REFLECTION	1			// --- The point is an arrival point of the path, namely, the receiver center.
#define DETECTION	2			// --- At that point, a reflection occurs, and then the point is on the surface.
#define DIFFRACTION 3			// --- At that point, a diffraction occurs, and then the point is on the surface edge.

/****************************/
/* STRUCT PATH INTERSECTION */
/****************************/
// --- Contains the different points composing a path. Each point is identified by the type (according to the same convention for inter_point) and the 
//     coordinates in pos. The other fields are used by the electromagnetic tool.
struct path_inters {
    
    int type;
    vec3 pos;
    
    int ID;			// --- Denotes the intersected triangle index in the case of reflection and the edge index in the case of diffraction
    float u, v;		// --- Honogeneous coordinates of the intersection point in the case of reflection.

};

/***************/
/* STRUCT PATH */
/***************/
// --- Such structure is partially filled by host_path_finder and device_path_finder which update only the fields code, rx, tx, orig and dir_v.
//     Subsequently, the last CPU tracing will find the rays and fill len_dist and intersections. This holds with one exception. For direct rays and
//     rays with single reflection (computed with image theory) which are created on the CPU, the data structure is completely filled. In this case,
//     the rays are not retraced.
struct path {

    uint4 code;				// --- Code which uniquely identifies each path.
    int tx;					// --- Transmitter index (contained also in code).
    int rx;					// --- Receiver index (contained also in code).
    
    vec3 orig;				// --- Ray origin (corresponds to the transmitter position).
    vec3 dir_v;				// --- Ray direction at the transmitter position.
    
    float len;				// --- Path length
    float dist;				// --- Distance of the last path ray from the receiver center.
     
    bool good;				// --- True if the path hits the receiver.
    
    std::vector<path_inters> intersections;
							// --- List of all the points (path_inters) of the path comprised between the origin (tx) and the end (the closest point to the
							//     receiver center.

    int get_num_refl() const;
    int get_num_diff() const;        
};

struct compare_ray_path { bool operator() (const path &p1,const path &p2); };

struct is_ray_path_equal { bool operator() (const path &p1, const path &p2); };

struct is_ray_path_not_good { bool operator() (const path &p); };

std::ostream &operator<<(std::ostream &os, const path &p);

#endif 
