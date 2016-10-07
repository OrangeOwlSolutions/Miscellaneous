#include <thrust/transform.h> 

#include "device_find_det_intersection.h"
#include "ray_primitive_test.h"
 
/*****************************************************************************************/
/* THRUST FUNCTOR FOR THE DETERMINATION OF THE INTERSECTION BETWEEN A RAY AND A DETECTOR */
/*****************************************************************************************/
struct find_det_intersection{

	RX rx;
    
    // --- Constructor
	find_det_intersection(const RX &rx_) : rx(rx_) {}
    
    //__device__ float operator()(R0_R1_tuple r, const geointer &ginter) const {
    __device__ float operator()(R0_R1_tuple r, geointer ginter) const {
    
        const float4 R0 = thrust::get<0>(r);  // --- [orig.x, orig.y, orig.z, k1]
        const float4 R1 = thrust::get<1>(r);  // --- [dir_v.x, dir_v.y, dir_v.z, k2]
        
        const vec3  orig  = make_vec3(R0.x, R0.y, R0.z);
        const vec3  dir_v = make_vec3(R1.x, R1.y, R1.z);

        const vec3  center = rx.pos;
        const float radius = rx.radius;
        
        // --- Finds a "median" distance between the ray origin and the detector
		float t = ray_sphere_median(orig, dir_v, center, radius); 
        
		//printf("Distance between detector at %f and ray origin %f\n", rx.pos.x, t);
		
		if(t > 0 && (t < ginter.t || ginter.t <0)) { // --- If the ray hits the detector and no object hides it
            //printf("Origin = %f %f %f; Direction = %f %f %f; Center = %f %f %f; Distance = %f\n", orig.x, orig.y, orig.z, dir_v.x, dir_v.y, dir_v.z,
				        //                                                                          center.x, center.y, center.z, t);
			return dist(orig + t * dir_v, center);
        }
        
        return -1.0f;
    }
    
};
 
/******************************************************/
/* FIND THE INTERSECTION BETWEEN A RAY AND A DETECTOR */
/******************************************************/
void device_find_det_inters(device_ray_list::R0_R1_iterator begin, device_ray_list::R0_R1_iterator end, thrust::device_vector<float>::iterator dinter_begin,
							thrust::device_vector<geointer>::iterator geointer_begin, const RX &rx) {
    
	thrust::transform(begin, end, geointer_begin, dinter_begin, find_det_intersection(rx));
}
