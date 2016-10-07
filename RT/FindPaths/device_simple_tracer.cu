#include <stdio.h>
 
#include "device_simple_tracer.h"
#include "ray_primitive_test.h"

/********************************/
/* SIMPLE TRACER THRUST FUNCTOR */
/********************************/
struct device_compute_geointers_simple {    
    
    // --- Iterators     
    thrust::device_vector<int4>::const_iterator		faces_begin;
    thrust::device_vector<int4>::const_iterator     faces_end;
    thrust::device_vector<float4>::const_iterator   verts_begin;
    
    // --- Constructor
	device_compute_geointers_simple(thrust::device_vector<int4>::const_iterator faces_begin_, thrust::device_vector<int4>::const_iterator faces_end_,
									thrust::device_vector<float4>::const_iterator verts_begin_) : faces_begin(faces_begin_), faces_end(faces_end_),
									verts_begin(verts_begin_) {}
    
    __device__ geointer operator ()(const R0_R1_tuple &r) { 
        
		const float4 R0 = thrust::get<0>(r);  // --- [orig.x, orig.y, orig.z, rho1]
        const float4 R1 = thrust::get<1>(r);  // --- [beta.x, beta.y, beta.z, rho2]
        
        geointer G0 = make_geointer(-1, INF_R, 0, 0);
        
        bool found = false;

		// --- Loop over the primitives
        for(thrust::device_vector<int4>::const_iterator obj_it = faces_begin; obj_it != faces_end; ++obj_it) {

            float3 tuv; // --- Output of the ray-primitive intersection code 
            
			int4 obj = *obj_it;
            
			const vec3 orig = make_vec3(R0.x, R0.y, R0.z);			// --- Origin of the ray
            const vec3 dir	= make_vec3(R1.x, R1.y, R1.z);			// --- Direction of the ray 
            
            // --- Get the vertex indices of the triangles
            const int v1_idx = obj.x;
            const int v2_idx = obj.y;
            const int v3_idx = obj.z;
                    
            const float4 data_v1 = verts_begin[v1_idx];
            const float4 data_v2 = verts_begin[v2_idx];
            const float4 data_v3 = verts_begin[v3_idx];
                    
            // --- Positions of the vertices
            const vec3 v1 = make_vec3(data_v1.x, data_v1.y, data_v1.z);
            const vec3 v2 = make_vec3(data_v2.x, data_v2.y, data_v2.z); 
            const vec3 v3 = make_vec3(data_v3.x, data_v3.y, data_v3.z);

			// --- Finds the intersection between the ray and the current primitive
            tuv = ray_triangle_nearest_with_uv(orig, dir, v1, v2, v3); 
                        
            if (tuv.x > 0 && tuv.x < G0.t) {
                
				// --- Compute object id
                int ID  = obj_it-faces_begin;                 
                G0 = make_geointer(ID, tuv.x, tuv.y, tuv.z);

				//printf("tuv %f %f %f G0.t %f\n", tuv.x, tuv.y, tuv.z, G0.t);
                found = true;             
            }
                
        }
        
        if(!found) G0.t = -1.0;  // --- Sets the t-componet in G0 to make invalid the intersection

		return G0;
    }
};



void device_simple_trace(
    device_ray_list::R0_R1_iterator begin,
    device_ray_list::R0_R1_iterator end,
    thrust::device_vector<int4>::const_iterator     faces_begin,
    thrust::device_vector<int4>::const_iterator     faces_end,
    thrust::device_vector<float4>::const_iterator   verts_begin,    
    thrust::device_vector<geointer>::iterator geointer_begin
    ){
    
    thrust::transform(
        begin,
        end,
        geointer_begin,
        device_compute_geointers_simple(
           faces_begin,
           faces_end,
           verts_begin
        )
    );
    
    
}
