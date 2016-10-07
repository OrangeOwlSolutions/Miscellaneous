#include "host_bvh_tracer.h"
#include "ray_primitive_test.h"
#include "aabb.h"
#include "bvh_node.h"

// --- Constructor
host_compute_geo_inters_with_bvh::host_compute_geo_inters_with_bvh(thrust::host_vector<int4>::const_iterator faces_begin_, 
		thrust::host_vector<float4>::const_iterator verts_begin_, thrust::host_vector<int>::const_iterator rmap_begin_, 
		thrust::host_vector<int2>::const_iterator N0_begin_, thrust::host_vector<float4>::const_iterator N1_begin_, 
		thrust::host_vector<float4>::const_iterator N2_begin_, thrust::host_vector<float4>::const_iterator N3_begin_, aabb scn_aabb_) :
        faces_begin(faces_begin_), verts_begin(verts_begin_), rmap_begin(rmap_begin_), N0_begin(N0_begin_), N1_begin(N1_begin_), N2_begin(N2_begin_),
        N3_begin(N3_begin_), scn_aabb(scn_aabb_) { }
    
    
// --- Functions to extract near or far distance 
inline float host_compute_geo_inters_with_bvh::near_fun(const float2 &near_far) { return near_far.x; }
inline float host_compute_geo_inters_with_bvh::far_fun (const float2 &near_far) { return near_far.y; } 
    
geointer  host_compute_geo_inters_with_bvh::operator()(const R0_R1_tuple &r) { 
	const float4 R0  = thrust::get<0>(r);				// --- [orig.x, orig.y, orig.z, rho1]
    const float4 R1  = thrust::get<1>(r);				// --- [beta.x, beta.y, beta.z, rho2]
        
    const vec3  orig = make_vec3(R0.x,R0.y,R0.z);		// --- Origin of the ray
    const vec3  dir  = make_vec3(R1.x,R1.y,R1.z);		// --- Direction of the ray 
    
    return this->operator()(orig,dir);
}
        
    
geointer host_compute_geo_inters_with_bvh::operator()(const vec3 &orig, const vec3 &dir) {
    
	const int STACK_SIZE = 32;           // --- Local stack max size for the persistent threads approach
        
	// --- Finds near and far intersections between a ray and an AABB
    const float2 scn_near_far = ray_aabb_near_far(orig, dir, scn_aabb.min_vec, scn_aabb.max_vec);

	// --- Initialize geometrical intersection with the primitive
    geointer G0 = make_geointer(-1, far_fun(scn_near_far) * (1 + EPS_R), 0, 0);
        
    bool found = false;		// --- Found an intersection ?
        
    int traversal_stack[STACK_SIZE];
        
    int  stack_top = 0;		// --- Current position in traversal stack.
    int current_node_idx;	// --- Current internal node.
 
    if(near_fun(scn_near_far) <= far_fun(scn_near_far) ){ // --- If the ray hit the scene bb 

		// --- Persistent thread stack: start from the root
        traversal_stack[stack_top++] = 0; // --- Push the root
            
        while(stack_top >0) {

			current_node_idx = traversal_stack[--stack_top]; // --- Pop
                
			while(!is_leaf(current_node_idx)) { // --- The current_node_idx will be negative for a leaf
                    
				// --- Determine where to go next

                // --- Fetch the AABB
                const int2   N0 = *(N0_begin + current_node_idx); // --- c0 and c1 indices (indices of the two child nodes)
                const float4 N1 = *(N1_begin + current_node_idx); // --- x, y and z coordinates of the extremals of the child AABBs
                const float4 N2 = *(N2_begin + current_node_idx);
                const float4 N3 = *(N3_begin + current_node_idx);
                    
                // --- Finds the coordinates of the "minimum point" and "maximum point" of the first and second child AABBs
                const vec3 c0_min_vec = get_c0_min_vec_from_N1_N3(N1, N3);
                const vec3 c1_min_vec = get_c1_min_vec_from_N2_N3(N2, N3);
                const vec3 c0_max_vec = get_c0_max_vec_from_N1_N3(N1, N3);
                const vec3 c1_max_vec = get_c1_max_vec_from_N2_N3(N2, N3);
                    
				// --- Intersections to the two children
                const float2 near_far_c0 = ray_aabb_near_far(orig, dir, c0_min_vec, c0_max_vec);
                const float2 near_far_c1 = ray_aabb_near_far(orig, dir, c1_min_vec, c1_max_vec);
                                                                        
                const bool traverse_c0 = (near_fun(near_far_c0) <= far_fun(near_far_c0)) && (far_fun(near_far_c0) >= 0) && (near_fun(near_far_c0) <= G0.t);
                const bool traverse_c1 = (near_fun(near_far_c1) <= far_fun(near_far_c1)) && (far_fun(near_far_c1) >= 0) && (near_fun(near_far_c1) <= G0.t);             
                    
                current_node_idx = get_c0_idx(N0);	// --- c0 index stored as int
                int c1_node_idx  = get_c1_idx(N0);	// --- c1 index stored as int
                    
                // --- Intersection with at least one child => go there.
                if(traverse_c0 != traverse_c1) {
					// --- Try c1
                    if (traverse_c1) current_node_idx = c1_node_idx;
                } else {

					if (!traverse_c0) {
						// --- Neither child was intersected 

						// --- Stack is empty => end
                        if(stack_top == 0) goto end_loop;
                            
                        // --- Stack is not empty => pop
						current_node_idx = traversal_stack[--stack_top];
                        
					} else {
						// --- Both children were intersected => push the farthest one.
                        if (near_fun(near_far_c1) < near_fun(near_far_c0)) {
							traversal_stack[stack_top++] = current_node_idx;                            
                            current_node_idx = c1_node_idx;
                        } else traversal_stack[stack_top++] = c1_node_idx;
                    }
                }       
                    
            } // --- End inner while
                
                
            // --- Current node is a leaf, remember to bit negate it
            // --- Fetch the start and end of the primitive list
            const int2 N0 = *(N0_begin + (~current_node_idx));	// --- N0.x is the start index; N0.y is the post end index
                
            for(thrust::host_vector<int>::const_iterator rmap_it = rmap_begin + N0.x; rmap_it != rmap_begin + N0.y; ++rmap_it) {

				int4 face = *(faces_begin + *rmap_it);

				// --- Get the vertex indices of the triangles
				const uint v1_idx = face.x;
                const uint v2_idx = face.y;
                const uint v3_idx = face.z;
                    
                const float4 data_v1 = verts_begin[v1_idx];
                const float4 data_v2 = verts_begin[v2_idx];
                const float4 data_v3 = verts_begin[v3_idx];
                    
                // --- Vertex positions
                const vec3 v1 = make_vec3(data_v1.x,data_v1.y,data_v1.z);
                const vec3 v2 = make_vec3(data_v2.x,data_v2.y,data_v2.z); 
                const vec3 v3 = make_vec3(data_v3.x,data_v3.y,data_v3.z);
					
				// --- Finds the intersection between the ray and the current primitive
				float3 tuv = ray_triangle_nearest_with_uv(orig, dir, v1, v2, v3);
                    
                if (tuv.x > 0 && tuv.x < G0.t) {
					int ID  = *rmap_it;                 
                    G0 = make_geointer(ID, tuv.x, tuv.y, tuv.z);
                    found = true; 
                }
                        
            } // --- End for
                
        } // --- End outer while
        
    } // --- End if
end_loop:        
    if(!found) G0.t = -1.0;  // --- Sets the t-componet in G0 to make invalid the intersection
        
    return G0;
        
}

/*******************/
/* BVH RAY TRACING */
/*******************/
void host_bvh_trace(host_ray_list::R0_R1_iterator begin, host_ray_list::R0_R1_iterator end, thrust::host_vector<int4>::const_iterator faces_begin, 
					thrust::host_vector<float4>::const_iterator verts_begin, thrust::host_vector<int>::const_iterator rmap_begin,
					thrust::host_vector<int2>::const_iterator N0_begin, thrust::host_vector<float4>::const_iterator N1_begin,
					thrust::host_vector<float4>::const_iterator N2_begin, thrust::host_vector<float4>::const_iterator N3_begin,
					thrust::host_vector<geointer>::iterator geointer_begin, const aabb &scn_aabb) {
    
    thrust::transform(begin, end, geointer_begin, host_compute_geo_inters_with_bvh(faces_begin, verts_begin, rmap_begin, N0_begin, N1_begin, N2_begin, 
		             N3_begin, scn_aabb));
    
}
