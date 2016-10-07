#include <cassert>

#include "trace_path.h" 
#include "host_simple_tracer.h" 
#include "host_bvh_tracer.h"
#include "ray_primitive_test.h"
#include "crc24.h"

/****************************/
/* TRACE PATH HOST FUNCTION */
/****************************/
void trace_path(std::vector<path> &global_path_list, const nvmesh &mesh, const std::vector<TX> &txs, const std::vector<RX> &rxs) {
    
    host_compute_geo_inters_with_bvh tracer(mesh.faces.begin(), mesh.verts.begin(), mesh.remapping_table.begin(), mesh.N0.begin(), mesh.N1.begin(),
											mesh.N2.begin(), mesh.N3.begin(), mesh.bb);
    
    for(size_t i = 0; i < global_path_list.size(); i++) {
        
        if (global_path_list[i].intersections.size() == 0) {
            
            vec3 orig	= global_path_list[i].orig;					// --- Ray origin
            vec3 dir_v	= global_path_list[i].dir_v;				// --- Ray direction
            const uint4 final_code = global_path_list[i].code;
                       
            uint4 current_code = make_R2(1, 0, 0, 0,  global_path_list[i].tx);
            
            RX rx = rxs[global_path_list[i].rx];

            {   // --- Add the starting point
                path_inters pinter;
                pinter.type = 0;
                pinter.pos  = orig;
                global_path_list[i].intersections.push_back(pinter);
            }
            
            bool done = false;
            
            do { 
 
                if(current_code.x == final_code.x && current_code.z == final_code.z && current_code.w == final_code.w) {

                    const vec3  center = rx.pos;                    
                    const float radius = rx.radius;
        
                    current_code.y = global_path_list[i].rx;
                    
                    float t = ray_sphere_median(orig, dir_v, center, radius) ; 

                    //if(t < 0) { std::cerr  << "no detector found" << std::endl; break; }
                    if(t < 0) break;
                    
                    vec3 pos = orig + t * dir_v;
                    
                    {   // --- Add the intersection nearest to the detector
                        path_inters pinter;
                        pinter.type = 2;
                        pinter.pos  = pos;                                       
                        global_path_list[i].intersections.push_back(pinter);
                        global_path_list[i].dist = dist(center, pos);                        
                    }
                    
                    done = true;
                    
                } else {
                    
					geointer G0		= tracer(orig, dir_v);                    
                    const int4 face = mesh.faces[G0.ID];

                    //if (G0.t <= 0) {
                    //    std::cerr  << "no object is intersected" << std::endl;
                    //    break;
                    //}
                    if (G0.t <= 0) break;

                    vec3 new_orig = orig + G0.t * dir_v;
                    
                    vec3 new_dir_v = orig + G0.t * dir_v;

                    // --- Unit normal at Qr
                    vec3 normal_v = mesh.get_normal(G0.ID, G0.u, G0.v);
                    
                    const vec3 dir_r_v = normalize(dir_v - 2 * dot(normal_v, dir_v) * normal_v);

                    const float s_detach = 0.001f; // --- A small value
 
                    orig = new_orig + dir_r_v * s_detach;
                    dir_v = dir_r_v;
                    
                    current_code.z += 1;
                    
                    current_code.x = crc_24(current_code.x, face.w); 

                    {   // --- Add a new point
                        path_inters pinter;
                        pinter.type = 1;
                        pinter.pos  = orig;                    
                        pinter.ID	= G0.ID;
                        pinter.u	= G0.u;
                        pinter.v	= G0.v;
                        
                        global_path_list[i].intersections.push_back(pinter);
                    }
                    
                }
                
            } while(!done);
            
            global_path_list[i].good = done;

        }
    }
    
    
}
