#include <cstdlib>
#include <ctime>

#include <thrust/host_vector.h>

#include "host_path_finder.h"
#include "aabb.h"
#include "RXTX.h"
#include "host_bvh_tracer.h"
#include "crc24.h"
#include "ray_primitive_test.h"

/*********************/
/* PATHFINDER STRUCT */
/*********************/
struct host_pathfinder: public pathfinder {
	
	const thrust::host_vector<float4>   &verts;
	const thrust::host_vector<int4>     &faces;
	const thrust::host_vector<float4>   &normals_and_k1;
	const thrust::host_vector<float4>   &x1_and_k2;
	const thrust::host_vector<int>      &remapping_table;
	const thrust::host_vector<int2>     &N0;
	const thrust::host_vector<float4>   &N1;
	const thrust::host_vector<float4>   &N2;
	const thrust::host_vector<float4>   &N3;
	const aabb							bb;
	
	const std::vector<edge>             &edge_list;

	std::vector<TX>						*txs;					// --- Transmitter locations. Set in the host_reset function.					
	std::vector<RX>						*rxs;					// --- Receiver locations. Set in the host_reset function.
  
	std::vector<path>					path_list;

	// --- Pathfinder constructor
	host_pathfinder(const thrust::host_vector<float4> &verts_, const thrust::host_vector<int4> &faces_, const thrust::host_vector<float4> &normals_and_k1_,
				    const thrust::host_vector<float4> &x1_and_k2_, const thrust::host_vector<int> &remapping_table_, const thrust::host_vector<int2> &N0_,
					const thrust::host_vector<float4> &N1_, const thrust::host_vector<float4> &N2_, const thrust::host_vector<float4> &N3_, const aabb &bb_,
					const std::vector<edge> &edge_list_) : pathfinder(CPU), verts(verts_), faces(faces_), remapping_table(remapping_table_),
					normals_and_k1(normals_and_k1_), x1_and_k2(x1_and_k2_), N0(N0_), N1(N1_), N2(N2_), N3(N3_), bb(bb_), edge_list(edge_list_) {}
  
};

/****************************/
/* CREATES A NEW PATHFINDER */
/****************************/
pathfinder_p host_create_path_finder(const thrust::host_vector<float4> &verts, const thrust::host_vector<int4> &faces, 
	                                 const thrust::host_vector<float4> &normals_and_k1, const thrust::host_vector<float4> &x1_and_k2, 
									 const thrust::host_vector<int> &remapping_table, const thrust::host_vector<int2> &N0, 
									 const thrust::host_vector<float4> &N1, const thrust::host_vector<float4> &N2, const thrust::host_vector<float4> &N3, 
									 const aabb &bb, const std::vector<edge> &edge_list) {
    
	host_pathfinder * hpf  = new host_pathfinder(verts,faces,normals_and_k1,x1_and_k2, remapping_table,N0,N1,N2,N3,bb,edge_list);
    
    return hpf;
} 

/********************************/
/* HOST INITIALIZATION FUNCTION */
/********************************/
void host_reset(pathfinder_p  hpf_, std::vector<TX> *txs, std::vector<RX> *rxs) { 
	
	host_pathfinder *hpf = reinterpret_cast<host_pathfinder*>(hpf_);
    
	srand (time(NULL));
    
	hpf->txs = txs;							// --- Sets the transmitter locations
    hpf->rxs = rxs;							// --- Sets the receiver locations
    
	hpf->path_list.resize(0);
}

/******************************/
/* FIND RAY PATHS ON THE HOST */
/******************************/
void host_find_paths(pathfinder_p hpf_){
    
	host_pathfinder *hpf = reinterpret_cast<host_pathfinder*>(hpf_);
    
	// --- Searches for direct paths
    host_compute_geo_inters_with_bvh tracer(hpf -> faces.begin(), hpf -> verts.begin(), hpf -> remapping_table.begin(), hpf -> N0.begin(), 
		                                    hpf -> N1.begin(), hpf -> N2.begin(), hpf -> N3.begin(), hpf -> bb);
    
    

    
    
    //find direct path
    for(size_t itx =0; itx < hpf->txs->size(); itx++){
        
        const TX tx =  (*hpf->txs)[itx];        
        for(size_t irx =0; irx < hpf->rxs->size(); irx++){
            
            const RX rx =  (*hpf->rxs)[irx];
            
            const vec3 orig  = tx.pos;
            const vec3 dir_v = normalize(rx.pos - tx.pos);
            const float s = dist(rx.pos,tx.pos);
            geointer G0 = tracer(orig,dir_v);
            
            if(G0.t < 0 || G0.t >  s){
                path p;
        
                p.code =  make_R2(1, irx, 0, 0,  itx);
                p.orig  = orig;
                p.dir_v = dir_v;
                
                p.tx = itx;
                p.rx = irx;

                hpf->path_list.push_back(p);
            }
        }    
    }
    
    
    //find reflected paths        
    for(size_t itx =0; itx < hpf->txs->size(); itx++){        
        const TX tx =  (*hpf->txs)[itx];        
        //seleziona un triangolo a caso
        for(size_t i =0; i < hpf->faces.size(); i++){
            size_t fdx = i;//rand()%hpf->faces.size();
            int4 face  = hpf->faces[fdx];
            
            const float4 a  = hpf->verts[face.x];
            const float4 b  = hpf->verts[face.y];
            const float4 c  = hpf->verts[face.z];
            
            const vec3 va   = make_vec3(a.x,a.y,a.z);
            const vec3 vb   = make_vec3(b.x,b.y,b.z);
            const vec3 vc   = make_vec3(c.x,c.y,c.z);
            
            
            const vec3 Q = (va+vb+vc)*(1.0/3.0);
            const vec3 normal = normalize(cross(vb-va,vc-va));
            
            const vec3 RTX = tx.pos - 2*dot(tx.pos-Q,normal)*normal;


            
            for(size_t irx =0; irx < hpf->rxs->size(); irx++){            
                const RX rx =  (*hpf->rxs)[irx];
                const vec3 dir_v = normalize(rx.pos-RTX);
                
                
                float3 tuv = ray_triangle_nearest_with_uv(
                    RTX, dir_v,
                    va,
                    vb,
                    vc
                );
                
                
                if(tuv.x > 0){
                    const vec3 RP = RTX+ dir_v*tuv.x;
                    
                    const vec3 dir_d_v = normalize(RP - tx.pos);
                    const vec3 dir_r_v = normalize(RP - rx.pos);
                    
                    geointer G0 = tracer(tx.pos,dir_d_v);
                    geointer G1 = tracer(rx.pos,dir_r_v);
                    
                    
                    if((G0.t > 0 && G0.ID ==fdx) && (G1.t > 0 && G1.ID ==fdx)){
                        path p;
                
                        p.code =  make_R2(crc_24(1,face.w), irx, 0, 1,  itx); 
                        p.orig  = tx.pos;
                        p.dir_v = dir_d_v;
                        
                        
                        p.tx = itx;
                        p.rx = irx;

                        hpf->path_list.push_back(p);
                    }

                }
            }
        }
    }
    
    
    if(hpf->edge_list.size() > 0){
        //find diffracted rays
        for(size_t itx =0; itx < hpf->txs->size(); itx++){        
            const TX tx =  (*hpf->txs)[itx];        
            //seleziona uno spigolo a caso
            for(size_t i =0; i < hpf->edge_list.size(); i++){
                const edge &ed = hpf->edge_list[i];
                        
                for(size_t irx =0; irx < hpf->rxs->size(); irx++){            
                    const RX rx =  (*hpf->rxs)[irx];                           
                    
                    vec3 diff_point;
                    
                    if(get_diffraction_point(tx.pos,rx.pos,ed,diff_point)){
                        
                        
                        const vec3 dir_d_v = normalize(diff_point - tx.pos);
                        const vec3 dir_r_v = normalize(diff_point - rx.pos);
                        
                        const float dist_tx_dp = dist(diff_point ,tx.pos);
                        const float dist_rx_dp = dist(diff_point ,rx.pos);
                        
                        geointer G0 = tracer(tx.pos,dir_d_v);
                        geointer G1 = tracer(rx.pos,dir_r_v);
                        
                        
                        if((G0.t < 0  || G0.t >= dist_tx_dp) && (G1.t < 0  || G1.t >= dist_rx_dp)){
                            path p;
    
                            p.code =  make_R2(crc_24(1,ed.ID), 
                                              irx, 
                                              1, 
                                              0, 
                                              itx);                            
                            p.orig  = tx.pos;
                            p.dir_v = dir_d_v;                            
                            
                            p.tx = itx;
                            p.rx = irx;
                            
                            p.dist = 0;
                            p.len  = dist_tx_dp+dist_rx_dp;
                            
                            p.good = true;
                            
                            {   //add the first point
                                path_inters pinter;
                                
                                pinter.type = EMISSION;
                                pinter.pos = tx.pos;                               
                                
                                p.intersections.push_back(pinter);
                            }

                            {   //add the diffraction
                                path_inters pinter;
                                
                                pinter.type = DIFFRACTION;
                                pinter.pos = diff_point;                               
                                pinter.ID  = i;
                                p.intersections.push_back(pinter);
                            }

                            {   //add the first point
                                path_inters pinter;
                                
                                pinter.type = DETECTION;
                                pinter.pos = rx.pos;                               
                                
                                p.intersections.push_back(pinter);
                            }
                            
                            hpf->path_list.push_back(p);
                            
                        }                    
                    }
                }
            }
        }
    }
}

/**********************************/
/* HOST APPEND NEW PATHS FUNCTION */
/**********************************/
void host_append_new_paths(pathfinder_p hpf_, std::vector<path> &global_path_list) {    
    
	host_pathfinder *hpf = reinterpret_cast<host_pathfinder*>(hpf_);

    global_path_list.insert(global_path_list.end(), hpf->path_list.begin(), hpf->path_list.end());
}

void host_destroy_path_finder(pathfinder_p hpf) { delete hpf; }


