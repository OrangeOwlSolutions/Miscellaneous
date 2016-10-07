//#include <cassert>
//
//#include <curand_kernel.h>
//
//#include <thrust/device_vector.h>
//#include <thrust/transform.h>
//#include <thrust/functional.h>
//#include <thrust/remove.h>
//#include <thrust/sequence.h>
//#include <thrust/reduce.h>
//#include <thrust/iterator/counting_iterator.h>
//#include <thrust/sort.h>
//#include <thrust/unique.h>
//#include <thrust/count.h>
//#include <thrust/functional.h>
//
//#include "device_path_finder.h"
//#include "aabb.h"
//#include "ray_list.h"
//#include "device_simple_tracer.h"
//#include "reflect_rays_functor.h"
//#include "RXTX.h"
//#include "device_find_det_intersection.h"
//#include "device_bvh_tracer.h"
//#include "Utilities.cuh"
//
//const size_t max_num_rays	= 512 * 512;
//const size_t max_num_paths	= 10000;
//
//#define  MAX_REFLECTION  2					// --- Actually, the maximum number of allowed reflections plus 1
//
////#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
////inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
////{
////   if (code != cudaSuccess) 
////   {
////      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
////      if (abort) exit(code);
////   }
////}
//
///*********************/
///* PATHFINDER STRUCT */
///*********************/
//struct device_pathfinder: public pathfinder {
//
//	thrust::device_vector<float4>		verts;
//	thrust::device_vector<int4>			faces;
//	thrust::device_vector<float4>		normals_and_k1;
//	thrust::device_vector<float4>		x1_and_k2; 
//	thrust::device_vector<int>			remapping_table;
//	thrust::device_vector<int2>			N0;
//	thrust::device_vector<float4>		N1;
//	thrust::device_vector<float4>		N2;
//	thrust::device_vector<float4>		N3;
//	aabb								bb;
//	
//	std::vector<TX>						*txs;				// --- Transmitter locations. Set in the device_reset function.
//	std::vector<RX>						*rxs;				// --- Receiver locations. Set in the device_reset function.
//
//	size_t								device_idx;			// --- Device identification
//
//	std::vector<path>					path_list;
//
//	device_ray_list						global_ray_list;	
//	device_ray_list						starting_ray_list;	// --- Starting positions of the rays
//	device_ray_list						saved_ray_list;		// --- Saved ray list
//  
//	thrust::device_vector<geointer>		ginters;			
//	thrust::device_vector<float>		dinters;			// --- Distance of each ray from the origin
//  
//	thrust::device_vector<curandState>	rnd_states;
//
//	thrust::device_vector<TX>			dtxs;				// --- Detectors (used by the generator of new rays) ???
//  
//	size_t								max_cnt;
//  
//	// --- Pathfinder constructor
//	device_pathfinder(const size_t device_idx_, const thrust::host_vector<float4> &verts_, const thrust::host_vector<int4> &faces_,
//					  const thrust::host_vector<float4> &normals_and_k1_, const thrust::host_vector<float4> &x1_and_k2_, 
//					  const thrust::host_vector<int> &remapping_table_, const thrust::host_vector<int2> &N0_, const thrust::host_vector<float4> &N1_,
//					  const thrust::host_vector<float4> &N2_, const thrust::host_vector<float4> &N3_, const aabb &bb_, const size_t max_cnt_) : 
//					  pathfinder(GPU), device_idx(device_idx_), verts(verts_), faces(faces_), remapping_table(remapping_table_), 
//					  normals_and_k1(normals_and_k1_), x1_and_k2(x1_and_k2_), N0(N0_), N1(N1_), N2(N2_), N3(N3_), bb(bb_),max_cnt(max_cnt_) {}
//  
//	// --- Pathfinder destructor
//	~device_pathfinder() {
//		verts.clear();				verts.shrink_to_fit();
//		faces.clear();				faces.shrink_to_fit();
//		normals_and_k1.clear();		normals_and_k1.shrink_to_fit();
//		x1_and_k2.clear();			x1_and_k2.shrink_to_fit();
//		remapping_table.clear();	remapping_table.shrink_to_fit();
//		N0.clear();					N0.shrink_to_fit();
//		N1.clear();					N1.shrink_to_fit();
//		N2.clear();					N2.shrink_to_fit();
//		N3.clear();					N3.shrink_to_fit();
//
//		ginters.clear();			ginters.shrink_to_fit();			
//		dinters.clear();			dinters.shrink_to_fit();
//  
//		rnd_states.clear();			rnd_states.shrink_to_fit();
//
//		dtxs.clear();				dtxs.shrink_to_fit();
//
//		printf("Destructor *****************\n");
//		
//		global_ray_list.resize(0);
//		starting_ray_list.resize(0);
//		saved_ray_list.resize(0);
//
//	}
//};
//
///****************************/
///* CREATES A NEW PATHFINDER */
///****************************/
//pathfinder_p device_create_path_finder(const size_t device_idx, const thrust::host_vector<float4> &verts, const thrust::host_vector<int4> &faces,
//									   const thrust::host_vector<float4> &normals_and_k1, const thrust::host_vector<float4> &x1_and_k2, 
//									   const thrust::host_vector<int> &remapping_table, const thrust::host_vector<int2> &N0, 
//									   const thrust::host_vector<float4> &N1, const thrust::host_vector<float4> &N2, const thrust::host_vector<float4> &N3,
//									   const aabb &bb, const size_t max_cnt) {
//
//	gpuErrchk(cudaSetDevice(device_idx));
//    device_pathfinder *dpf  = new device_pathfinder(device_idx, verts, faces, normals_and_k1, x1_and_k2, remapping_table, N0, N1, N2, N3, bb, max_cnt);
//    
//    return dpf;
//}
//
///*******************************/
///* IS VALID RAY THRUST FUNCTOR */
///*******************************/
//// --- Length must be larger than 0 and the number of reflections less than the maximum number of reflections
//struct is_valid_ray { __device__ bool operator()(const geointer &g, const uint4 code) const { return g.t >= 0 && code.z < MAX_REFLECTION; } };
//
///**************************************************/
///* THRUST FUNCTOR TO REMOVE NON-INTERSECTING RAYS */
///**************************************************/
//static size_t device_clean_rays(pathfinder_p dpf_) {
//
//    device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
//    
//	gpuErrchk(cudaSetDevice(dpf -> device_idx));
//    
//    // --- map of the rays to remove
//    thrust::device_vector<bool> map(dpf->global_ray_list.size());
//    
//	// --- Sets map with 1 if valid, 0 if non valid
//	//thrust::transform(dpf -> ginters.begin(), dpf -> ginters.end(), map.begin(), is_valid_intersection());
//    thrust::transform(dpf -> ginters.begin(), dpf -> ginters.end(), dpf -> global_ray_list.R2.begin(), map.begin(), is_valid_ray());
//    
//    // --- Removes the rays from the global list
//	thrust::remove_if(dpf -> global_ray_list.begin(), dpf -> global_ray_list.end(), map.begin(), thrust::logical_not<bool>());
//
//    // --- Removes the rays from the starting list
//    thrust::remove_if(dpf -> starting_ray_list.begin(), dpf -> starting_ray_list.end(), map.begin(), thrust::logical_not<bool>());
//
//    // --- Removes the intersections
//	size_t new_size = thrust::remove_if(dpf -> ginters.begin(), dpf -> ginters.end(), map.begin(), thrust::logical_not<bool>()) - dpf->ginters.begin();                          
//                                     
//    return new_size;
//}
//
///*****************************/
///* GENERATE NEW RAYS FUNCTOR */
///*****************************/
//struct device_generate_new_rays{
//    
//    thrust::device_vector<TX>::iterator tx_begin;
//    size_t num_tx;
//    curandState * state_p;
//    
//    device_generate_new_rays(thrust::device_vector<TX>::iterator tx_begin_, size_t num_tx_, curandState * state_p_) : 
//		tx_begin(tx_begin_), num_tx(num_tx_), state_p(state_p_) {}
//    
//    
//    // --- Generates random ray directions
//    __device__ vec3 generate_random_unit_dir(curandState *local_state) {
//        
//		float x1, x2, d;
//        do {
//			// --- x1 and x2 are uniformly distributed in (-1, 1). Generates until the norm of (x1, x2) is larger than 1.
//			x1 = 2 * (curand_uniform(local_state) - 0.5);
//            x2 = 2 * (curand_uniform(local_state) - 0.5);
//            d = x1 * x1 + x2 * x2;
//        } while(d >= 1.0f);
//        
//        vec3 dir_v = normalize(make_vec3(2 * x1 * sqrtf(1 - d), 2 * x2 * sqrtf(1 - d), 1 - 2 * d));
//        
//        return dir_v;
//    }
//
//    // --- operator()
//	__device__ ray_tuple operator()(int idx) {
//        
//        curandState local_state = state_p[idx];
//        
//        uint tx_code	= idx % num_tx;
//        TX tx			= (*(tx_begin + tx_code));
//        
//        vec3 pos		= tx.pos;
//        vec3 dir_v		= generate_random_unit_dir(&local_state);
//        
//        
//        state_p[idx] = local_state;
//        
//        float4 R0 = make_R0_from_orig_and_len(pos, 0);
//        float4 R1 = make_R1_from_dir_v_and_dist(dir_v, 0);
//        uint4  R2 = {1, 0, 0, tx_code};
//        
//        return thrust::make_tuple(R0,R1,R2);
//    }
//
//};
//
///***********************************/
///* GENERATE NEW RAYS HOST FUNCTION */
///***********************************/
//// --- Generating new rays means randomly denerating their directions and initializing R0, R1 and R2.
//static void generate_new_rays(pathfinder_p dpf_, device_ray_list::ray_list_iterator gbegin, device_ray_list::ray_list_iterator gend, 
//	                          device_ray_list::ray_list_iterator sbegin) {
//    
//	device_pathfinder *dpf = reinterpret_cast<device_pathfinder *>(dpf_);
//    
//	cudaSetDevice(dpf->device_idx);
//    
//    size_t num_new_rays = gend - gbegin;
//    //std::cerr << "generate new rays "<< num_new_rays << std::endl;
//    
//    size_t offset = gbegin - dpf->global_ray_list.begin();
//    
//    curandState * state_p =  thrust::raw_pointer_cast(dpf->rnd_states.data())+ offset;
//    
//    thrust::transform(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(0) + num_new_rays, gbegin, device_generate_new_rays(
//													     dpf->dtxs.begin(), dpf->dtxs.size(), state_p));    
//     
//    thrust::copy(gbegin, gend, sbegin); 
//
//}
//
///**************************/
///* SAVE PATHS TO THE HOST */
///**************************/
//void save_paths(pathfinder_p dpf_) {
//    
//	device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
//    
//	cudaSetDevice(dpf->device_idx);
//    
//    host_ray_list h_saved_ray_list = dpf -> saved_ray_list;
//    
//    for (size_t i = 0; i < h_saved_ray_list.size(); i++) {
//        
//		const float4 R0 = h_saved_ray_list.R0[i];
//        const float4 R1 = h_saved_ray_list.R1[i];
//        
//        path p;
//        
//        p.code	= h_saved_ray_list.R2[i];
//        p.orig	= make_vec3(R0.x, R0.y, R0.z);
//        p.dir_v = make_vec3(R1.x, R1.y, R1.z);
//        p.tx	= p.code.w;
//        p.rx	= p.code.y;
//
//        dpf -> path_list.push_back(p);
//        
//    }
//
//}
//
///**********************************************************************/
///* THRUST FUNCTOR SETTING THE NUMBER OF THE HIT DETECTOR FOR THE RAYS */
///**********************************************************************/
//struct set_rx_detector{
//    
//	uint rx_code;
//    
//    set_rx_detector(uint rx_code_) : rx_code(rx_code_) {}
//    
//    __device__ inline uint4 operator()(uint4 v) const { v.y = rx_code; return v; }
//};
//
///****************************************************/
///* THRUST FUNCTOR SETTING THE RAY-DETECTOR DISTANCE */
///****************************************************/
//struct set_detector_distance {
//        
//	__device__ inline float4 operator()(float4 R1, float distance) const { 
//		R1.w = distance;
//        return R1;
//    }
//};
//
///***************************************/
///* THRUST FUNCTOR TO COMPARE RAY PATHS */
///***************************************/
//struct device_compare_ray_path {
//    
//	__device__ bool operator() (const ray_tuple &p1,const  ray_tuple &p2) {
//    
//		const float4 p1_R1  = thrust::get<1>(p1);
//		const float4 p2_R1  = thrust::get<1>(p2);
//    
//		const uint4 p1_R2	= thrust::get<2>(p1);
//		const uint4 p2_R2	= thrust::get<2>(p2);
//      
//		// --- Codes for rays #1 and #2
//		const unsigned int code1[4] = {p1_R2.x, p1_R2.y, p1_R2.z, p1_R2.w};
//		const unsigned int code2[4] = {p2_R2.x, p2_R2.y, p2_R2.z, p2_R2.w};    
//    
//		// --- If the rays are different, returns
//		for (int i = 0 ; i <4; i++) { if(code1[i] != code2[i]) { return code1[i] < code2[i]; }}
//    
//		// --- Here the rays are "the same", so compare the distance from detector
//		return  p1_R1.w < p2_R1.w;
//
//	}
//  
//};
//
///*************************************************/
///* THRUST FUNCTOR TO CHECK IF TWO RAYS ARE EQUAL */
///*************************************************/
//// --- Si può inglobare nel funtore di sopra ???
//struct device_is_ray_path_equal {
//  
//	__device__  bool operator() (const ray_tuple &p1, const ray_tuple &p2) {
//      
//		const uint4 p1_R2  = thrust::get<2>(p1);
//		const uint4 p2_R2  = thrust::get<2>(p2);
//      
//		return p1_R2.x == p2_R2.x && p1_R2.y == p2_R2.y && p1_R2.z == p2_R2.z && p1_R2.w == p2_R2.w;      
//	}
//  
//};
//
///********************************/
///* FIND RAY PATHS ON THE DEVICE */
///********************************/
//void device_find_paths(pathfinder_p dpf_) {
//    
//	device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
//    
//	cudaSetDevice(dpf->device_idx);
//    
//	bool done = false;
//    
//	size_t num_active_rays = dpf -> global_ray_list.size();
//    
//    //std::cout << "num_active_rays = " << num_active_rays << "\n";
//	
//	size_t cnt = 0;
//    
//	do {
//    
//        // --- The number of active rays must be the maximum possible
//        assert(dpf -> global_ray_list.size() == num_active_rays);
//        
//        // --- Trace the rays 
//        //device_simple_trace(dpf->global_ray_list.begin_R0_R1(), dpf->global_ray_list.end_R0_R1(), dpf->faces.begin(), dpf->faces.end(), 
//		      //                dpf->verts.begin(), dpf->ginters.begin());
//		device_bvh_trace(dpf -> global_ray_list.begin_R0_R1(), dpf -> global_ray_list.end_R0_R1(), dpf -> faces.begin(), dpf -> verts.begin(),
//						 dpf -> remapping_table.begin(), dpf -> N0.begin(), dpf -> N1.begin(), dpf -> N2.begin(), dpf -> N3.begin(), dpf -> ginters.begin(),
//						 dpf->bb);
//		
//		//for (int k = 0; k < dpf -> ginters.size(); k++) {
//		//	geointer test = dpf -> ginters[k];
//		//	if (test.ID != -1) std::cout << "Intersection found with triangle # " << test.ID << " with distance " << test.t << "\n";
//		//}
//		//gpuErrchk(cudaPeekAtLastError());
//		//gpuErrchk(cudaDeviceSynchronize());
//        
//        // --- Finds the intersections with the detectors. Initially, the intersections involve direct ray paths from tx to rx, while in the 
//		//     subsequent iterations the intersections involve reflectedy ray paths from the intersection point to rx.
//        for (size_t i = 0; i < dpf -> rxs -> size(); i++) {
//            
//			RX rx = (*dpf -> rxs)[i];
//            
//			using namespace thrust::placeholders;
//            
//            // --- Finds the intersections between the rays and the detectors
//			device_find_det_inters(dpf -> global_ray_list.begin_R0_R1(), dpf -> global_ray_list.begin_R0_R1() + num_active_rays, dpf -> dinters.begin(), 
//				                   dpf -> ginters.begin(), rx);
//	            
//			//for (int k = 0; k < dpf -> ginters.size(); k++) {
//			//	float test = dpf -> dinters[k];
//			//	if (test > 0.f) std::cout << "Distance from detector = " << test << " with distance " << test << "\n";
//			//}
//            //std::cout << "num_active_rays " << num_active_rays << "\n";
//			
//			// --- Counts the number of paths intersecting the current receiver
//			size_t num_new_path = thrust::count_if(dpf -> dinters.begin(), dpf -> dinters.end(), _1 > 0);
//            
//            //std::cerr << "num_new_path" << num_new_path << std::endl;
//                        
//            size_t num_saved_rays = dpf -> saved_ray_list.size();
//            //std::cerr << "num_saved_rays" << num_saved_rays << std::endl;
//
//            // --- If the number of paths exceeds the maximum number of handable paths, then the ray list is reduced.
//			if(num_new_path > 0 && num_new_path + num_saved_rays >= max_num_paths) {
//                
//                // --- Sorts the rays according to their path length to the detector
//				thrust::sort(dpf -> saved_ray_list.begin(), dpf -> saved_ray_list.end(), device_compare_ray_path());
//                
//                // --- Removes duplicates and computes new size list
//				size_t new_size = thrust::unique(dpf -> saved_ray_list.begin(), dpf -> saved_ray_list.end(), device_is_ray_path_equal()) - 
//					                             dpf -> saved_ray_list.begin();
//        
//                // --- Resizes the ray list
//				dpf -> saved_ray_list.resize(new_size);
//                num_saved_rays = new_size;        
//            }
//            //std::cerr << "num_saved_rays after" << num_saved_rays << std::endl;
//           
//            // --- If the number of paths DOES NOT exceed the maximum number of handable paths, then the ray list is reduced.
//            if (num_new_path + num_saved_rays < max_num_paths) {          
//                
//                if (num_new_path > 0) {
//					
//					// --- Resizes the ray list
//					dpf -> saved_ray_list.resize(num_saved_rays + num_new_path);
//                    
//                    // --- Copies R0 and R1 from starting ray list to saved ray list for only the rays hitting the detector
//					thrust::copy_if(dpf -> starting_ray_list.begin_R0_R1(), dpf -> starting_ray_list.end_R0_R1(), dpf -> dinters.begin(),
//									dpf -> saved_ray_list.begin_R0_R1() + num_saved_rays, _1 > 0);
//                    
//                    // --- Copies R2 from global ray list to saved ray list for only the rays hitting the detector
//                    thrust::copy_if(dpf -> global_ray_list.R2.begin(), dpf -> global_ray_list.R2.end(), dpf -> dinters.begin(), 
//						            dpf -> saved_ray_list.R2.begin() + num_saved_rays, _1 > 0);
//                    
//                    // --- Sets the number of hit detector
//					thrust::transform(dpf -> saved_ray_list.R2.begin() + num_saved_rays, dpf -> saved_ray_list.R2.end(), 
//						              dpf -> saved_ray_list.R2.begin() + num_saved_rays, set_rx_detector(i));
//                    
//                    // --- Removes the rays that do not hit the detector
//					thrust::remove_if(dpf -> dinters.begin(), dpf -> dinters.end(), _1 <= 0) ;
//                    
//                    // --- Set the distance to the detector for each ray
//					thrust::transform(dpf -> saved_ray_list.R1.begin() + num_saved_rays, dpf -> saved_ray_list.R1.end(), dpf -> dinters.begin(),              
//									  dpf -> saved_ray_list.R1.begin() + num_saved_rays, set_detector_distance());
//                }
//                
//            } else std::cerr << "buffer full" << std::endl;
//            
//        } 
//        
//        // --- Removes the rays that do not have intersections
//        num_active_rays = device_clean_rays(dpf) ;
//
//        // --- Reflect rays
//        thrust::transform(dpf -> global_ray_list.begin(), dpf -> global_ray_list.begin() + num_active_rays, dpf -> ginters.begin(), 
//						  dpf -> global_ray_list.begin(), make_compute_reflected_ray(dpf -> verts.begin(), dpf -> faces.begin(), 
//						  dpf -> normals_and_k1.begin()));
//        
//        // --- Adds new rays from the transmitters
//        size_t num_available_slots = dpf->global_ray_list.size()-num_active_rays;
//        generate_new_rays(dpf, dpf -> global_ray_list.begin() + num_active_rays, dpf -> global_ray_list.end(),		
//						  dpf -> starting_ray_list.begin() + num_active_rays);
//        num_active_rays = dpf -> global_ray_list.size();
//
//        if(cnt == dpf->max_cnt) done = true;
//
//		cnt++;
//    } while(!done);
//
//    save_paths(dpf);
//}
//
///**********************/
///* SETUP_RND FUNCTION */
///**********************/
//struct setup_rnd{
//    
//    curandState * state_p;
//    
//    setup_rnd(curandState *state_p_) : state_p(state_p_) {}
//        
//    __device__ void operator ()(int id) { curand_init(1234, id, 0, state_p + id); }
//};
//
///**********************************************/
///* INITIALIZE RANDOM STATE FOR RAY GENERATION */
///**********************************************/
//static void init_random_state(pathfinder_p dpf_) {
//    
//	device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
//    
//	cudaSetDevice(dpf->device_idx);
//    
//	dpf->rnd_states.resize(max_num_rays);
//    
//    thrust::for_each(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(0)+max_num_rays, setup_rnd(thrust::raw_pointer_cast(dpf -> 
//		             rnd_states.data())));    
//}
//
///**********************************/
///* DEVICE INITIALIZATION FUNCTION */
///**********************************/
//void device_reset(pathfinder_p dpf_, std::vector<TX> *txs, std::vector<RX> *rxs) {
//    
//	device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
//    
//	cudaSetDevice(dpf -> device_idx);							// --- Sets the device id
//    
//    dpf -> txs = txs;											// --- Sets the transmitter locations
//    dpf -> rxs = rxs;											// --- Sets the receiver locations
//    
//    dpf -> global_ray_list.resize(max_num_rays);				// --- Initializes the device ray list with the maximum number of rays
//    dpf -> starting_ray_list.resize(max_num_rays);    
//    dpf -> saved_ray_list.resize(0);
//    dpf -> saved_ray_list.reserve(max_num_paths);
//    //std::cerr << "resizinf " << std::endl;
//
//	dpf -> ginters.resize(max_num_rays);
//    dpf -> dinters.resize(max_num_rays);
//    
//    init_random_state(dpf);
//    //std::cerr << "init random " << std::endl;
//
//	dpf -> dtxs.resize(txs -> size());							// --- Copy the detectors to the device
//    thrust::copy(txs -> begin(), txs -> end(), dpf -> dtxs.begin());
//	//std::cerr << "copy " << std::endl;
//    
//	// --- Generate rays
//	generate_new_rays(dpf, dpf -> global_ray_list.begin(), dpf -> global_ray_list.end(), dpf -> starting_ray_list.begin());
//    //std::cerr << "generating new rays " << std::endl;
//
//    dpf -> path_list.resize(0);
//    
//    //for(size_t i =0; i < txs->size();i++) {
//    //    std::cerr << "i:  " << (*txs)[i].pos.x << " " <<  (*txs)[i].pos.y << " " << (*txs)[i].pos.z << std::endl;
//    //}
//    
//}
//
///************************************/
///* DEVICE APPEND NEW PATHS FUNCTION */
///************************************/
//void device_append_new_paths(pathfinder_p dpf_, std::vector<path> &global_path_list) {
//    
//	device_pathfinder *dpf = reinterpret_cast<device_pathfinder *>(dpf_);
//    
//	cudaSetDevice(dpf->device_idx);
//
//    global_path_list.insert(global_path_list.end(), dpf->path_list.begin(), dpf->path_list.end());
//}
//
//
//
//void device_destroy_path_finder(pathfinder_p dpf_){
//    device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
//    cudaSetDevice(dpf->device_idx);    
//    delete dpf;
//}
//
//
#include <cassert>

#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/count.h>
#include <thrust/functional.h>

#include "device_path_finder.h"
#include "aabb.h"
#include "device_simple_tracer.h"
#include "reflect_rays_functor.h"
#include "RXTX.h"
#include "device_find_det_intersection.h"
#include "device_bvh_tracer.h"
#include "Utilities.cuh"

const size_t max_num_rays	= 512 * 512;
const size_t max_num_paths	= 10000;

#define  MAX_REFLECTION  2					// --- Actually, the maximum number of allowed reflections plus 1

device_pathfinder::device_pathfinder(const size_t device_idx_, const thrust::host_vector<float4> &verts_, const thrust::host_vector<int4> &faces_,
  				   const thrust::host_vector<float4> &normals_and_k1_, const thrust::host_vector<float4> &x1_and_k2_, 
					  const thrust::host_vector<int> &remapping_table_, const thrust::host_vector<int2> &N0_, const thrust::host_vector<float4> &N1_,
					  const thrust::host_vector<float4> &N2_, const thrust::host_vector<float4> &N3_, const aabb &bb_, const size_t max_cnt_) : 
					  pathfinder(GPU), device_idx(device_idx_), verts(verts_), faces(faces_), remapping_table(remapping_table_), 
					  normals_and_k1(normals_and_k1_), x1_and_k2(x1_and_k2_), N0(N0_), N1(N1_), N2(N2_), N3(N3_), bb(bb_),max_cnt(max_cnt_) {}

/***********************/
/* CREATE A PATHFINDER */
/***********************/
pathfinder_p device_create_path_finder(const size_t device_idx, const thrust::host_vector<float4> &verts, const thrust::host_vector<int4> &faces,
									   const thrust::host_vector<float4> &normals_and_k1, const thrust::host_vector<float4> &x1_and_k2, 
									   const thrust::host_vector<int> &remapping_table, const thrust::host_vector<int2> &N0, 
									   const thrust::host_vector<float4> &N1, const thrust::host_vector<float4> &N2, const thrust::host_vector<float4> &N3,
									   const aabb &bb, const size_t max_cnt) {

	gpuErrchk(cudaSetDevice(device_idx));
    device_pathfinder * dpf  = new device_pathfinder(device_idx, verts, faces, normals_and_k1, x1_and_k2, remapping_table, N0, N1, N2, N3, bb, max_cnt);
    
    return dpf;
}

/*******************************/
/* IS VALID RAY THRUST FUNCTOR */
/*******************************/
// --- Length must be larger than 0 and the number of reflections less than the maximum number of reflections
struct is_valid_ray { __device__ bool operator()(const geointer &g, const uint4 code) const { return g.t >= 0 && code.z < MAX_REFLECTION; } };

/**************************************************/
/* THRUST FUNCTOR TO REMOVE NON-INTERSECTING RAYS */
/**************************************************/
static size_t device_clean_rays(pathfinder_p dpf_) {

    device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf -> device_idx));
    
    // --- map of the rays to remove
    thrust::device_vector<bool> map(dpf->global_ray_list.size());
    
	// --- Sets map with 1 if valid, 0 if non valid
	//thrust::transform(dpf -> ginters.begin(), dpf -> ginters.end(), map.begin(), is_valid_intersection());
    thrust::transform(dpf -> ginters.begin(), dpf -> ginters.end(), dpf -> global_ray_list.R2.begin(), map.begin(), is_valid_ray());
    
    // --- Removes the rays from the global list
	thrust::remove_if(dpf -> global_ray_list.begin(), dpf -> global_ray_list.end(), map.begin(), thrust::logical_not<bool>());

    // --- Removes the rays from the starting list
    thrust::remove_if(dpf -> starting_ray_list.begin(), dpf -> starting_ray_list.end(), map.begin(), thrust::logical_not<bool>());

    // --- Removes the intersections
	size_t new_size = thrust::remove_if(dpf -> ginters.begin(), dpf -> ginters.end(), map.begin(), thrust::logical_not<bool>()) - dpf->ginters.begin();                          
                                     
    return new_size;
}

/*****************************/
/* GENERATE NEW RAYS FUNCTOR */
/*****************************/
struct device_generate_new_rays{
    
    thrust::device_vector<TX>::iterator tx_begin;
    size_t num_tx;
    curandState * state_p;
    
    device_generate_new_rays(thrust::device_vector<TX>::iterator tx_begin_, size_t num_tx_, curandState * state_p_) : 
		tx_begin(tx_begin_), num_tx(num_tx_), state_p(state_p_) {}
    
    
    // --- Generates random ray directions
    __device__ vec3 generate_random_unit_dir(curandState *local_state) {
        
		float x1, x2, d;
        do {
			// --- x1 and x2 are uniformly distributed in (-1, 1). Generates until the norm of (x1, x2) is larger than 1.
			x1 = 2 * (curand_uniform(local_state) - 0.5);
            x2 = 2 * (curand_uniform(local_state) - 0.5);
            d = x1 * x1 + x2 * x2;
        } while(d >= 1.0f);
        
        vec3 dir_v = normalize(make_vec3(2 * x1 * sqrtf(1 - d), 2 * x2 * sqrtf(1 - d), 1 - 2 * d));
        
        return dir_v;
    }

    // --- operator()
	__device__ ray_tuple operator()(int idx) {
        
        curandState local_state = state_p[idx];
        
        uint tx_code	= idx % num_tx;
        TX tx			= (*(tx_begin + tx_code));
        
        vec3 pos		= tx.pos;
        vec3 dir_v		= generate_random_unit_dir(&local_state);
        
        
        state_p[idx] = local_state;
        
        float4 R0 = make_R0_from_orig_and_len(pos, 0);
        float4 R1 = make_R1_from_dir_v_and_dist(dir_v, 0);
        uint4  R2 = {1, 0, 0, tx_code};
        
        return thrust::make_tuple(R0,R1,R2);
    }

};

/***********************************/
/* GENERATE NEW RAYS HOST FUNCTION */
/***********************************/
// --- Generating new rays means randomly denerating their directions and initializing R0, R1 and R2.
static void generate_new_rays(pathfinder_p dpf_, device_ray_list::ray_list_iterator gbegin, device_ray_list::ray_list_iterator gend, 
	                          device_ray_list::ray_list_iterator sbegin) {
    
	device_pathfinder *dpf = reinterpret_cast<device_pathfinder *>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf->device_idx));
    
    size_t num_new_rays = gend - gbegin;
    //std::cerr << "generate new rays "<< num_new_rays << std::endl;
    
    size_t offset = gbegin - dpf->global_ray_list.begin();
    
    curandState * state_p =  thrust::raw_pointer_cast(dpf->rnd_states.data())+ offset;
    
    thrust::transform(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(0) + num_new_rays, gbegin, device_generate_new_rays(
													     dpf->dtxs.begin(), dpf->dtxs.size(), state_p));    
     
    thrust::copy(gbegin, gend, sbegin); 

}

/**************************/
/* SAVE PATHS TO THE HOST */
/**************************/
void save_paths(pathfinder_p dpf_) {
    
	device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf->device_idx));
    
    host_ray_list h_saved_ray_list = dpf -> saved_ray_list;
    
    for (size_t i = 0; i < h_saved_ray_list.size(); i++) {
        
		const float4 R0 = h_saved_ray_list.R0[i];
        const float4 R1 = h_saved_ray_list.R1[i];
        
        path p;
        
        p.code	= h_saved_ray_list.R2[i];
        p.orig	= make_vec3(R0.x, R0.y, R0.z);
        p.dir_v = make_vec3(R1.x, R1.y, R1.z);
        p.tx	= p.code.w;
        p.rx	= p.code.y;

        dpf -> path_list.push_back(p);
        
    }

}

/**********************************************************************/
/* THRUST FUNCTOR SETTING THE NUMBER OF THE HIT DETECTOR FOR THE RAYS */
/**********************************************************************/
struct set_rx_detector{
    
	uint rx_code;
    
    set_rx_detector(uint rx_code_) : rx_code(rx_code_) {}
    
    __device__ inline uint4 operator()(uint4 v) const { v.y = rx_code; return v; }
};

/****************************************************/
/* THRUST FUNCTOR SETTING THE RAY-DETECTOR DISTANCE */
/****************************************************/
struct set_detector_distance {
        
	__device__ inline float4 operator()(float4 R1, float distance) const { 
		R1.w = distance;
        return R1;
    }
};

/***************************************/
/* THRUST FUNCTOR TO COMPARE RAY PATHS */
/***************************************/
struct device_compare_ray_path {
    
	__device__ bool operator() (const ray_tuple &p1,const  ray_tuple &p2) {
    
		const float4 p1_R1  = thrust::get<1>(p1);
		const float4 p2_R1  = thrust::get<1>(p2);
    
		const uint4 p1_R2	= thrust::get<2>(p1);
		const uint4 p2_R2	= thrust::get<2>(p2);
      
		// --- Codes for rays #1 and #2
		const unsigned int code1[4] = {p1_R2.x, p1_R2.y, p1_R2.z, p1_R2.w};
		const unsigned int code2[4] = {p2_R2.x, p2_R2.y, p2_R2.z, p2_R2.w};    
    
		// --- If the rays are different, returns
		for (int i = 0 ; i <4; i++) { if(code1[i] != code2[i]) { return code1[i] < code2[i]; }}
    
		// --- Here the rays are "the same", so compare the distance from detector
		return  p1_R1.w < p2_R1.w;

	}
  
};

/*************************************************/
/* THRUST FUNCTOR TO CHECK IF TWO RAYS ARE EQUAL */
/*************************************************/
// --- Si può inglobare nel funtore di sopra ???
struct device_is_ray_path_equal {
  
	__device__  bool operator() (const ray_tuple &p1, const ray_tuple &p2) {
      
		const uint4 p1_R2  = thrust::get<2>(p1);
		const uint4 p2_R2  = thrust::get<2>(p2);
      
		return p1_R2.x == p2_R2.x && p1_R2.y == p2_R2.y && p1_R2.z == p2_R2.z && p1_R2.w == p2_R2.w;      
	}
  
};

/********************************/
/* FIND RAY PATHS ON THE DEVICE */
/********************************/
void device_find_paths(pathfinder_p dpf_) {
    
	device_pathfinder *dpf = reinterpret_cast<device_pathfinder *>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf -> device_idx));
    
	bool done = false;
    
	size_t num_active_rays = dpf -> global_ray_list.size();
    
    //std::cout << "num_active_rays = " << num_active_rays << "\n";
	
	size_t cnt = 0;
    
	do {
    
        // --- The number of active rays must be the maximum possible
        assert(dpf -> global_ray_list.size() == num_active_rays);
        
        // --- Trace the rays 
        //device_simple_trace(dpf->global_ray_list.begin_R0_R1(), dpf->global_ray_list.end_R0_R1(), dpf->faces.begin(), dpf->faces.end(), 
		      //                dpf->verts.begin(), dpf->ginters.begin());
		device_bvh_trace(dpf -> global_ray_list.begin_R0_R1(), dpf -> global_ray_list.end_R0_R1(), dpf -> faces.begin(), dpf -> verts.begin(),
						 dpf -> remapping_table.begin(), dpf -> N0.begin(), dpf -> N1.begin(), dpf -> N2.begin(), dpf -> N3.begin(), dpf -> ginters.begin(),
						 dpf->bb);
		
		//for (int k = 0; k < dpf -> ginters.size(); k++) {
		//	geointer test = dpf -> ginters[k];
		//	if (test.ID != -1) std::cout << "Intersection found with triangle # " << test.ID << " with distance " << test.t << "\n";
		//}
		//gpuErrchk(cudaPeekAtLastError());
		//gpuErrchk(cudaDeviceSynchronize());
        
        // --- Finds the intersections with the detectors. Initially, the intersections involve direct ray paths from tx to rx, while in the 
		//     subsequent iterations the intersections involve reflectedy ray paths from the intersection point to rx.
        for (size_t i = 0; i < dpf -> rxs -> size(); i++) {
            
			RX rx = (*dpf -> rxs)[i];
            
			using namespace thrust::placeholders;
            
            // --- Finds the intersections between the rays and the detectors
			device_find_det_inters(dpf -> global_ray_list.begin_R0_R1(), dpf -> global_ray_list.begin_R0_R1() + num_active_rays, dpf -> dinters.begin(), 
				                   dpf -> ginters.begin(), rx);
	            
			//for (int k = 0; k < dpf -> ginters.size(); k++) {
			//	float test = dpf -> dinters[k];
			//	if (test > 0.f) std::cout << "Distance from detector = " << test << " with distance " << test << "\n";
			//}
            //std::cout << "num_active_rays " << num_active_rays << "\n";
			
			// --- Counts the number of paths intersecting the current receiver
			size_t num_new_path = thrust::count_if(dpf -> dinters.begin(), dpf -> dinters.end(), _1 > 0);
            
            //std::cerr << "num_new_path" << num_new_path << std::endl;
                        
            size_t num_saved_rays = dpf -> saved_ray_list.size();
            //std::cerr << "num_saved_rays" << num_saved_rays << std::endl;

            // --- If the number of paths exceeds the maximum number of handable paths, then the ray list is reduced.
			if(num_new_path > 0 && num_new_path + num_saved_rays >= max_num_paths) {
                
                // --- Sorts the rays according to their path length to the detector
				thrust::sort(dpf -> saved_ray_list.begin(), dpf -> saved_ray_list.end(), device_compare_ray_path());
                
                // --- Removes duplicates and computes new size list
				size_t new_size = thrust::unique(dpf -> saved_ray_list.begin(), dpf -> saved_ray_list.end(), device_is_ray_path_equal()) - 
					                             dpf -> saved_ray_list.begin();
        
                // --- Resizes the ray list
				dpf -> saved_ray_list.resize(new_size);
                num_saved_rays = new_size;        
            }
            //std::cerr << "num_saved_rays after" << num_saved_rays << std::endl;
           
            // --- If the number of paths DOES NOT exceed the maximum number of handable paths, then the ray list is reduced.
            if (num_new_path + num_saved_rays < max_num_paths) {          
                
                if (num_new_path > 0) {
					
					// --- Resizes the ray list
					dpf -> saved_ray_list.resize(num_saved_rays + num_new_path);
                    
                    // --- Copies R0 and R1 from starting ray list to saved ray list for only the rays hitting the detector
					thrust::copy_if(dpf -> starting_ray_list.begin_R0_R1(), dpf -> starting_ray_list.end_R0_R1(), dpf -> dinters.begin(),
									dpf -> saved_ray_list.begin_R0_R1() + num_saved_rays, _1 > 0);
                    
                    // --- Copies R2 from global ray list to saved ray list for only the rays hitting the detector
                    thrust::copy_if(dpf -> global_ray_list.R2.begin(), dpf -> global_ray_list.R2.end(), dpf -> dinters.begin(), 
						            dpf -> saved_ray_list.R2.begin() + num_saved_rays, _1 > 0);
                    
                    // --- Sets the number of hit detector
					thrust::transform(dpf -> saved_ray_list.R2.begin() + num_saved_rays, dpf -> saved_ray_list.R2.end(), 
						              dpf -> saved_ray_list.R2.begin() + num_saved_rays, set_rx_detector(i));
                    
                    // --- Removes the rays that do not hit the detector
					thrust::remove_if(dpf -> dinters.begin(), dpf -> dinters.end(), _1 <= 0) ;
                    
                    // --- Set the distance to the detector for each ray
					thrust::transform(dpf -> saved_ray_list.R1.begin() + num_saved_rays, dpf -> saved_ray_list.R1.end(), dpf -> dinters.begin(),              
									  dpf -> saved_ray_list.R1.begin() + num_saved_rays, set_detector_distance());
                }
                
            } else std::cerr << "buffer full" << std::endl;
            
        } 
        
        // --- Removes the rays that do not have intersections
        num_active_rays = device_clean_rays(dpf) ;

        // --- Reflect rays
        thrust::transform(dpf -> global_ray_list.begin(), dpf -> global_ray_list.begin() + num_active_rays, dpf -> ginters.begin(), 
						  dpf -> global_ray_list.begin(), make_compute_reflected_ray(dpf -> verts.begin(), dpf -> faces.begin(), 
						  dpf -> normals_and_k1.begin()));
        
        // --- Adds new rays from the transmitters
        size_t num_available_slots = dpf->global_ray_list.size()-num_active_rays;
        generate_new_rays(dpf, dpf -> global_ray_list.begin() + num_active_rays, dpf -> global_ray_list.end(),		
						  dpf -> starting_ray_list.begin() + num_active_rays);
        num_active_rays = dpf -> global_ray_list.size();

        if(cnt == dpf->max_cnt) done = true;

		cnt++;
    } while(!done);

    save_paths(dpf);
}

/**********************/
/* SETUP_RND FUNCTION */
/**********************/
struct setup_rnd{
    
    curandState * state_p;
    
    setup_rnd(curandState *state_p_) : state_p(state_p_) {}
        
    __device__ void operator ()(int id) { curand_init(1234, id, 0, state_p + id); }
};

/**********************************************/
/* INITIALIZE RANDOM STATE FOR RAY GENERATION */
/**********************************************/
static void init_random_state(pathfinder_p dpf_) {
    
	device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf->device_idx));
    
	dpf->rnd_states.resize(max_num_rays);
    
    thrust::for_each(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(0)+max_num_rays, setup_rnd(thrust::raw_pointer_cast(dpf -> 
		             rnd_states.data())));    
}

/**********************************/
/* DEVICE INITIALIZATION FUNCTION */
/**********************************/
void device_reset(pathfinder_p dpf_, std::vector<TX> *txs, std::vector<RX> *rxs) {
    
	device_pathfinder *dpf = reinterpret_cast<device_pathfinder*>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf -> device_idx));				// --- Sets the device id
    
    dpf -> txs = txs;											// --- Sets the transmitter locations
    dpf -> rxs = rxs;											// --- Sets the receiver locations
    
    dpf -> global_ray_list.resize(max_num_rays);				// --- Initializes the device ray list with the maximum number of rays
    dpf -> starting_ray_list.resize(max_num_rays);    
    dpf -> saved_ray_list.resize(0);
    dpf -> saved_ray_list.reserve(max_num_paths);
    //std::cerr << "resizinf " << std::endl;

	dpf -> ginters.resize(max_num_rays);
    dpf -> dinters.resize(max_num_rays);
    
    init_random_state(dpf);
    //std::cerr << "init random " << std::endl;

	dpf -> dtxs.resize(txs -> size());							// --- Copy the detectors to the device
    thrust::copy(txs -> begin(), txs -> end(), dpf -> dtxs.begin());
	//std::cerr << "copy " << std::endl;
    
	// --- Generate rays
	generate_new_rays(dpf, dpf -> global_ray_list.begin(), dpf -> global_ray_list.end(), dpf -> starting_ray_list.begin());
    //std::cerr << "generating new rays " << std::endl;

    dpf -> path_list.resize(0);
    
    //for(size_t i =0; i < txs->size();i++) {
    //    std::cerr << "i:  " << (*txs)[i].pos.x << " " <<  (*txs)[i].pos.y << " " << (*txs)[i].pos.z << std::endl;
    //}
    
}

/************************************/
/* DEVICE APPEND NEW PATHS FUNCTION */
/************************************/
void device_append_new_paths(pathfinder_p dpf_, std::vector<path> &global_path_list) {
    
	device_pathfinder *dpf = reinterpret_cast<device_pathfinder *>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf->device_idx));

    global_path_list.insert(global_path_list.end(), dpf->path_list.begin(), dpf->path_list.end());
}

/*****************************/
/* DEVICE DESTROY PATHFINDER */
/*****************************/
void device_destroy_path_finder(pathfinder_p dpf_) {
    
	device_pathfinder *dpf = reinterpret_cast<device_pathfinder *>(dpf_);
    
	gpuErrchk(cudaSetDevice(dpf->device_idx));    
    
	delete dpf;
}


