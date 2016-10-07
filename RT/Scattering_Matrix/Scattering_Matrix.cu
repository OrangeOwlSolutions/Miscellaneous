//#include <vector>
//
//#include <omp.h>
//
//#include <algorithm>
//
//#include "host_path_finder.h"
//#include "device_path_finder.h"
//#include "Scattering_Matrix.cuh"
//#include "trace_path.h"
//#include "GLFunctions.h"
//#include "cfloat.h"
//#include "transport_ray_field.h"
//#include "Utilities.cuh"
//#include "device_path_finder.h"
//
//void computeScatteringMatrix(int argc, char** argv, const int num_CPU, const int num_GPU, const size_t max_cnt, 
//	                         std::vector<pathfinder_p> &pathfinder_list, nvmesh &mesh, std::vector<edge> &edge_list, 
//							 std::vector<path> &global_path_list, const int num_TX, const int num_RX, std::vector<TX> &txs, 
//							 std::vector<RX> &rxs, const double lambda, std::vector<cfloat> &H, bool final_iter) {
//
//	int num_coprocs = num_CPU + num_GPU;
//
//	size_t freeMemory, totalMemory;
//	printf("num_coprocs %i num_CPU %i num_GPU %i\n", num_coprocs, num_CPU, num_GPU);
//	//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//	//printf("Prima di device_create_path_finder Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//	for (size_t i = 0; i < num_GPU; i++) 
//		pathfinder_list.push_back(device_create_path_finder(i, mesh.verts, mesh.faces, mesh.normals_and_k1, mesh.x1_and_k2, mesh.remapping_table, 
//															mesh.N0, mesh.N1, mesh.N2, mesh.N3, mesh.bb, max_cnt));
//	//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//	//printf("Dopo di device_create_path_finder Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
// //   system("pause");
//
//    for (size_t i = 0; i < num_CPU; i++)
//        pathfinder_list.push_back(host_create_path_finder(mesh.verts, mesh.faces, mesh.normals_and_k1, mesh.x1_and_k2, mesh.remapping_table, mesh.N0,
//														  mesh.N1, mesh.N2, mesh.N3, mesh.bb, edge_list));
//	
//	//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//	//printf("Prima di device_reset Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//	for (size_t i = 0; i < pathfinder_list.size(); i++) {
//		switch ((pathfinder_list[i]) -> type) {
//			case GPU : { device_reset(pathfinder_list[i], &txs, &rxs); } break;
//            case CPU : { host_reset  (pathfinder_list[i], &txs, &rxs); } break;
//        }
//	}
//	//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//	//printf("Dopo di device_reset Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
// //   system("pause");
//        
//	int nProcessors = omp_get_max_threads();
//	std::cout<<"Maximum number of CPU threads = " << nProcessors << std::endl;
//
//	int cpu_thread_id;
//	omp_set_num_threads(num_coprocs);  // --- Create as many CPU threads as there are CPU threads and CUDA devices
//	#pragma omp parallel private(cpu_thread_id)
//	{                        
//		cpu_thread_id = omp_get_thread_num();
//		//printf("CPU THREAD ID %i\n", cpu_thread_id);
//        if ((pathfinder_list[cpu_thread_id]) -> type == GPU) {
//			//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//			//printf("Prima di device_find_paths Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//			device_find_paths(pathfinder_list[cpu_thread_id]);
//			//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//			//printf("Dopo di device_find_paths Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//			//system("pause");
//            //std::cerr << "GPU: " << cpu_thread_id << std::endl;
//        } else {
//            host_find_paths(pathfinder_list[cpu_thread_id]);
//            //std::cerr << "CPU: " << cpu_thread_id << std::endl;
//        }
//	}
//
//	//printf("Dopo la creazione della lista di pathfinder\n");
//	for(size_t i = 0; i < pathfinder_list.size(); i++) {
//		
//		switch((pathfinder_list[i]) -> type) {
//                case GPU : { 
//					//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//					//printf("Prima di device_append_new_paths( Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//					device_append_new_paths(pathfinder_list[i], global_path_list);  break;
//					//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//					//printf("Dopo di device_append_new_paths( Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//				 //   system("pause");
//				}
//				case CPU : { host_append_new_paths  (pathfinder_list[i], global_path_list); } break;
//        }
//	}
//
//	//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//	//printf("Prima di trace_path Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//	trace_path(global_path_list, mesh, txs, rxs);        
//	//gpuErrchk(cudaMemGetInfo(&freeMemory, &totalMemory));
//	//printf("Dopo di trace_path Total memory: %i Free memory: %i\n", totalMemory / 1024 / 1024, freeMemory / 1024 / 1024);
//	//system("pause");
//
//	size_t new_size = std::remove_if(global_path_list.begin(), global_path_list.end(), is_ray_path_not_good()) - global_path_list.begin();
//    global_path_list.resize(new_size);
//
//	std::sort(global_path_list.begin(), global_path_list.end(), compare_ray_path());
//
//	new_size = std::unique(global_path_list.begin(), global_path_list.end(), is_ray_path_equal()) - global_path_list.begin();
//	global_path_list.resize(new_size);
//
//	//std::cout << "Total number of paths = " << global_path_list.size()  << std::endl;
//	//for(size_t i = 0; i < global_path_list.size(); i++){
//	//	std::cout << global_path_list[i] << std::endl;
//	//}
//        
//	// --- Visualization tool
//	if (final_iter) Visualize(&argc, argv);
//	//Visualize(&argc, argv);
//
//	size_t num_paths = global_path_list.size();
//	for (size_t i = 0; i < num_paths; i++) {
//		
//		size_t num_points = global_path_list[i].intersections.size();
//
//		size_t tx = global_path_list[i].tx;
//		size_t rx = global_path_list[i].rx;
//            
//        std::vector<vec3> positions;
//            
//        cfloat & Htr = H[tx + rx * num_TX];
//               
//        for (size_t j = 0; j < num_points; j++) {
//			vec3 pos = global_path_list[i].intersections[j].pos;
//            positions.push_back(pos);
//        }
//           
//		Htr = Htr + transport_ray_field(positions, txs[tx], rxs[rx], lambda);
//
//		//cfloat test = transport_ray_field(positions, txs[tx], rxs[rx], lambda);
//		//std::cout << Htr.x << " " << Htr.y << "\n";
//		//std::cout << lambda << " " << tx << " " << rx << " " << H[tx + rx * num_TX].x << " " << H[tx + rx * num_TX].y << " " << test.x << " " << test.y << "\n";
//		//for (int p = 0; p < positions.size(); p++) std::cout << positions[p].x << " " << positions[p].y << " " << positions[p].z << "\n";
//
//	}
//
//	global_path_list.resize(0);
//
//	for (size_t i = 0; i < num_GPU; i++) 
//		delete pathfinder_list[i];
//
//	//std::cout << "#num TX  num RX" << std::endl;
//	//std::cout << num_TX << " " << num_RX << std::endl;
//	//for (size_t r = 0; r < num_RX; r++) {
//	//	for (size_t c = 0; c < num_TX; c++) {
//	//		std::cout << REAL(H[c + r * num_TX]) << " " << IMAG(H[c + r * num_TX]) << " ";  
// //       }
// //       std::cout << std::endl;
// //   }
//
//	omp_set_num_threads(1);
//	gpuErrchk(cudaSetDevice(0));
//							 
//}

#include <vector>

#include <omp.h>

#include <algorithm>

#include "host_path_finder.h"
#include "device_path_finder.h"
#include "Scattering_Matrix.cuh"
#include "trace_path.h"
#include "GLFunctions.h"
#include "cfloat.h"
#include "transport_ray_field.h"
#include "device_path_finder.h"
#include "Utilities.cuh"

void computeScatteringMatrix(int argc, char** argv, const int numCPU, const int numGPU, const size_t max_cnt, 
	                         std::vector<pathfinder_p> &pathfinder_list, nvmesh &mesh, std::vector<edge> &edge_list, 
							 std::vector<path> &global_path_list, const int num_TX, const int num_RX, std::vector<TX> &txs, 
							 std::vector<RX> &rxs, const double lambda, std::vector<cfloat> &H, bool final_iter) {

	int numProcessors = numCPU + numGPU;

	for (size_t i = 0; i < pathfinder_list.size(); i++) {
		switch ((pathfinder_list[i]) -> type) {
			case GPU : { device_reset(pathfinder_list[i], &txs, &rxs); } break;
            case CPU : { host_reset  (pathfinder_list[i], &txs, &rxs); } break;
        }
	}
        
	int maxNumProcessors = omp_get_max_threads();
	std::cout << "Maximum number of CPU threads = " << maxNumProcessors << std::endl;

	omp_set_num_threads(numProcessors);  // --- Create as many CPU threads as there are CPU threads and CUDA devices
	#pragma omp parallel 
	{                        
		unsigned int cpu_thread_id = omp_get_thread_num();
		if ((pathfinder_list[cpu_thread_id]) -> type == GPU) {
			device_find_paths(pathfinder_list[cpu_thread_id]);
			printf("CPU thread ID %i is of GPU type\n", cpu_thread_id);
        } else {
            host_find_paths(pathfinder_list[cpu_thread_id]);
			printf("CPU thread ID %i is of CPU type\n", cpu_thread_id);
        }
	}

	for(size_t i = 0; i < pathfinder_list.size(); i++) {
		
		switch((pathfinder_list[i]) -> type) {
                case GPU : { device_append_new_paths(pathfinder_list[i], global_path_list); } break;
                case CPU : { host_append_new_paths  (pathfinder_list[i], global_path_list); } break;
        }
	}

	trace_path(global_path_list, mesh, txs, rxs);        

	size_t new_size = std::remove_if(global_path_list.begin(), global_path_list.end(), is_ray_path_not_good()) - global_path_list.begin();
    global_path_list.resize(new_size);

	std::sort(global_path_list.begin(), global_path_list.end(), compare_ray_path());

	new_size = std::unique(global_path_list.begin(), global_path_list.end(), is_ray_path_equal()) - global_path_list.begin();
	global_path_list.resize(new_size);

	// --- Visualization tool
	if (final_iter) Visualize(&argc, argv);
	//Visualize(&argc, argv);

	size_t num_paths = global_path_list.size();
	for (size_t i = 0; i < num_paths; i++) {
		
		size_t num_points = global_path_list[i].intersections.size();
			
		size_t tx = global_path_list[i].tx;
		size_t rx = global_path_list[i].rx;
            
        std::vector<vec3> positions;
            
        cfloat & Htr = H[tx + rx * num_TX];
               
        for (size_t j = 0; j < num_points; j++) {
			vec3 pos = global_path_list[i].intersections[j].pos;
            positions.push_back(pos);
        }
           
		//Htr = Htr + transport_ray_field(positions, txs[tx], rxs[rx], lambda);
		Htr = Htr + transport_ray_field(global_path_list[i], txs[tx], rxs[rx], mesh, edge_list, lambda);

	}

	omp_set_num_threads(1);
	gpuErrchk(cudaSetDevice(0));
							 
}
