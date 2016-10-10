// --- Controllare almostEqual se si può fare una versione float

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <iomanip>
#include <utility>
#include <map>
#include <set>
#include <cassert>
#include "vec3.h"
#include "vec4.h"
#include "aabb.h"
#include "bvh.h"
#include "normals_and_curvatures.h"

#include "Mesh.h"

/********/
/* MAIN */
/********/
int main() {
	
	/**************/
	/* INPUT MESH */
	/**************/
	std::vector<mesh> globalInputMesh;					// --- Global PLY, SLT or OBJ BVHmesh to be loaded

	// --- Versione per lavoro ACES
	//mesh partialInputMesh;								// --- Partial PLY, STL or OBJ BVHmesh to be loaded
	////partialInputMesh.filename		= "C:\\Users\\angelo\\Documents\\Project\\Breglia\\cilindro.ply";		
	////partialInputMesh.filename		= "C:\\Users\\angelo\\Documents\\Project\\Breglia\\BVH_Constructor_2\\BVH_Constructor\\cilindro.ply";		
	//partialInputMesh.filename		= "D:\\Project\\Breglia\\BVH_Constructor_3\\BVH_Constructor\\cylinder_rectangular_MIMO.ply";		
	////partialInputMesh.filename		= "C:\\Users\\angelo\\Documents\\Project\\Breglia\\BVH_Constructor_2\\BVH_Constructor\\cube.ply";	
	////partialInputMesh.filename		= "C:\\Users\\angelo\\Documents\\Project\\Breglia\\BVH_Constructor_2\\BVH_Constructor\\two_spheres.ply";	
	//partialInputMesh.format			= PLY;
	//globalInputMesh.push_back(partialInputMesh);			
	////partialInputMesh.filename		= "C:\\Users\\angelo\\Documents\\Project\\Breglia\\BVH_Constructor_2\\BVH_Constructor\\plate_v2.ply";		
	//partialInputMesh.filename		= "D:\\Project\\Breglia\\BVH_Constructor_3\\BVH_Constructor\\plate_MIMO.ply";		
	//globalInputMesh.push_back(partialInputMesh);			

	//// --- Loading the input meshes
	//if(!loadPly(globalInputMesh[0])) { std::cerr << "Error opening file: " << globalInputMesh[0].filename << std::endl; return -1; }
	//if(!loadPly(globalInputMesh[1])) { std::cerr << "Error opening file: " << globalInputMesh[1].filename << std::endl; return -1; }

	mesh partialInputMesh;								// --- Partial PLY, STL or OBJ BVHmesh to be loaded
	partialInputMesh.filename = "D:\\Project\\Breglia\\BVH_Constructor_3\\BVH_Constructor\\nave.ply";
	partialInputMesh.format = PLY;
	globalInputMesh.push_back(partialInputMesh);

	// --- Loading the input meshes
	if (!loadPly(globalInputMesh[0])) { std::cerr << "Error opening file: " << globalInputMesh[0].filename << std::endl; return -1; }

	/***************/
	/* OUTPUT MESH */
	/***************/
	BVHmesh globalOutputMesh;
	
    bool assume_flat			= false;					// --- Assume flat shape flag
    bool write_normal_curvature = true;						// --- Write normal curvature to file flag

	OUTPUT_FORMAT oformat		= NBIN;

	std::ofstream outputFile;
	//outputFile.open("cilindro.nbin", std::ios::binary);
	outputFile.open("nave.nbin", std::ios::binary);
	//outputFile.open("cube.nbin", std::ios::binary);
	//outputFile.open("two_spheres.nbin", std::ios::binary);

	if (write_normal_curvature && !assume_flat) {
        
		for (size_t i = 0; i < globalInputMesh.size() ; i++) {
            
			// --- Compute normals if not available
			if (!globalInputMesh[i].normalsAvailable) compute_normals(globalInputMesh[i].vertices, globalInputMesh[i].faces, globalInputMesh[i].normals);
            
			// --- Compute curvatures if not available
            if (!globalInputMesh[i].curvatureAvailable) compute_curvatures(globalInputMesh[i].vertices, globalInputMesh[i].faces, 
				                                                            globalInputMesh[i].normals, globalInputMesh[i].principalDirections, 
																			globalInputMesh[i].firstCurvatures,    globalInputMesh[i].secondCurvatures);             
        }

    }
	
	/**************************************************/
	/* MERGE THE INPUT MESHES IN A UNIQUE OUTPUT MESH */
	/**************************************************/
	printf("Test\n");
	{	size_t numVertices  = 0;							// --- Total number of BVHmesh vertices
		size_t numFaces		= 0;							// --- Total number of BVHmesh faces
		
		// --- Count the total number of vertices and faces
		for (size_t i = 0; i < globalInputMesh.size(); i++) {
			numVertices += globalInputMesh[i].vertices.size();
			numFaces	+= globalInputMesh[i].faces.size();			
		}
		printf("Total number of vertices = %i; Total number of faces = %i\n", numVertices, numFaces);
		
		globalOutputMesh.vertices.resize(numVertices);
		globalOutputMesh.faces.resize(numFaces);
                
        if(write_normal_curvature && !assume_flat) {
            globalOutputMesh.normals.resize  (numVertices);
            globalOutputMesh.principalDirections.resize (numVertices);
            globalOutputMesh.firstCurvatures.resize(numVertices);
            globalOutputMesh.secondCurvatures.resize(numVertices);
        }
        
		// --- Fill the BVHmesh
		size_t globalVertexIndices = 0;
		size_t globalFaceIndices = 0;
		for (size_t i = 0; i < globalInputMesh.size(); i++) {
    		
			int face_offset = globalVertexIndices;
			for(size_t j = 0; j < globalInputMesh[i].vertices.size(); j++) {
				
				const vec3f v   = globalInputMesh[i].vertices[j];

				globalOutputMesh.vertices[globalVertexIndices] = make_vec3d(v.x, v.y, v.z);
                
                if (write_normal_curvature && !assume_flat) {

                    const vec3f n   = globalInputMesh[i].normals[j];
                    const vec3f x1  = globalInputMesh[i].principalDirections[j];
                    const float firstCurvatures  = globalInputMesh[i].firstCurvatures[j];
                    const float secondCurvatures  = globalInputMesh[i].secondCurvatures[j];
                    
					globalOutputMesh.normals[globalVertexIndices]  = make_vec3d(n.x, n.y, n.z);
                    globalOutputMesh.principalDirections[globalVertexIndices] = make_vec3d(x1.x, x1.y, x1.z);
                    
                    globalOutputMesh.firstCurvatures[globalVertexIndices] = firstCurvatures;
                    globalOutputMesh.secondCurvatures[globalVertexIndices] = secondCurvatures;
                
                }
				globalVertexIndices++;
			}
			
			for (size_t j = 0; j < globalInputMesh[i].faces.size(); j++) {
				
				const vec3i face = globalInputMesh[i].faces[j];

				globalOutputMesh.faces[globalFaceIndices] = make_vec4i(face.a + face_offset, face.b + face_offset, face.c + face_offset, i);
				//printf("%i %i %i %i\n", face.a + face_offset, face.b + face_offset, face.c + face_offset, i);
				globalFaceIndices++; 
			
			}			
		}
	}
	
	globalOutputMesh.bb = compute_aabb(globalOutputMesh.vertices, globalOutputMesh.faces);
	
	std::cerr << "bounding box: " << globalOutputMesh.bb << std::endl;
	
	bvh_builder bbuild(bvh_platform(), globalOutputMesh);
	
	bbuild.build();
	
    if (assume_flat && write_normal_curvature) recompute_flat_verts_normals_and_curvature(globalOutputMesh);
    
	switch (oformat) {
		case SSV : {
			writeTrianglesSSVMollerTrumbore(std::cout, globalOutputMesh);
			bbuild.write_bvh_ssv(std::cout,oformat);
            if (write_normal_curvature) write_normcurv_ssv_mt(std::cout, globalOutputMesh);
            
		} break;
		case BIN : case NBIN : {
			writeTrianglesBinaryMollerTrumbore(outputFile, globalOutputMesh);
			bbuild.write_bvh_ssv(outputFile, oformat);
            if (write_normal_curvature) write_normcurv_bin_mt(outputFile,globalOutputMesh);
		} break;
	};
    
	return 0;
}











