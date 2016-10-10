#ifndef __MESH_H__
#define __MESH_H__

#include "vec3.h"
#include "vec4.h"
#include "aabb.h"
#include <vector>

/***********************************/
/* STANDARD (STL, PLY OR OBJ) MESH */
/***********************************/
enum INPUT_FORMAT {STL, PLY, OBJ};

struct mesh {
	
	std::string			filename;
	
	std::vector<vec3f>	vertices; 
	std::vector<vec3i>	faces; 
    std::vector<vec3f>	normals; 
    std::vector<vec3f>	principalDirections; 
    std::vector<float>	firstCurvatures; 
    std::vector<float>	secondCurvatures;
    
	INPUT_FORMAT		format;	
    
    bool				normalsAvailable;
    bool				curvatureAvailable;
};

/************/
/* BVH MESH */
/************/
struct BVHmesh {
	
	std::vector<vec3d>	vertices;
	std::vector<vec4i>	faces;	
    std::vector<vec3d>	normals; 
    std::vector<vec3d>	principalDirections; 
    std::vector<double> firstCurvatures; 
    std::vector<double> secondCurvatures; 
	
	aabb bb;
};

// --- Mesh loaders
bool loadPly(mesh &);
bool load_obj(mesh &);

// --- Mesh writers
void writeTrianglesBinaryMollerTrumbore(std::ostream &, const BVHmesh &);
void writeTrianglesSSVMollerTrumbore   (std::ostream &, const BVHmesh &);

// --- Line handlers
bool get_valid_line(std::istream &, std::string &);
std::string trim(const std::string &, const std::string &whitespace = " \t");

#endif
