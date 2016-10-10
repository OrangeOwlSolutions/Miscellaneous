#ifndef __NORMALS_AND_CURVATURES_H__
#define __NORMALS_AND_CURVATURES_H__

#include <vector>

#include "vec3.h"
#include "Mesh.h"

void compute_normals(
    const std::vector<vec3f> &vertices, //list of vertex
    const std::vector<vec3i> &faces, //list of faces
    std::vector<vec3f> &normals //list of normals              
);


void compute_curvatures(
    const std::vector<vec3f> &vertices, //list of vertex
    const std::vector<vec3i> &faces, //list of faces
    const std::vector<vec3f> &normals, //list of normals              
    std::vector<vec3f> &principalDirections, //list of normals              
    std::vector<float> &k1list, //list of normals              
    std::vector<float> &k2list //list of normals              
);


void recompute_flat_verts_normals_and_curvature(BVHmesh &globalOutputMesh);

void write_normcurv_ssv_mt(std::ostream &os,const  BVHmesh &partialInputMesh);

void write_normcurv_bin_mt(std::ostream &os,const  BVHmesh &partialInputMesh);


#endif
