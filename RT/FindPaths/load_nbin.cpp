////#include <thrust\host_vector.h>
//#include "load_nbin.h"
//
//#ifdef _WIN32 || _WIN64
//typedef unsigned int uint;
//#endif
//
//#include <fstream>
//#include <sstream>
//#include <cassert>
//#include <string>
//
//#include "utils.h"
//
///********************/
///* GET MESH NORMALS */
///********************/
//vec3 nvmesh::get_normal(const int ID, float u, float  v) const {
//    
//    // --- Get the vertex indices of the triangle
//    const int4 face   = faces[ID];
//    const uint v1_idx = face.x;
//    const uint v2_idx = face.y;
//    const uint v3_idx = face.z;
//       
//    // --- Get the normals
//    const float4 data_l2_v1 = normals_and_k1[v1_idx];
//    const float4 data_l2_v2 = normals_and_k1[v2_idx];
//    const float4 data_l2_v3 = normals_and_k1[v3_idx];
//    
//    // --- Normals in each vertex
//    const vec3 normal_v1  = make_vec3(data_l2_v1.x, data_l2_v1.y, data_l2_v1.z);
//    const vec3 normal_v2  = make_vec3(data_l2_v2.x, data_l2_v2.y, data_l2_v2.z);
//    const vec3 normal_v3  = make_vec3(data_l2_v3.x, data_l2_v3.y, data_l2_v3.z);                
//    
//    // --- Perform interpolation
//    const vec3 normal_v = normalize(normal_v1 + u * (normal_v2 - normal_v1) + v * (normal_v3 - normal_v1));
//    
//    return normal_v;
//                                                    
//}
//
///**************/
///* PRINT MESH */
///**************/
//std::ostream &operator << (std::ostream &os, const nvmesh &m){
//    
//    os << "#aabb" << std::endl;
//    os << m.bb <<  std::endl;
//        
//    os << "#vertices" << std::endl;
//    for(size_t i = 0; i < m.verts.size(); i++){
//        os << m.verts[i].x << " " << m.verts[i].y << " " << m.verts[i].z << " " << m.verts[i].w << std::endl;
//    }
//    
//    os << "#faces" << std::endl;
//    for(size_t i = 0; i < m.faces.size(); i++){
//        os << m.faces[i].x << " " << m.faces[i].y << " " << m.faces[i].z << " " << m.faces[i].w << std::endl;
//    }
//    
//    os << "#remapping table" << std::endl;
//    for(size_t i = 0; i < m.remapping_table.size(); i++){
//        os << m.remapping_table[i] << std::endl;
//    }
//    
//    os << "#normals and k1" << std::endl;
//    for(size_t i = 0; i < m.normals_and_k1.size(); i++){
//        os << m.normals_and_k1[i].x << " " << m.normals_and_k1[i].y << " " << m.normals_and_k1[i].z << " " << m.normals_and_k1[i].w << std::endl;
//    }
//
//    os << "#x1 and k2" << std::endl;
//    for(size_t i = 0; i < m.x1_and_k2.size(); i++){
//        os << m.x1_and_k2[i].x << " " << m.x1_and_k2[i].y << " " << m.x1_and_k2[i].z << " " << m.x1_and_k2[i].w << std::endl;
//    }
//
//    os << "#N0" << std::endl;
//    for(size_t i = 0; i < m.N0.size(); i++){
//        os << m.N0[i].x << " " << m.N0[i].y  << std::endl;
//    }
//    
//    os << "#N1" << std::endl;
//    for(size_t i = 0; i < m.N1.size(); i++){
//        os << m.N1[i].x << " " << m.N1[i].y  << " " << m.N1[i].z << " " << m.N1[i].w << std::endl;        
//    }
//    
//    os << "#N2" << std::endl;
//    for(size_t i = 0; i < m.N2.size(); i++){
//        os << m.N2[i].x << " " << m.N2[i].y  << " " << m.N2[i].z << " " << m.N2[i].w << std::endl;        
//    }
//    
//    os << "#N3" << std::endl;
//    for(size_t i = 0; i < m.N3.size(); i++){
//        os << m.N3[i].x << " " << m.N3[i].y  << " " << m.N3[i].z << " " << m.N3[i].w << std::endl;        
//    }
//    
//    return os;
//}
//
///**********************/
///* LOAD MESH FUNCTION */
///**********************/
//bool load_model_nbin(const std::string filename, nvmesh &m) {
//
//    using namespace nv_types;    
//    std::string line;
//
//	size_t num_verts			= 0;
//    size_t num_normals			= 0;
//    size_t num_faces			= 0;
//    size_t remapping_table_size = 0;
//    size_t node_list_size		= 0;
//    
//	std::ifstream ifile(filename.c_str(), std::ios::binary);
//        
//    if (!ifile.is_open()) return false;
//    
//	// --- Read the number of vertices and faces
//	{   
//        get_valid_line(ifile, line);
//        std::stringstream ss(line);
//        ss >> num_verts >> num_faces;       
//    }
//
//	printf("Number of vertices = %i; Number of faces = %i\n", num_verts, num_faces);
//
//	// --- Load vertices and faces
//	{   
//        nv_float3 *verts = new nv_float3[num_verts];
//        
//        nv_int4   *faces = new nv_int4[num_faces];
//        
//        ifile.read(reinterpret_cast<char*>(verts), sizeof(nv_float3) * num_verts);
//        ifile.read(reinterpret_cast<char*>(faces), sizeof(nv_int4) * num_faces);
//
//        printf("Loading vertices...\n");
//		m.verts.resize(num_verts);
//        for(size_t i =0; i < num_verts;i++){
//            nv_float3 v = verts[i];
//            m.verts[i] = make_float4(v.x, v.y, v.z, 1.0);
//        }
//		printf("Done.\n");
//
//        printf("Loading faces...\n");
//        m.faces.resize(num_faces);
//        for(size_t i =0; i < num_faces;i++){
//            nv_int4 f = faces[i];
//            m.faces[i] = make_int4(f.x, f.y, f.z, f.w);
//        }
//		printf("Done.\n");
//
//        delete [] verts;
//        delete [] faces;
//               
//        // --- Skip the new line character
//        int  newline = ifile.get();
//        assert(newline == '\n');
//    }
//    
//	// --- Read the sizes of remapping table and nonde list
//	{
//        get_valid_line(ifile,line);
//        std::stringstream ss(line);
//        ss >> remapping_table_size >> node_list_size;
//		printf("Size of the remapping table = %i; Size of the node list = %i\n", remapping_table_size, node_list_size);
//    }
//    
//    // --- Load the remapping table
//	{
//
//        printf("Loading the remapping table...\n");
//		m.remapping_table.resize(remapping_table_size);
//		ifile.read(reinterpret_cast<char*>(thrust::raw_pointer_cast(m.remapping_table.data())), sizeof(int) * remapping_table_size);
//        printf("Done.\n");
//        
//	}
//    
//    // --- Load N0
//	{    
//        printf("Loading inner/leaf node information...\n");
//		nv_int2 *N0 = new nv_int2[node_list_size];
//        ifile.read(reinterpret_cast<char*>(N0),sizeof(nv_int2)*node_list_size);
//
//        m.N0.resize(node_list_size);
//        for(size_t i = 0; i < node_list_size; i++){
//            const nv_types::nv_int2 n = N0[i];
//            m.N0[i] = make_int2(n.x,n.y);
//        }    
//        delete [] N0;
//        printf("Done.\n");
//    }
//    
//    // --- Load N1
//    {
//        printf("Loading x/y bounds for left child...\n");
//        nv_float4 *N1 = new nv_float4[node_list_size]; 
//        ifile.read(reinterpret_cast<char*>(N1),sizeof(nv_float4)*node_list_size);
//        
//        m.N1.resize(node_list_size);
//        for(size_t i = 0; i < node_list_size; i++){
//            const nv_types::nv_float4 n = N1[i];
//            m.N1[i] = make_float4(n.x,n.y,n.z,n.w);
//        }    
//        delete [] N1;
//        printf("Done.\n");
//    }
//    
//    // --- Load N2
//	{   
//        printf("Loading x/y bounds for right child...\n");
//        nv_float4 *N2 = new nv_float4[node_list_size]; 
//        ifile.read(reinterpret_cast<char*>(N2),sizeof(nv_float4)*node_list_size);
//        m.N2.resize(node_list_size);
//        for(size_t i = 0; i < node_list_size; i++){
//            const nv_types::nv_float4 n = N2[i];
//            m.N2[i] = make_float4(n.x,n.y,n.z,n.w);
//        }    
//        delete [] N2;
//        printf("Done.\n");
//    }
//    
//    // --- Load N3
//    {   
//        printf("Loading z bounds for left and right child...\n");
//        nv_float4 *N3 = new nv_float4[node_list_size];
//        ifile.read(reinterpret_cast<char*>(N3),sizeof(nv_float4)*node_list_size);
//
//        m.N3.resize(node_list_size);
//        for(size_t i = 0; i < node_list_size; i++){
//            const nv_types::nv_float4 n = N3[i];
//            m.N3[i] = make_float4(n.x,n.y,n.z,n.w);
//        }    
//        delete [] N3;
//        printf("Done.\n");
//    }
//    
//    // --- Skip the new line character
//    int  newline = ifile.get();
//    assert(newline == '\n');
//    
//    {   //read the number of vertex and faces
//        get_valid_line(ifile,line);
//        std::stringstream ss(line);
//        ss >> num_normals ;       
//		printf("Number of normals = %i\n", num_normals);
//	}
//
//	assert(num_normals == num_verts);
//    
//    // --- Load normals and k1
//	{   
//        printf("Loading normals and principal curvatures...\n");
//        nv_float4 *normals_and_k1 = new nv_float4[num_normals];
//
//        ifile.read(reinterpret_cast<char*>(normals_and_k1),sizeof(nv_float4)*num_normals);
//
//        m.normals_and_k1.resize(num_normals);
//        for(size_t i =0; i < num_normals;i++){
//            nv_float4 nk1 = normals_and_k1[i];
//            m.normals_and_k1[i] = make_float4(nk1.x,nk1.y,nk1.z,nk1.w);
//        }
//
//        delete [] normals_and_k1;
//        printf("Done.\n");
//    }
//
//    // --- Load x1 and k2
//	{   
//        printf("Loading principal directions at the vertices and second principal curvatures...\n");
//        nv_float4 *x1_and_k2 = new nv_float4[num_normals];
//
//        ifile.read(reinterpret_cast<char*>(x1_and_k2),sizeof(nv_float4)*num_normals);
//
//        m.x1_and_k2.resize(num_normals);
//        for(size_t i =0; i < num_normals; i++){
//            nv_float4 x1k2 = x1_and_k2[i];
//            m.x1_and_k2[i] = make_float4(x1k2.x,x1k2.y,x1k2.z,x1k2.w);
//        }
//
//        delete [] x1_and_k2;
//        printf("Done.\n");
//    }
//    
//	printf("Computing the minimum edge size for visualization...\n");
// //   //m.min_len = norm(m.verts[m.faces[0].x] - m.verts[m.faces[0].y]);
//	//printf("m.faces[0].x %i m.verts.size %i \n", m.faces[0].x, m.verts.size());
//	////printf("m.verts[m.faces[0].x].x %f m.verts[m.faces[0].x].y %f m.verts[m.faces[0].x].z %f\n", m.verts[m.faces[0].x].x, m.verts[m.faces[0].x].y, m.verts[m.faces[0].x].z);
//	
//	
//	
//    //model.min_len = norm(model.vlist[model.flist[0].a]- model.vlist[model.flist[0].b]);
//    ////compute the min edge_size
//    //for(size_t i =0; i< num_faces; i++){        
//    //    
//    //    float len_a = norm(model.vlist[model.flist[i].a]- model.vlist[model.flist[i].b]);
//    //    float len_b = norm(model.vlist[model.flist[i].b]- model.vlist[model.flist[i].c]);
//    //    float len_c = norm(model.vlist[model.flist[i].c]- model.vlist[model.flist[i].a]);;
//    //    
//    //    model.min_len = fminf(model.min_len,fminf(len_a,fminf(len_b,len_c)));
//    //    
//    //}
//	
//	
//	printf("%i %i\n", m.faces.size(), m.verts.size());
//	printf("%i\n", m.faces[0].x);
//	
//	vec3 v1 = make_vec3(m.verts[m.faces[0].x].x, m.verts[m.faces[0].x].y, m.verts[m.faces[0].x].z);
//	printf("Here13.a\n");
//    vec3 v2 = make_vec3(m.verts[m.faces[0].y].x, m.verts[m.faces[0].y].y, m.verts[m.faces[0].y].z);
//	printf("Here13.b\n");
//    m.min_len = norm(v1 - v2);
//	//compute the min edge_size
//    for(size_t i = 0; i< num_faces; i++){        
//        
//        v1 = make_vec3(m.verts[m.faces[i].x].x, m.verts[m.faces[i].x].y, m.verts[m.faces[i].x].z);
//		v2 = make_vec3(m.verts[m.faces[i].y].x, m.verts[m.faces[i].y].y, m.verts[m.faces[i].y].z);
//		float len_a = norm(v1 - v2);
//		v1 = make_vec3(m.verts[m.faces[i].y].x, m.verts[m.faces[i].y].y, m.verts[m.faces[i].y].z);
//		v2 = make_vec3(m.verts[m.faces[i].z].x, m.verts[m.faces[i].z].y, m.verts[m.faces[i].z].z);
//		float len_b = norm(v1 - v2);
//		v1 = make_vec3(m.verts[m.faces[i].z].x, m.verts[m.faces[i].z].y, m.verts[m.faces[i].z].z);
//		v2 = make_vec3(m.verts[m.faces[i].x].x, m.verts[m.faces[i].x].y, m.verts[m.faces[i].x].z);
//        float len_c = norm(v1 - v2);
//        
//#ifdef __linux__
//		model.min_len = fminf(model.min_len,fminf(len_a,fminf(len_b,len_c)));
//#elif _WIN32 || _WIN64
//		//model.min_len = min(model.min_len, min(len_a, min(len_b, len_c)));
//		m.min_len = minimum(m.min_len, minimum(len_a, minimum(len_b, len_c)));
//#endif
//    }
//    printf("Done.\n");
//    
//    for(size_t i = 0; i < m.verts.size(); i++){
//        m.bb.grow(make_vec3(m.verts[i].x,m.verts[i].y,m.verts[i].z));        
//    }
//   
//}
//
//#include <thrust\host_vector.h>
#include "load_nbin.h"

#ifdef _WIN32 || _WIN64
typedef unsigned int uint;
#endif

#include <fstream>
#include <sstream>
#include <cassert>
#include <string>

#include "utils.h"

/********************/
/* GET MESH NORMALS */
/********************/
vec3 nvmesh::get_normal(const int ID, float u, float  v) const {
    
    // --- Get the vertex indices of the triangle
    const int4 face   = faces[ID];
    const uint v1_idx = face.x;
    const uint v2_idx = face.y;
    const uint v3_idx = face.z;
       
    // --- Get the normals
    const float4 data_l2_v1 = normals_and_k1[v1_idx];
    const float4 data_l2_v2 = normals_and_k1[v2_idx];
    const float4 data_l2_v3 = normals_and_k1[v3_idx];
    
    // --- Normals in each vertex
    const vec3 normal_v1  = make_vec3(data_l2_v1.x, data_l2_v1.y, data_l2_v1.z);
    const vec3 normal_v2  = make_vec3(data_l2_v2.x, data_l2_v2.y, data_l2_v2.z);
    const vec3 normal_v3  = make_vec3(data_l2_v3.x, data_l2_v3.y, data_l2_v3.z);                
    
    // --- Perform interpolation
    const vec3 normal_v = normalize(normal_v1 + u * (normal_v2 - normal_v1) + v * (normal_v3 - normal_v1));
    
    return normal_v;
                                                    
}

/******************************/
/* GET NORMALS AND CURVATURES */
/******************************/
void  nvmesh::get_normal_and_curvature(const int ID, float u, float v, vec3 &normal_v, vec3 &U1_v, vec3 &U2_v, float &a1, float &a2) const {
    
    // --- Gather the data 
    
	// --- Get the vertex indices of the triangle
    const int4 face   = faces[ID];
    const uint v1_idx = face.x;
    const uint v2_idx = face.y;
    const uint v3_idx = face.z;
    
    // --- Get the normals
    const float4 data_l2_v1 = normals_and_k1[v1_idx];
    const float4 data_l2_v2 = normals_and_k1[v2_idx];
    const float4 data_l2_v3 = normals_and_k1[v3_idx];

    const float4 data_l3_v1 = x1_and_k2[v1_idx];
    const float4 data_l3_v2 = x1_and_k2[v2_idx];
    const float4 data_l3_v3 = x1_and_k2[v3_idx];
    
    // --- Normals in each vertex
    const vec3 normal_v1  = make_vec3(data_l2_v1.x, data_l2_v1.y, data_l2_v1.z);
    const vec3 normal_v2  = make_vec3(data_l2_v2.x, data_l2_v2.y, data_l2_v2.z);
    const vec3 normal_v3  = make_vec3(data_l2_v3.x, data_l2_v3.y, data_l2_v3.z);                

    const vec3 x1_v1  = make_vec3(data_l3_v1.x, data_l3_v1.y, data_l3_v1.z);
    const vec3 x1_v2  = make_vec3(data_l3_v2.x, data_l3_v2.y, data_l3_v2.z);
    const vec3 x1_v3  = make_vec3(data_l3_v3.x, data_l3_v3.y, data_l3_v3.z);                

    const float k1_v1   =  data_l2_v1.w;
    const float k1_v2   =  data_l2_v2.w;
    const float k1_v3   =  data_l2_v3.w;
    
    const float k2_v1   =  data_l3_v1.w;
    const float k2_v2   =  data_l3_v2.w;
    const float k2_v3   =  data_l3_v3.w;
    
    // --- Perform interpolation
    normal_v = normalize(normal_v1 + u * (normal_v2 - normal_v1) + v * (normal_v3 - normal_v1));
    
    U1_v = normalize(x1_v1 + u * (x1_v2 - x1_v1) + v * (x1_v3 - x1_v1));
    
    const float k1 = k1_v1 + u * (k1_v2 - k1_v1) + v * (k1_v3 - k1_v1);
    
    const float k2 = k2_v1 + u * (k2_v2 - k2_v1) + v * (k2_v3 - k2_v1);
                                                     
    a1= 1 / k1;
    
    a2= 1 / k2;
    
    U2_v = cross(normal_v, U1_v);
}

/**************/
/* PRINT MESH */
/**************/
std::ostream &operator << (std::ostream &os, const nvmesh &m){
    
    os << "#aabb" << std::endl;
    os << m.bb <<  std::endl;
        
    os << "#vertices" << std::endl;
    for(size_t i = 0; i < m.verts.size(); i++){
        os << m.verts[i].x << " " << m.verts[i].y << " " << m.verts[i].z << " " << m.verts[i].w << std::endl;
    }
    
    os << "#faces" << std::endl;
    for(size_t i = 0; i < m.faces.size(); i++){
        os << m.faces[i].x << " " << m.faces[i].y << " " << m.faces[i].z << " " << m.faces[i].w << std::endl;
    }
    
    os << "#remapping table" << std::endl;
    for(size_t i = 0; i < m.remapping_table.size(); i++){
        os << m.remapping_table[i] << std::endl;
    }
    
    os << "#normals and k1" << std::endl;
    for(size_t i = 0; i < m.normals_and_k1.size(); i++){
        os << m.normals_and_k1[i].x << " " << m.normals_and_k1[i].y << " " << m.normals_and_k1[i].z << " " << m.normals_and_k1[i].w << std::endl;
    }

    os << "#x1 and k2" << std::endl;
    for(size_t i = 0; i < m.x1_and_k2.size(); i++){
        os << m.x1_and_k2[i].x << " " << m.x1_and_k2[i].y << " " << m.x1_and_k2[i].z << " " << m.x1_and_k2[i].w << std::endl;
    }

    os << "#N0" << std::endl;
    for(size_t i = 0; i < m.N0.size(); i++){
        os << m.N0[i].x << " " << m.N0[i].y  << std::endl;
    }
    
    os << "#N1" << std::endl;
    for(size_t i = 0; i < m.N1.size(); i++){
        os << m.N1[i].x << " " << m.N1[i].y  << " " << m.N1[i].z << " " << m.N1[i].w << std::endl;        
    }
    
    os << "#N2" << std::endl;
    for(size_t i = 0; i < m.N2.size(); i++){
        os << m.N2[i].x << " " << m.N2[i].y  << " " << m.N2[i].z << " " << m.N2[i].w << std::endl;        
    }
    
    os << "#N3" << std::endl;
    for(size_t i = 0; i < m.N3.size(); i++){
        os << m.N3[i].x << " " << m.N3[i].y  << " " << m.N3[i].z << " " << m.N3[i].w << std::endl;        
    }
    
    return os;
}

/**********************/
/* LOAD MESH FUNCTION */
/**********************/
bool load_model_nbin(const std::string filename, nvmesh &m) {

    using namespace nv_types;    
    std::string line;

	size_t num_verts			= 0;
    size_t num_normals			= 0;
    size_t num_faces			= 0;
    size_t remapping_table_size = 0;
    size_t node_list_size		= 0;
    
	std::ifstream ifile(filename.c_str(), std::ios::binary);
        
    if (!ifile.is_open()) return false;
    
	// --- Read the number of vertices and faces
	{   
        get_valid_line(ifile, line);
        std::stringstream ss(line);
        ss >> num_verts >> num_faces;       
    }

	printf("Number of vertices = %i; Number of faces = %i\n", num_verts, num_faces);

	// --- Load vertices and faces
	{   
        nv_float3 *verts = new nv_float3[num_verts];
        
        nv_int4   *faces = new nv_int4[num_faces];
        
        ifile.read(reinterpret_cast<char*>(verts), sizeof(nv_float3) * num_verts);
        ifile.read(reinterpret_cast<char*>(faces), sizeof(nv_int4) * num_faces);

        printf("Loading vertices...\n");
		m.verts.resize(num_verts);
        for(size_t i =0; i < num_verts;i++){
            nv_float3 v = verts[i];
			printf("Vert %i out of %i\n", i, num_verts);
            m.verts[i] = make_float4(v.x, v.y, v.z, 1.0);
        }
		printf("Done.\n");

        printf("Loading faces...\n");
        m.faces.resize(num_faces);
        for(size_t i =0; i < num_faces;i++){
            nv_int4 f = faces[i];
            m.faces[i] = make_int4(f.x, f.y, f.z, f.w);
        }
		printf("Done.\n");

        delete [] verts;
        delete [] faces;
               
        // --- Skip the new line character
        int  newline = ifile.get();
        assert(newline == '\n');
    }
    
	// --- Read the sizes of remapping table and nonde list
	{
        get_valid_line(ifile,line);
        std::stringstream ss(line);
        ss >> remapping_table_size >> node_list_size;
		printf("Size of the remapping table = %i; Size of the node list = %i\n", remapping_table_size, node_list_size);
    }
    
    // --- Load the remapping table
	{

        printf("Loading the remapping table...\n");
		m.remapping_table.resize(remapping_table_size);
		ifile.read(reinterpret_cast<char*>(thrust::raw_pointer_cast(m.remapping_table.data())), sizeof(int) * remapping_table_size);
        printf("Done.\n");
        
	}
    
    // --- Load N0
	{    
        printf("Loading inner/leaf node information...\n");
		nv_int2 *N0 = new nv_int2[node_list_size];
        ifile.read(reinterpret_cast<char*>(N0),sizeof(nv_int2)*node_list_size);

        m.N0.resize(node_list_size);
        for(size_t i = 0; i < node_list_size; i++){
            const nv_types::nv_int2 n = N0[i];
            m.N0[i] = make_int2(n.x,n.y);
        }    
        delete [] N0;
        printf("Done.\n");
    }
    
    // --- Load N1
    {
        printf("Loading x/y bounds for left child...\n");
        nv_float4 *N1 = new nv_float4[node_list_size]; 
        ifile.read(reinterpret_cast<char*>(N1),sizeof(nv_float4)*node_list_size);
        
        m.N1.resize(node_list_size);
        for(size_t i = 0; i < node_list_size; i++){
            const nv_types::nv_float4 n = N1[i];
            m.N1[i] = make_float4(n.x,n.y,n.z,n.w);
        }    
        delete [] N1;
        printf("Done.\n");
    }
    
    // --- Load N2
	{   
        printf("Loading x/y bounds for right child...\n");
        nv_float4 *N2 = new nv_float4[node_list_size]; 
        ifile.read(reinterpret_cast<char*>(N2),sizeof(nv_float4)*node_list_size);
        m.N2.resize(node_list_size);
        for(size_t i = 0; i < node_list_size; i++){
            const nv_types::nv_float4 n = N2[i];
            m.N2[i] = make_float4(n.x,n.y,n.z,n.w);
        }    
        delete [] N2;
        printf("Done.\n");
    }
    
    // --- Load N3
    {   
        printf("Loading z bounds for left and right child...\n");
        nv_float4 *N3 = new nv_float4[node_list_size];
        ifile.read(reinterpret_cast<char*>(N3),sizeof(nv_float4)*node_list_size);

        m.N3.resize(node_list_size);
        for(size_t i = 0; i < node_list_size; i++){
            const nv_types::nv_float4 n = N3[i];
            m.N3[i] = make_float4(n.x,n.y,n.z,n.w);
        }    
        delete [] N3;
        printf("Done.\n");
    }
    
    // --- Skip the new line character
    int  newline = ifile.get();
    assert(newline == '\n');
    
    {   //read the number of vertex and faces
        get_valid_line(ifile,line);
        std::stringstream ss(line);
        ss >> num_normals ;       
		printf("Number of normals = %i\n", num_normals);
	}

	assert(num_normals == num_verts);
    
    // --- Load normals and k1
	{   
        printf("Loading normals and principal curvatures...\n");
        nv_float4 *normals_and_k1 = new nv_float4[num_normals];

        ifile.read(reinterpret_cast<char*>(normals_and_k1),sizeof(nv_float4)*num_normals);

        m.normals_and_k1.resize(num_normals);
        for(size_t i =0; i < num_normals;i++){
            nv_float4 nk1 = normals_and_k1[i];
            m.normals_and_k1[i] = make_float4(nk1.x,nk1.y,nk1.z,nk1.w);
        }

        delete [] normals_and_k1;
        printf("Done.\n");
    }

    // --- Load x1 and k2
	{   
        printf("Loading principal directions at the vertices and second principal curvatures...\n");
        nv_float4 *x1_and_k2 = new nv_float4[num_normals];

        ifile.read(reinterpret_cast<char*>(x1_and_k2),sizeof(nv_float4)*num_normals);

        m.x1_and_k2.resize(num_normals);
        for(size_t i =0; i < num_normals; i++){
            nv_float4 x1k2 = x1_and_k2[i];
            m.x1_and_k2[i] = make_float4(x1k2.x,x1k2.y,x1k2.z,x1k2.w);
        }

        delete [] x1_and_k2;
        printf("Done.\n");
    }
    
	printf("Computing the minimum edge size for visualization...\n");
 	vec3 v1 = make_vec3(m.verts[m.faces[0].x].x, m.verts[m.faces[0].x].y, m.verts[m.faces[0].x].z);
    vec3 v2 = make_vec3(m.verts[m.faces[0].y].x, m.verts[m.faces[0].y].y, m.verts[m.faces[0].y].z);
    m.min_len = norm(v1 - v2);
	//compute the min edge_size
    for(size_t i = 0; i< num_faces; i++){        
        
        v1 = make_vec3(m.verts[m.faces[i].x].x, m.verts[m.faces[i].x].y, m.verts[m.faces[i].x].z);
		v2 = make_vec3(m.verts[m.faces[i].y].x, m.verts[m.faces[i].y].y, m.verts[m.faces[i].y].z);
		float len_a = norm(v1 - v2);
		v1 = make_vec3(m.verts[m.faces[i].y].x, m.verts[m.faces[i].y].y, m.verts[m.faces[i].y].z);
		v2 = make_vec3(m.verts[m.faces[i].z].x, m.verts[m.faces[i].z].y, m.verts[m.faces[i].z].z);
		float len_b = norm(v1 - v2);
		v1 = make_vec3(m.verts[m.faces[i].z].x, m.verts[m.faces[i].z].y, m.verts[m.faces[i].z].z);
		v2 = make_vec3(m.verts[m.faces[i].x].x, m.verts[m.faces[i].x].y, m.verts[m.faces[i].x].z);
        float len_c = norm(v1 - v2);
        
#ifdef __linux__
		model.min_len = fminf(model.min_len,fminf(len_a,fminf(len_b,len_c)));
#elif _WIN32 || _WIN64
		//model.min_len = min(model.min_len, min(len_a, min(len_b, len_c)));
		m.min_len = minimum(m.min_len, minimum(len_a, minimum(len_b, len_c)));
#endif
    }
    printf("Done.\n");
    
    for(size_t i = 0; i < m.verts.size(); i++){
        m.bb.grow(make_vec3(m.verts[i].x,m.verts[i].y,m.verts[i].z));        
    }
   
}

