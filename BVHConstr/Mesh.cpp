#include <fstream> 
#include <sstream>
#include <cassert>
#include <iomanip>

#include <thrust\host_vector.h>

#include "Mesh.h"
#include "packed_struct.h"

/*****************/
/* LOAD PLY MESH */
/*****************/
bool loadPly(mesh &plyMesh) {
						
	enum DATA_FORMAT {ASCII, BINARY_LITTLE_ENDIAN, BINARY_BIG_ENDIAN} format;
    
	std::string			line;
    std::string			token;
    size_t				numVertices = 0;
    size_t				numFaces = 0;
    plyMesh.normalsAvailable = false;
    plyMesh.curvatureAvailable = false;
    
    std::ifstream ifile(plyMesh.filename.c_str());
		
	if (!ifile.is_open()) return false;
	
	{   // --- Read the magic number
	    get_valid_line(ifile, line);
	    std::stringstream ss(line);
	    ss >> token;
	    assert(token == "ply");
	}
	
	{	// --- Read the format
		get_valid_line(ifile, line);
	    std::stringstream ss(line);
		ss >> token;
	    assert(token == "format");
	    ss >> token;
	    if (token == "ascii") format = ASCII;
		else if (token == "binary_little_endian") {
			format = BINARY_LITTLE_ENDIAN;
			assert(false);
		} else if(token == "binary_little_endian") {
			format = BINARY_BIG_ENDIAN;
			assert(false);
		} else assert(false);
	}
	
	// --- Read the header
	do {
    	get_valid_line(ifile, line);
    	std::stringstream ss(line);
    	ss >> token ; 
    	
    	if(token == "element") {
            ss >> token;									// --- read the next token
            if(token == "vertex")  ss >> numVertices;			// --- read the number of vertices
            else if (token == "face") ss >> numFaces;		// --- read the number of faces
    	} else if(token == "property") {
            ss >> token;									// --- read the next token
            ss >> token;
            if (token == "nx" || token == "ny" || token == "nz") plyMesh.normalsAvailable = true;
            else if (token == "x1x" || token == "x1y" || token == "x1z" || token == "firstCurvatures"  || token == "secondCurvatures")
					plyMesh.curvatureAvailable = true;
        }
	} while (token != "end_header");
	
    if (plyMesh.curvatureAvailable == true && plyMesh.normalsAvailable == false) assert(false);
    
	plyMesh.vertices.resize(0);
	plyMesh.faces.resize(0);
	plyMesh.normals.resize(0);
    plyMesh.principalDirections.resize(0);
    plyMesh.firstCurvatures.resize(0);
    plyMesh.secondCurvatures.resize(0);
    
	plyMesh.vertices.reserve(numVertices);
	plyMesh.faces.reserve(numFaces);
    if (plyMesh.normalsAvailable) plyMesh.normals.reserve(numVertices);
    
    if (plyMesh.curvatureAvailable) {
        plyMesh.principalDirections.reserve(numVertices);
        plyMesh.firstCurvatures.reserve(numVertices);
        plyMesh.secondCurvatures.reserve(numVertices);
    }
    
	if (format == ASCII) {
		
		// --- Read the vertices
		for (size_t i = 0; i< numVertices; i++) {
			
			vec3f v;
            
			get_valid_line(ifile, line);
			std::stringstream ss(line);
			ss >> v.x >> v.y >> v.z;
			plyMesh.vertices.push_back(v);
            
            if (plyMesh.normalsAvailable) {
                vec3f n;
                ss >> n.x >> n.y >> n.z;
                plyMesh.normals.push_back(n);                
            }

            if (plyMesh.curvatureAvailable) {
                vec3f x1;
                float firstCurvatures;
                float secondCurvatures;
                
                ss >> x1.x >> x1.y >> x1.z >> firstCurvatures >> secondCurvatures;
                plyMesh.principalDirections.push_back(x1);
                plyMesh.firstCurvatures.push_back(firstCurvatures);
                plyMesh.secondCurvatures.push_back(secondCurvatures);
            }
            
		}
		assert(plyMesh.vertices.size()==numVertices);
		
		// --- Read the faces
		for (size_t i = 0; i < numFaces; i++) {
			get_valid_line(ifile, line);
			std::stringstream ss(line);
			int list_len,a,b,c;
			ss >> list_len >> a >> b >> c;
			assert(list_len == 3);
			plyMesh.faces.push_back(make_vec3i(a, b, c));	    
		}
		assert(plyMesh.faces.size()==numFaces);
	}

	return true;
}

/*****************/
/* LOAD OBJ MESH */
/*****************/
bool load_obj(mesh &plyMesh) {
						
    std::string line;
    std::string prefix;
    
    plyMesh.normalsAvailable	= false;
    plyMesh.curvatureAvailable	= false;
    
    std::ifstream ifile(plyMesh.filename.c_str());
		
	if (!ifile.is_open()) return false;
	
	plyMesh.vertices.resize(0);
	plyMesh.faces.resize(0);
	
	while (get_valid_line(ifile, line)) {
		
		std::stringstream ss(line);
		ss >> prefix;
		
		if (prefix == "v") {
			float x, y, z;
			ss >> x >> y >> z;   
			plyMesh.vertices.push_back(make_vec3f(x, y, z));	         
		} else if (prefix == "f") {
			int a, b, c;
			ss >> a >> b >> c;
			plyMesh.faces.push_back(make_vec3i(a - 1, b - 1, c - 1));	  
		}
	}
	return true;
}

/******************/
/* GET VALID LINE */
/******************/
bool get_valid_line(std::istream &ifile, std::string &line) {
    do {
        std::getline(ifile, line);
        line = trim(line); // --- trim leading and trailing white spaces
    } while (ifile.good() && (line[0]=='#' || line.size()==0));
    
    return ifile.good();
}

/*****************/
/* TRIM FUNCTION */
/*****************/
std::string trim(const std::string& str, const std::string& whitespace) {
    
	const size_t strBegin = str.find_first_not_of(whitespace);
    
	if (strBegin == std::string::npos) return ""; // --- No content

    const size_t strEnd		= str.find_last_not_of(whitespace);
    const size_t strRange	= strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

/*********************************************************************/
/* WRITE THE OUTPUT MESH IN THE MOLLER-TRUMBORE FORMAT - NBIN OUTPUT */
/*********************************************************************/
void writeTrianglesBinaryMollerTrumbore(std::ostream &os, const BVHmesh &outputMesh) {
	
	os << "#Number of vertices and faces: " << std::endl;
    os << outputMesh.vertices.size() << " " << outputMesh.faces.size() << std::endl;
	
	// --- Create a packed list of vertices
    {
        thrust::host_vector<vertex3f> packed_verts(outputMesh.vertices.begin(),outputMesh.vertices.end());
        vertex3f * packed_verts_ptr = thrust::raw_pointer_cast(packed_verts.data());
        os.write(reinterpret_cast<char*>(packed_verts_ptr), sizeof(packed_verts[0])*packed_verts.size());
    }
    
	// --- Create a packed list of faces
    {
        thrust::host_vector<face4i> packed_faces(outputMesh.faces.begin(), outputMesh.faces.end());
        face4i *packed_faces_ptr = thrust::raw_pointer_cast(packed_faces.data());
        os.write(reinterpret_cast<char*>(packed_faces_ptr), sizeof(packed_faces[0])*packed_faces.size());
    }
	os << std::endl;
	
}

/*****************************************************************/
/* WRITE THE BVH MESH IN THE MOLLER-TRUMBORE FORMAT - SSV OUTPUT */
/*****************************************************************/
void writeTrianglesSSVMollerTrumbore(std::ostream &os, const BVHmesh &outputMesh) {
	
	os << "#number of vertices and faces: " << std::endl;
	os << outputMesh.vertices.size() << " " << outputMesh.faces.size() << std::endl;
	
	os << "##Start vertices list: " << std::endl;
	os << "# x y z" << std::endl;
	for(size_t i = 0; i < outputMesh.vertices.size(); i++) {
		os << std::setprecision(15) << outputMesh.vertices[i].x << " " << outputMesh.vertices[i].y << " "  << outputMesh.vertices[i].z << " ";
		os<< std::endl;
	}	
	os << "##End vertex list." << std::endl;

	os << "##Start faces list: " << std::endl;
	os << "# a b c mesh_index" << std::endl;
	for (size_t i = 0; i < outputMesh.faces.size(); i++) {
		os << outputMesh.faces[i].a << " " << outputMesh.faces[i].b << " "  << outputMesh.faces[i].c << " " << outputMesh.faces[i].d << std::endl;
	}	
	os << "##End face list." << std::endl;
	os << std::endl;
	
}
