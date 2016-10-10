#include <iomanip>
#include <map>
#include <set>

#include "Eigen\Dense"
#include "normals_and_curvatures.h"
#include "Mesh.h"
#include "packed_struct.h"
#include "thrust/host_vector.h"

template<typename DerivedA>
Eigen::Matrix< typename DerivedA::Scalar, 
               DerivedA::ColsAtCompileTime, 
               DerivedA::RowsAtCompileTime> 
PseudoInverse( const Eigen::MatrixBase<DerivedA> & a,                         
               double epsilon = std::numeric_limits<typename DerivedA::Scalar>::epsilon())
{
 

    if(a.rows()<a.cols()){
        std::cerr << " a.rows() " << a.rows() << std::endl;
        std::cerr << " a.rows() " << a.rows() << std::endl;
    }
    
    assert(a.rows()>=a.cols());
    
    typedef Eigen::Matrix<typename DerivedA::Scalar,
                       DerivedA::RowsAtCompileTime,
                       DerivedA::ColsAtCompileTime> InputType;
 
    typedef Eigen::Matrix<typename DerivedA::Scalar,
                        DerivedA::ColsAtCompileTime,
                        DerivedA::RowsAtCompileTime> ReturnType;
 
    Eigen::JacobiSVD<InputType> svd = a.jacobiSvd(Eigen::ComputeFullU |Eigen::ComputeFullV);

    double tolerance = epsilon * std::max(a.cols(),a.rows()) * svd.singularValues().array().abs().maxCoeff();


    ReturnType sigma = ReturnType::Zero(a.cols(), a.rows());


    sigma.block(0, 0, a.cols(), a.cols()) = (svd.singularValues().array().abs()>tolerance).
                                           select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal();
                                           

    return svd.matrixV() * sigma * svd.matrixU().adjoint();

    
}

void find_best_local_reference(const vec3f& np,
                               vec3f &localX,
                               vec3f &localY,
                               vec3f &localZ)
{

                        
    vec3f axis[3] = {make_vec3f(1.0,0.0,0.0),
                     make_vec3f(0.0,1.0,0.0),
                     make_vec3f(0.0,0.0,1.0)};
    
    int best_axis = 0;
    float min_abs_dot = fabsf(dot(np,axis[best_axis]));
   
    for(int i = 1; i <= 2; i++){  
        float abs_dot = fabsf(dot(np,axis[i]));
        if(abs_dot < min_abs_dot){        
            min_abs_dot = abs_dot;
            best_axis = i;
        }
    }

    localX = normalize(cross(np,axis[best_axis]));
    localY = normalize(cross(np,localX));
    localZ = normalize(localZ);

}



vec3f  transform_to_local_reference(
    const vec3f& p,
    const vec3f& orig,
    const vec3f& localX,
    const vec3f& localY,
    const vec3f& localZ)    
{
    
    vec3f new_p;
    vec3f dp = p-orig;
    new_p.x = dot(dp,localX);
    new_p.y = dot(dp,localY);
    new_p.z = dot(dp,localZ);
    
    return new_p;
} 


bool QSA(const vec3f  &p, 
         const  vec3f &np,
         const std::vector<vec3f > &point_list,
         float &firstCurvatures, vec3f & X1,
         float &secondCurvatures, vec3f & X2
        )
{



    size_t num_points = point_list.size();
    
    vec3f localX;
    vec3f localY;
    vec3f localZ;
    
    //find a local orthonormal coordinate system with z-axis along Np    
    find_best_local_reference(np,localX,localY,localZ); //localZ is equal to NP;
    
    //local coordinates of vertices         
  
    Eigen::Matrix<double,Eigen::Dynamic,3> M(num_points,3);
    Eigen::Matrix<double,Eigen::Dynamic,1> Z_vec(num_points);

        
    for(size_t i = 0; i< num_points; i++){

        vec3f new_p = transform_to_local_reference(point_list[i],p,localX,localY,localZ);
        double xi = new_p.x;
        double yi = new_p.y;
        Z_vec[i] = new_p.z;
        
        M(i,0) = 0.5*(xi*xi);
        M(i,1) = (xi*yi);
        M(i,2) = 0.5*(yi*yi);
    }   
    
    Eigen::Matrix<double,3,1> K = PseudoInverse(M)*Z_vec;
    
    Eigen::Matrix<double,2,2> W;
    W(0,0) =  K[0];
    W(0,1) =  K[1];
    W(1,0) =  K[1] ;
    W(1,1) =  K[2];
    
    Eigen::EigenSolver<Eigen::Matrix<double,2,2> > es(W);
    
    Eigen::Matrix<double,2,2> D  = es.pseudoEigenvalueMatrix();
    Eigen::Matrix<double,2,2> V  = es.pseudoEigenvectors(); 
        
    firstCurvatures = -D(0,0);
    secondCurvatures = -D(1,1);
    
    
    X1 = normalize(V(0,0)*localX+V(1,0)*localY);
    X2 = normalize(V(0,1)*localX+V(1,1)*localY);
    
    return true;
}

void compute_curvatures(
    const std::vector<vec3f> &vertices, //list of vertex
    const std::vector<vec3i> &faces, //list of faces
    const std::vector<vec3f> &normals, //list of normals              
    std::vector<vec3f> &principalDirections, //list of normals              
    std::vector<float> &k1list, //list of normals              
    std::vector<float> &k2list //list of normals              
){
    
    principalDirections.clear();
    k1list.clear();
    k2list.clear();   
    
    principalDirections.reserve(vertices.size());
    k1list.reserve(vertices.size());
    k2list.reserve(vertices.size());
    
    
    std::vector< 
        std::vector<int> 
    >vert_to_face_map(vertices.size());
    
    
    //for each vertex make a list of faces which share that particular vertex 
    //and compute the normal of each face
    for(size_t face = 0; face < faces.size(); face++){
       vert_to_face_map[faces[face].a].push_back(face);
       vert_to_face_map[faces[face].b].push_back(face);
       vert_to_face_map[faces[face].c].push_back(face);
    }
    
    //compute the curvatures
    {
     
        std::set<int> adjacent_verts;
        std::vector<vec3f> point_list;
        //compute principal direction Vectors and curvatures for each vertex   
        for (size_t v =0; v< vertices.size(); v++){
            adjacent_verts.clear();
            point_list.clear();
            //make a list of adjacent vertices
            for(size_t face =0; face<vert_to_face_map[v].size(); face++)
            {
                adjacent_verts.insert(faces[vert_to_face_map[v][face]].a);
                adjacent_verts.insert(faces[vert_to_face_map[v][face]].b);
                adjacent_verts.insert(faces[vert_to_face_map[v][face]].c);
            }
            //erase the vertex v itself
            adjacent_verts.erase(v);        
            
            {//add other vertices
                std::set<int> copy_of_adjacent_verts (adjacent_verts.begin(),adjacent_verts.end());
                for (std::set<int>::iterator it=copy_of_adjacent_verts.begin(); 
                     it!=copy_of_adjacent_verts.end(); ++it){
                     
                    for(size_t face =0; face<vert_to_face_map[*it].size(); face++){
                        adjacent_verts.insert(faces[vert_to_face_map[*it][face]].a);
                        adjacent_verts.insert(faces[vert_to_face_map[*it][face]].b);
                        adjacent_verts.insert(faces[vert_to_face_map[*it][face]].c);
                    }
                }        
                //erase the vertex at v
                adjacent_verts.erase(v);
            }
            
            for (std::set<int>::iterator it=adjacent_verts.begin(); it!=adjacent_verts.end(); ++it){
                point_list.push_back(vertices[*it]);
            }
            
            float firstCurvatures,secondCurvatures;
            vec3f X1;
            vec3f X2;
            
            if(point_list.size() < 3){
                std::cerr << "v = " << v  << std::endl;
            }
            assert(point_list.size()>=3);
            
            //Compute QSA
            QSA(vertices[v],normals[v], point_list, firstCurvatures, X1, secondCurvatures, X2);
                                                
                        
            principalDirections.push_back(X1);
            k1list.push_back(firstCurvatures);
            k2list.push_back(secondCurvatures);
            
        }
    
    }
    
}

void compute_normals(
    const std::vector<vec3f> &vertices, //list of vertex
    const std::vector<vec3i> &faces, //list of faces
    std::vector<vec3f> &normals //list of normals              
     ){

    normals.clear();

    
    normals.reserve(vertices.size());

    
    std::vector< 
        std::vector<int> 
    >vert_to_face_map(vertices.size());
    
    std::vector<vec3f> face_normal_list;
    face_normal_list.reserve(faces.size()); 
    
    
    //preallocate memory space for each vertex 
    for (size_t  v = 0; v< vertices.size(); v++){
        //we expect on average that each vertex touches six faces 
        vert_to_face_map[v].reserve(6);
    }
    
    //for each vertex make a list of faces which share that particular vertex 
    //and compute the normal of each face
    for(size_t face = 0; face < faces.size(); face++){
       vert_to_face_map[faces[face].a].push_back(face);
       vert_to_face_map[faces[face].b].push_back(face);
       vert_to_face_map[faces[face].c].push_back(face);
       
       //compute normal for this face
       const vec3f &A = vertices[faces[face].a];
       const vec3f &B = vertices[faces[face].b];
       const vec3f &C = vertices[faces[face].c];
       
       face_normal_list.push_back(cross((B-A),C-A));
    }
    
    //from here 
    //we know for each vertex a list of faces that share that particular vertex
    //and the normal of each face

    //Compute  normal for each vertex by summing the normals of faces that use 
    // that vertex and normalizing.
    for (size_t v = 0; v< vertices.size(); v++){
    
        vec3f N = make_vec3f(0,0,0);
        for(size_t face = 0; face < vert_to_face_map[v].size(); face++){
            N = N + face_normal_list[vert_to_face_map[v][face]];
        }        

        normals.push_back(normalize(N));
        
    }   
    
    
    
}


void recompute_flat_verts_normals_and_curvature(BVHmesh &globalOutputMesh){
    
    
    {   //recompute faces and verts
        size_t new_num_verts = 3*globalOutputMesh.faces.size();


        std::vector<vec3d> new_verts_list;
        std::vector<vec4i> new_faces_list;
        
        new_verts_list.reserve(new_num_verts);
        new_faces_list.reserve(globalOutputMesh.faces.size());
        
        for(size_t f = 0; f < globalOutputMesh.faces.size(); f++){
            vec4i  old_faces_face = globalOutputMesh.faces[f];
            
            int ia = new_verts_list.size();        
            new_verts_list.push_back(globalOutputMesh.vertices[old_faces_face.a]);
            int ib = new_verts_list.size();        
            new_verts_list.push_back(globalOutputMesh.vertices[old_faces_face.b]);
            int ic = new_verts_list.size();        
            new_verts_list.push_back(globalOutputMesh.vertices[old_faces_face.c]);
            
            vec4i new_face = make_vec4i(ia,ib,ic,old_faces_face.d);       
            
            new_faces_list.push_back(new_face);    
            
        }
        
        globalOutputMesh.vertices.swap(new_verts_list);
        globalOutputMesh.faces.swap(new_faces_list);
    }
    
    
    //compute nornals and principal curvature vectors
    
    globalOutputMesh.normals.resize(globalOutputMesh.vertices.size());
    globalOutputMesh.principalDirections.resize(globalOutputMesh.vertices.size());
    globalOutputMesh.firstCurvatures.resize(globalOutputMesh.vertices.size());
    globalOutputMesh.secondCurvatures.resize(globalOutputMesh.vertices.size());
    
    for(size_t f = 0; f < globalOutputMesh.faces.size(); f++){
        
        const vec3d A = globalOutputMesh.vertices[globalOutputMesh.faces[f].a];
        const vec3d B = globalOutputMesh.vertices[globalOutputMesh.faces[f].b];
        const vec3d C = globalOutputMesh.vertices[globalOutputMesh.faces[f].c];
        //compute faces normal
        const vec3d n_v =  normalize(cross(B-A,C-A));
        //compute the principal curvature vector
        const vec3d X1_v = normalize(B-A);
        
        //here we assume flat triangles
        const double firstCurvatures = 0;
        const double secondCurvatures = 0;

        globalOutputMesh.normals[globalOutputMesh.faces[f].a] = n_v;
        globalOutputMesh.normals[globalOutputMesh.faces[f].b] = n_v;
        globalOutputMesh.normals[globalOutputMesh.faces[f].c] = n_v;
        
        globalOutputMesh.principalDirections[globalOutputMesh.faces[f].a] = X1_v;
        globalOutputMesh.principalDirections[globalOutputMesh.faces[f].b] = X1_v;
        globalOutputMesh.principalDirections[globalOutputMesh.faces[f].c] = X1_v;
        

        globalOutputMesh.firstCurvatures[globalOutputMesh.faces[f].a] = firstCurvatures;
        globalOutputMesh.firstCurvatures[globalOutputMesh.faces[f].b] = firstCurvatures;
        globalOutputMesh.firstCurvatures[globalOutputMesh.faces[f].c] = firstCurvatures;

        
        globalOutputMesh.secondCurvatures[globalOutputMesh.faces[f].a] = secondCurvatures;
        globalOutputMesh.secondCurvatures[globalOutputMesh.faces[f].b] = secondCurvatures;
        globalOutputMesh.secondCurvatures[globalOutputMesh.faces[f].c] = secondCurvatures;        
        
    }
    
    
}




void write_normcurv_ssv_mt(std::ostream &os,const  BVHmesh &partialInputMesh){
    os << "#number of vertices"<< std::endl;
    os << partialInputMesh.vertices.size() <<std::endl;
    os << "##start normal list" << std::endl;
    os << "# x y z" << std::endl;
    for(size_t i = 0; i < partialInputMesh.normals.size(); i++){
        os << std::setprecision(15) << partialInputMesh.normals[i].x << " " << partialInputMesh.normals[i].y << " "  << partialInputMesh.normals[i].z << " " << std::endl;
    }   
    os << "##end normals list" << std::endl;
    os << "##start first principal direction  list" << std::endl;
    for(size_t i = 0; i < partialInputMesh.principalDirections.size(); i++){
        os << std::setprecision(15) << partialInputMesh.principalDirections[i].x << " " << partialInputMesh.principalDirections[i].y << " "  << partialInputMesh.principalDirections[i].z << " " << std::endl;
    }
    os << "##end normals list" << std::endl;
    
    os << "##start curvature list" << std::endl;
    for(size_t i = 0; i < partialInputMesh.principalDirections.size(); i++){
        os << std::setprecision(15) << partialInputMesh.firstCurvatures[i] << " " << partialInputMesh.secondCurvatures[i] << std::endl;
    }  
    os << "##end curvature list" << std::endl;

    os << std::endl;
    
}


//write the BVHmesh in Moller–Trumbore format
void write_normcurv_bin_mt(std::ostream &os,const  BVHmesh &partialInputMesh){
    os << "#number of vertices"<< std::endl;
    os << partialInputMesh.vertices.size() <<std::endl;
    //create a packed list of normals
    {
        thrust::host_vector<vertex4f> packed_normals_k1(partialInputMesh.normals.size());
        for(size_t i =0; i < partialInputMesh.normals.size(); i++){
            packed_normals_k1[i].x = partialInputMesh.normals[i].x;
            packed_normals_k1[i].y = partialInputMesh.normals[i].y;            
            packed_normals_k1[i].z = partialInputMesh.normals[i].z;
            packed_normals_k1[i].w = partialInputMesh.firstCurvatures[i];
        }
            
        vertex4f * packed_normals_k1_ptr = thrust::raw_pointer_cast(packed_normals_k1.data());
        os.write(reinterpret_cast<char*>(packed_normals_k1_ptr), sizeof(packed_normals_k1[0])*packed_normals_k1.size());
    }

    {
        thrust::host_vector<vertex4f> packed_x1_k2(partialInputMesh.principalDirections.size());
        for(size_t i =0; i < partialInputMesh.principalDirections.size(); i++){
            packed_x1_k2[i].x = partialInputMesh.principalDirections[i].x;
            packed_x1_k2[i].y = partialInputMesh.principalDirections[i].y;            
            packed_x1_k2[i].z = partialInputMesh.principalDirections[i].z;
            packed_x1_k2[i].w = partialInputMesh.secondCurvatures[i];
        }
            
        vertex4f * packed_x1_k2_ptr = thrust::raw_pointer_cast(packed_x1_k2.data());
        os.write(reinterpret_cast<char*>(packed_x1_k2_ptr), sizeof(packed_x1_k2[0])*packed_x1_k2.size());
    }
    

    
    
    os << std::endl;
    
}

