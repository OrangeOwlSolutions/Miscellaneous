#include "bvh.h"
#include <thrust/transform_reduce.h>
#include <thrust/sequence.h>
#include <thrust/tuple.h>
#include <thrust/sort.h>
#include <stack>
#include <cassert>
#include <iomanip>

#include <stdint.h>

inner_node::inner_node(const aabb &bb_,
					   bvh_node * node1_,
	                   bvh_node * node2_):
	                   bvh_node(bb_){
	nodes[0] = node1_;
	nodes[1] = node2_;
}

bool inner_node::is_leaf() const{
	return false;
}

size_t inner_node::get_ntris() const{
	return nodes[0]->get_ntris()+nodes[1]->get_ntris();
}

void inner_node::append_tris(thrust::host_vector<size_t> & grmap){
	nodes[0]->append_tris(grmap);
	nodes[1]->append_tris(grmap);
}

size_t inner_node::count(){
	return 1+nodes[0]->count()+nodes[1]->count();
}



void inner_node::print(int depth) const{
	for(int i =0; i < depth; i++)std::cout << "\t" ;
	std::cout << "inner node: " << this <<std::endl;
	for(int i =0; i < depth; i++)std::cout << "\t" ;
	std::cout << "aabb: "<< bb << std::endl;
	for(int i =0; i < depth; i++)std::cout << "\t" ;	
	std::cout << "n0: "<< nodes[0] <<  " n1:" << nodes[1] << std::endl;
	
	nodes[0]->print(depth + 1);
	nodes[1]->print(depth + 1);
}


leaf_node::leaf_node(const aabb &bb_):
bvh_node(bb_){}


bool leaf_node::is_leaf() const{
	return true;
}

size_t leaf_node::get_ntris() const{
	return tri_list.size();
}

size_t leaf_node::count(){
	return 1;
}


void leaf_node::append_tris(thrust::host_vector<size_t> & grmap){
	start_idx = grmap.size();
	grmap.insert(grmap.end(),tri_list.begin(),tri_list.end());
	end_idx = grmap.size();
}


void leaf_node::print(int depth) const{
	for(int i =0; i < depth; i++)std::cout << "\t" ;	
	std::cout << "leaf node: " << this <<std::endl;
	for(int i =0; i < depth; i++)std::cout << "\t" ;
	std::cout << "aabb: "<< bb << std::endl;
	for(int i =0; i < depth; i++)std::cout << "\t" ;
	std::cout << "num objs " << tri_list.size() << std::endl; 
}


static aabb get_tri_aabb(const vec3d &va,const vec3d &vb,const vec3d &vc){
	const vec3d lowerLeftPoint = cwisemin(cwisemin(va,vb),vc);
	const vec3d upperRightPoint = cwisemax(cwisemax(va,vb),vc);    
	return aabb(lowerLeftPoint,upperRightPoint);
}


aabb compute_aabb(const std::vector<vec3d> &vertices,
				  const std::vector<vec4i> &faces){

	
    aabb bb;
    for(size_t i = 0; i< faces.size(); i++){
		const vec3d va = vertices[faces[i].a];
		const vec3d vb = vertices[faces[i].b];
		const vec3d vc = vertices[faces[i].c];	
	
        bb.grow(get_tri_aabb(va,vb,vc));
    }

	return bb;
}


namespace detail{
	
	struct compute_tri_aabb{
		const std::vector<vec3d> &vertices;
		const std::vector<vec4i> &faces;
		
		compute_tri_aabb(const std::vector<vec3d> &vlist_,
						 const std::vector<vec4i> &flist_):
			vertices(vlist_),faces(flist_){}			 
		
		aabb operator()(size_t i){
			const vec3d va = vertices[faces[i].a];
			const vec3d vb = vertices[faces[i].b];
			const vec3d vc = vertices[faces[i].c];	
			const vec3d lowerLeftPoint = cwisemin(cwisemin(va,vb),vc);
			const vec3d upperRightPoint = cwisemax(cwisemax(va,vb),vc);  
			
			return aabb(lowerLeftPoint,upperRightPoint);
		 }
		
	};

	
	struct get_tri_min_max{
		//const std::vector<vec3d> &vertices;
		//const std::vector<vec4i> &faces;
		const thrust::device_vector<aabb>::iterator tri_bounds_begin;
		
		//get_tri_min_max(//const std::vector<vec3d> &vlist_,
						//const std::vector<vec4i> &flist_):
			//vertices(vlist_),faces(flist_){}			 

		get_tri_min_max(const thrust::device_vector<aabb>::iterator tri_bounds_begin_
						):tri_bounds_begin(tri_bounds_begin_){}
		
		thrust::tuple<vec3d,vec3d>
		 operator()(size_t i){
			const aabb  bb = tri_bounds_begin[i];
			return thrust::make_tuple(bb.lowerLeftPoint, bb.upperRightPoint);
		 }
		
	};
	
	struct min_max{
		thrust::tuple<vec3d,vec3d>
			operator()(thrust::tuple<vec3d,vec3d> a,
				       thrust::tuple<vec3d,vec3d> b ){
				
				const vec3d lowerLeftPoint = cwisemin(
											   thrust::get<0>(a),
											   thrust::get<0>(b)
											  );
											  
				const vec3d upperRightPoint = cwisemax(
											   thrust::get<1>(a),
											   thrust::get<1>(b)
											  );  
			
				return thrust::make_tuple(lowerLeftPoint,upperRightPoint); 
			}
	};
	
}



bvh_builder::bvh_builder(bvh_platform platform_,
			const BVHmesh & m_):platform(platform_),
							 partialInputMesh(m_),
							 compare_obj_centroid(*this){
	root = NULL;
	num_tris = partialInputMesh.faces.size();
}


double bvh_builder::get_sah_cost(const aabb &bb, 
                                 const aabb &lbb,
                                 const aabb &rbb, 
                                 size_t NL, size_t NR) const{
    double lsa = lbb.get_surface_area();
    double rsa = rbb.get_surface_area();
    double sa  = bb.get_surface_area() ;                        
                               
    return platform.KT + platform.KI*(lsa/sa*NL + rsa/sa*NR);
}


void bvh_builder::build(){
	thrust::device_vector<size_t> remap_table(num_tris);
	r_bounds.resize(num_tris);
	l_bounds.resize(num_tris);
    tri_bounds.resize(num_tris);
       
	thrust::sequence(remap_table.begin(),remap_table.end());
	thrust::transform(remap_table.begin(),
					  remap_table.end(),
					  tri_bounds.begin(),
				  	  detail::compute_tri_aabb(partialInputMesh.vertices,partialInputMesh.faces)
	);
	
	
	std::cerr << "start building  bvh" << std::endl; 
	
	root = rec_build(remap_table.begin(),
					 remap_table.end(),
					0);
			  
}

bool bvh_builder::
     cmp_obj_centroid::operator()(size_t obj1_idx, size_t obj2_idx){ 
     const vec3d  obj1_centroid = parent_ref.get_object_aabb(obj1_idx).get_center();
     const vec3d obj2_centroid = parent_ref.get_object_aabb(obj2_idx).get_center();     
     bool val(false);
     
     //const double * obj1_centroid_ptr = reinterpret_cast<const double*>(&obj1_centroid);
     //const double * obj2_centroid_ptr = reinterpret_cast<const double*>(&obj2_centroid);
     
     //return obj1_centroid_ptr[axis] < obj2_centroid_ptr[axis];
     
     switch(axis){
        case 0:{
            val = obj1_centroid.x < obj2_centroid.x;
        }break;        
        case 1:{
            val = obj1_centroid.y < obj2_centroid.y;
        }break;        
        case 2:{
            val = obj1_centroid.z < obj2_centroid.z;
        }break;        
     }   
     return val;
}


bvh_node* bvh_builder::rec_build(thrust::device_vector<size_t>::iterator remap_begin,
								 thrust::device_vector<size_t>::iterator remap_end,
								 int depth){
	
	size_t tri_list_size = remap_end-remap_begin;
	const aabb bb = get_aabb_from(remap_begin,remap_end);
	
	
	if (depth > platform.max_depth || tri_list_size < 2){//preemptive leaf
		
		leaf_node * lnode = new leaf_node(bb);
		//allocate the space to store the triangle list in this leaf
		lnode->tri_list.resize(tri_list_size);
		thrust::copy(remap_begin,remap_end,lnode->tri_list.begin());
		
		return lnode;
    };
	
	//split
	double min_cost;
	thrust::device_vector<size_t>::iterator  best_obj_split;
	
	split(remap_begin,remap_end,depth,best_obj_split, min_cost);
	
	
	if(platform.KI* tri_list_size < min_cost){
		leaf_node * lnode = new leaf_node(bb);
		//allocate the space to store the triangle list in this leaf
		lnode->tri_list.resize(tri_list_size);
		thrust::copy(remap_begin,remap_end,lnode->tri_list.begin());
		
		return lnode;
	}
	
	//determine if it is needed to allocate a new vector containing the remap_table
	//I think that it is possible to not allocate any vector 
	//but it is necessary  a vector that contain the common tris.
	
	
	bvh_node * node1 = rec_build(remap_begin,best_obj_split,depth + 1);
	if(depth == 0){
		std::cerr << "50%"  << std::endl; 
	}
	
	bvh_node * node2 = rec_build(best_obj_split,remap_end, depth + 1);
	if(depth == 0){
		std::cerr << "100%"  << std::endl; 
	}
	 
	
	return new inner_node(bb,node1,node2);;
}


void bvh_builder::split(const thrust::device_vector<size_t>::iterator remap_begin, 
					    const thrust::device_vector<size_t>::iterator remap_end, 
					    const int depth,
					    thrust::device_vector<size_t>::iterator &obj_split,
					    double &min_cost){
						   
    assert(remap_end -remap_begin >=2);
    int best_axis ;
    min_cost = std::numeric_limits<double>::infinity(); 
    size_t best_obj_split = 1;
    
    size_t num_objs = remap_end -remap_begin;
    
    for(int axis =0; axis < 3; axis++){//consider all 3 axis in turn
        compare_obj_centroid.axis = axis;
        thrust::sort(remap_begin,
                     remap_end,
                     compare_obj_centroid
                     );
        
                     
        //re-init rbound
        {
            aabb rbb;
            for(size_t j = num_objs; j > 0;  ){
				j--;
                rbb.grow(get_object_aabb(remap_begin[j]));
                r_bounds[j] = rbb; 
            }
        }
        
        //re-init lbound
        {
            aabb lbb;
            for(size_t j = 0 ; j< num_objs;  j++){
                lbb.grow(get_object_aabb(remap_begin[j]));
                l_bounds[j] = lbb;                
            }
        }
        
        
        // note r_bounds[0] is the aabb 
        const aabb bb = r_bounds[0];
       
       
        for(size_t i =0; i < num_objs-1; i++){
            size_t NL =  i+1;
            size_t NR =  num_objs - NL; 
            double cost = get_sah_cost(bb,l_bounds[i],r_bounds[i+1],NL,NR);
            if(cost < min_cost){
                best_axis =  axis;
                best_obj_split = i+1;
                min_cost = cost;  
            }
        }
                     
    }
    
    //sort objs again along the best axis you found
    compare_obj_centroid.axis = best_axis;
    thrust::sort(remap_begin,
                 remap_end,
                 compare_obj_centroid
                 );
    
    obj_split =  remap_begin+ best_obj_split;
    
    assert(obj_split>remap_begin);
    
			   
}





aabb bvh_builder::get_aabb_from(
					thrust::device_vector<size_t>::iterator remap_begin,
					thrust::device_vector<size_t>::iterator remap_end
					){
	
	
	//detail::get_tri_min_max taabb(partialInputMesh.vertices,partialInputMesh.faces);
	detail::get_tri_min_max taabb(tri_bounds.begin());
	
	thrust::tuple<vec3d,vec3d>  min_max_v = 
			thrust::transform_reduce(
									 remap_begin,
									 remap_end,
									 taabb,
									 taabb(*remap_begin),
									 detail::min_max()									 
									);


	return aabb(thrust::get<0>(min_max_v),thrust::get<1>(min_max_v));
}


aabb bvh_builder::get_object_aabb(size_t i)  const{
	
	/*const vec3d va = partialInputMesh.vertices[partialInputMesh.faces[i].a];
	const vec3d vb = partialInputMesh.vertices[partialInputMesh.faces[i].b];
	const vec3d vc = partialInputMesh.vertices[partialInputMesh.faces[i].c];	
	const vec3d lowerLeftPoint = cwisemin(cwisemin(va,vb),vc);
	const vec3d upperRightPoint = cwisemax(cwisemax(va,vb),vc);  
	
	return aabb(lowerLeftPoint,upperRightPoint);
*/

	return tri_bounds[i];
}


void bvh_builder::build_remap_table(){
	//count actual the number of tris
	size_t ntris = root->get_ntris();
	grmap.reserve(ntris);
	//fill the grmap with tris and assign to leaves the start and final index
	root->append_tris(grmap);	
}


#include "packed.h"
struct large_node {
//struct __attribute__((__packed__)) large_node {
	
	unsigned char type; //0 inner node,  1 leaf node
	int node0_idx_or_start_idx;	
	int node1_idx_or_end_idx;
	
	double x_min, y_min, z_min;
	double x_max, y_max, z_max;
	
	void operator=(const bvh_node *node){
		if(!node->is_leaf()){			
			type = 0;
			const inner_node* inode  = static_cast<const inner_node*>(node);
			node0_idx_or_start_idx = inode->nodes_idx[0];
			node1_idx_or_end_idx = inode->nodes_idx[1];
			
		}else{			
			type = 1;
			const  leaf_node* lnode  = static_cast<const leaf_node*>(node);
			node0_idx_or_start_idx = lnode->start_idx;
			node1_idx_or_end_idx = lnode->end_idx;
		}
		
		x_min  = node->bb.lowerLeftPoint.x;
		y_min  = node->bb.lowerLeftPoint.y;
		z_min  = node->bb.lowerLeftPoint.z;
		
		x_max  = node->bb.upperRightPoint.x;
		y_max  = node->bb.upperRightPoint.y;
		z_max  = node->bb.upperRightPoint.z;		
		
	}
	
} PACKED;
#include "endpacked.h"

namespace nvcuda_types{

	#include "packed.h"
    //struct __attribute__((__packed__)) nv_int2 {
    struct nv_int2 {
        int32_t x, y;
        //int x, y;
    //};
    } PACKED;
	#include "endpacked.h"

	#include "packed.h"
    //struct __attribute__((__packed__)) nv_float4{
    struct nv_float4 {
        float x, y, z, w;
    //};
    } PACKED;
	#include "endpacked.h"

}


void bvh_builder::write_bvh_ssv(std::ostream &os, OUTPUT_FORMAT oformat){
	
	build_remap_table();
	//count the number of nodes		
	size_t n_nodes = root->count();  
			
	os << "# remapping table size      size of node list " << std::endl;
	os << grmap.size() << " " << n_nodes<< std::endl;
	
	if(oformat == BIN){
		size_t * grmap_ptr = thrust::raw_pointer_cast(grmap.data());
		os.write(reinterpret_cast<char*>(grmap_ptr), sizeof(grmap[0])*grmap.size());
		
	}else if(oformat == NBIN){
		//convert the size_t value of remapping table to int32_t 
		thrust::host_vector<int32_t> igrmap(grmap.begin(),grmap.end());
		int32_t * igrmap_ptr = thrust::raw_pointer_cast(igrmap.data());
		//thrust::host_vector<int> igrmap(grmap.begin(),grmap.end());
		//int * igrmap_ptr = thrust::raw_pointer_cast(igrmap.data());
		os.write(reinterpret_cast<char*>(igrmap_ptr), sizeof(igrmap[0])*igrmap.size());
	}else if(oformat == SSV){
		for (size_t i = 0; i < grmap.size(); i++)	
			os << grmap[i] << " ";
		os << std::endl; 
			
	}
	
	std::stack<bvh_node *> node_stack;
	//create a list of nodes as they has been visited 
	std::vector<bvh_node *> node_list; 
	node_list.reserve(n_nodes);
	
	size_t idx = 0;
	node_stack.push(root);
	
	//TODO modifica qui in modo da ottenere alrti attraversamenti
	while(!node_stack.empty()){
		bvh_node * node  = node_stack.top();
		node->idx = idx;
		node_stack.pop();
		node_list.push_back(node);
		
		if(!node->is_leaf()){
			inner_node* inode  = static_cast<inner_node*>(node);

			node_stack.push(inode->nodes[0]);
			node_stack.push(inode->nodes[1]);
		}
		idx++;
	}	
	assert(node_stack.empty());
	
	
	//update idx in the node structure
	node_stack.push(root);
	while(!node_stack.empty()){
		bvh_node * node  = node_stack.top();
		node_stack.pop();
		if(!node->is_leaf()){
			inner_node* inode  = static_cast<inner_node*>(node);
			inode->nodes_idx[0] = inode->nodes[0]->idx;
			inode->nodes_idx[1] = inode->nodes[1]->idx;
		
			node_stack.push(inode->nodes[0]);
			node_stack.push(inode->nodes[1]);
		}		
	}
	

	//write nodes in the order which the indices are assigned to.
	if(oformat == SSV){
		os << "# inner node format. the first number is the type of node 0 for inner node " << std::endl; 
		os << "# 0 node1_idx node2_idx aabb min (x y z) aabb max (x y z)" << std::endl; 
		os << "# leaf node format. the first number is the type of node  1 for leaf node " << std::endl; 
		os << "# 1 start_idx end_idx aabb min (x y z) aabb max (x y z)" << std::endl; 
	
	
		for(size_t i = 0; i < node_list.size(); i++ ){
			bvh_node * node = node_list[i];
			
			if(!node->is_leaf()){
				inner_node* inode  = static_cast<inner_node*>(node);
				os <<  std::setprecision(15)  
				   <<  0 << " " <<  inode->nodes_idx[0] << " " << inode->nodes_idx[1]  << " " 
				   <<  inode->bb.lowerLeftPoint.x << " " << inode->bb.lowerLeftPoint.y << " " << inode->bb.lowerLeftPoint.z << " " 
				   <<  inode->bb.upperRightPoint.x << " " << inode->bb.upperRightPoint.y << " " << inode->bb.upperRightPoint.z << std::endl; 			   
			}else{
				leaf_node* lnode  = static_cast<leaf_node*>(node);
				 os <<  std::setprecision(15)  
					<<  1 << " " <<  lnode->start_idx << " " << lnode->end_idx   << " " 
					<<  lnode->bb.lowerLeftPoint.x << " " << lnode->bb.lowerLeftPoint.y << " " << lnode->bb.lowerLeftPoint.z << " " 
					<<  lnode->bb.upperRightPoint.x << " " << lnode->bb.upperRightPoint.y << " " << lnode->bb.upperRightPoint.z << std::endl; 			   
			}
		}
	}else if(oformat == BIN){//large node binary format
		thrust::host_vector<large_node>  packed_nodes(node_list.begin(), node_list.end());
		large_node * packed_nodes_ptr = thrust::raw_pointer_cast(packed_nodes.data());
		os.write(reinterpret_cast<char*>(packed_nodes_ptr), sizeof(packed_nodes[0])*packed_nodes.size());
		
	}else if(oformat == NBIN){
		
		//  NVIDIA NODE FORMAT
		//  N0:int2  [c0.inner or ~c0.leaf or begin, c1.inner or ~c1.leaf or end]
		//  N1:real4 [c0.lo.x, c0.hi.x, c0.lo.y, c0.hi.y]
		//  N2:real4 [c1.lo.x, c1.hi.x, c1.lo.y, c1.hi.y]
		//  N3:real4 [c0.lo.z, c0.hi.z, c1.lo.z, c1.hi.z]
		
		
		//TODO Inserire un controllo per l'overflow dei nodi				
		thrust::host_vector<nvcuda_types::nv_int2>   N0(n_nodes);
		thrust::host_vector<nvcuda_types::nv_float4> N1(n_nodes);
		thrust::host_vector<nvcuda_types::nv_float4> N2(n_nodes);
		thrust::host_vector<nvcuda_types::nv_float4> N3(n_nodes);
		
		for(size_t i = 0; i < n_nodes;  i++){
			if(!node_list[i]->is_leaf()){
				const inner_node* inode  = static_cast<const inner_node*>(node_list[i]);
				const int node0_idx = inode->nodes_idx[0];
				const int node1_idx = inode->nodes_idx[1];
				
				N0[i].x = (!node_list[node0_idx]->is_leaf())?node0_idx : (~node0_idx);
				N0[i].y = (!node_list[node1_idx]->is_leaf())?node1_idx : (~node1_idx);
			
				//child nodes 
				bvh_node *cnodes[2]  = {node_list[node0_idx],node_list[node1_idx]};
				
				N1[i].x  = cnodes[0]->bb.lowerLeftPoint.x; 
				N1[i].y  = cnodes[0]->bb.upperRightPoint.x;
				N1[i].z  = cnodes[0]->bb.lowerLeftPoint.y;
				N1[i].w  = cnodes[0]->bb.upperRightPoint.y;

				N2[i].x  = cnodes[1]->bb.lowerLeftPoint.x; 
				N2[i].y  = cnodes[1]->bb.upperRightPoint.x;
				N2[i].z  = cnodes[1]->bb.lowerLeftPoint.y;
				N2[i].w  = cnodes[1]->bb.upperRightPoint.y;

				N3[i].x  = cnodes[0]->bb.lowerLeftPoint.z; 
				N3[i].y  = cnodes[0]->bb.upperRightPoint.z;
				N3[i].z  = cnodes[1]->bb.lowerLeftPoint.z;
				N3[i].w  = cnodes[1]->bb.upperRightPoint.z;
				
			}else {
				const  leaf_node* lnode  = static_cast<const leaf_node*>(node_list[i]);
				int start_idx = lnode->start_idx;
				int end_idx   = lnode->end_idx;
				
				N0[i].x  = start_idx;
				N0[i].y  = end_idx;				
			}
	
			
		}
		nvcuda_types::nv_int2 * N0_ptr = thrust::raw_pointer_cast(N0.data());
		os.write(reinterpret_cast<char*>(N0_ptr), sizeof(N0[0])*N0.size());

		nvcuda_types::nv_float4 * N1_ptr = thrust::raw_pointer_cast(N1.data());
		os.write(reinterpret_cast<char*>(N1_ptr), sizeof(N1[0])*N1.size());

		nvcuda_types::nv_float4 * N2_ptr = thrust::raw_pointer_cast(N2.data());
		os.write(reinterpret_cast<char*>(N2_ptr), sizeof(N2[0])*N2.size());

		nvcuda_types::nv_float4 * N3_ptr = thrust::raw_pointer_cast(N3.data());
		os.write(reinterpret_cast<char*>(N3_ptr), sizeof(N3[0])*N3.size());		
		
		
	}	
	
	os << std::endl;
}



