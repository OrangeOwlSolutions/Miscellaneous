#ifndef __BVH_H__
#define __BVH_H__

#include "aabb.h"
#include "Mesh.h"
#include <vector>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

enum OUTPUT_FORMAT {SSV,BIN ,NBIN};

/***********************/
/* BVH PLATFORM STRUCT */
/***********************/
struct bvh_platform {
    
	float KT;						// --- Triangle intersection cost
    float KI;						// --- AABB intersection cost
    int max_depth;					// --- Maximum tree depth
    bool enable_spatial_split;
    
    bvh_platform(float KT_ = 15, float KI_ = 20, int max_depth_ = 16, bool enable_spatial_split_ = false):
			KT(KT_), KI(KI_), max_depth(max_depth_), enable_spatial_split(enable_spatial_split_) {};
};


struct bvh_node{
	size_t idx;
	aabb bb;
	bvh_node(const aabb &bb_):idx(0),bb(bb_){}
	virtual bool is_leaf() const  = 0;
	virtual void print(int depth =0 ) const = 0;
	virtual size_t get_ntris() const = 0;
	virtual void append_tris(thrust::host_vector<size_t> & grmap) = 0;	
	virtual size_t count() = 0;
	
};

struct inner_node: public bvh_node{
	bvh_node * nodes[2];
	size_t   nodes_idx[2];
	inner_node(const aabb &bb,
			   bvh_node * node1,
	           bvh_node * node2);
	
	virtual bool is_leaf() const;
	
	virtual void print(int depth) const ;

	virtual size_t get_ntris() const;

	virtual void append_tris(thrust::host_vector<size_t> & grmap);	
	
	virtual size_t count();
	
};

struct leaf_node: public bvh_node{
	thrust::host_vector<size_t> tri_list;
	size_t start_idx, end_idx;
	
	
	leaf_node(const aabb &bb);
	
	virtual bool is_leaf() const;
	
	virtual void print(int depth) const ;

	virtual size_t get_ntris() const;

	virtual void append_tris(thrust::host_vector<size_t> & grmap);	

	virtual size_t count();
};



struct bvh_builder{
	bvh_platform platform;
	const BVHmesh & partialInputMesh;
	bvh_node *root;
	size_t num_tris;
	thrust::host_vector<size_t> grmap;//global remapping_table
	thrust::device_vector<aabb>   r_bounds;
    thrust::device_vector<aabb>   l_bounds;
	thrust::device_vector<aabb>   tri_bounds;
	
	bvh_builder(bvh_platform platform,
				const BVHmesh & partialInputMesh);
				
	void build();
	
	bvh_node* rec_build(thrust::device_vector<size_t>::iterator remap_begin,
					   thrust::device_vector<size_t>::iterator remap_end,
					   int depth);

	void write_bvh_ssv(std::ostream &os, OUTPUT_FORMAT oformat );


private:

    struct cmp_obj_centroid: public thrust::binary_function<
			thrust::device_vector<size_t>::iterator,
			thrust::device_vector<size_t>::iterator,bool>{
        
        const bvh_builder &parent_ref;
        int axis;    
        cmp_obj_centroid(const bvh_builder &p_ref):parent_ref(p_ref),axis(0){}
    
        bool operator()(size_t obj1_idx, 
						size_t obj2_idx);        
    }compare_obj_centroid;


	void build_remap_table();

	aabb get_aabb_from(thrust::device_vector<size_t>::iterator remap_begin,
					   thrust::device_vector<size_t>::iterator remap_end);
				
    aabb get_object_aabb(size_t i)  const ;		
					   
	double get_sah_cost(const aabb &bb, 
                                 const aabb &lbb,
                                 const aabb &rbb, 
                                 size_t NL, size_t NR) const;
                                 
	void split(const thrust::device_vector<size_t>::iterator obj_begin, 
			   const thrust::device_vector<size_t>::iterator obj_end, 
			   const int depth,
			   thrust::device_vector<size_t>::iterator &best_obj_split,
			   double &min_cost);
			   
			   
			   
};


aabb compute_aabb(const std::vector<vec3d> &vertices,
				  const std::vector<vec4i> &faces);



	

#endif

