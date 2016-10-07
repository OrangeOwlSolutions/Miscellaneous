//#include "ray_path.h"
//
///****************************/
///* DISPLAY PATH FUNCTION << */
///****************************/
//std::ostream &operator << (std::ostream &os, const path &p) {
//    
//	os << "code.x = " << p.code.x << " code.y = " << p.code.y << " code.z = " << p.code.z << " code.w = " << p.code.w << " p.tx = " << p.tx << " p.rx = " << p.rx << " nr intersections = " << p.intersections.size() << std::endl;
//    
//    for(size_t i = 0; i < p.intersections.size(); i++) {
//        
//		path_inters pi = p.intersections[i];
//        
//		if (pi.type == REFLECTION) { os << " type = " << pi.type << " pi.pos.x = " << pi.pos.x << " pi.pos.y = " << pi.pos.y << " pi.pos.z = " << pi.pos.z << " pi.ID = " << pi.ID << " pi.u = " << pi.u << " pi.v = " << pi.v; }
//		else if (pi.type == DIFFRACTION) { os << " type = " << pi.type << " pi.pos.x = " << pi.pos.x << " pi.pos.y = " << pi.pos.y << " pi.pos.z = " << pi.pos.z << " pi.ID = " << pi.ID; }
//		else { os << " type = " << pi.type << " pi.pos.x = " << pi.pos.x << " pi.pos.y = " << pi.pos.y << " pi.pos.z = " << pi.pos.z; }
//        if(i != p.intersections.size() - 1) os << std::endl;  
//    }   
//    
//    return os;
//} 
//
///*****************************/
///* COMPARE RAY PATH OPERATOR */
///*****************************/
//bool compare_ray_path::operator() (const path &p1, const path &p2) {
//
//    unsigned int code1[4] = {p1.code.x, p1.code.y, p1.code.z, p1.code.w};
//    unsigned int code2[4] = {p2.code.x, p2.code.y, p2.code.z, p2.code.w};    
//    
//    for(int i = 0 ; i < 4; i++) if (code1[i] != code2[i]) return code1[i] < code2[i];
//    
//    // --- Here the rays are "the same", so compare the distance from detector
//	return  p1.dist < p2.dist;
//    
//}
//
///*********************************/
///* IS RAY PATH NOT GOOD OPERATOR */
///*********************************/
//bool is_ray_path_not_good::operator() (const path &p) { return !p.good; }
//
///******************************/
///* IS RAY PATH EQUAL OPERATOR */
///******************************/
//bool is_ray_path_equal::operator() (const path &p1, const path &p2) { 
//	return p1.code.x == p2.code.x && p1.code.y == p2.code.y && p1.code.z == p2.code.z && p1.code.w == p2.code.w;
//}
//
#include "ray_path.h"

/**********************************/
/* GET NUMBER OF REFLECTED FIELDS */
/**********************************/
//int path::get_num_refl() const { return code[2] & 0xffff; }
int path::get_num_refl() const { return code.z & 0xffff; }

/***********************************/
/* GET NUMBER OF DIFFRACTED FIELDS */
/***********************************/
//int path::get_num_diff() const { return code[2] >> 16; }
int path::get_num_diff() const { return code.z >> 16; }

/****************************/
/* DISPLAY PATH FUNCTION << */
/****************************/
std::ostream &operator << (std::ostream &os, const path &p) {
    
	os << "code.x = " << p.code.x << " code.y = " << p.code.y << " code.z = " << p.code.z << " code.w = " << p.code.w << " p.tx = " << p.tx << " p.rx = " << p.rx << " nr intersections = " << p.intersections.size() << std::endl;
    
    for(size_t i = 0; i < p.intersections.size(); i++) {
        
		path_inters pi = p.intersections[i];
        
		if (pi.type == REFLECTION) { os << " type = " << pi.type << " pi.pos.x = " << pi.pos.x << " pi.pos.y = " << pi.pos.y << " pi.pos.z = " << pi.pos.z << " pi.ID = " << pi.ID << " pi.u = " << pi.u << " pi.v = " << pi.v; }
		else if (pi.type == DIFFRACTION) { os << " type = " << pi.type << " pi.pos.x = " << pi.pos.x << " pi.pos.y = " << pi.pos.y << " pi.pos.z = " << pi.pos.z << " pi.ID = " << pi.ID; }
		else { os << " type = " << pi.type << " pi.pos.x = " << pi.pos.x << " pi.pos.y = " << pi.pos.y << " pi.pos.z = " << pi.pos.z; }
        if(i != p.intersections.size() - 1) os << std::endl;  
    }   
    
    return os;
} 

/*****************************/
/* COMPARE RAY PATH OPERATOR */
/*****************************/
bool compare_ray_path::operator() (const path &p1, const path &p2) {

    unsigned int code1[4] = {p1.code.x, p1.code.y, p1.code.z, p1.code.w};
    unsigned int code2[4] = {p2.code.x, p2.code.y, p2.code.z, p2.code.w};    
    
    for(int i = 0 ; i < 4; i++) if (code1[i] != code2[i]) return code1[i] < code2[i];
    
    // --- Here the rays are "the same", so compare the distance from detector
	return  p1.dist < p2.dist;
    
}

/*********************************/
/* IS RAY PATH NOT GOOD OPERATOR */
/*********************************/
bool is_ray_path_not_good::operator() (const path &p) { return !p.good; }

/******************************/
/* IS RAY PATH EQUAL OPERATOR */
/******************************/
bool is_ray_path_equal::operator() (const path &p1, const path &p2) { 
	return p1.code.x == p2.code.x && p1.code.y == p2.code.y && p1.code.z == p2.code.z && p1.code.w == p2.code.w;
}

