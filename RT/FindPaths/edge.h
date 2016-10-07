//#ifndef __EDGE_H__
//#define __EDGE_H__
//
//#include <iostream>
//#include <vector>
//#include <string>
//
//#include "vec3.h"
//
///***************/
///* EDGE STRUCT */
///***************/
//// --- Needed for diffraction. Here, each edge is represented by two points.
//struct edge{
//    vec3 a;
//    vec3 b;
//    unsigned int ID;
//};
//
//void load_edge_list(const std::string &filename, const unsigned int first_ID, std::vector<edge> &edge_list);
//
//std::ostream &operator << (std::ostream & os, const edge& e);
//
//bool get_diffraction_point(const vec3 &tx_pos,const vec3 &rx_pos, const edge &ed, vec3 &diff_point);
//
//#endif
#ifndef __EDGE_H__
#define __EDGE_H__

#include <iostream>
#include <vector>
#include <string>

#include "vec3.h"

/***************/
/* EDGE STRUCT */
/***************/
// --- Needed for diffraction. Here, each edge is represented by two points.
struct edge {
    vec3 a;
    vec3 b;
    unsigned int ID;

    vec3   no;
    vec3   nn;
    double n;
};

void load_edge_list(const std::string &filename, const unsigned int first_ID, std::vector<edge> &edge_list);

std::ostream &operator << (std::ostream & os, const edge& e);

bool get_diffraction_point(const vec3 &tx_pos,const vec3 &rx_pos, const edge &ed, vec3 &diff_point);

#endif
