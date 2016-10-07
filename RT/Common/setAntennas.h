#ifndef __SETANTENNAS_H__
#define __SETANTENNAS_H__

#include <vector>

#include "vec3.h"
#include "RXTX.h"

/***********/
/* SEGMENT */
/***********/
struct segment{
    vec3 orig;		// --- Segment middle point
    vec3 dir_v;		// --- Segment direction (unit vector)
    float L;		// --- Segment half-length
    segment(vec3 orig_, vec3 dir_v_, float L_): orig(orig_), dir_v(normalize(dir_v_)), L(L_) {}
        
    segment() {};    
};

extern "C" void setAntennas(segment *, segment *, vec3 origintx, vec3 originrx, std::vector<vec3> &, std::vector<vec3> &, TX *, RX *, 
	                        std::vector<std::string> &, std::vector<std::string> &, vec3, vec3, float, 
							std::string, std::string, const int, const int);

void set_pattern(const std::vector<std::string> &, std::string *, unsigned int);

extern "C" void set_azr(const std::vector<vec3> &, vec3 *, unsigned int);

#endif
