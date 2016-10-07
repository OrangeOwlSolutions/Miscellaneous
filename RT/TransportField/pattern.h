#ifndef __PATTERN_H__
#define __PATTERN_H__

#include <string.h>
#include <vector>
#include <iostream>

#include "cfloat.h"

struct pattern{

    std::string name;

    cfloat Zin;

    size_t num_phi;
    size_t num_theta;
    
	std::vector<double> phi;
    std::vector<double> theta;   
    
	std::vector<cfloat> h_theta;
    std::vector<cfloat> h_phi;    
        
    pattern(const std::string &name);
    
};

std::ostream & operator << (std::ostream & os, const pattern &p );

#endif
