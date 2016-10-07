#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>

#include "interpolate_pattern.h"
#include "pattern.h"
#include "utils.h"
#include "Utilities.cuh"

pattern::pattern(const std::string &name_) { 

    name = name_;
    std::string line;
    std::ifstream ifs(name.c_str());

    { // --- Read the input impedence 
        get_valid_line(ifs, line);
        std::stringstream ss(line);    
        ss >> REAL(Zin) >> IMAG(Zin);
    }
	
    { // --- Read the number of azimuth points 
        get_valid_line(ifs, line);
        std::stringstream ss(line);    
        ss >> num_theta >> num_phi;
    }
    
    for (size_t i = 0; i < num_phi * num_theta; i++) {
        
		get_valid_line(ifs, line);
        std::stringstream ss(line);    
        
        double theta_, phi_;
        cfloat h_theta_, h_phi_;
        
        ss >> theta_ >> phi_;
        ss >> REAL(h_theta_) >> IMAG(h_theta_) >> REAL(h_phi_) >> IMAG(h_phi_); 
        
        theta.push_back(deg2rad(theta_));
        phi.push_back(deg2rad(phi_));
        h_theta.push_back(h_theta_);
        h_phi.push_back(h_phi_);
    
    }

}


std::ostream & operator << (std::ostream & os, const pattern &p) {
    
	os << REAL(p.Zin) << " " << IMAG(p.Zin) << std::endl;
    
    return os;
}
