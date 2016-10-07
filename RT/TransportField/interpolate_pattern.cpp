#include <cassert>

#include "interpolate_pattern.h"

/**************************/
/* LAGRANGE INTERPOLATION */
/**************************/
cfloat u(cfloat a, cfloat b, float s){
    return a + (b-a)*s;
}

/***********************/
/* INTERPOLATE PATTERN */
/***********************/
void interpolate_pattern(const double phi, const double theta, const pattern *patt, cfloat &h_theta, cfloat &h_phi, vec3 &theta_v, vec3 &phi_v) {
    
	// --- Unit vectors i_theta and i_phi in the transmitting antenna reference system
	theta_v = normalize(make_vec3(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta)));
    phi_v	= normalize(make_vec3(            -sin(phi),              cos(phi),          0));

    size_t theta_idx =0;
    size_t phi_idx =0;
    
    // --- Finds the closest "tabulated phi" to the required value of phi
	for (size_t i = 0; i < patt->num_phi; i++) {
        if (patt -> phi[i * patt -> num_theta] > phi) {
            phi_idx = i - 1;
            break;
        }
    }    
    
    assert(phi_idx >= 0 && phi_idx < patt -> num_phi - 1);
    
    // --- Finds the closest "tabulated theta" to the required value of theta
    for (size_t i = 0; i < patt -> num_theta; i++) {
        if (patt -> theta[i] > theta) {
            theta_idx = i - 1;
            break;
        }
    }
    
	assert(theta_idx >= 0 && theta_idx < patt -> num_theta - 1);

	double phi_u	= (phi   - patt -> phi[phi_idx*patt -> num_theta]) / (patt -> phi[(phi_idx + 1) * patt -> num_theta]- patt->phi[(phi_idx) * patt -> num_theta]);
    double theta_u	= (theta - patt -> theta[theta_idx]              ) / (patt -> theta[theta_idx + 1]                  - patt -> theta[theta_idx]);
    
    cfloat a_theta, b_theta, c_theta, d_theta;
    cfloat a_phi, b_phi, c_phi, d_phi;
    
    a_theta = patt -> h_theta[theta_idx +      phi_idx      * patt -> num_theta];
    b_theta = patt -> h_theta[theta_idx + 1 +  phi_idx      * patt -> num_theta];
    c_theta = patt -> h_theta[theta_idx +     (phi_idx + 1) * patt -> num_theta];
    d_theta = patt -> h_theta[theta_idx + 1 + (phi_idx + 1) * patt -> num_theta];
    
    a_phi   = patt -> h_phi  [theta_idx   +    phi_idx      * patt -> num_theta];
    b_phi   = patt -> h_phi  [theta_idx + 1 +  phi_idx      * patt -> num_theta];
    c_phi   = patt -> h_phi  [theta_idx   +   (phi_idx + 1) * patt -> num_theta];
    d_phi   = patt -> h_phi  [theta_idx + 1 + (phi_idx + 1) * patt -> num_theta];    
    
    h_theta = u(u(a_theta, b_theta, theta_u), u(c_theta, d_theta, theta_u), phi_u);
    h_phi	= u(u(a_phi, b_phi, theta_u), u(c_phi, d_phi, theta_u), phi_u); 
 
}
 
