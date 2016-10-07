//#include "transport_ray_field.h"
//#include "complex_vec3.h"
//#include "interpolate_pattern.h"
//
///***************************************/
///* COMPUTE ORTHONORMAL REFERENCE FRAME */
///***************************************/
//void compute_ortonormal_reference_frame(vec3 &x_v, vec3 &y_v, const vec3 &z_v) {
//
//	vec3 ax[3] = { make_vec3(1, 0, 0), make_vec3(0, 1, 0), make_vec3(0, 0, 1) };
//
//    float min_dot = fabsf(dot(z_v, ax[0]));
//    
//	int min_axis = 0;
//    
//	for(int i = 1; i < 3; i++) {
//        
//		float dot_prod = fabsf(dot(z_v, ax[i])); 
//        
//		if(dot_prod < min_dot){
//            min_dot = dot_prod;
//            min_axis = i;
//        }
//    }
//    
//	// --- ax[min_axis] is the axis "mostly orthogonal" to z_v. x_v and y_v are constructed to be a plane "mostly orthogonal" to z_v.
//	y_v = normalize(cross(z_v, ax[min_axis]));
//    x_v = normalize(cross(y_v, z_v));    
//
//} 
//
///**********************/
///* GENERATE FIRST RAY */
///**********************/
//void gen_first_ray(const TX &tx, const vec3 &si_v, cfloat &Ax1m, cfloat &Ax2m, vec3 &x1m, vec3 &x2m, double &rho_m, double lambda) {
//    
//    
//    cfloat J = make_cfloat(0, 1);					// --- Imaginary unit
//    
//	// --- Components of the unit vector si_v pointing from/to the first two path points in the reference system of the transmitting antenna
//	vec3 local_si_v;
//    local_si_v.x = dot(si_v, tx.x_v);
//    local_si_v.y = dot(si_v, tx.y_v);
//    local_si_v.z = dot(si_v, tx.z_v);
//
//	//std::cout << tx.z_v.x << " " << tx.z_v.y << " " << tx.z_v.z << "\n";
//
//	// --- Azimuth and zenith angles of the second path point in the transmitting antenna reference system
//	double azimuth   = atan2(local_si_v.y, local_si_v.x);
//    double zenith	 = acos(local_si_v.z);
//
//	//std::cout << azimuth << " " << zenith << "\n";
//	
//	// --- Components of the antenna effective length
//	cfloat h_theta;
//    cfloat h_phi;
//
//    // --- Calculates the antenna effective length in the direction of the second point, when the origin is at the first point
//    interpolate_pattern(azimuth, zenith, tx.patt, h_theta, h_phi, x1m, x2m);
//    
//	//std::cout << h_theta.x << " " << h_theta.y << "\n";
//	//std::cout << h_phi.x << " " << h_phi.y << "\n";
//
//	cfloat Z0  = make_cfloat(50.0, 0);						// --- Characteristic impedance of the line      
//    cfloat Zin = tx.patt->Zin;								// --- Input impedance of the antenna
//    cfloat Gamma = (Zin - Z0) / (Zin + Z0);					// --- Reflection coefficient at the antenna terminals
//    
//    x1m = x1m.x*tx.x_v + x1m.y*tx.y_v+x1m.z*tx.z_v;
//    x2m = x2m.x*tx.x_v + x2m.y*tx.y_v+x2m.z*tx.z_v;
//
//    cfloat I =  2.0 / (Zin + Z0);							// --- Input current
//    
//	Ax1m = h_theta * J * ETA_0 * I / (2 * lambda);			// --- theta component of the field impinging on the second path point
//    Ax2m = h_phi   * J * ETA_0 * I / (2 * lambda);			// --- phi   component of the field impinging on the second path point
//
//	//std::cout << Ax2m.x << " " << Ax2m.y << "\n";
//
//	rho_m = 0;
//}
//
///*******************/
///* FIELD TRANSPORT */
///*******************/
//extern "C" cfloat transport_ray_field(std::vector<vec3> positions, TX &tx, RX &rx, double lambda){
//
//    // --- Positions is a vector containing the positions of the emission/reflection/diffraction/detection points of the generic path
//	size_t num_points = positions.size();
//                
//    cfloat J = make_cfloat(0, 1);				// --- Immaginary unit
//    
//    cfloat H = make_cfloat(0, 0);				// --- Output initialization
//    
//    cfloat Ax1m, Ax1p;
//    cfloat Ax2m, Ax2p;
//       
//    // --- x1m and x2m are the unit vectors i_theta and i_phi in the reference system defined by the tx antenna
//	vec3 x1m, x1p;
//    vec3 x2m, x2p;
//    
//    double rho_m, rho_p;
// 
//	{	// --- Generate the first ray
//		double s0 = dist(positions[1], positions[0]);					// --- Distance between the first two path points
//        vec3 s0_v = normalize(positions[1] - positions[0]);				// --- Unit vector pointing from the first to the second path points
//        gen_first_ray(tx, s0_v, Ax1m, Ax2m, x1m, x2m, rho_m, lambda);
//	    //std::cout << s0 << " " << s0_v.x << " " << s0_v.y << " " << s0_v.z << "\n";
//	    //std::cout << lambda << " " << Ax1m.x << " " << Ax1m.y << " " << Ax2m.x << " " << Ax2m.y << " " << "\n";
//    }
//    
//	
//	// --- Transport the field
//    for (size_t i = 0; i < num_points-1; i++) {
//        
//		double s  = dist(positions[i+1], positions[i]);					// --- Distance between two consecutive path points
//        vec3 si_v = normalize(positions[i + 1] - positions[i]);			// --- Unit vector pointing from one path point to the next 
//																		//     (incident propagation unit vector)
//        
//        // --- Calculates the propagation factor
//        cfloat F;        
//        if (rho_m == 0.0) F = (1 / s) * expf(-J * 2 * PI_R / lambda * s);
//		else			  F = (rho_m / (rho_m + s)) * expf(-J * 2 * PI_R / 0.125 * s);
//        
//        // --- Propagates the field
//        Ax1p = Ax1m * F;
//        Ax2p = Ax2m * F;
//        
//		x1p = x1m;
//        x2p = x2m;
//        rho_p = rho_m + s;												// --- Updates the path length
//        
//        // --- Performs reflections
//        if (i < num_points-2){
//            
//            vec3 sr_v = normalize(positions[i + 2] - positions[i + 1]);	// --- Reflection propagation unit vector
//                
//            vec3 normal_v = make_vec3(0, 0, 1);
//
//            vec3 en, epi, epr;
//
//            // --- Computes the reflected field
//            if (fabs(dot(si_v, normal_v)) < 0.99998) { // --- Oblique incidence
//                  en  = normalize(cross(normal_v, si_v));
//                  epi = normalize(cross(si_v, en));
//                  epr = normalize(cross(sr_v, en));
//            } else {
//                  compute_ortonormal_reference_frame(en, epi, normal_v); 
//                  epr = -epi;
//            }
//            
//            const cfloat E0in = Ax1p * dot(x1p, en ) + Ax2p * dot(x2p, en );		// --- Normal component of the incident field
//            const cfloat E0ip = Ax1p * dot(x1p, epi) + Ax2p * dot(x2p, epi);		// --- Parallel component of the incient field
//
//            // --- Use PEC
//            const cfloat E0rn = -E0in;												// --- 
//            const cfloat E0rp =  E0ip;
//
//            rho_m = rho_p;
//            x1m   = en;
//            x2m   = epr;
//            Ax1m  = E0rn;
//            Ax2m  = E0rp;
//        }
//    }
//    
//    // --- Calculates the scattering parameter
//    {
//        cvec3 Ei = Ax1p * x1p + Ax2p * x2p;
//
//        vec3 sr_v = normalize(positions[num_points-2] - positions[num_points-1]);
//        
//        vec3 local_sr_v;
//        local_sr_v.x = dot(sr_v, rx.x_v);
//        local_sr_v.y = dot(sr_v, rx.y_v);
//        local_sr_v.z = dot(sr_v, rx.z_v);
//
//        double azimuth = atan2(local_sr_v.y, local_sr_v.x);
//        double zenith  = acos (local_sr_v.z);
//        
//        cfloat h_theta_r,h_phi_r;
//        vec3   theta_v, phi_v;
//
//        interpolate_pattern(azimuth, zenith, rx.patt, h_theta_r, h_phi_r, theta_v, phi_v);          
//        
//        theta_v = theta_v.x * rx.x_v + theta_v.y * rx.y_v + theta_v.z * rx.z_v;
//        phi_v   = phi_v.x   * rx.x_v + phi_v.y   * rx.y_v + phi_v.z   * rx.z_v;
//
//        cvec3 hr = h_theta_r * theta_v + h_phi_r * phi_v;
//        
//        cfloat Zin = rx.patt -> Zin;
//        cfloat Z0  = make_cfloat(50.0, 0);
//        
//        H = edot(Ei, hr) * (Z0 / (Z0 + Zin));
//    }
//    
//    return H;
//}
#include "transport_ray_field.h"
#include "complex_vec3.h"
#include "interpolate_pattern.h"

/***************************************/
/* COMPUTE ORTHONORMAL REFERENCE FRAME */
/***************************************/
void compute_ortonormal_reference_frame(vec3 &x_v, vec3 &y_v, const vec3 &z_v) {

	vec3 ax[3] = { make_vec3(1, 0, 0), make_vec3(0, 1, 0), make_vec3(0, 0, 1) };

    float min_dot = fabsf(dot(z_v, ax[0]));
    
	int min_axis = 0;
    
	for(int i = 1; i < 3; i++) {
        
		float dot_prod = fabsf(dot(z_v, ax[i])); 
        
		if(dot_prod < min_dot){
            min_dot = dot_prod;
            min_axis = i;
        }
    }
    
	// --- ax[min_axis] is the axis "mostly orthogonal" to z_v. x_v and y_v are constructed to be a plane "mostly orthogonal" to z_v.
	y_v = normalize(cross(z_v, ax[min_axis]));
    x_v = normalize(cross(y_v, z_v));    

} 

/**********************/
/* GENERATE FIRST RAY */
/**********************/
void gen_first_ray(const TX &tx, const vec3 &si_v, cfloat &Ax1m, cfloat &Ax2m, vec3 &x1m, vec3 &x2m, double lambda) {
    
    cfloat J = make_cfloat(0, 1);					// --- Imaginary unit

	// --- Components of the unit vector si_v pointing from/to the first two path points in the reference system of the transmitting antenna
	vec3 local_si_v;
    local_si_v.x = dot(si_v, tx.x_v);
    local_si_v.y = dot(si_v, tx.y_v);
    local_si_v.z = dot(si_v, tx.z_v);

	// --- Azimuth and zenith angles of the second path point in the transmitting antenna reference system
    double azimuth = atan2(local_si_v.y, local_si_v.x);
    double zenith  = acos(local_si_v.z);

	// --- Components of the antenna effective length
    cfloat h_theta;
    cfloat h_phi;

	// --- Calculates the antenna effective length in the direction of the second point, when the origin is at the first point
    interpolate_pattern(azimuth, zenith, tx.patt, h_theta, h_phi, x1m, x2m);
    
	cfloat Z0  = make_cfloat(50.0, 0);						// --- Characteristic impedance of the line     
    cfloat Zin = tx.patt->Zin;								// --- Input impedance of the antenna
    cfloat Gamma = (Zin - Z0) / (Zin + Z0);					// --- Reflection coefficient at the antenna terminals
    
    x1m = x1m.x * tx.x_v + x1m.y * tx.y_v + x1m.z * tx.z_v;
    x2m = x2m.x * tx.x_v + x2m.y * tx.y_v + x2m.z * tx.z_v;

    cfloat I =  2.0 / (Zin + Z0);							// --- Input current
    
    Ax1m = h_theta * J * ETA_0 * I / (2 * lambda); 			// --- theta component of the field impinging on the second path point
    Ax2m = h_phi   * J * ETA_0 * I / (2 * lambda);			// --- phi   component of the field impinging on the second path point

}

//void gen_first_ray(const TX &tx, const vec3 &si_v, cfloat &Ax1m, cfloat &Ax2m, vec3 &x1m, vec3 &x2m, double &rho_m, double lambda) {
//    
//    
//    cfloat J = make_cfloat(0, 1);					// --- Imaginary unit
//    
//	// --- Components of the unit vector si_v pointing from/to the first two path points in the reference system of the transmitting antenna
//	vec3 local_si_v;
//    local_si_v.x = dot(si_v, tx.x_v);
//    local_si_v.y = dot(si_v, tx.y_v);
//    local_si_v.z = dot(si_v, tx.z_v);
//
//	//std::cout << tx.z_v.x << " " << tx.z_v.y << " " << tx.z_v.z << "\n";
//
//	// --- Azimuth and zenith angles of the second path point in the transmitting antenna reference system
//	double azimuth   = atan2(local_si_v.y, local_si_v.x);
//    double zenith	 = acos(local_si_v.z);
//
//	//std::cout << azimuth << " " << zenith << "\n";
//	
//	// --- Components of the antenna effective length
//	cfloat h_theta;
//    cfloat h_phi;
//
//    // --- Calculates the antenna effective length in the direction of the second point, when the origin is at the first point
//    interpolate_pattern(azimuth, zenith, tx.patt, h_theta, h_phi, x1m, x2m);
//    
//	//std::cout << h_theta.x << " " << h_theta.y << "\n";
//	//std::cout << h_phi.x << " " << h_phi.y << "\n";
//
//	cfloat Z0  = make_cfloat(50.0, 0);						// --- Characteristic impedance of the line      
//    cfloat Zin = tx.patt->Zin;								// --- Input impedance of the antenna
//    cfloat Gamma = (Zin - Z0) / (Zin + Z0);					// --- Reflection coefficient at the antenna terminals
//    
//    x1m = x1m.x*tx.x_v + x1m.y*tx.y_v+x1m.z*tx.z_v;
//    x2m = x2m.x*tx.x_v + x2m.y*tx.y_v+x2m.z*tx.z_v;
//
//    cfloat I =  2.0 / (Zin + Z0);							// --- Input current
//    
//    Ax1m = h_theta * J * ETA_0 * I / (2 * lambda);			// --- theta component of the field impinging on the second path point
//    Ax2m = h_phi   * J * ETA_0 * I / (2 * lambda);			// --- phi   component of the field impinging on the second path point
//
//	//std::cout << Ax2m.x << " " << Ax2m.y << "\n";
//
//	rho_m = 0;
//}

/**************************************/
/* COMPUTING THE SCATTERING PARAMETER */
/**************************************/
static cfloat compute_scattering_parameter(const cvec3 &Ei, const vec3 &si_v, RX &rx) {
        
    vec3 local_si_v;
    local_si_v.x = dot(si_v, rx.x_v);
    local_si_v.y = dot(si_v, rx.y_v);
    local_si_v.z = dot(si_v, rx.z_v);

    double azimuth = atan2(local_si_v.y, local_si_v.x);
    double zenith  = acos (local_si_v.z);
    
    cfloat h_theta_r, h_phi_r;
    vec3   theta_v, phi_v;

    interpolate_pattern(azimuth, zenith, rx.patt, h_theta_r, h_phi_r, theta_v, phi_v);          

    theta_v = theta_v.x * rx.x_v + theta_v.y * rx.y_v + theta_v.z * rx.z_v;
    phi_v	= phi_v.x   * rx.x_v + phi_v.y   * rx.y_v + phi_v.z   * rx.z_v;
    
    cvec3 hr = h_theta_r * theta_v + h_phi_r * phi_v;
    
    cfloat Zin = rx.patt->Zin;
    cfloat Z0  = make_cfloat(50.0, 0);
         
    return edot(Ei, hr) * (Z0 / (Z0 + Zin));
}

/********************************************/
/* FIELD TRANSPORT FOR NON-DIFFRACTIVE RAYS */
/********************************************/
cfloat transport_ray_field_generic_reflection(const path &p, TX &tx, RX &rx, const nvmesh &mesh, double lambda) {

    // --- p.intersections contains the positions of the emission/reflection/diffraction/detection points of the generic path
    size_t num_points = p.intersections.size();
                
    cfloat J = make_cfloat(0, 1);				// --- Immaginary unit
    
    cfloat H = make_cfloat(0, 0);				// --- Output initialization
    
    cfloat Ax1m, Ax2m, Ax1p, Ax2p;
       
    // --- x1m and x2m are the unit vectors i_theta and i_phi in the reference system defined by the tx antenna
	vec3 x1m, x2m, x1p, x2p;
    
    double rho1_m, rho2_m, rho1_p, rho2_p;
   
    {	// --- Generate the first ray
        double s0 = dist(p.intersections[1].pos, p.intersections[0].pos);			// --- Distance between the first two path points
        vec3 s0_v = normalize(p.intersections[1].pos - p.intersections[0].pos);		// --- Unit vector pointing from the first to the second path points
        gen_first_ray(tx, s0_v, Ax1m, Ax2m, x1m, x2m, lambda);
        
        rho1_m = 0.;
        rho2_m = 0.;
    }
    
	// --- Transport the field
    for(size_t i = 0; i < num_points - 1; i++) {
        
		double s = dist(p.intersections[i + 1].pos, p.intersections[i].pos);			// --- Distance between two consecutive path points
        vec3 si_v = normalize(p.intersections[i + 1].pos - p.intersections[i].pos);		// --- Unit vector pointing from one path point to the next 
																						//     (incident propagation unit vector)
        
        // --- Calculates the propagation factor
        cfloat F1, F2;
        
        if (rho1_m == 0.0)		 F1 = sqrtf(make_cfloat(1. / s, 0.));
		else if (rho1_m < INF_R) F1 = sqrtf(make_cfloat(rho1_m / (rho1_m + s), 0.));
        else					 F1 = make_cfloat(1., 0.);

        if (rho2_m == 0.0)		 F2 = sqrtf(make_cfloat(1. / s, 0.));
        else if (rho2_m < INF_R) F2 = sqrtf(make_cfloat(rho2_m / (rho2_m + s), 0.));
        else					 F2 = make_cfloat(1., 0.);
        
        F1 = F1;
        F2 = F2;
        cfloat EXPF = expf(-J * 2. * PI_R / lambda * s);
        
        // --- Propagates the field
        Ax1p = Ax1m * F1 * F2 * EXPF;
        Ax2p = Ax2m * F1 * F2 * EXPF;
        
        x1p = x1m;
        x2p = x2m;
        rho1_p = rho1_m + s;
        rho2_p = rho2_m + s;										// --- Updates the path length
        

        // --- Performs reflections
        if (i < num_points - 2) {
            
            vec3 sr_v = normalize(p.intersections[i + 2].pos - p.intersections[i + 1].pos);	// --- Reflection propagation unit vector			
            
            vec3 normal_v, U1_v, U2_v;
            
			float a1, a2;
            
			{ // --- Get the normal
                assert(p.intersections[i+1].type == REFLECTION);
                int ID	= p.intersections[i + 1].ID;
                float u = p.intersections[i + 1].u;
                float v = p.intersections[i + 1].v;
                mesh.get_normal_and_curvature(ID, u, v, normal_v, U1_v, U2_v, a1, a2);
            }

            const float cos_theta_i = -dot(normal_v, si_v);

            const float Gamma11 =   dot(x1p, U1_v);
            const float Gamma12 =   dot(x1p, U2_v);
            const float Gamma21 =   dot(x2p, U1_v);
            const float Gamma22 =   dot(x2p, U2_v);
            
            const float Gamma11_sqr = Gamma11 * Gamma11; 
            const float Gamma12_sqr = Gamma12 * Gamma12;
            const float Gamma21_sqr = Gamma21 * Gamma21;
            const float Gamma22_sqr = Gamma22 * Gamma22;                        

            const float det_Gamma		= Gamma11 * Gamma22 - Gamma12 * Gamma21;
            const float det_Gamma_sqr	= det_Gamma * det_Gamma;
            
            const float tmp_a = cos_theta_i / det_Gamma_sqr;
            
            const float tmp_b = (Gamma22_sqr + Gamma12_sqr) / a1 + (Gamma21_sqr + Gamma11_sqr) / a2;

            const float tmp_c = pow(1 / rho1_p - 1 / rho2_p, 2);
            
			const float tmp_d = tmp_a * tmp_a;
            
            const float tmp_e = (Gamma22_sqr - Gamma12_sqr) / a1 + (Gamma21_sqr - Gamma11_sqr) / a2;
                                
            const float tmp_f = tmp_a * tmp_b;

            const float tmp_g = tmp_c + tmp_c * 4 * tmp_a * tmp_e + 4 * tmp_d * (tmp_b * tmp_b - 4 * det_Gamma_sqr / (a1 * a2));                           

            // --- Note: tmp_g should be always non-negative. Unfortunately, this is not guaranteed with floating point arithmetic especially 
			//     when tmp_p is close to zero. To prevent getting nan from the sqrt we use the absolute value of tmp_g.
            const float tmp_g_abs = fabsf(tmp_g);
            
            const float tmp_h = 0.5f * sqrtf(tmp_g_abs);
            const float r_f1  = tmp_f + tmp_h;
            const float r_f2  = tmp_f - tmp_h;                     
            
            const float c1    =  0.5f * (1 / rho1_p + 1 / rho2_p) + r_f1;
            const float c2    =  0.5f * (1 / rho1_p + 1 / rho2_p) + r_f2;        

            float rho1r = 1.0 / c1;
            float rho2r = 1.0 / c2;
            
            const float Q12   = -2.0f * tmp_a * (Gamma22 * Gamma12 / a1 +Gamma11 * Gamma21 / a2);
            const float Q22   = 1.0f / rho2_p + 2.0f * tmp_a * (Gamma12_sqr / a1 + Gamma11_sqr / a2);
            
            const vec3 sigma1_v = x1p - 2 * dot(normal_v, x1p) * normal_v;
            const vec3 sigma2_v = x2p - 2 * dot(normal_v, x2p) * normal_v;        

            vec3 X1r_v;
            vec3 X2r_v;
            if (norm((Q22 - 1 / rho1r) * sigma1_v - Q12 * sigma2_v) < EPS_R) compute_ortonormal_reference_frame(X1r_v,X2r_v,sr_v);
            else {    
                X1r_v  = normalize((Q22 - 1 / rho1r) * sigma1_v - Q12 * sigma2_v);
                X2r_v  = normalize(cross(sr_v, X1r_v));
            }
            
            vec3 en, epi, epr;

			// --- Computes the reflected field
            if (fabs(dot(si_v, normal_v)) < 0.99998) { // --- Oblique incidence
                  en  = normalize(cross(normal_v, si_v));
                  epi = normalize(cross(si_v, en));
                  epr = normalize(cross(sr_v, en));
            } else {
                  compute_ortonormal_reference_frame(en, epi, normal_v); 
                  epr = -epi;
            }
            
            const cfloat E0in = Ax1p * dot(x1p, en)  + Ax2p * dot(x2p, en);		// --- Normal component of the incident field
            const cfloat E0ip = Ax1p * dot(x1p, epi) + Ax2p * dot(x2p, epi);	// --- Parallel component of the incient field

            
            // --- Use PEC
            const cfloat E0rn = -E0in;
            const cfloat E0rp =  E0ip;
            cvec3 E0r =  E0rn * en + E0rp * epr;
            
            cfloat AX1r = edot(E0r, X1r_v);
            cfloat AX2r = edot(E0r, X2r_v); 
            
            rho1_m = rho1r;
            rho2_m = rho2r;
            
            x1m    = X1r_v;
            x2m    = X2r_v;
            
            Ax1m   = AX1r;
            Ax2m   = AX2r;
        }
    }
    
    // --- Calculates the scattering parameter
    {
        cvec3 Ei = Ax1p * x1p + Ax2p * x2p;
        vec3 sr_v = normalize(p.intersections[num_points - 2].pos - p.intersections[num_points - 1].pos);

        H = compute_scattering_parameter(Ei, sr_v, rx);
    }
    return H;
}

/*******************/
/* SIGNUM FUNCTION */
static inline float sgn(float a) { return (a >= 0.0f) ? (1.0f) : (-1.0f); }
/*******************/

/*******************************/
/* COMPUTE TRANSITION FUNCTION */
/*******************************/
cfloat transition_function(float x_) {
    
    cfloat F;
    cfloat J = make_cfloat(0, 1);
    
    float x;
    // --- Table computed with Wolfram alpha
    //http://www.wolframalpha.com/input/?i=N%5BTable%5B%282+sqrt%28-1%29+sqrt%28x%29+exp%28sqrt%28-1%29+x%29+sqrt%28pi%2F2%29+%28%281%2F2-FresnelC%28sqrt%282%2Fpi%29+*sqrt%28x%29%29%29-sqrt%28-1%29+%281%2F2-fresnelS%28sqrt%282%2Fpi%29+*sqrt%28x%29%29%29%29%29+%2C+{x%2C+0.1%2C+6.1%2C+0.2}%5D%5D
    
    cfloat Fxi [] = { 
        make_cfloat(0.368104, 0.234453), 
        make_cfloat(0.571713, 0.272992), 
        make_cfloat(0.676763, 0.268233), 
        make_cfloat(0.74395,  0.254857),
        make_cfloat(0.79096,  0.239719), 
        make_cfloat(0.825639, 0.224869), 
        make_cfloat(0.852159, 0.210971), 
        make_cfloat(0.872989, 0.198208), 
        make_cfloat(0.889692, 0.186577), 
        make_cfloat(0.903313, 0.176005), 
        make_cfloat(0.914577, 0.166396), 
        make_cfloat(0.924004, 0.157651), 
        make_cfloat(0.931974, 0.149677),
        make_cfloat(0.938774, 0.142388), 
        make_cfloat(0.944621, 0.13571), 
        make_cfloat(0.949686, 0.129575), 
        make_cfloat(0.9541,   0.123925),
        make_cfloat(0.957971, 0.11871), 
        make_cfloat(0.961383, 0.113883), 
        make_cfloat(0.964405, 0.109407), 
        make_cfloat(0.967094, 0.105245), 
        make_cfloat(0.969497, 0.101369), 
        make_cfloat(0.971651, 0.0977505), 
        make_cfloat(0.973591, 0.0943665), 
        make_cfloat(0.975342, 0.0911958), 
        make_cfloat(0.976929, 0.0882199), 
        make_cfloat(0.978372, 0.085422), 
        make_cfloat(0.979686, 0.0827873), 
        make_cfloat(0.980886, 0.0803025), 
        make_cfloat(0.981986, 0.0779556), 
        make_cfloat(0.982995, 0.0757358)        
    };
    
    float xi[] = {
        0.1, 0.3, 0.5, 0.7, 0.9,
        1.1, 1.3, 1.5, 1.7, 1.9,
        2.1, 2.3, 2.5, 2.7, 2.9,
        3.1, 3.3, 3.5, 3.7, 3.9,
        4.1, 4.3, 4.5, 4.7, 4.9,
        5.1, 5.3, 5.5, 5.7, 5.9,
        6.1
    };
    
    x = fabs(x_);
    
    if(x <= 0.1 ) F = (sqrtf(PI_R * x) - 2 * x * expf(J * PI_R / 4) - 2 * x * x * expf(-J * PI_R / 4) / 3) * expf(J * (x + PI_R / 4));
    else if (x >= 6.1) F = 1.0f + J * (1 / (2 * x)) - (3 / (4 * x * x)) - J * (15 / (8 *x * x * x)) + (75 / (16 * x * x * x * x));
    else {
		size_t i = 0;
        for (; i < sizeof(xi) / sizeof(float); i++) if (x <= xi[i]) break;
        
        i--;
        
        F = Fxi[i] + (Fxi[i + 1] - Fxi[i]) * ((x - xi[i]) / (xi[i + 1] - xi[i]));
    }
    
    if (x_ < 0) F = conjugate(F);
 
    return F;
}

/**************************/
/* COTANGENT - FLOAT CASE */
/**************************/
static float fcot(float x){ return tan(PI_R / 2 - x); }

/*******************/
/* A_FUNC FUNCTION */
/*******************/
static float a_func(float x) { return 2 * (cos(x / 2) * cos(x / 2)); }

/*******************/
/* A_FUNC FUNCTION */
/*******************/
static float ap_func(float n, float x) { return 2 * (cos((2 * PI_R * n - x) / 2) * cos((2 * PI_R * n - x) / 2)); }

/************************************/
/* COMPUTE DIFFRACTION COEFFICIENTS */
/************************************/
static void compute_diffraction_coefficients(float k, float L, float phi, float phip, float sin_beta0, float n, cfloat & Ds, cfloat & Dh) {
    
	const cfloat J = make_cfloat(0, 1);
    
    cfloat F = expf(-J * PI_R / 4) / (-2 * n * sqrt(2 * PI_R * k) * sin_beta0);
    
    cfloat D1 = F * fcot((PI_R + (phi - phip)) / (2 * n)) * transition_function(k * L * a_func(phi - phip));
    cfloat D2 = F * fcot((PI_R - (phi - phip)) / (2 * n)) * transition_function(k * L * a_func(phi - phip));
    cfloat D3 = F * fcot((PI_R + (phi + phip)) / (2 * n)) * transition_function(k * L * ap_func(n, phi + phip));
    cfloat D4 = F * fcot((PI_R - (phi + phip)) / (2 * n)) * transition_function(k * L * a_func(phi + phip));
    
    Ds = D1 + D2 - (D3 + D4);
    Dh = D1 + D2 + (D3 + D4);    
}

/****************************************/
/* FIELD TRANSPORT FOR DIFFRACTIVE RAYS */
/****************************************/
cfloat transport_ray_field_direct_diffraction(const path &p, TX &tx, RX &rx, const nvmesh &mesh, const std::vector<edge> &edge_list, double lambda) {
    
    cfloat J	= make_cfloat(0., 1.);
    float k		= 2 * PI_R / lambda;
    
    size_t num_points = p.intersections.size();    
    assert(num_points == 3);
    assert(p.intersections[1].type == DIFFRACTION);
    
    edge ed = edge_list[p.intersections[1].ID];

	const vec3 e_v  = normalize(ed.b - ed.a);
    const vec3 no_v = normalize(ed.no);
    const vec3 nn_v = normalize(ed.nn);
    
    const vec3 T  =  p.intersections[0].pos;
    const vec3 Qe =  p.intersections[1].pos;
    const vec3 R  =  p.intersections[2].pos;
    
    const vec3 sp_v = normalize(Qe - T);
    const float sp  = dist(Qe, T);
    
    const vec3 s_v  = normalize(R - Qe);
    const float s   = dist(R, Qe);
    
    const float sin_beta0 = norm(cross(s_v, e_v));
     
    const vec3 phi_p_v = -normalize(cross(e_v, sp_v));
    const vec3 phi_v   =  normalize(cross(e_v, s_v ));
    
    const vec3 beta0_p_v = cross (phi_p_v, sp_v);
    const vec3 beta0_v   = cross (phi_v, s_v);
    
    const vec3 to_v = cross(no_v, e_v);
    
    const vec3 spt_v = normalize(sp_v -dot(sp_v, e_v) * e_v);
    const vec3 st_v  = normalize(s_v  -dot(s_v,  e_v) * e_v);
    
    const float phip = PI_R - (PI_R - acos(-dot(spt_v, to_v)) ) * sgn(-dot(spt_v, no_v));
    const float phi  = PI_R - (PI_R - acos( dot(st_v,  to_v)) ) * sgn(-dot(st_v,  no_v));
    
    cvec3 Ei;
    
    {	// --- Generate the first ray
        
        cfloat Ax1, Ax2;
       
        vec3 x1, x2;
        gen_first_ray(tx, -sp_v, Ax1, Ax2, x1, x2, lambda);
        
        Ei = (Ax1 * x1 + Ax2 * x2) * (1 / sp) * expf(-J * 2 * PI_R / lambda * sp);
    }
    
    cfloat Ei_beta0p_v = edot(Ei, beta0_p_v);
    cfloat Ei_phi_p_v  = edot(Ei, phi_p_v);
    
    float L = (s * sp) / (s + sp) * (sin_beta0 * sin_beta0);

    cfloat Ds, Dh;
    
    compute_diffraction_coefficients(k, L, phi, phip, sin_beta0, ed.n, Ds, Dh);
    
    const cfloat SF = sqrtf(sp / (s * (sp + s))) * expf(-J * k * s);

    const cfloat Ed_beta0_v = Ei_beta0p_v * (-Ds * SF);
    const cfloat Ed_phi_v   = Ei_phi_p_v  * (-Dh * SF);
    
    const cvec3 Ed = Ed_beta0_v * beta0_v + Ed_phi_v * phi_v;
    
    return compute_scattering_parameter(Ed, -s_v, rx);
}

/*******************/
/* FIELD TRANSPORT */
/*******************/
extern "C" cfloat transport_ray_field(const path &p, TX &tx, RX &rx, const nvmesh &mesh, const std::vector<edge> &edge_list, double lambda) {
    
    // --- Determining the ray type
    if (p.get_num_diff() > 0) 
        if (p.get_num_diff() == 1 && p.intersections.size() == 3) return transport_ray_field_direct_diffraction(p, tx, rx, mesh, edge_list, lambda);
    else return transport_ray_field_generic_reflection(p, tx, rx, mesh, lambda);
}

//extern "C" cfloat transport_ray_field(std::vector<vec3> positions, TX &tx, RX &rx, double lambda){
//
//    // --- Positions is a vector containing the positions of the emission/reflection/diffraction/detection points of the generic path
//	size_t num_points = positions.size();
//                
//    cfloat J = make_cfloat(0, 1);				// --- Immaginary unit
//    
//    cfloat H = make_cfloat(0, 0);				// --- Output initialization
//    
//    cfloat Ax1m, Ax1p;
//    cfloat Ax2m, Ax2p;
//       
//    // --- x1m and x2m are the unit vectors i_theta and i_phi in the reference system defined by the tx antenna
//	vec3 x1m, x1p;
//    vec3 x2m, x2p;
//    
//    double rho_m, rho_p;
// 
//    {	// --- Generate the first ray
//        double s0 = dist(positions[1], positions[0]);					// --- Distance between the first two path points
//        vec3 s0_v = normalize(positions[1] - positions[0]);				// --- Unit vector pointing from the first to the second path points
//        gen_first_ray(tx, s0_v, Ax1m, Ax2m, x1m, x2m, rho_m, lambda);
//	    //std::cout << s0 << " " << s0_v.x << " " << s0_v.y << " " << s0_v.z << "\n";
//	    //std::cout << lambda << " " << Ax1m.x << " " << Ax1m.y << " " << Ax2m.x << " " << Ax2m.y << " " << "\n";
//    }
//    
//	
//	// --- Transport the field
//    for (size_t i = 0; i < num_points-1; i++) {
//        
//		double s  = dist(positions[i+1], positions[i]);					// --- Distance between two consecutive path points
//        vec3 si_v = normalize(positions[i + 1] - positions[i]);			// --- Unit vector pointing from one path point to the next 
//																		//     (incident propagation unit vector)
//        
//        // --- Calculates the propagation factor
//        cfloat F;        
//        if (rho_m == 0.0) F = (1 / s) * expf(-J * 2 * PI_R / lambda * s);
//		else			  F = (rho_m / (rho_m + s)) * expf(-J * 2 * PI_R / 0.125 * s);
//        
//        // --- Propagates the field
//        Ax1p = Ax1m * F;
//        Ax2p = Ax2m * F;
//        
//		x1p = x1m;
//        x2p = x2m;
//        rho_p = rho_m + s;												// --- Updates the path length
//        
//        // --- Performs reflections
//        if (i < num_points-2){
//            
//            vec3 sr_v = normalize(positions[i + 2] - positions[i + 1]);	// --- Reflection propagation unit vector
//                
//            vec3 normal_v = make_vec3(0, 0, 1);
//
//            vec3 en, epi, epr;
//
//            // --- Computes the reflected field
//            if (fabs(dot(si_v, normal_v)) < 0.99998) { // --- Oblique incidence
//                  en  = normalize(cross(normal_v, si_v));
//                  epi = normalize(cross(si_v, en));
//                  epr = normalize(cross(sr_v, en));
//            } else {
//                  compute_ortonormal_reference_frame(en, epi, normal_v); 
//                  epr = -epi;
//            }
//            
//            const cfloat E0in = Ax1p * dot(x1p, en ) + Ax2p * dot(x2p, en );		// --- Normal component of the incident field
//            const cfloat E0ip = Ax1p * dot(x1p, epi) + Ax2p * dot(x2p, epi);		// --- Parallel component of the incient field
//
//            // --- Use PEC
//            const cfloat E0rn = -E0in;												// --- 
//            const cfloat E0rp =  E0ip;
//
//            rho_m = rho_p;
//            x1m   = en;
//            x2m   = epr;
//            Ax1m  = E0rn;
//            Ax2m  = E0rp;
//        }
//    }
//    
//    // --- Calculates the scattering parameter
//    {
//        cvec3 Ei = Ax1p * x1p + Ax2p * x2p;
//
//        vec3 sr_v = normalize(positions[num_points-2] - positions[num_points-1]);
//        
//        vec3 local_sr_v;
//        local_sr_v.x = dot(sr_v, rx.x_v);
//        local_sr_v.y = dot(sr_v, rx.y_v);
//        local_sr_v.z = dot(sr_v, rx.z_v);
//
//        double azimuth = atan2(local_sr_v.y, local_sr_v.x);
//        double zenith  = acos (local_sr_v.z);
//        
//        cfloat h_theta_r,h_phi_r;
//        vec3   theta_v, phi_v;
//
//        interpolate_pattern(azimuth, zenith, rx.patt, h_theta_r, h_phi_r, theta_v, phi_v);          
//        
//        theta_v = theta_v.x * rx.x_v + theta_v.y * rx.y_v + theta_v.z * rx.z_v;
//        phi_v   = phi_v.x   * rx.x_v + phi_v.y   * rx.y_v + phi_v.z   * rx.z_v;
//
//        cvec3 hr = h_theta_r * theta_v + h_phi_r * phi_v;
//        
//        cfloat Zin = rx.patt -> Zin;
//        cfloat Z0  = make_cfloat(50.0, 0);
//        
//        H = edot(Ei, hr) * (Z0 / (Z0 + Zin));
//    }
//    
//    return H;
//}
