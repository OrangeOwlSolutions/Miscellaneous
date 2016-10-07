#include <iomanip>
#include <time.h>

#include <cuda_runtime.h>

//#include <Eigen/Dense>
//#include <boost/math/special_functions.hpp>

#include <thrust\host_vector.h>
#include <thrust\device_vector.h>
#include <thrust\extrema.h>
#include <thrust\random\uniform_real_distribution.h>
#include <thrust\sort.h>
#include <thrust/tuple.h>
#include <thrust\random\linear_congruential_engine.h>

#include "projectOnCurve.h"
#include "setAntennas.h"

#include "optimizerAncillary.cuh"
#include "Matlab_like.cuh"
#include "Polynomials.cuh"
#include "Utilities.cuh"

thrust::minstd_rand rng(time(NULL));

// --- Constructor
group::group(uint NTX_, uint NTXCoeff_, uint NRX_ , uint NRXCoeff_ , uint NL_, real Rdet_, segment TXLine_, segment RXLine_):
	NTX(NTX_), NRX(NRX_), NTXCoeff(NTXCoeff_), NRXCoeff(NRXCoeff_), TXLine(TXLine_), RXLine(RXLine_), Nconf(NL_), Rdet(Rdet_) {
        
	// --- Best group fitness initialization
    best_fitness = -INF_R;
    best_fitness_index = 0;
        
    TXcoeff.resize(NTXCoeff);
    RXcoeff.resize(Nconf * NRXCoeff);
    h_rx.resize(NRX);
    h_tx.resize(NTX);
    pattern_rx.resize(NRX);
    pattern_tx.resize(NTX);          
}
    
// --- Default constructor
group::group(){};
    
// --- Returns the pointer to the i-th receiving configuration of the group
real * group::get_RXconfiguration(uint i) { assert(i < Nconf); return &(RXcoeff[i * NRXCoeff]); }

const real * group::get_RXconfiguration(uint i) const { assert(i < Nconf); return &(RXcoeff[i*NRXCoeff]); }

real * group::get_TXconfiguration() { return &(TXcoeff[0]); }

const real * group::get_TXconfiguration() const { return &(TXcoeff[0]); }
    
// --- Returns the pointer to the tx and rx orientations of the element
vec3 * group::get_azr_rx() { return &(h_rx[0]); }

const vec3 * group::get_azr_rx() const { return &(h_rx[0]); }
    
vec3 * group::get_azr_tx() { return &(h_tx[0]); }

const vec3 * group::get_azr_tx() const { return &(h_tx[0]); }
    
// --- Returns the pointer to the filenames of the files containing the tx and rx patterns
std::string * group::get_pattern_rx() { return &(pattern_rx[0]); }

const std::string * group::get_pattern_rx() const { return &(pattern_rx[0]); }
     
std::string * group::get_pattern_tx() { return &(pattern_tx[0]); }

const std::string * group::get_pattern_tx() const { return &(pattern_tx[0]); }      

/*******************************************/
/* RANDOM POSITION COEFFICIENTS GENERATION */
/*******************************************/
void gen_random_coeff(float *conf, uint NX, uint NXCoeff,std::vector<double> &LegendreCoeff, double &step, std::vector<double> &xsi, 
	                  float subrange, float offset) {
    
	assert(subrange <= 2);
    assert(subrange > 0);        
    assert(offset <=  1 - subrange / 2);
    assert(offset >= -1 + subrange / 2);

    thrust::uniform_real_distribution<double> dist(-subrange / 2, subrange / 2);

	// --- Element positions
    std::vector<double> spos(NX);
    // --- Positions initialization. By default, the positions are randomly initialized in (-1, 1) around the offset
    for(uint i = 0; i < NX; i++) spos[i] = dist(rng) + offset;
	// --- Positions sorting
    thrust::sort(spos.begin(), spos.end());

	// --- Expansion coefficients of the positions
    std::vector<double> coeff(NXCoeff); 
    
    // --- From element locations to Legendre coefficients
    project_on_curve(spos, coeff, LegendreCoeff, step, xsi);
    
    // --- Copy the result to output
    for(uint i = 0; i < NXCoeff; i++) conf[i] = coeff[i];

}

/************************/
/* GROUP INITIALIZATION */
/************************/
extern "C" void init_groups(std::vector<group> &pop, uint NTX, uint NTXCoeff, uint NRX, uint NRXCoeff, uint Nconf, float Rdet, const segment &TXline, 
	             const segment &RXline, std::vector<vec3> &azr_tx, std::vector<vec3> &azr_rx, std::vector<std::string> &pattern_tx, 
				 std::vector<std::string> &pattern_rx, double &step_tx, double &step_rx, std::vector<double> &xsi_tx, std::vector<double> &xsi_rx,               
                 std::vector<double> &LegendreCoeffTX, std::vector<double> &LegendreCoeffRX) {
                 
    // --- Creation of the groups. Each group is a member of the population.
	for (uint i = 0; i < pop.size(); i++) {
		pop[i] = group(NTX, NTXCoeff, NRX, NRXCoeff, Nconf, Rdet, TXline, RXline);
        //printf("ORIG %f %f %f\n", RXline.orig.x, RXline.orig.y, RXline.orig.z);
        //printf("ORIG %f %f %f\n", pop[i].RXLine.orig.x, pop[i].RXLine.orig.y, pop[i].RXLine.orig.z);
	}

	// --- Sets the tx and rx antenna orientations of each group
    for(uint i = 0; i < pop.size(); i++) {
        group &g = pop[i];
        set_azr(azr_tx, g.get_azr_tx(), NTX);
        set_azr(azr_rx, g.get_azr_rx(), NRX);
    }
    
    // --- Sets the tx and rx antenna patterns of each group
    for(uint i = 0; i < pop.size(); i++) {
		group &g = pop[i];
        set_pattern(pattern_tx, g.get_pattern_tx(), NTX);
        set_pattern(pattern_rx, g.get_pattern_rx(), NRX);
    }
   
    // --- Random generation of the tx and rx antenna position coefficients
	for(uint i = 0; i < pop.size(); i++) {
		group &g = pop[i];
        // --- Random generation of the tx antenna position coefficients
        gen_random_coeff(g.get_TXconfiguration(), g.NTX, g.NTXCoeff, LegendreCoeffTX, step_tx, xsi_tx);
        // --- Random generation of the rx antenna position coefficients
		for (uint j = 0; j < g.Nconf; j++) {
			gen_random_coeff(g.get_RXconfiguration(j), g.NRX, g.NRXCoeff, LegendreCoeffRX, step_rx, xsi_rx);       
		}
	}
   
}

/***********************************/
/* CHILD POPULATION INITIALIZATION */
/***********************************/
void init_child(std::vector<group> &pop, uint NTX, uint NTXCoeff, uint NRX, uint NRXCoeff, uint Nconf, float radius, const segment &TXline, 
			    const segment &RXline, std::vector<vec3> &azr_tx, std::vector<vec3> &azr_rx, std::vector<std::string> &pattern_tx,
                std::vector<std::string> &pattern_rx, double &step_tx, double &step_rx, std::vector<double> &xsi_tx, std::vector<double> &xsi_rx,               
                std::vector<double> &LegendreCoeffTX, std::vector<double> &LegendreCoeffRX) {
                 
    // --- Random creation of groups
    for (uint i = 0; i < pop.size(); i++) pop[i] = group(NTX, NTXCoeff, NRX, NRXCoeff, Nconf, radius, TXline, RXline);
    
    // --- Set the TX and RX antenna orientations
    for (uint i = 0; i < pop.size(); i++) {
        group & g = pop[i];
        set_azr(azr_tx, g.get_azr_tx(), NTX);
        set_azr(azr_rx, g.get_azr_rx(), NRX);
    }
    
    // --- Set the TX and RX antenna patterns
    for(uint i = 0; i < pop.size(); i++) {
        group & g = pop[i];
        set_pattern(pattern_tx, g.get_pattern_tx(), NTX);
        set_pattern(pattern_rx, g.get_pattern_rx(), NRX);
    }
}

/******************************************/
/* FIND THE BEST ELEMENT AMONG ALL GROUPS */
/******************************************/
unsigned int find_best_among_all_groups(std::vector<group> pop, double &maximum) {

    unsigned int i			= 0;
    unsigned int ind_max	= 0;
    maximum					= pop[0].best_fitness;
    
	for (i = 1; i < pop.size(); i++) {
        
		if(pop[i].best_fitness > maximum) {
            maximum = pop[i].best_fitness;
            ind_max = i;
        }
    }
    
	return ind_max;  
}

/*****************************************************/
/* FIND THE BEST ELEMENT IN THE GROUP - HOST VERSION */
/*****************************************************/
unsigned int find_best_configuration_in_the_group_h(thrust::host_vector<double>::iterator begin, thrust::host_vector<double>::iterator end) {
    
	return thrust::max_element(begin, end) - begin;

}

/*******************************************************/
/* FIND THE BEST ELEMENT IN THE GROUP - DEVICE VERSION */
/*******************************************************/
unsigned int find_best_configuration_in_the_group_d(thrust::device_vector<double>::iterator begin, thrust::device_vector<double>::iterator end) {
    
	return thrust::max_element(begin, end) - begin;

}

/******************/
/* CLAMP FUNCTION */
/******************/
double clamp(double x, double low, double high) { return  ( x > high ) ? high : ((x < low) ? low : x); } 

/*********************************************************/
/* BUILD POSITIONS FUNCTION WITHOUT POSITION CONSTRAINTS */
/*********************************************************/
extern "C" void build_positions(const std::vector<group> &pop, thrust::host_vector<vec3> &pos_tx, thrust::host_vector<vec3> &pos_rx, 
								const std::vector<double> &LegendreCoeffTX,	const std::vector<double> &LegendreCoeffRX) {
 
	// --- Calculate transmitter and receiver positions starting from the Legendre coefficients.
	//     LegendreCoeffTX and LegendreCoeffRX are row-major.
	double sum_tx = 0.;
    double sum_rx = 0.;
        
	for (uint g = 0; g < pop.size(); g++) {
        // --- Get the coefficients of the TX configuration
		const real *TX_conf = pop[g].get_TXconfiguration();
        for (uint i = 0; i < pop[g].NTX; i++) {
			sum_tx = 0.;
			for (uint j = 0; j < pop[g].NTXCoeff; j++) sum_tx += TX_conf[j] * LegendreCoeffTX[j + i * pop[g].NTXCoeff];
		    pos_tx[i + g * pop[g].NTX] = clamp(sum_tx, -1, 1) * pop[g].TXLine.dir_v * pop[g].TXLine.L + pop[g].TXLine.orig;
			//h_tx[i + g*pop[g].NTX] = pop[g].h_tx;
		}
		        
        for (uint l = 0; l < pop[g].Nconf; l++) {
			const real *RX_conf = pop[g].get_RXconfiguration(l);
        	for (uint i = 0; i < pop[g].NRX; i++) {
				sum_rx = 0;
			    for (uint j = 0; j < pop[g].NRXCoeff; j++) sum_rx+= RX_conf[j] * LegendreCoeffRX[j + i * pop[g].NRXCoeff];
			    pos_rx[i + l * pop[g].NRX + g * pop[g].Nconf * pop[g].NRX] = clamp(sum_rx, -1, 1) * pop[g].RXLine.dir_v * pop[g].RXLine.L + 
					                        pop[g].RXLine.orig;
			    //h_rx[i + l*pop[g].NRX + g*pop[g].Nconf*pop[g].NRX] = pop[g].h_rx;
		    }
		}
	}
}

/***********************************/
/* DERIVATIVE LEGENDRE POLYNOMIALS */
/***********************************/
double diffLegendrePoly(unsigned int n , double x) {

    // --- Set initial points 
    double diff_0	= 0.0;
    double diff_1	= 1.0;
    double val_im1	= 0.;
    double val		= 0.;
    
    if (n == 0) return diff_0;
    else if(n == 1) return diff_1;
    else {
        for(unsigned int k = 1; k <= n; k++) {
            //val = k * boost::math::legendre_p(k - 1, x) + x * val_im1; 
            val = k * LegendreN(k - 1, x) + x * val_im1; 
            val_im1 = val;
        } 
        return val;
    }
}

/***********************************************/
/* FIND MINIMUM AND MAXIMUM OF A REAL FUNCTION */
/***********************************************/
inline void find_min_max(thrust::device_vector<double> &d_arrval, const int num_points, double *minimum, double *maximum) {
    
	thrust::pair<thrust::device_vector<double>::iterator,thrust::device_vector<double>::iterator> tuple;
    tuple = thrust::minmax_element(d_arrval.begin(), d_arrval.end());
    *minimum = *(tuple.first);
    *maximum = *tuple.second;

}

/*************************************************/
/* FIND MINIMUM AND MAXIMUM INTERELEMENT SPACING */
/*************************************************/
thrust::tuple<double, double> find_min_max_legendre(const float *conf, const unsigned int num_points, const unsigned int NXCoeff) {

	// --- Construct an array of num_points points in the interval [-1, 1]
    //Eigen::VectorXd xsi(num_points);
    //xsi.setLinSpaced(num_points, -0.99999999999, 0.99999999999);    
    double *xsi = h_linspace(-0.99999999999, 0.99999999999, num_points);
    
	// --- Array of the sampled position function derivative
	thrust::host_vector<double> arrval(num_points);

    double temp = 0.0f;    
    
    // --- Construct the position function derivative
	for (unsigned int i = 0 ; i < num_points; i++) {
		temp = 0.0;
        //for (unsigned int j = 0; j < NXCoeff; j++) temp += conf[j] * diffLegendrePoly(j + 1 , xsi(i));
        for (unsigned int j = 0; j < NXCoeff; j++) temp += conf[j] * diffLegendrePoly(j + 1 , xsi[i]);
        arrval[i] = temp;   
    }

	thrust::device_vector<double> d_arrval(num_points);
	d_arrval = arrval;

	double minimum, maximum;
	find_min_max(d_arrval, num_points, &minimum, &maximum);

	// --- Find minimum and maximum interelement spacings
	return thrust::make_tuple(minimum, maximum);   
}

/******************************************/
/* FIND MINIMUM AND MAXIMUM ARRAY EXTENTS */
/******************************************/
thrust::tuple<double, double> find_min_max_func(const real *conf, const unsigned int num_points, const double a, const double b, const unsigned int NXCoeff) {

	// --- Construct an array of num_points points in the interval [-1, 1]
    //Eigen::VectorXd xsi(num_points);
    //xsi.setLinSpaced(num_points, -0.99999999999, 0.99999999999); 
    double *xsi = h_linspace(-0.99999999999, 0.99999999999, num_points);
    
    // --- Array of the sampled position function 
    //double *arrval = (double*) malloc(num_points * sizeof(double));
	thrust::host_vector<double> arrval(num_points);

	//double *d_arrval; gpuErrchk(cudaMalloc(&d_arrval, num_points * sizeof(double)));
   
    double temp = 0.0;    
    
    // --- Construct the array positions
    for(unsigned int i = 0 ; i < num_points; i++){
		temp = 0.0;
        for(unsigned int j = 0; j < NXCoeff; j++){
            ////temp += fabs(a*( conf[j]*boost::math::legendre_p(j + 1,xsi(i))) + b*xsi(i));
            ////temp += (a*(conf[j]*boost::math::legendre_p(j + 1,xsi(i))));
            //temp += (a*(conf[j]*boost::math::legendre_p(j + 1, xsi[i])));
            temp += (a * (conf[j] * LegendreN(j + 1, xsi[i])));
        } 
        arrval[i] = fabs(temp + b * xsi[i]);   
    }        
    
    //gpuErrchk(cudaMemcpy(d_arrval, arrval, num_points * sizeof(double), cudaMemcpyHostToDevice));
	thrust::device_vector<double> d_arrval(num_points);
	d_arrval = arrval;

	double minimum, maximum;
	find_min_max(d_arrval, num_points, &minimum, &maximum);

	// --- Find minimum and maximum interelement spacings
	return thrust::make_tuple(minimum, maximum);   

	//thrust::device_ptr<double> arrval_ptr(d_arrval);                
    
    //thrust::device_ptr<double> maximum = thrust::max_element(arrval_ptr, arrval_ptr + num_points); 
    
    //return *maximum;
}

/********************************************************************************/
/* CHECK IF MINIMUM AND MAXIMUM INTERELEMENT SPACING ARE WITHIN THE CONSTRAINTS */
/********************************************************************************/
void debug_der_func(const float *conf, const unsigned int num_points, const double a, const double b, const double max_der, const double min_der, 
					const unsigned int NXCoeff) {
                
	// --- Construct an array of num_points points in the interval [-1, 1]
	//Eigen::VectorXd xsi(num_points);
 //   xsi.setLinSpaced(num_points,-0.99999999999, 0.99999999999);
    double *xsi = h_linspace(-0.99999999999, 0.99999999999, num_points);
    
    // --- Array of the sampled position function derivative
    double *arrval = (double*) malloc(num_points * sizeof(double));
    double *d_arrval; gpuErrchk(cudaMalloc(&d_arrval, num_points * sizeof(double)));
   
    double temp = 0.0;
    
    for (unsigned int i = 0 ; i < num_points; i++) {
		temp = 0.0;
        //for (unsigned int j = 0; j < NXCoeff; j++) temp += (a * (conf[j] * diffLegendrePoly(j + 1 , xsi(i))));
        for (unsigned int j = 0; j < NXCoeff; j++) temp += (a * (conf[j] * diffLegendrePoly(j + 1 , xsi[i])));
        arrval[i] = temp + b;   
    }     
    
    for (unsigned int i = 0 ; i < num_points; i++) {
    
        if ((arrval[i] > max_der + max_der/10)) std::cerr << "Out of range max value = " << std::setprecision(10) << arrval[i] << std::endl;
        
        if ((arrval[i] < min_der - min_der/10)) std::cerr << "Out of range min value = " << std::setprecision(10) << arrval[i] << std::endl;

	}
    
	gpuErrchk(cudaFree(d_arrval));
}

/***************************************************/
/* CHECK IF ARRAY EXTENT IS WITHIN THE CONSTRAINTS */
/***************************************************/
void debug_func(const float *conf, const unsigned int num_points, const double a, const double b, const double max_func, const unsigned int NXCoeff) {
                
	// --- Construct an array of num_points points in the interval [-1, 1]
    //Eigen::VectorXd xsi(num_points);
    //xsi.setLinSpaced(num_points, -0.99999999999, 0.99999999999);
    double *xsi = h_linspace(-0.99999999999, 0.99999999999, num_points);
    
    // --- Array of the sampled position function derivative
    double *arrval = (double*) malloc(num_points * sizeof(double));
    double *d_arrval; gpuErrchk(cudaMalloc(&d_arrval, num_points * sizeof(double)));
   
    double temp = 0.0;
    
    for(unsigned int i = 0 ; i < num_points; i++) {
		temp = 0.0;
        for (unsigned int j = 0; j < NXCoeff; j++) {
            ////temp += (a*( conf[j]*boost::math::legendre_p(j + 1,xsi(i)) ) + b*xsi(i))/max_func;
            ////temp += (a * (conf[j] * boost::math::legendre_p(j + 1, xsi(i))));
            //temp += (a * (conf[j] * boost::math::legendre_p(j + 1, xsi[i])));
            temp += (a * (conf[j] * LegendreN(j + 1, xsi[i])));
        } 
        arrval[i] = (temp  + b * xsi[i]) / max_func;   
    }     
    
    
    for (unsigned int i = 0 ; i < num_points; i++) {
    
        if((arrval[i] >= 1) || (arrval[i] <= -1) ) std::cerr << "Out of range value = " << std::setprecision(10) << arrval[i] << std::endl;
        
        assert(arrval[i] <=  1);
        assert(arrval[i] >= -1);
    }

	gpuErrchk(cudaFree(d_arrval));
}

/******************************************************/
/* BUILD POSITIONS FUNCTION WITH POSITION CONSTRAINTS */
/******************************************************/
extern "C" void build_positions_with_constraints(const std::vector<group> &pop, std::vector<vec3> &pos_tx, std::vector<vec3> &pos_rx,
												 const std::vector<double> &xsi_tx, const std::vector<double> &xsi_rx, const double lambda, 
												 const double step_tx, const double step_rx, const std::vector<double> &LegendreCoeffTX, 
												 const std::vector<double> &LegendreCoeffRX) {
                                      
	// --- Minimum derivative values
	const double d_min_tx = (lambda / 2. + 4. * lambda / 20.) / step_tx; 
    const double d_min_rx = (lambda / 2. + 4. * lambda / 20.) / step_rx;
        
	std::cerr << "d_min_tx = " << d_min_tx << std::endl;
    std::cerr << "d_min_rx = " << d_min_rx << std::endl;        
        
	// --- Maximum derivative values
	const double d_max_tx = 0.9; 
    const double d_max_rx = 0.9;
	
	// --- Number of function sampling points
	const uint num_points = 10000;				
    
	double minimum, maximum, test_minimum;
    double minimum_func, maximum_func;
    double a, b, det;
        
	// --- Calculate transmitter and receiver positions starting from the Legendre coefficients.
	//     LegendreCoeffTX and LegendreCoeffRX are row-major.
	double sum_tx = 0.;
	double sum_rx = 0.;
        
	for (uint g = 0; g < pop.size(); g++) {
        
		// --- Get the coefficients of the TX configuration
        const real *TX_conf = pop[g].get_TXconfiguration();
        
		// --- Find minimum and maximum interlement spacings
		thrust::tie(minimum, maximum) = find_min_max_legendre(TX_conf, num_points, pop[g].NTXCoeff);
        // --- Check if minimum and maximum interlement spacings are within the constraints
		assert(maximum >= minimum);               
		det = (minimum - maximum);
		a = (d_min_tx - d_max_tx) / (det + 1E-16);
        b = (d_max_tx * minimum - d_min_tx * maximum) / (det + 1E-16);
		debug_der_func(TX_conf, num_points, a, b, d_max_tx, d_min_tx, pop[g].NTXCoeff); 
                 
        // --- Find max|f(x)| to normalize the function
		thrust::tie(minimum_func, maximum_func) = find_min_max_func(TX_conf, num_points, a, b, pop[g].NTXCoeff); 
		//debug_func(TX_conf, num_points, a, b, maximum_func, pop[g].NTXCoeff);  
                
		// --- Construct the transmitter positions
		for (uint i = 0; i < pop[g].NTX; i++) {
			sum_tx = 0;
			for(uint j = 0; j < pop[g].NTXCoeff; j++) sum_tx+= TX_conf[j] * LegendreCoeffTX[j + i * pop[g].NTXCoeff];
            sum_tx = (a * sum_tx + b * xsi_tx[i]) / maximum_func;
			pos_tx[i + g*pop[g].NTX] = sum_tx * pop[g].TXLine.dir_v * pop[g].TXLine.L + pop[g].TXLine.orig;
		}
		      
		// --- Construct the receiver positions
        for (uint l = 0; l < pop[g].Nconf; l++) {
			const real *RX_conf = pop[g].get_RXconfiguration(l);
            thrust::tie(minimum, maximum) = find_min_max_legendre(RX_conf, num_points, pop[g].NRXCoeff);
            assert(maximum >= minimum);
            det = (minimum - maximum);
            a = (d_min_rx - d_max_rx) / (det + 1E-16);
            b = (d_max_rx * minimum - d_min_rx * maximum) / (det + 1E-16);  
            debug_der_func(RX_conf, num_points, a, b, d_max_rx, d_min_rx, pop[g].NRXCoeff);
            thrust::tie(minimum_func, maximum_func) = find_min_max_func(RX_conf, num_points, a, b, pop[g].NRXCoeff);   
            //debug_func(RX_conf, num_points, a, b, maximum_func, pop[g].NRXCoeff);                 
                    
        	for (uint i = 0; i < pop[g].NRX; i++) {
				sum_rx = 0.;
		        for (uint j = 0; j < pop[g].NRXCoeff; j++) sum_rx+= RX_conf[j] * LegendreCoeffRX[j + i * pop[g].NRXCoeff];
                sum_rx = (a * sum_rx + b * xsi_rx[i]) / maximum_func;
                pos_rx[i + l * pop[g].NRX + g * pop[g].Nconf * pop[g].NRX] = sum_rx * pop[g].RXLine.dir_v * pop[g].RXLine.L + pop[g].RXLine.orig;		        	   
			}
        }
	}
	std::cerr << "Building up positions finished..." << std::endl;
}

/***********************************************************/
/* BUILD BEST POSITIONS FUNCTION WITH POSITION CONSTRAINTS */
/***********************************************************/
void build_best_positions_with_constraints(const group &g, const std::vector<double> &LegendreCoeffTX, const std::vector<double> &LegendreCoeffRX, 
										   thrust::host_vector<vec3> &pos_tx, thrust::host_vector<vec3> &pos_rx, thrust::host_vector<vec3> &h_tx,
                                           thrust::host_vector<vec3> &h_rx, const std::vector<double> &xsi_tx, const std::vector<double> &xsi_rx, 
										   const double lambda, const double step_tx, const double step_rx, unsigned int index_max_element) {
         
	// --- Minimum distance between antennas
	const double d_min_rx = (lambda / 2. + 4. * lambda / 20.) / step_rx; 
    const double d_min_tx = (lambda / 2. + 4. * lambda / 20.) / step_tx;;        

	// --- Maximum distance between antennas
	const double d_max_tx = 0.9; 
    const double d_max_rx = 0.9;
                                                   
	// --- Number of function sampling points
	const uint num_points = 10000;				

	double minimum, maximum, test_minimum;
    double minimum_func, maximum_func;
    double a, b, det;                                          
                          
	// --- Calculate transmitter and receiver positions starting from the Legendre coefficients.
	//     LegendreCoeffTX and LegendreCoeffRX are row-major.
	double sum_tx = 0.;
	double sum_rx = 0.;
       
    // --- Get the coefficients of the best TX configuration    
    const real *TX_conf = g.get_TXconfiguration();
    
	// --- Find minimum and maximum interlement spacings
	thrust::tie(minimum, maximum) = find_min_max_legendre(TX_conf, num_points, g.NTXCoeff);

	det = (minimum - maximum);
    a = (d_min_tx - d_max_tx) / (det + 1E-16);
    b = (d_max_tx * minimum - d_min_tx * maximum) / (det + 1E-16);

	// --- Find max|f(x)| to normalize the function
	thrust::tie(minimum_func, maximum_func) = find_min_max_func(TX_conf, num_points, a, b, g.NTXCoeff);       
        
	// --- Construct the best transmitter positions
	for (uint i = 0; i < g.NTX; i++) {
		sum_tx = 0;
	    for(uint j = 0; j < g.NTXCoeff; j++) sum_tx+= TX_conf[j] * LegendreCoeffTX[j + i * g.NTXCoeff];
	    sum_tx = (a * sum_tx + b * xsi_tx[i]) / maximum_func;
		pos_tx[i] = sum_tx * g.TXLine.dir_v * g.TXLine.L + g.TXLine.orig;
    }
		        
	// --- Construct the receiver positions
    const real *RX_conf = g.get_RXconfiguration(index_max_element);
	thrust::tie(minimum, maximum) = find_min_max_legendre(RX_conf, num_points, g.NRXCoeff);
	det = (minimum - maximum);
	a = (d_min_rx - d_max_rx) / (det + 1E-16);
	b = (d_max_rx*minimum - d_min_rx*maximum) / (det + 1E-16);  
	thrust::tie(minimum_func, maximum_func) = find_min_max_func(RX_conf, num_points, a, b, g.NRXCoeff);         
                    
	for (uint i = 0; i < g.NRX; i++) {
		sum_rx = 0.;
        for (uint j = 0; j < g.NRXCoeff; j++) sum_rx+= RX_conf[j] * LegendreCoeffRX[j + i * g.NRXCoeff];
    	sum_rx = (a*sum_rx + b*xsi_rx[i])/maximum_func;
    	pos_rx[i] = sum_rx*g.RXLine.dir_v*g.RXLine.L + g.RXLine.orig;
	}
}
