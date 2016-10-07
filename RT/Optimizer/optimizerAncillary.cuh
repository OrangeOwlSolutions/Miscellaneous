#ifndef __OPTIMIZERANCILLARY_H__
#define __OPTIMIZERANCILLARY_H__

#include <iostream>
#include <vector>
#include <assert.h>

#include "setAntennas.h"
#include "vec3.h"

#ifndef uint
#define uint unsigned int
#endif

#define real float

/*********/
/* GROUP */
/*********/
struct group{    
    uint NTX;										// --- Number of transmitting antennas
    uint NTXCoeff;									// --- Number of expansion coefficients for the transmitting antennas

    uint NRX;										// --- Number of receiving antennas
    uint NRXCoeff;									// --- Number of expansion coefficients for the transmitting antennas
    uint Nconf;										// --- Number of elements per group

    std::vector<real> TXcoeff;						// --- Legendre coefficients for the transmitting antenna locations (row-major ordering) [NTXCoeff]
    std::vector<real> RXcoeff;						// --- Legendre coefficients for the receiving antenna locations	(row-major ordering) [NRXCoeff x Nconf]
    
    real best_fitness;								// --- Best fitness within the group
    uint best_fitness_index;						// --- Index of the best fitness configuration within the group
    segment TXLine, RXLine;							// --- Transmitting and receiving antenna segments
    std::vector<vec3> h_tx;							// --- [Azimuth Zenith Roll] for the transmitting antennas
    std::vector<vec3> h_rx;							// --- [Azimuth Zenith Roll] for the receiving antennas
    std::vector<std::string> pattern_rx;			// --- Patterns for the transmitting antennas (filname)
    std::vector<std::string> pattern_tx;			// --- Patterns for the receiving antennas (filename)
    real Rdet;										// --- Radius of the detector
    
    group(uint, uint, uint, uint, uint, real, segment, segment);
    
    group();
    
    // --- Returns the pointer to the i-th receiving configuration of the group
    real *get_RXconfiguration(uint i);

    const real *get_RXconfiguration(uint i) const;

    real *get_TXconfiguration();

    const real *get_TXconfiguration() const;
    
    // --- Returns the pointer to the tx and rx orientations of the element
    vec3 *get_azr_rx();

    const vec3 *get_azr_rx() const;
    
    vec3 *get_azr_tx();

    const vec3 *get_azr_tx() const;
    
    // --- Returns the pointer to the filenames of the files containing the tx and rx patterns
    std::string *get_pattern_rx();

    const std::string *get_pattern_rx() const;
     
    std::string *get_pattern_tx();

    const std::string *get_pattern_tx() const;

};

void gen_random_coeff(float *, uint, uint,std::vector<double> &, double &, std::vector<double> &, float subrange = 2, float offset = 0);

extern "C" void init_groups(std::vector<group> &, uint, uint, uint, uint, uint, float, const segment &, const segment &, std::vector<vec3> &, 
	             std::vector<vec3> &, std::vector<std::string> &, std::vector<std::string> &, double &, double &, std::vector<double> &, 
				 std::vector<double> &, std::vector<double> &, std::vector<double> &);

void init_child(std::vector<group> &, uint, uint, uint, uint, uint, float, const segment &, const segment &, std::vector<vec3> &, std::vector<vec3> &, 
				std::vector<std::string> &, std::vector<std::string> &pattern_rx, double &, double &, std::vector<double> &, std::vector<double> &,               
                std::vector<double> &, std::vector<double> &);

unsigned int find_best_among_all_groups(std::vector<group>, double &);

unsigned int find_best_configuration_in_the_group_d(thrust::device_vector<double>::iterator, thrust::device_vector<double>::iterator);
unsigned int find_best_configuration_in_the_group_h(thrust::host_vector<double>::iterator, thrust::host_vector<double>::iterator);

extern "C" void build_positions(const std::vector<group> &, thrust::host_vector<vec3> &, thrust::host_vector<vec3> &, const std::vector<double> &, 
								const std::vector<double> &);

extern "C" void build_positions_with_constraints(const std::vector<group> &, std::vector<vec3> &, std::vector<vec3> &,
												 const std::vector<double> &, const std::vector<double> &, const double, const double, const double,                 
												 const std::vector<double> &, const std::vector<double> &);
void build_best_positions_with_constraints(const group &, const std::vector<double> &, const std::vector<double> &, thrust::host_vector<vec3> &, 
										   thrust::host_vector<vec3> &, thrust::host_vector<vec3> &, thrust::host_vector<vec3> &, 
										   const std::vector<double> &, const std::vector<double> &, const double, const double, const double, 
										   unsigned int);
#endif
