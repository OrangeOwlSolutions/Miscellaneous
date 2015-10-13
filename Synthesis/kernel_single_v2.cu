#include <stdio.h>
#include <math.h>
#include <iostream>

#include <thrust\device_vector.h>
#include <thrust\transform_reduce.h>
#include <thrust\reduce.h>
#include <thrust\tuple.h>

#include "BBComplex.h"
#include "InputOutput.cuh"
#include "Utilities.cuh"
#include "Matlab_like.cuh"
#include "Polynomials.cuh"
#include "Synthesis.cuh"
#include "NFFT2_2D.cuh"
#include "NDFT2_2D.cuh"
#include "PSO.cuh"

#define pi 3.141592653589793238463

// --- Algorithm parameters
#define freq				((14.25)*(1e9))					// Operating frequency 
//#define lambda				((3e8)/(freq))                  // Wavelength
extern const float	lambda_f	= (float) 3e8 / freq;		// Wavelength				
extern const double	lambda		= (double)3e8 / freq;		// Wavelength				
extern const float  beta_f		= 2. * pi / lambda;			// Wavenumber
extern const double beta		= 2. * pi / lambda;			// Wavenumber

extern const int M_x = 44;									// Number of reflectarray elements along the x-axis
extern const int M_y = 44;									// Number of reflectarray elements along the y-axis

#define dx					((0.5)*(lambda))                // dist elem x (per array)
#define dy					((0.5)*(lambda))                // dist elem y (per array)

#define aap					(((M_x)-(1))*((dx)/(2)))        // Reflectarray semi-dimension along the x-axis
#define bap					(((M_y)-(1)))*((dy)/(2))        // Reflectarray semi-dimension along the y-axis

extern const float  mfact_f = 12.f;                              // Feed pattern: cos^mfact(theta)
extern const double mfact   = 12.;                               // Feed pattern: cos^mfact(theta)

// ??? INUTILI ???
#define dmin				((0.51)*(lambda))               // Minimum allowed inter-element spacing
#define dmax				((0.7 )*(lambda))               // Maximum allowed inter-element spacing

#define dmin_x				((0.51)*(lambda))               // Minimum allowed inter-element spacing along the x-axis
#define dmin_y				((0.51)*(lambda))               // Minimum allowed inter-element spacing along the y-axis
#define dmax_x				((0.7 )*(lambda))               // Maximum allowed inter-element spacing along the x-axis
#define dmax_y				((0.7 )*(lambda))               // Maximum allowed inter-element spacing along the y-axis

#define z0					((2)*(0.8)*(sqrt((aap)*(aap)+(bap)*(bap))))
															// Focal length of the reflectarray surface

extern const float	feed_center_x_f	= 0.f;								
extern const float	feed_center_y_f	= 1.15f * bap;
extern const float	feed_center_z_f  = -z0;
extern const double feed_center_x	 = 0.;								
extern const double feed_center_y	 = 1.15 * bap;
extern const double feed_center_z	 = -z0;

extern const float  alfa_f = -atan(feed_center_y_f / feed_center_z_f);// Feed illumination angle
extern const double alfa   = -atan(feed_center_y   / feed_center_z);  // Feed illumination angle

extern const int Num_unknowns_x = 5;								// Number of unknowns for the element positions along the x-axis
extern const int Num_unknowns_y = 5;								// Number of unknowns for the element positions along the y-axis

extern const int Num_unknowns_phases = 6;							// Number of unknowns for the phase representation

#define chi_u_prime			4                               // Spectral oversampling factor along u
#define chi_v_prime			4                               // Spectral oversampling factor along v

//#define a_prime				((chi_u_prime)*(aap))			
//#define b_prime				((chi_v_prime)*(bap))
extern const float  a_prime_f	= (float)((chi_u_prime)*(aap));
extern const float  b_prime_f	= (float)((chi_v_prime)*(bap));
extern const double a_prime		= ((chi_u_prime)*(aap));
extern const double b_prime		= ((chi_v_prime)*(bap));

#define u_max				((beta)/(2.))					// Maximum value of the spectral region along the u axis
#define u_min				(-(beta)/(2.))					// Minimum value of the spectral region along the u axis
#define v_max				((beta)/(2.))					// Maximum value of the spectral region along the v axis
#define v_min				(-(beta)/(2.))					// Minimum value of the spectral region along the v axis

extern int	Nu;
extern int	Nv;

extern cublasHandle_t cublasHandleSynthesis; 

extern float  *d_U_discrete_f			= NULL;
extern float  *d_V_discrete_f			= NULL;
extern float  *d_Filter_f				= NULL;
extern float  *d_LEG_f					= NULL;
extern float  *d_ZERNIKE_f				= NULL;
extern float  *d_Internal_Coverage_f	= NULL;
extern float  *d_External_Coverage_f	= NULL;
extern double *d_U_discrete				= NULL;
extern double *d_V_discrete				= NULL;
extern double *d_Filter					= NULL;
extern double *d_LEG					= NULL;
extern double *d_ZERNIKE				= NULL;
extern double *d_Internal_Coverage		= NULL;
extern double *d_External_Coverage		= NULL;

#define DEBUG

/********/
/* MAIN */
/********/
int main()
{
	cublasSafeCall(cublasCreate(&cublasHandleSynthesis));
	
	// --- Defining spectral quantities
	thrust::pair<thrust::pair<float *, float *>, float *> d_SpectralTuple = defineSpectralQuantities((float)u_max, (float)v_max, (float)a_prime_f, (float)b_prime_f, (float)beta, &Nu, &Nv);
	thrust::pair<float *, float *> d_UV_discrete = d_SpectralTuple.first;
	d_U_discrete_f = d_UV_discrete.first;
	d_V_discrete_f = d_UV_discrete.second;
	d_Filter_f	   = d_SpectralTuple.second;

	// --- Generating the (csi, eta) grid and the Legendre polynomials
	thrust::pair<thrust::pair<float *, float *>, float *> d_LegendreTuple = generateLegendreFactorized<float>(Num_unknowns_x, Num_unknowns_y, M_x, M_y);
	thrust::pair<float *, float *> d_CSI_ETA = d_LegendreTuple.first;
	float *d_CSI = d_CSI_ETA.first;
	float *d_ETA = d_CSI_ETA.second;
	d_LEG_f = d_LegendreTuple.second;
	
	// --- Generating the Zernike polynomials
	d_ZERNIKE_f = generateZernikep(d_CSI, d_ETA, Num_unknowns_phases, M_x, M_y);
	
	// --- Loading the masks
	d_Internal_Coverage_f = loadGPUrealtxt("C:\\Users\\angelo\\Documents\\CEM\\ParticleSwarm\\ParticleSwarmSynthesis\\ParticleSwarmSynthesisMatlab\\Internal_Coverage.txt", d_Internal_Coverage_f, (2 * Nu) * (2 * Nv));
	d_External_Coverage_f = loadGPUrealtxt("C:\\Users\\angelo\\Documents\\CEM\\ParticleSwarm\\ParticleSwarmSynthesis\\ParticleSwarmSynthesisMatlab\\External_Coverage.txt", d_External_Coverage_f, (2 * Nu) * (2 * Nv));

	/**********************/
	/* PSO INITIALIZATION */
	/**********************/
	h_PSO_Initialize();
	h_PSO_Optimize();

	return 0;
}
