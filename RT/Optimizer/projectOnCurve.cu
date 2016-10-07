#include <Eigen\Dense>
//#include "C:\Users\angelo\Documents\Project\ThirdPartyLibraries\Eigen\Eigen\Dense"
#include <limits>
#include <vector>
//#include <boost/math/special_functions.hpp>
#include "Polynomials.cuh"
//#include "debug.h"
#include <fstream>
#include <iomanip>
#include <math.h>

#include "setAntennas.h"

typedef unsigned int uint;

/************************************/
/* MATRIX PSEUDOINVERSE CALCULATION */
/************************************/
template<typename DerivedA> Eigen::Matrix< typename DerivedA::Scalar, DerivedA::ColsAtCompileTime, DerivedA::RowsAtCompileTime> 
			PseudoInverse(const Eigen::MatrixBase<DerivedA> & a, double epsilon = std::numeric_limits<typename DerivedA::Scalar>::epsilon())
{
    assert(a.rows()>=a.cols());
    
	typedef Eigen::Matrix<typename DerivedA::Scalar, DerivedA::RowsAtCompileTime, DerivedA::ColsAtCompileTime> InputType;
 
    typedef Eigen::Matrix<typename DerivedA::Scalar, DerivedA::ColsAtCompileTime, DerivedA::RowsAtCompileTime> ReturnType;
 
    Eigen::JacobiSVD<InputType> svd = a.jacobiSvd(Eigen::ComputeFullU |Eigen::ComputeFullV);

    double tolerance = epsilon * std::max(a.cols(),a.rows()) * svd.singularValues().array().abs().maxCoeff();

    ReturnType sigma = ReturnType::Zero(a.cols(), a.rows());

    sigma.block(0, 0, a.cols(), a.cols()) = (svd.singularValues().array().abs()>tolerance).
                                             select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal();

    return svd.matrixV() * sigma * svd.matrixU().adjoint();
}

/***************************************************/
/* FROM ELEMENT POSITIONS TO LEGENDRE COEFFICIENTS */
/***************************************************/
// --- Projects the tx-rx antenna positions on Legendre polynomials
void project_on_curve(const std::vector<double> &spos_,  std::vector<double> &coeff_, std::vector<double> &LegendreCoeff_, 
                      double &step_, std::vector<double> &xsi_){
                      
	const uint Npos		= spos_.size();
    const uint Ncoeff	= coeff_.size();

	Eigen::MatrixXd		spos(Npos,1);
    Eigen::MatrixXd		coeff(Ncoeff,1);

    for(uint i = 0; i < Npos; i++) spos(i) = spos_[i];

	Eigen::MatrixXd LegendreCoeff(Npos, Ncoeff);			// --- Legendre matrix coefficients
    
    Eigen::VectorXd		xsi(Npos);
    xsi.setLinSpaced(-1, 1);								// --- Sampling grid of the Legendre coefficients
    
    for(uint j = 0; j < Npos; j++) xsi_[j] = xsi(j);
    
    step_ = fabs(xsi(1) - xsi(0));
  
	// --- Computes the Legendre coefficients matrix in row-major order
	for(uint r = 0; r < Npos; r++ ) {
        for(uint c = 0; c < Ncoeff; c++ ) {
			//printf("LegendreCoeff(r,c) = boost::math::legendre_p(c + 1, xsi(r));\n");
			//printf("LegendreCoeff(r,c) = %f\n", LegendreCoeff(r,c));
			////printf("boost::math::legendre_p(c + 1, xsi(r)) %f %f\n", boost::math::legendre_p(c + 1, xsi(r)), LegendreN(c + 1, xsi(r)));
			//printf("boost::math::legendre_p(c + 1, xsi(r)) %f\n", LegendreN(c + 1, xsi(r)));
   //         //LegendreCoeff(r,c) = boost::math::legendre_p(c + 1, xsi(r)); 
            LegendreCoeff(r,c) = LegendreN(c + 1, xsi(r)); 
			//printf("LegendreCoeff_[r * Ncoeff + c] = LegendreCoeff(r,c); \n");
            LegendreCoeff_[r * Ncoeff + c] = LegendreCoeff(r,c); 
        }
    }
    
	coeff = PseudoInverse(LegendreCoeff) * spos;

    // --- Copy the result to output
    for(uint i = 0; i < Ncoeff; i++) coeff_[i] = coeff(i);
}
