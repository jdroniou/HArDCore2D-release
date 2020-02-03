// Class to provide a non-linear function

//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
* This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/


#ifndef _TEST_CASE_NONLINEARITY_HPP
#define _TEST_CASE_NONLINEARITY_HPP

#include <functional>
#include <string>

namespace HArDCore2D {}

using namespace HArDCore2D;

/*!
* @defgroup TestCases
*	@brief Defines test cases (exact solutions, source terms, boundary condition, etc.)
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

// @addtogroup TestCases
//@{

/// The TestCaseNonLinearity class provides definition of a nonlinear function, and related functions
class TestCaseNonLinearity {

public:
  using nonlinearity_function_type = std::function<double(double,std::string)>;    ///< type for nonlinear function

  /// Initialise data
  TestCaseNonLinearity(
    const int iTCNL,  ///< The id of the test case
    const double m=1.0    ///< Optional, power of the porous medium non-linearity
  );

  /// Nonlinearity
  double nonlinearity(
    const double s,
    const std::string type      ///< "type" = "fct" for the function itself, "nlin" for the nonlinear part gamma such that fct(u) = gamma(u) u, "der" for derivative, "hess" for 2nd derivative
  );

private:
  // id of test case
  const int m_iTCNL;
  // PME nonlinearity
  const int m_m;

	// Sign function
  double sign(const double s) const {
  	double val = 0;
  	if (s>0){
  		val = 1;
  	}else if (s<0){
  		val = -1;
  	}
  	return val;	
  };

};


//@}

#endif //_TEST_CASE_NONLINEARITY_HPP
