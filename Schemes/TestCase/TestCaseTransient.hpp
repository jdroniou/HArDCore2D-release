// Class to provide various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*	This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/


#ifndef _TEST_CASE_TRANSIENT_HPP
#define _TEST_CASE_TRANSIENT_HPP

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "basis.hpp"

#include <Eigen/Dense>

namespace HArDCore2D {		// Forward declaration
	class Cell;
}

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

/// The TestCaseTransient class provides definition of test cases
class TestCaseTransient {

public:
  /// Initialise data
  TestCaseTransient(
    const std::vector<int> iTC,  ///< The vector id of the test case: (id of solution, id of diffusion)
    const double mPME  ///< PME power, for the Barenblatt solution
  );

	/// Returns the exact solution at time t and point p
	std::function<double(const double &, const VectorRd &)> solution();

	/// Returns the gradient of the exact solution at time t and point p
  /// In case of discontinuity, we need to know the cell we're in to select the correct formula
  /// Note: the gradient for the Barenblatt solution is NOT computed (left at 0) as it is not useful to compute the source term
	std::function<VectorRd(const double&, const VectorRd&, const Cell*)> grad_solution();

	/// Returns the Hessian of the exact solution at time t and point p
  /// Note: the Hessian for the Barenblatt solution is NOT computed (left at 0) as it is not useful to compute the source term
	std::function<Eigen::Matrix2d(const double&, const VectorRd&, const Cell*)> hess_solution();

	/// Returns the time derivative of the exact solution at time t and point p
	std::function<double(const double &, const VectorRd &)> delt_solution();

	/// Returns the diffusion matrix at the point x, y (not depending on time for the moment)
	Eigen::Matrix2d diff(
		const double x,
		const double y,
		const Cell* cell			///< In case of discontinuity, we need to know the cell we're in to select the correct formula
	);

	/// Returns the divergence by row of the diffusion matrix at the point x, y
	Eigen::Vector2d div_diff(
		const double x,
		const double y,
		const Cell* cell			///< In case of discontinuity, we need to know the cell we're in to select the correct formula
	);

	/// Returns div(diff \nabla) of the exact solution at the point t, x, y
	double div_diff_grad(
    const double t,
		const double x,
		const double y,
		const Cell* cell			///< In case of discontinuity, we need to know the cell we're in to select the correct formula
	);

  /// Returns the value of the parameter lambda
	inline double get_lambda();

	/// Check if the provided test cases are valid (within range, and combination of solution/diffusion valid)
	void validate();

	/// Returns the degree of the diffusion tensor (useful to set up quadrature rules of proper degree)
	inline size_t get_deg_diff();

private:
  // Parameters: id of test case, pi
  const std::vector<int> m_iTC;
	const double pi = acos(-1);
	const double sqrt2 = std::pow(2,0.5);
	const double gamma = 1.0/3.0;
	const double eps = 1e-5;		// constant for rotating diffusion case

	size_t _deg_diff;
  // Generic parameter (e.g. for changing anisotropy etc.)
  double _lambda;
  
  // Coefficients for Barenblatt: initial time, C_B, gamma
  double m_t0 = 0.1;
  double m_CB = 0.005;
  double m_mPME;
  double m_alpha;
  double m_rho;
  double m_gamma;

};

inline size_t TestCaseTransient::get_deg_diff() { return _deg_diff; };
inline double TestCaseTransient::get_lambda() { return _lambda; };

//@}

#endif //_TEST_CASE_TRANSIENT_HPP
