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


#ifndef _TEST_CASE_HPP
#define _TEST_CASE_HPP

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace HArDCore2D {		// Forward declaration
	class Cell;
}

using namespace HArDCore2D;

/*!
* @defgroup TestCases
*	@brief Defines test cases (exact solutions, source terms etc.)
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

// @addtogroup TestCases
//@{

/// The TestCase class provides definition of test cases
class TestCase {

public:
  /// Initialise data
  TestCase(
    const std::vector<int> iTC  ///< The vector id of the test case: (id of solution, id of diffusion)
  );

	/// Returns the exact solution at the points x, y
	double sol(
		const double x,
		const double y
	);

	/// Returns the gradient of the exact solution at the points x, y
	Eigen::Vector2d grad_sol(
		const double x,
		const double y,
		const Cell* cell			///< In case of discontinuity, we need to know the cell we're in to select the correct formula
	);

	/// Returns the Hessian of the exact solution at the points x, y
	Eigen::Matrix2d hess_sol(
		const double x,
		const double y,
		const Cell* cell			///< In case of discontinuity, we need to know the cell we're in to select the correct formula
	);


	/// Returns the diffusion matrix at the points x, y
	Eigen::Matrix2d diff(
		const double x,
		const double y,
		const Cell* cell			///< In case of discontinuity, we need to know the cell we're in to select the correct formula
	);

	/// Returns the divergence by row of the diffusion matrix at the points x, y
	Eigen::Vector2d div_diff(
		const double x,
		const double y,
		const Cell* cell			///< In case of discontinuity, we need to know the cell we're in to select the correct formula
	);

	/// Returns the source term at the points x, y
	double source(
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
  const std::vector<int> iTC;
	const double pi = acos(-1);
	const double eps = 1e-5;		// constant for rotating diffusion case

	size_t _deg_diff;
  // Generic parameter (e.g. for changing anisotropy etc.)
  double _lambda;


};

inline size_t TestCase::get_deg_diff() { return _deg_diff; };
inline double TestCase::get_lambda() { return _lambda; };

//@}

#endif //_TEST_CASE_HPP
