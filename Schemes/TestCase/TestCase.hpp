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

#include "basis.hpp"

namespace HArDCore2D
{ // Forward declaration
class Cell;
}

using namespace HArDCore2D;

/*!
* @defgroup TestCases
*	@brief Defines test cases (exact solutions, source terms, boundary condition, etc.)
*/

// @addtogroup TestCases
//@{

// Types
template <typename T>
using CellFType = std::function<T(const VectorRd &, const Cell *)>; ///< type for function of a point, that could be discontinuous between cells. T is the type of value of the function

/// Structure to store a solution and its derivatives
struct Solution
  {
    // Constructor
    Solution(
      const FType<double> value,       /// Solution
      const CellFType<VectorRd> gradient,  /// Gradient of the solution
      const CellFType<MatrixRd> hessian    /// Hessian of the solution
    )
    : value(value),
      gradient(gradient),
      hessian(hessian)
    {
      // Do nothing
    }
    
    // Members
    const FType<double> value;
    const CellFType<VectorRd> gradient;
    const CellFType<MatrixRd> hessian;
  };

/// Structure to store a diffusion tensor and its row-wise divergence
struct Diffusion
  {
    // Constructor
    Diffusion(
      const CellFType<MatrixRd> value,       /// Diffusion tensor
      const CellFType<VectorRd> divergence,  /// Divergence, taken row by row
      const size_t degree                    /// Polynomial degree of the tensor
    )
    : value(value),
      divergence(divergence),
      degree(degree)
    {
      // Do nothing
    }
    
    // Members
    const CellFType<MatrixRd> value;
    const CellFType<VectorRd> divergence;
    const size_t degree;
  };

/// Structure to store an advection velocity and its divergence, together with various flags.
struct Advection
  {
    // Constructor
    Advection(
      const CellFType<VectorRd> value,      /// Velocity
      const CellFType<double> divergence,   /// Divergence of velocity
      const bool is_zero,          /// True if velocity is zero
      const bool is_divergence_zero,        /// True if velocity has zero divergence
      const bool is_divergence_constant    /// True if divergence is constant
    )
    : value(value),
      divergence(divergence),
      is_zero(is_zero),
      is_divergence_zero(is_divergence_zero),
      is_divergence_constant(is_divergence_constant)
    {
      // Do nothing
    }
    
    // Default constructor
    Advection()
      : value([](const VectorRd x, const Cell *cell) -> VectorRd { return VectorRd::Zero();}),
        divergence([](const VectorRd x, const Cell *cell) -> double { return 0;}),
        is_zero(true),
        is_divergence_zero(true),
        is_divergence_constant(true)
    {
      // Do nothing
    }
    
    // Members
    CellFType<VectorRd> value;
    CellFType<double> divergence;
    bool is_zero;
    bool is_divergence_zero;
    bool is_divergence_constant;
  };

/// Structure to store a reaction term, together with various flags.
struct Reaction
  {
    // Constructor
    Reaction(
      const CellFType<double> value,   /// Reaction coefficient
      const bool is_zero,                /// True if coefficient is zero
      const bool is_constant             /// True if coefficient is constant
    )
    : value(value),
      is_zero(is_zero),
      is_constant(is_constant)
    {
      // Do nothing
    }
    
    // Default constructor
    Reaction()
      : value([](const VectorRd x, const Cell *cell) -> double { return 0.;}),
        is_zero(true),
        is_constant(true)
    {
      // Do nothing
    }
    
    // Members
    CellFType<double> value;
    bool is_zero;
    bool is_constant;
  };
  
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The TestCase class provides definition of test cases
class TestCase
{

public:
     /// Initialise data
     TestCase(
         std::vector<int> iTC ///< The vector id of the test case: (id of solution, id of diffusion)
     );

     /// Returns the solution
     inline Solution get_solution();

     /// Returns the diffusion
     inline Diffusion get_diffusion();

     /// Returns the advection
     inline Advection get_advection();

     /// Returns the reaction
     inline Reaction get_reaction();

     /// Returns the diffusion source term
     CellFType<double> diff_source();

     /// Returns the diffusion-advection-reaction source term
     CellFType<double> diff_advec_reac_source();

     /// Returns the value of the parameter lambda
     inline double get_lambda();

     /// Check if the provided test cases are valid (within range, and combination of solution/diffusion valid)
     void validate();



private:
     // Parameters: id of test case, pi
     std::vector<int> m_iTC;
     const double pi = acos(-1);
     const double sqrt2 = std::pow(2, 0.5);
     const double eps = 1e-5; // constant for rotating diffusion case

     // Generic parameter (e.g. for changing anisotropy etc.)
     double _lambda;
     
     // Advection, reaction, etc and functions to create them
     Solution m_sol;
     Diffusion m_diff;
     Advection m_advec;
     Reaction m_reac;

     Solution create_solution(const size_t n);
     Diffusion create_diffusion(const size_t n);
     Advection create_advection(const size_t n);
     Reaction create_reaction(const size_t n);

};

inline double TestCase::get_lambda() { return _lambda; };
inline Solution TestCase::get_solution() { return m_sol; };
inline Diffusion TestCase::get_diffusion() { return m_diff; };
inline Advection TestCase::get_advection() { return m_advec; };
inline Reaction TestCase::get_reaction() { return m_reac; };

//@}

#endif //_TEST_CASE_HPP
