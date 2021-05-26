// Structure to provide weight functions for integrals: provides value, and polynomial degree (cell by cell)

/*
 *
 *      This library was developed around HHO methods, although some parts of it have a more
 * general purpose. If you use this code or part of it in a scientific publication, 
 * please mention the following book as a reference for the underlying principles
 * of HHO schemes:
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
 *  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
 *  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
 *  url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 */

#ifndef INTEGRALWEIGHT_HPP
#define INTEGRALWEIGHT_HPP

#include <cell.hpp>

namespace HArDCore2D
{

  /*!
   *	\addtogroup Common
   * @{
   */
   
    /// Structure for weights (scalar, at the moment) in integral.
    /** Each weight is represented as a piecewise function defined cell-by-cell, together with its local polynomial degree to determine the offset for quadrature rules (degree=0 means that the weight is constant and some calculations are easier) */
    struct IntegralWeight
    {
       // Generic constructor
        IntegralWeight(
              const std::function<double (const Cell & T, const Eigen::Vector2d & x)> _value, ///< Value of weight
              std::function<size_t (const Cell & T)> _deg ///< Local degree of weight
              )
        : value(_value),
          deg(_deg)
        {
          // Do nothing
        }

       // Easy constructor for constant weights
       IntegralWeight(double val)
        : IntegralWeight( [val](const Cell &T, const Eigen::Vector2d &x)->double {return val;}, [](const Cell &T)->size_t {return 0;} )
       {
        // Do nothing
       }
        
      std::function<double (const Cell & T, const Eigen::Vector2d & x)> value;
      std::function<size_t (const Cell & T)> deg;
    };

  //@}

} // end of namespace HArDCore2D

#endif
