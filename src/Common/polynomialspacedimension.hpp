// Core data structures and methods required to implement the discrete de Rham sequence in 2D
//
// Provides:
//  - Dimension of ull and partial polynomial spaces on the element and edges
//
// Authors: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr), Jerome Droniou (jerome.droniou@monash.edu)
//

/*
 *
 * This library was developed around HHO methods, although some parts of it have a more
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

/*
 * The DDR sequence has been designed in
 *
 *  Fully discrete polynomial de Rham sequences of arbitrary degree on polygons and polyhedra.
 *   D. A. Di Pietro, J. Droniou, and F. Rapetti, 33p, 2019. url: https://arxiv.org/abs/1911.03616.
 *
 * If you use this code in a scientific publication, please mention the above article.
 *
 */
#ifndef POLYNOMIALSPACEDIMENSION_HPP
#define POLYNOMIALSPACEDIMENSION_HPP

#include <cell.hpp>
#include <edge.hpp>

namespace HArDCore2D
{

  /*!	
 * @defgroup Common 
 * @brief Various general functions and classes
 */

  /*!
   *	\addtogroup Common
   * @{
   */


  /// Basis dimensions for various polynomial spaces on edges/faces/elements (when relevant): Pk, Gk, Rk and complements.
  template<typename GeometricSupport>
  struct PolynomialSpaceDimension
  {
    // Only specializations are relevant
  };
  
  template<>
  struct PolynomialSpaceDimension<Cell>
  {
    /// Dimension of Pk(T)
    static size_t Poly(int k)
    {
      return (k >= 0 ? (k + 1) * (k + 2) / 2 : 0);
    }    
    /// Dimension of Gk(T)
    static size_t Goly(int k)
    {
      return (k >= 0 ? PolynomialSpaceDimension<Cell>::Poly(k + 1) - 1 : 0);
    }
    /// Dimension of Gck(T)
    static size_t GolyCompl(int k)
    {
      return 2 * PolynomialSpaceDimension<Cell>::Poly(k) - PolynomialSpaceDimension<Cell>::Goly(k);
    }

    /// Dimension of Rk(T)
    static size_t Roly(int k)
    {
      return (k >= 0 ? PolynomialSpaceDimension<Cell>::Poly(k + 1) - 1 : 0);
    }
    /// Dimension of Rck(T)
    static size_t RolyCompl(int k)
    {
      return 2 * PolynomialSpaceDimension<Cell>::Poly(k) - PolynomialSpaceDimension<Cell>::Roly(k);
    }

    /// Dimension of Hk(T)
    static size_t Holy(int k)
    {
      return (k+4) * (k+3) / 2 - 3;
    }

    /// Dimension of Hck(T)
    static size_t HolyCompl(int k)
    {
      return (k+1) * k;
    }
  };

  template<>
  struct PolynomialSpaceDimension<Edge>
  {
    /// Dimension of Pk(E)
    static size_t Poly(int k)
    {
      return (k >= 0 ? k + 1 : 0);
    }
  };

  //@}
  
} // namespace HArDCore2D
#endif
