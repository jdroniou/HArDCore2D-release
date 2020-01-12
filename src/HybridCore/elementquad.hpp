// Class to create and store values of cell and edge basis functions on quadrature points
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*  This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/



#ifndef ELEMENTQUAD_HPP
#define ELEMENTQUAD_HPP

#include <Eigen/Dense>
#include <hybridcore.hpp>
#include <quadraturerule.hpp>


/*!  
* @defgroup HybridCore 
* @brief Classes providing cell and edge quadrature rules, and values of basis functions at the nodes
*/

namespace HArDCore2D {

/*!
*  \addtogroup HybridCore
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/** The ElementQuad class creates cell and edge quadrature rules, and vectors of values of basis functions and
*  gradients at these points
**/

class ElementQuad {

public:
  ///@brief Class constructor: loads the quadrature rules and values of basis functions/gradients at these points
  ElementQuad(
            const HybridCore& hho, ///< A reference to the hybridcore instance
            const size_t iT, ///< Number of cell
            const size_t doeT, ///< The degree of exactness for cell quadratures 
            const size_t doeF ///< The degree of exactness of edge quadratures
             ); 

  inline QuadratureRule get_quadT() const;  ///< Returns quadrature rules in cell
  inline std::vector<Eigen::ArrayXd> get_phiT_quadT() const;  ///< Returns values of cell basis functions at cell quadrature rules in cell
  inline std::vector<Eigen::ArrayXXd> get_dphiT_quadT() const;  ///< Returns values of gradients of cell basis functions at cell quadrature rules in cell

  inline QuadratureRule get_quadF(size_t ilF) const;  ///< Returns quadrature rules on edge with local number ilF
  inline std::vector<Eigen::ArrayXd> get_phiT_quadF(size_t ilF) const;  ///< Returns values of cell basis functions at cell quadrature rules on edge with local number ilF
  inline std::vector<Eigen::ArrayXd> get_phiF_quadF(size_t ilF) const;  ///< Returns values of edges basis functions at cell quadrature rules on edge with local number ilF
  inline std::vector<Eigen::ArrayXXd> get_dphiT_quadF(size_t ilF) const;  ///< Returns values of gradients of cell basis functions at cell quadrature rules on edge with local number ilF


private:
  /// Mesh, cell, degrees
  const HybridCore& _hho;  // reference to the hybridcore instance
  const size_t _iT; // cell number
  const size_t _doeT; // degree of exactness of cell quadrature rules
  const size_t _doeF; // degree of exactness of edge quadrature rules

  /// Quadrature and values of basis functions in cells
  QuadratureRule _quadT;
  std::vector<Eigen::ArrayXd> _phiT_quadT;
  std::vector<Eigen::ArrayXXd> _dphiT_quadT;

  /// Quadratures and values of basis functions on edges
  std::vector<QuadratureRule> _quadF;
  std::vector<std::vector<Eigen::ArrayXd>> _phiT_quadF;
  std::vector<std::vector<Eigen::ArrayXd>> _phiF_quadF;
  std::vector<std::vector<Eigen::ArrayXXd>> _dphiT_quadF;

};



// --------------------------------------------------------------------------------------------------
// ------- Functions that return class elements


QuadratureRule ElementQuad::get_quadT() const { return _quadT; }
std::vector<Eigen::ArrayXd> ElementQuad::get_phiT_quadT() const { return _phiT_quadT; }
std::vector<Eigen::ArrayXXd> ElementQuad::get_dphiT_quadT() const { return _dphiT_quadT; }
QuadratureRule ElementQuad::get_quadF(size_t ilF) const { return _quadF[ilF]; }
std::vector<Eigen::ArrayXd> ElementQuad::get_phiT_quadF(size_t ilF) const { return _phiT_quadF[ilF]; }
std::vector<Eigen::ArrayXd> ElementQuad::get_phiF_quadF(size_t ilF) const { return _phiF_quadF[ilF]; }
std::vector<Eigen::ArrayXXd> ElementQuad::get_dphiT_quadF(size_t ilF) const { return _dphiT_quadF[ilF]; }

//@}

}  // end of namespace HArDCore2D

#endif /* ELEMENTQUAD_HPP */
