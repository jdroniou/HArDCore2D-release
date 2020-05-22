// Class to create and store values of cell and edge basis functions on quadrature points
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
 *  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
 *  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
 *  url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 */

#ifndef ELEMENTQUAD_HPP
#define ELEMENTQUAD_HPP

#include <iostream>
#include <Eigen/Dense>
#include <hybridcore.hpp>
#include <quadraturerule.hpp>


/*!	
 * @defgroup HybridCore 
 * @brief Classes providing cell and edge quadrature rules, and values of basis functions at the nodes
 */

namespace HArDCore2D {

  /*!
   *	\addtogroup HybridCore
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
		const size_t iT, ///< Global index of cell
		const size_t doeT, ///< The degree of exactness for cell quadratures 
		const size_t doeF ///< The degree of exactness of edge quadratures
		); 

	  /// Returns quadrature rules in cell
    inline const QuadratureRule & get_quadT() const 
      { 
        return m_quadT; 
      };

    /// Returns quadrature rules on edge with local number ilF
    inline const QuadratureRule & get_quadF(size_t ilF) const
      {
        return m_quadF[ilF]; 
      };

    /// Returns values of cell basis functions at cell quadrature nodes
    inline const boost::multi_array<double, 2> & get_phiT_quadT() const
    { 
      return m_phiT_quadT; 
    };

    /// Returns values of gradients of cell basis functions at cell quadrature nodes
    inline const boost::multi_array<VectorRd, 2> & get_dphiT_quadT() const
    {
      return m_dphiT_quadT;
    };

    /// Returns values of cell basis functions at edge quadrature nodes, for edge with local number ilF
    inline const boost::multi_array<double, 2> & get_phiT_quadF(size_t ilF) const
    {
      return m_phiT_quadF[ilF];
    };

    /// Returns values of edge basis functions at edge quadrature nodes, for edge with local number ilF
    inline const boost::multi_array<double, 2> get_phiF_quadF(size_t ilF) const
    { 
      return m_phiF_quadF[ilF];       
    };

    /// Returns values of gradients of cell basis functions at edge quadrature nodes, for edge with local number ilF
    inline const boost::multi_array<VectorRd, 2> & get_dphiT_quadF(size_t ilF) const
    {
      return m_dphiT_quadF[ilF];
    };


    /// Builds on the fly the values of vector cell basis functions at cell quadrature nodes. The vector basis is obtained by tensorization of the scalar one: \f$(\phi_1,0), (\phi_2,0), ..., (\phi_N,0), (0,\phi_1), (0,\phi_2) ... (0,\phi_N)\f$.
    boost::multi_array<VectorRd, 2> get_vec_phiT_quadT(
      size_t degree   /// maximal degree of basis functions that is required
    ) const; ///< @returns vec_phiT_quadT such that, for r=0,..,dim-1 and i=0,..,DimPoly<Cell>(degree)-1, vec_phiT_quadT[r*DimPoly<Cell>(degree)-1 + i] has, on its row r, the values of the i-th scalar basis function at the quadrature nodes, and 0 on its other rows. 

    /// Builds on the fly the values of vector cell basis functions at edge quadrature nodes. The vector basis is obtained by tensorization of the scalar one as in get_vec_phiT_quadT.
    boost::multi_array<VectorRd, 2> get_vec_phiT_quadF(
        size_t ilF,   /// local number of edge
        size_t degree   /// maximum degree of basis functions required
    ) const; ///< @returns vec_phiT_quadT such that, for r=0,..,dim-1 and i=0,..,DimPoly<Cell>(degree)-1, vec_phiT_quadT[r*DimPoly<Cell>(degree)-1 + i] has, on its row r, the values of the i-th scalar basis function at the quadrature nodes, and 0 on its other rows. 

    /// Builds on the fly the values of vector edge basis functions at edge quadrature nodes. The vector basis is obtained by tensorization of the scalar one as in get_vec_phiT_quadT.
    boost::multi_array<VectorRd, 2> get_vec_phiF_quadF(
        size_t ilF,   /// local number of edge
        size_t degree   /// required degree of basis function
     ) const; ///< @returns vec_phiT_quadT such that, for r=0,..,dim-1 and i=0,..,DimPoly<Cell>(degree)-1, vec_phiT_quadT[r*DimPoly<Cell>(degree)-1 + i] has, on its row r, the values of the i-th scalar basis function at the quadrature nodes, and 0 on its other rows. 


  private:
    /// Mesh, cell, degrees
    const HybridCore& m_hcore;  // reference to the hybridcore instance
    const size_t m_iT; // cell number
    const size_t m_doeT; // degree of exactness of cell quadrature rules
    const size_t m_doeF; // degree of exactness of edge quadrature rules

    /// Quadrature and values of basis functions in cells
    const Cell& m_T;
    QuadratureRule m_quadT;

    boost::multi_array<double, 2> m_phiT_quadT;
    boost::multi_array<VectorRd, 2> m_dphiT_quadT;
    boost::multi_array<VectorRd, 2> m_vec_phiT_quadT;

    /// Quadratures and values of basis functions on edges
    std::vector<QuadratureRule> m_quadF;

    std::vector<boost::multi_array<double, 2>> m_phiT_quadF;
    std::vector<boost::multi_array<double, 2>> m_phiF_quadF;
    std::vector<boost::multi_array<VectorRd, 2>> m_dphiT_quadF;
    std::vector<boost::multi_array<VectorRd, 2>> m_vec_phiT_quadF;
    std::vector<boost::multi_array<VectorRd, 2>> m_vec_phiF_quadF;

  };


  //@}

}  // end of namespace HArDCore2D

#endif /* ELEMENTQUAD_HPP */
