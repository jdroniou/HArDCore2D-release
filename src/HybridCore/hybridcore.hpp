// Data structures and methods to implement hybrid schemes in 2D, with polynomial unknowns
// in the cells and on the edges, such as Hybrid High-order (HHO) schemes.
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

#ifndef HYBRIDCORE_HPP
#define HYBRIDCORE_HPP

#include <cassert>
#include <cmath>

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <mesh.hpp>
#include <cell.hpp>
#include <edge.hpp>
#include <quadraturerule.hpp>
#include <basis.hpp>
#include <iostream>
#include <parallel_for.hpp>

/*!	
 * @defgroup HybridCore 
 * @brief Classes providing tools to implement schemes having polynomial unknowns in the cells and on the edges
 */

namespace HArDCore2D {


  /*!
   *	\addtogroup HybridCore
   * @{
   */

  // -----------------------------------------
  //          Free Functions
  // -----------------------------------------
  
  // Dimension of polynomial spaces -- only specialisations are relevant
  template<typename GeometricSupport>
  inline size_t const DimPoly(int m);

  /// Compute the size of the basis of 2-variate polynomials up to degree m
  template<>
  inline const size_t DimPoly<Cell>(const int m)   /**< Polynomial degree */
    {
      return (m >= 0 ? (m + 1) * (m + 2) / 2 : 0);
    };

  /// Compute the size of the basis of 1-variate polynomials up to degree m
  template<>
  inline const size_t DimPoly<Edge>(const int m)   /**< Polynomial degree */
    {
      return (m >= 0 ? m + 1 : 0);
    }

  //-----------------------------------------------------------------------------
  //                          UVector class definition
  //-----------------------------------------------------------------------------
  /** The UVector class describes a vector of unknowns for discrete functions defined by polynomials in the cells
and polynomials on the edges (polynomials on the edges are not taken into account yet, but this can easily be changed).
    The class describes the respective degrees of cell and edges polynomials considered. In the values stored in an element of this class, the coefficients on the cell basis functions come first (in the order of the cells), and then all the edge basis functions (in the order of the edges) **/
  class UVector
  {
    public:
      UVector(
        const Eigen::VectorXd values,   ///< values of the vector
        const Mesh& mesh,               ///< reference to the mesh
        const int cell_deg,             ///< polynomial degrees in cell
        const size_t edge_deg           ///< polynomial degrees on edge
      );

    /// Return the values as an Eigen vector
    inline Eigen::VectorXd & asVectorXd() const {
      return m_values;
    }

    /// Return the cell degree
    inline const int get_cell_deg() const{
      return m_cell_deg;
    }

    /// Return the edge degree
    inline const size_t get_edge_deg() const{
      return m_edge_deg;
    }
  
    /// Number of dofs in each cell
    inline const size_t n_cell_dofs() const{
      return DimPoly<Cell>( std::max(m_cell_deg,0) );
    }

    /// Number of dofs on each edge
    inline const size_t n_edge_dofs() const{
      return DimPoly<Edge>( m_edge_deg );
    }

    /// Total number of cell dofs (in the vector, this is where the edge dofs start)
    inline const size_t n_total_cell_dofs() const{
      return m_mesh.n_cells() * n_cell_dofs();
    }

    /// Extract the restriction of the unknowns corresponding to cell iT and its edges
    Eigen::VectorXd restr(size_t iT) const;

    /// Overloads the addition: adds the coefficients
    UVector operator+(const UVector& b){
      assert(m_cell_deg == b.get_cell_deg() || m_edge_deg == b.get_edge_deg() );
      return UVector(m_values + b.asVectorXd(), m_mesh, m_cell_deg, m_edge_deg);
    }

    /// Overloads the subtraction: subtracts the coefficients
    UVector operator-(const UVector& b){
      assert(m_cell_deg == b.get_cell_deg() || m_edge_deg == b.get_edge_deg() );
      return UVector(m_values - b.asVectorXd(), m_mesh, m_cell_deg, m_edge_deg);
    }

    /// Overloads the (): returns the corresponding coefficient
    double operator()(size_t index) const{
      return m_values(index);
    }

    private:
      mutable Eigen::VectorXd m_values;
      const Mesh& m_mesh;
      const int m_cell_deg;
      const size_t m_edge_deg;

  };


  // ----------------------------------------------------------------------------
  //                            HybridCore class definition
  // ----------------------------------------------------------------------------

  /** The HybridCore class provides an interface for generating polynomial basis functions on cell and edges,
  *  interpolation of continuous functions, discrete norms of vectors of coefficients, and methods to evaluate
  *  discrete functions (given by vectors of coefficients) in the cells, on the edges, or at vertices (averaged of
  *  cell or edge values)
  *
  *   The current implementation has the following behaviours/expectations:
  *     - Edge polynomials must be at least of degree 0, and edge basis functions are always generated
  *     - Cell polynomials could be of degree -1, or 0+. In the former case, basis functions of degree 0 are generated,
  *     but a function is provided to compute weights to express the cell values in terms of linearly exact averages of
  *     edge values. This function is used, e.g., when interpolating a continuous function.
  *
  **/

  class HybridCore 
  {

  public:
    ///@brief Class constructor: initialises the data structure with the given mesh, and desired polynomial degrees of the basis functions.
    /** The orthonormalisation comes at a cost in terms of manipulation of the basis functions. 
  This should only be used when the polynomial degree is large and/or the cell is distorted. 
  However, in these cases, it can make a huge difference on the observed convergence rate. */
    HybridCore(
	       const Mesh* mesh_ptr, ///< A pointer to the loaded mesh
	       const int cell_deg, ///< The degree of the cell polynomials 
	       const size_t edge_deg, ///< The degree of the edge polynomials 
	       const bool use_threads = true, ///< Optional argument to indicate if threads should be used
	       std::ostream & output = std::cout ///< Optional argument for specifying outputs of messages.
	       ); 

    //---- Basis functions types -----//
    typedef Family<MonomialScalarBasisCell> PolyCellBasisType; ///< type for cell basis
    typedef Family<MonomialScalarBasisEdge> PolyEdgeBasisType; ///< type for edge basis

    //-------- Getters --------//
    /// Returns a pointer to the mesh
    inline const Mesh* get_mesh() const {return m_mesh;}

    /// Return the degree of cell polynomials
    inline const int CellDegree() const { return m_cell_deg; };
    inline const int CellDegreePos() const { return m_cell_deg_pos; };
    /// Return the degree of edge polynomials
    inline const size_t EdgeDegree() const { return m_edge_deg; };
 
    /// Return cell basis for element with global index iT
    inline const PolyCellBasisType & CellBasis(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_cell_basis[iT] );
      return *m_cell_basis[iT].get();
    }

    /// Return edge basis for edge with global index iE
    inline const PolyEdgeBasisType & EdgeBasis(size_t iE) const
    {
      // Make sure that the basis has been created
      assert( m_edge_basis[iE] );
      return *m_edge_basis[iE].get();
    }

    //-------- Functions ------------//
    /// Compute L2 norm of a discrete function (using cell values)
    double L2norm(
      const UVector &Xh     ///< Vector of unknowns
    ) const;

    /// Compute discrete H1 norm of a discrete function
    double H1norm(
      const UVector &Xh     ///< Vector of unknowns
      ) const;

    /// Compute the interpolant in the discrete space of a continuous function
    template<typename ContinuousFunction>
    UVector interpolate(
				const ContinuousFunction& f, 	///< function to interpolate
        const int deg_cell,    ///< degree of cell polynomials for interpolation
        const size_t deg_edge,    ///< degree of edge polynomials for interpolation
				size_t doe	///< degree of exactness of the quadrature rules to compute the interpolate
				) const;  /**< @returns vector Xh of coefficients on the basis functions, as described in class UVector. **/

	  /// Computes the weights to get cell values from edge values when l=-1
    Eigen::VectorXd compute_weights(size_t iT) const;
	
    /// Evaluates a discrete function in the cell iT at point x
    double evaluate_in_cell(const UVector Xh, size_t iT, VectorRd x) const; 
    /// Evaluates a discrete function on the edge iE at point x
    double evaluate_in_edge(const UVector Xh, size_t iE, VectorRd x) const;

    /// From a hybrid function, computes a vector of values at the vertices of the mesh
    Eigen::VectorXd VertexValues(
				 const UVector Xh,   ///< vector of discrete unknowns on cell and edge polynomials
				 const std::string from_dofs ///< Type of unknowns to use: "cell" or "edge"
				 );

  private:
    // Mesh
    const Mesh* m_mesh;  // Pointer to mesh data
    // Degree of the cell polynomials (can be -1).
    const int m_cell_deg;
    const size_t m_cell_deg_pos;
    // Degree of the edge polynomials
    const int m_edge_deg;
    // Using threads
    const bool m_use_threads;
    // Output stream
    std::ostream & m_output;

    // Cell and edges bases
    std::vector<std::unique_ptr<PolyCellBasisType>> m_cell_basis;
    std::vector<std::unique_ptr<PolyEdgeBasisType>> m_edge_basis;

    // Creates the cell and edge bases
    PolyCellBasisType _construct_cell_basis(size_t iT);
    PolyEdgeBasisType _construct_edge_basis(size_t iF);

    // offset for quadrature rules, should be 0 except for testing purposes
    int _offset_doe;	

  };



  // -------------------------------------------------------------
  // ------- Interpolates continuous function on discrete space

  template<typename ContinuousFunction>
  UVector HybridCore::interpolate(const ContinuousFunction& f, int cell_deg, size_t edge_deg, size_t doe) const {
        
    size_t cell_deg_pos = std::max(cell_deg, 0);
    assert(cell_deg <= CellDegree() || edge_deg <= EdgeDegree() || cell_deg > -1);
    
    // Sizes of local interpolation
    size_t n_cell_dofs = DimPoly<Cell>(cell_deg_pos);
    size_t n_total_cell_dofs = m_mesh->n_cells() * n_cell_dofs;
    size_t n_edge_dofs = DimPoly<Edge>(edge_deg);
    size_t n_total_edge_dofs = m_mesh->n_edges() * n_edge_dofs;
    Eigen::VectorXd XTF = Eigen::VectorXd::Zero( n_total_cell_dofs + n_total_edge_dofs );

    // Edge projections
    std::function<void(size_t, size_t)> construct_all_edge_projections
      = [&](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          Edge E = *m_mesh->edge(iE);

          // Compute the L2 projection on edge
          QuadratureRule quadE = generate_quadrature_rule(E, doe);
          boost::multi_array<double, 2> phiE_quadE = evaluate_quad<Function>::compute(EdgeBasis(iE), quadE);
          Eigen::VectorXd UF = l2_projection<PolyEdgeBasisType>(f, EdgeBasis(iE), quadE, phiE_quadE);

          // Fill in the complete vector of DOFs
          XTF.segment(n_total_cell_dofs + iE * n_edge_dofs, n_edge_dofs) = UF;
        }
      };
    parallel_for(m_mesh->n_edges(), construct_all_edge_projections, m_use_threads);

    // Cell projections
    std::function<void(size_t, size_t)> construct_all_cell_projections
      = [&](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          Cell* cell = m_mesh->cell(iT);

          // L2 projection in cell
          // Mass matrix in cell
          QuadratureRule quadT = generate_quadrature_rule(*cell, doe);
          boost::multi_array<double, 2> phiT_quadT = evaluate_quad<Function>::compute(CellBasis(iT), quadT);
          Eigen::MatrixXd MT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, n_cell_dofs, n_cell_dofs, "sym");

          // Vector of integrals of f against basis functions
          Eigen::VectorXd bT = Eigen::VectorXd::Zero(n_cell_dofs);
          for (size_t i = 0; i < n_cell_dofs; i++) {
          	for (size_t iqn = 0; iqn < quadT.size(); iqn++){
          	  bT(i) += quadT[iqn].w * phiT_quadT[i][iqn] * f(quadT[iqn].vector());
          	}
          }

          // Vector of coefficients (on cell basis functions) of the L2(T) projection of f
          Eigen::VectorXd UT = MT.ldlt().solve(bT);

          // Fill in the complete vector of DOFs
          size_t offset_T = iT * n_cell_dofs;
          XTF.segment(offset_T, n_cell_dofs) = UT;

          // Special case of L=-1, we replace the previously computed cell value with the proper average of edge values
          if ( cell_deg==-1 ){
            size_t nedges = cell->n_edges();
	          // Weights corresponding to the values in cells/on edges
	          Eigen::VectorXd barycoefT = compute_weights(iT);
	          // We transform these weights so that they apply to the coefficients on the basis functions
	          VectorRd xT = cell->center_mass();
	          double phiT_cst = CellBasis(iT).function(0, xT);
	          for (size_t ilE = 0; ilE < nedges; ilE++){
	            VectorRd xE = cell->edge(ilE)->center_mass();
	            size_t iE = cell->edge(ilE)->global_index();
	            double phiE_cst = EdgeBasis(iE).function(0, xE);
	            barycoefT(ilE) *= phiE_cst / phiT_cst;
	          }

	          XTF(iT) = 0;
	          for (size_t ilE = 0; ilE < nedges; ilE++) {
	            size_t iE = cell->edge(ilE)->global_index();
	            XTF(iT) += barycoefT(ilE) * XTF(n_total_cell_dofs + iE);
	          }
          }

        }
      };
    parallel_for(m_mesh->n_cells(), construct_all_cell_projections, m_use_threads);

    return UVector(XTF, *get_mesh(), cell_deg, edge_deg);
  }



  // --------------------------------------------------------------------------------------------------
  // ------- Functions that return class elements



  //@}

}  // end of namespace HArDCore2D

#endif /* HYBRIDCORE_HPP */
