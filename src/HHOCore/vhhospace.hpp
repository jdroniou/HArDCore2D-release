// Core data structures and methods required to implement the Hybrid High-Order in 2D, for vector-valued functions
//
// Provides:
//  - Polynomial spaces on the element and edges
//  - Interpolator of smooth functions
//  - Full gradient, potential and stabilisation bilinear form in the elements
//
// Author: Alessandra Crippa (alessandra.crippa@ifpen.fr)
//

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


#ifndef VHHOSPACE_HPP
#define VHHOSPACE_HPP

#include <iostream>
#include <globaldofspace.hpp>
#include <basis.hpp>
#include <polynomialspacedimension.hpp>


namespace HArDCore2D
{

  /*!
   *	\addtogroup HHOSpace
   * @{
   */


  //------------------------------------------------------------------------------

  // Type for selecting some cells (return 0 or 1 for each cell T)
  typedef std::function<bool(const Cell &)> CellSelection;
  static const CellSelection allcells = [](const Cell &)->bool {return true;};

  /// Class definition: polynomial bases and operators
  class VHHOSpace : public GlobalDOFSpace
  {
  public:
    // Types for element bases
    typedef Family<MonomialScalarBasisCell> PolyBasisCellType;
    typedef TensorizedVectorFamily<PolyBasisCellType, dimspace> PolydBasisCellType;
    typedef MatrixFamily<PolyBasisCellType, dimspace> PolydxdBasisCellType;

    // Types for edge bases
    typedef Family<MonomialScalarBasisEdge> PolyBasisEdgeType;
    typedef Family<TensorizedVectorFamily<PolyBasisEdgeType, dimspace>> PolydBasisEdgeType;

    // Types for functions to interpolate
    typedef std::function<VectorRd(const VectorRd &)> FunctionType;
    
    /// Structure to store element bases
    /** 'Poly': basis of polynomial space.\n
        'k' and 'kpo' (k+1) determines the degree.\n
        'd' for vector-valued, 'dxd' for matrix-valued
      */
    struct CellBases
    {
      /// Geometric support
      typedef Cell GeometricSupport;

      std::unique_ptr<PolyBasisCellType> Polykpo;
      std::unique_ptr<PolyBasisCellType> Polyk;
      std::unique_ptr<PolydBasisCellType> Polykpod;
      std::unique_ptr<PolydBasisCellType> Polykd;
      std::unique_ptr<PolydxdBasisCellType> Polykdxd;
    };

    /// Structure to store edge bases
    /** See CellBases for details */
    struct EdgeBases
    {
      /// Geometric support
      typedef Edge GeometricSupport;

      std::unique_ptr<PolyBasisEdgeType> Polyk;
      std::unique_ptr<PolydBasisEdgeType> Polykd;
    };
    
    /// A structure to store local operators (gradient, potential, stabilisation)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _gradient, ///< Gradient operator
                     const Eigen::MatrixXd & _potential, ///< HHO Potential operator
                     const Eigen::MatrixXd & _stabilisation ///< H1 Stabilisation bilinear form associated to the HHO potential
                     )
        : gradient(_gradient),
          potential(_potential),
          stabilisation(_stabilisation)
      {
        // Do nothing
      }
      
      Eigen::MatrixXd gradient;
      Eigen::MatrixXd potential;
      Eigen::MatrixXd stabilisation;
    };

    /// Constructor (with function to select cells in which boundary faces are used in the stabilisation)
    VHHOSpace(const Mesh & mesh, size_t K, const CellSelection & BoundaryStab, bool use_threads = true, std::ostream & output = std::cout);

    /// Overloaded constructor when the selection of boundary stabilisation is not entered (all boundary faces are then used)
    VHHOSpace(const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout)
            : VHHOSpace(mesh, K, allcells, use_threads, output) {};
    
    /// Return a const reference to the mesh
    const Mesh & mesh() const
    {
      return m_mesh;
    }

    /// Return the polynomial degree (common face and elements)
    const size_t & degree() const
    {
      return m_K;
    }

    /// Return the function to select the cells with boundary stabilisation
    const CellSelection & boundaryStab() const
    {
      return m_boundary_stab;
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const FunctionType & q, ///< The function to interpolate
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1 ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;

    /// Local Interpolator of a continuous function on the element iT
    Eigen::VectorXd local_interpolate(
            size_t iT, ///< The element of the local interpolator
            const FunctionType & q, ///< The function to interpolate
            const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
            const int doe_face = -1 ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
    ) const;

    /// Return cell bases for element iT
    inline const CellBases & cellBases(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_cell_bases[iT] );
      return *m_cell_bases[iT].get();
    }

    /// Return cell bases for cell T
    inline const CellBases & cellBases(const Cell & T) const
    {
      return cellBases(T.global_index());
    }
    
    /// Return edge bases for edge iE
    inline const EdgeBases & edgeBases(size_t iE) const
    {
      // Make sure that the basis has been created
      assert( m_edge_bases[iE] );
      return *m_edge_bases[iE].get();
    }

    /// Return cell bases for edge E
    inline const EdgeBases & edgeBases(const Edge & E) const
    {
      return edgeBases(E.global_index());
    }

    /// Return operators for the cell of index iT
    inline const LocalOperators & operators(size_t iT) const
    {
      assert( m_operators[iT] );
      return *m_operators[iT];
    }

    /// Return cell operators for cell T
    inline const LocalOperators & operators(const Cell & T) const
    {
      return operators(T.global_index());
    }
    
    /// Computes the discrete L2 (cell unknowns only) and H1 norms of a list of vectors
    std::vector<std::pair<double,double>> computeNorms(
                   const std::vector<Eigen::VectorXd> & list_dofs   ///< The list of vectors representing the dofs
                  ) const;

    /// Computes the values of the potential reconstruction at the mesh vertices
    /// and stores them in a way that is compatible with vtu_writer
    std::vector<Eigen::VectorXd> computeVertexValues(
                                                     const Eigen::VectorXd & u ///< DOFs in the discrete space
                                                     ) const;
   
  private:
    /// Compute the bases on an element T
    CellBases _construct_cell_bases(size_t iT);

    /// Compute the bases on a edge E
    EdgeBases _construct_edge_bases(size_t iE);
    
    /// Compute operators in an element T
    LocalOperators _compute_operators(size_t iT);
    
    // Pointer to the mesh
    const Mesh & m_mesh;
    // Degrees
    const size_t m_K;
    // Choice of boundary stabilisation
    CellSelection m_boundary_stab;
    // Parallel or not
    bool m_use_threads;
    // Output stream
    std::ostream & m_output;    
    
    // Cell bases
    std::vector<std::unique_ptr<CellBases> > m_cell_bases;
    // Edge bases
    std::vector<std::unique_ptr<EdgeBases> > m_edge_bases;

    // Local operators
    std::vector<std::unique_ptr<LocalOperators> > m_operators;
        
  };

  
} // end of namespace HArDCore2D

#endif // VHHOSPACE_HPP
