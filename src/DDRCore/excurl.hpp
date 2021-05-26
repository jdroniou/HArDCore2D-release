#ifndef EXCURL_HPP
#define EXCURL_HPP

#include <xcurl.hpp>

namespace HArDCore2D
{
  /*!
   *  \addtogroup DDRCore
   * @{
   */

  /// Extended XCurl space, with vector-valued polynomials on the edges.
  /** Following DDRSpace, the DOFs are organised this way: DOFs of edges (for each edge: normal components, then tangential components - both oriented according to the intrinsic normal and tangent to the edge), DOFs of cells. 
  The matrices m_reduction are used to locally recover DOFs for the XCurl space by removing the normal components to the edges. */

  /// Extended discrete Hcurl space: local operators, L2 product and global interpolator
  class EXCurl : public DDRSpace
  {
  public:
    typedef std::function<Eigen::Vector2d(const Eigen::Vector2d &)> FunctionType;
    typedef MatrixFamily<RestrictedBasis<DDRCore::PolyBasisCellType>, dimspace> Polyk2x2Type;
    typedef TensorizedVectorFamily<DDRCore::PolyBasisCellType, dimspace> Polykpo2Type;
    
    /// A structure to store the local HHO operators
    struct hhoLocalOperators
    {
      hhoLocalOperators(
                     const Eigen::MatrixXd & _gradient,     ///< gradient operator
                     const Eigen::MatrixXd & _symmetric_gradient,     ///< symmetric gradient operator
                     const Eigen::MatrixXd & _potential, ///< Potential operator
                     const Eigen::MatrixXd & _stabilisation ///< Stabilisation
                     )
        : gradient(_gradient),
          symmetric_gradient(_symmetric_gradient),
          potential(_potential),
          stabilisation(_stabilisation)
      {
        // Do nothing
      }

      Eigen::MatrixXd gradient;
      Eigen::MatrixXd symmetric_gradient;
      Eigen::MatrixXd potential;
      Eigen::MatrixXd stabilisation;
    };

    /// Constructor
    EXCurl(const DDRCore & ddr_core, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_ddr_core.mesh();
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_ddr_core.degree();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
                const FunctionType & v, ///< The function to interpolate
                const int deg_quad = -1 ///< The optional degre of quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
                ) const;

    /// Return cell potential for the cell of index iT
    // We don't store the potential and curl, as they are very cheap to construct from the (stored) operators in m_xcurl
    inline const Eigen::MatrixXd ddrPotential(size_t iT) const
    {
      return m_xcurl.cellOperators(iT).potential * m_reduction[iT];
    };

    /// Return cell (scalar) curl for the cell of index iT
    inline const Eigen::MatrixXd ddrCurl(size_t iT) const
    {
      return m_xcurl.cellOperators(iT).curl * m_reduction[iT];
    };

    /// Return cell potential for cell T
    inline const Eigen::MatrixXd ddrPotential(const Cell & T) const
    {
      return ddrPotential(T.global_index());
    };

    /// Return cell (scalar) curl for cell T
    inline const Eigen::MatrixXd ddrCurl(const Cell & T) const
    {
      return ddrCurl(T.global_index());
    };

    /// Return hho operators for the cell of index iT
    inline const hhoLocalOperators & hhoOperators(size_t iT) const
    {
      assert( m_hho_operators[iT] );
      return *m_hho_operators[iT];
    }

    /// Return cell operators for cell T
    inline const hhoLocalOperators & hhoOperators(const Cell & T) const
    {
      return * m_hho_operators[T.global_index()];  
    }

    /// Return ddrcore cell bases for the cell of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_ddr_core.cellBases(iT);
    }

    /// Return ddrcore cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_ddr_core.cellBases(T.global_index());
    }
        
    /// Return basis for Matricial P^k space
    inline const Polyk2x2Type & Polyk2x2(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_Polyk2x2[iT] );
      return *m_Polyk2x2[iT].get();
    }

    /// Return basis for the (P^{k+1})^2 space
    inline const Polykpo2Type & Polykpo2(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_Polykpo2[iT] );
      return *m_Polykpo2[iT].get();
    }

    /// Return ddrcore edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_ddr_core.edgeBases(iE);
    }

    /// Return ddrcore edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_ddr_core.edgeBases(E.global_index());
    }

    /// Compute the matrix of the (weighted) L2-product for the cell of index iT.
    // The mass matrix of P^k(T)^2 is the most expensive mass matrix in the calculation of this norm, which
    // is why there's the option of passing it as parameter if it's been already pre-computed when the norm is called.
    Eigen::MatrixXd computeL2Product(
                                     const size_t iT, ///< index of the cell
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk2_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to 1
                                     ) const;
                                     
    /// Compute the matrix of the (weighted) L2-product as 'computeL2Product', with application of the discrete gradient on the left/right/both sides (depending on argument "side").
    Eigen::MatrixXd computeL2ProductGradient(
                                     const size_t iT, ///< index of the cell
                                     const XGrad & x_grad, ///< instance of XGrad to access the full gradients
                                     const std::string & side, ///< which side (left, right, both) we apply the gradient to
                                     const double & penalty_factor = 1., ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & mass_Pk2_T = Eigen::MatrixXd::Zero(1,1), ///< if pre-computed, the mass matrix of (P^k(T))^3; if none is pre-computed, passing Eigen::MatrixXd::Zero(1,1) will force the calculation
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< weight function in the L2 product, defaults to constant 1.
                                     ) const;

    /// Compute the matrix of the L2 product, applying leftOp and rightOp to the variables. Probably not directly called, mostly invoked through the wrappers computeL2Product and computeL2ProductGradient
    Eigen::MatrixXd computeL2Product_with_Ops(
                                     const size_t iT, ///< index of the cell
                                     const std::vector<Eigen::MatrixXd> & leftOp, ///< edge and element operators to apply on the left
                                     const std::vector<Eigen::MatrixXd> & rightOp, ///< edge and element operators to apply on the right
                                     const double & penalty_factor, ///< pre-factor for stabilisation term
                                     const Eigen::MatrixXd & w_mass_Pk2_T, ///< mass matrix of (P^k(T))^3 weighted by weight
                                     const IntegralWeight & weight ///< weight function in the L2 product
                                     ) const;
                                

  private:
    hhoLocalOperators _compute_hho_operators(size_t iT);
      
    const DDRCore & m_ddr_core;
    const XCurl m_xcurl; // Underlying standard XCurl space
    std::vector<std::unique_ptr<Polyk2x2Type> > m_Polyk2x2;
    std::vector<std::unique_ptr<Polykpo2Type> > m_Polykpo2;
    std::vector<Eigen::MatrixXd> m_reduction; // reduction matrix to go from EXcurl DOFs to XCurl DOFs 
    bool m_use_threads;
    std::ostream & m_output;
    
    // Containers for local hho operators
    std::vector<std::unique_ptr<hhoLocalOperators> > m_hho_operators;
  };

} // end of namespace HArDCore2D
#endif
