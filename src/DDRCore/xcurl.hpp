#ifndef XCURL_HPP
#define XCURL_HPP

#include <ddrcore.hpp>
#include <ddrspace.hpp>
#include <integralweight.hpp>
#include <xgrad.hpp>

namespace HArDCore2D
{
  /*!
   *  \addtogroup DDRCore
   * @{
   */

  /// Discrete Hcurl space: local operators, L2 product and global interpolator
  class XCurl : public DDRSpace
  {
  public:
    typedef std::function<Eigen::Vector2d(const Eigen::Vector2d &)> FunctionType;

    /// A structure to store the local operators (curl and potential)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _curl,     ///< Curl operator
                     const Eigen::MatrixXd & _potential ///< Potential operator
                     )
        : curl(_curl),
          potential(_potential)
      {
        // Do nothing
      }

      Eigen::MatrixXd curl;
      Eigen::MatrixXd potential;
    };

    /// Constructor
    XCurl(const DDRCore & ddr_core, bool use_threads = true, std::ostream & output = std::cout);

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

    /// Return cell operators for the cell of index iT
    inline const LocalOperators & cellOperators(size_t iT) const
    {
      return *m_cell_operators[iT];
    }

    /// Return cell operators for cell T
    inline const LocalOperators & cellOperators(const Cell & T) const
    {
      return * m_cell_operators[T.global_index()];  
    }

    /// Return cell bases for the cell of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_ddr_core.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_ddr_core.cellBases(T.global_index());
    }
        
    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_ddr_core.edgeBases(iE);
    }

    /// Return edge bases for edge E
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
    LocalOperators _compute_cell_curl_potential(size_t iT);
    
    const DDRCore & m_ddr_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };

} // end of namespace HArDCore2D
#endif
