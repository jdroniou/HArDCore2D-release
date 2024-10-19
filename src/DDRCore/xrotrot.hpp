#ifndef XROTROT_HPP
#define XROTROT_HPP

#include <globaldofspace.hpp>
#include <ddrcore.hpp>
#include <integralweight.hpp>
#include <xgrad.hpp>
#include <xrot.hpp>
#include <GMpoly_cell.hpp>

namespace HArDCore2D
{
  /*!
   *  \addtogroup DDRCore
   * @{
   */

  /// Discrete HRotRot space: local operators, L2 product and global interpolator
  class XRotRot : public GlobalDOFSpace
  {
  public:
    typedef std::function<Eigen::Vector2d(const Eigen::Vector2d &)> FunctionType;
    typedef std::function<double(const Eigen::Vector2d &)> RotType;

    /// A structure to store the local operators (scalar rotor and potential)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _rotor,         ///< Rot operator
                     const Eigen::MatrixXd & _potential      ///< Potential operator
                     )
        : rotor(_rotor),
          potential(_potential)
      {
        // Do nothing
      }

      Eigen::MatrixXd rotor;
      Eigen::MatrixXd potential;
    };

    /// Constructor
    XRotRot(const DDRCore & ddr_core, bool use_threads = true, std::ostream & output = std::cout);

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
                                const RotType & rot_v,  ///< Rotor of the function to interpolate
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

    /// Return cell rotor for cell of index iT
    inline Eigen::MatrixXd cellRotor(size_t iT) const
    {
      const Cell & T = *mesh().cell(iT);
      return cellOperators(iT).rotor.block(m_rot_dofs.localOffset(T), 0, PolynomialSpaceDimension<Cell>::Poly(degree()), dimensionCell(iT));
    }

    /// Return cell rotor for cell T
    inline Eigen::MatrixXd cellRotor(const Cell & T) const
    {
      return cellRotor(T.global_index());
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
                                     size_t iT,                                         ///< Index of the cell
                                     const double & penalty_factor = 1.,                ///< Pre-factor for stabilisation term
                                     const IntegralWeight & weight = IntegralWeight(1.) ///< Weight function in the L2 product, defaults to 1
                                     ) const;
                                     
    /// Compute the matrix of the (weighted) L2-product as 'computeL2Product', with application of the discrete gradient on the left side
    Eigen::MatrixXd computeGradientPotentialL2Product(
                                                      size_t iT,                                         ///< Index of the cell
                                                      const XGrad * x_grad,                              ///< Instance of XGrad to access the full gradients
                                                      const double & penalty_factor = 1.,                ///< Pre-factor for stabilisation term
                                                      const IntegralWeight & weight = IntegralWeight(1.) ///< Weight function in the L2 product, defaults to constant 1.
                                                      ) const;

    /// Compute the matrix of the (weighted) gradient-gradient L2-product
    Eigen::MatrixXd computeGradientL2Product(
                                             size_t iT,                                         ///< Index of the cell
                                             const XGrad * x_grad,                              ///< Instance of XGrad to access the full gradients
                                             const double & penalty_factor = 1.,                ///< Pre-factor for stabilisation term
                                             const IntegralWeight & weight = IntegralWeight(1.) ///< Weight function in the L2 product, defaults to constant 1.
                                             ) const;

    
    /// Compute the matrix of the L2 product given operators filling functions
    template<typename LeftOperatorFillerType, typename RightOperatorFillerType>
    Eigen::MatrixXd computeL2ProductWithOperatorFillers(
                                                        size_t iT,                          ///< Index of the cell
                                                        const double & penalty_factor,      ///< Penalty factor for stabilisation term
                                                        const IntegralWeight & weight,      ///< Weight function in the L2 product
                                                        LeftOperatorFillerType fillLeftOp,  ///< Function to fill the left-hand side operator
                                                        RightOperatorFillerType fillRightOp ///< Function to fill the right-hand side operator
                                                        ) const;

    /// Compute the L2-norm of a vector of the space
    double computeL2Norm(const Eigen::VectorXd & v) const;

    /// Compute the L2-norm of the discrete gradient of a vector in XGrad
    double computeGradientL2Norm(
                                 const Eigen::VectorXd & v,
                                 const XGrad * x_grad
                                 ) const;
    
  private:
    LocalOperators _compute_cell_rotor_potential(size_t iT);
    double _compute_squared_l2_norm(size_t iT, const Eigen::VectorXd & vT) const;
    double _compute_squared_gradient_l2_norm(
                                             size_t iT,
                                             const XGrad * x_grad,
                                             const Eigen::VectorXd & vT
                                             ) const;
    
    const DDRCore & m_ddr_core;
    bool m_use_threads;
    std::ostream & m_output;

    // DOFs for the rot (these DOFs coincide with those of XRot)
    GlobalDOFSpace m_rot_dofs;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;

  };

  //------------------------------------------------------------------------------
  // Implementation of template functions
  //------------------------------------------------------------------------------
  
  template<typename LeftOperatorFillerType, typename RightOperatorFillerType>
  Eigen::MatrixXd XRotRot::computeL2ProductWithOperatorFillers(
                                                               size_t iT,
                                                               const double & penalty_factor,
                                                               const IntegralWeight & weight,
                                                               LeftOperatorFillerType fillLeftOp,
                                                               RightOperatorFillerType fillRightOp
                                                               ) const
  {
    const Cell & T = *mesh().cell(iT);

    // Create the weighted mass matrix, with simple product if weight is constant
    Eigen::MatrixXd w_mass_Pk2_T;
    if (weight.deg(T) == 0) { // Constant weight
      MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2); 
      w_mass_Pk2_T = weight.value(T, T.center_mass()) * GramMatrix(T, *cellBases(iT).Polyk2, int_mono_2kp2);
    } else { // Weight is not constant, we create a weighted mass matrix
      QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
      auto basis_Pk2_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2kpw_T);
      std::function<double(const Eigen::Vector2d &)> weight_T 
        = [&T, &weight](const Eigen::Vector2d &x)->double {
          return weight.value(T, x);
        };
      w_mass_Pk2_T = compute_weighted_gram_matrix(weight_T, basis_Pk2_T_quad, basis_Pk2_T_quad, quad_2kpw_T, "sym");
    }

    // Fill left-hand side and right-hand side operators
    std::vector<Eigen::MatrixXd> leftOp  = fillLeftOp(iT);
    std::vector<Eigen::MatrixXd> rightOp = fillRightOp(iT);

    // Compute the L2-product
    double hT = T.diam();

    size_t offset_PT = T.n_vertices() + 2 * T.n_edges();
    size_t offset_RT = offset_PT + 1;

    Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(leftOp[0].cols(), rightOp[0].cols());
  
    // Vertex penalty terms
    for (size_t iV = 0; iV < T.n_vertices(); iV++) {
      const Vertex & V = *T.vertex(iV);
      VectorRd xV = V.coords();
      Eigen::MatrixXd basis_PkT_xV = Eigen::MatrixXd::Zero(1, cellBases(iT).Polyk->dimension());
      for (size_t i = 0; i < cellBases(iT).Polyk->dimension(); i++) {
        basis_PkT_xV(0, i) = cellBases(iT).Polyk->function(i, xV);
      } // for
      Eigen::MatrixXd left_RT_xV = basis_PkT_xV * leftOp[offset_RT];
      Eigen::MatrixXd right_RT_xV = basis_PkT_xV * rightOp[offset_RT];

      double w_hT4_xV = weight.value(T, xV) * std::pow(hT, 4);
      L2P += w_hT4_xV * (left_RT_xV - leftOp[iV]).transpose() * (right_RT_xV - rightOp[iV]);
    } // for iV

    // Edge penalty terms
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
      const Edge & E = *T.edge(iE);
      VectorRd tE = E.tangent();

      size_t offset_PE = T.n_vertices() + 2 * iE;
      size_t offset_RE = offset_PE + 1;
        
      QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());    
    
      // Weight and scaling hE
      double max_weight_quad_E = weight.value(T, quad_2k_E[0].vector());
      // If the weight is not constant, we want to take the largest along the edge
      if (weight.deg(T)>0) {
        for (size_t iqn = 1; iqn < quad_2k_E.size(); iqn++) {
          max_weight_quad_E = std::max(max_weight_quad_E, weight.value(T, quad_2k_E[iqn].vector()));
        } // for
      } // if
      double w_hE = max_weight_quad_E * E.measure();
           
      // Penalty term h_E int_E (PT w . tE - w_E) * (PT v . tE - v_E)
      auto basis_Pk2_T_dot_tE_quad = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2k_E), tE);
      auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polyk, quad_2k_E);
      Eigen::MatrixXd gram_Pk2T_dot_tE_PkE = compute_gram_matrix(basis_Pk2_T_dot_tE_quad, basis_Pk_E_quad, quad_2k_E);
    
      L2P += w_hE * (
                     leftOp[offset_PT].transpose() * compute_gram_matrix(basis_Pk2_T_dot_tE_quad, quad_2k_E) * rightOp[offset_PT]
                     - leftOp[offset_PT].transpose() * gram_Pk2T_dot_tE_PkE * rightOp[offset_PE]
                     - leftOp[offset_PE].transpose() * gram_Pk2T_dot_tE_PkE.transpose() * rightOp[offset_PT]
                     + leftOp[offset_PE].transpose() * compute_gram_matrix(basis_Pk_E_quad, quad_2k_E) * rightOp[offset_PE]
                     );

      if(degree() >= 1) {
        double w_hE3 = max_weight_quad_E * pow(E.measure(), 3);
         
        // Penalty term h_E^3 \int_E (pi^{k-1}_E RT w - C_{w,E}) * (pi^{k-1}_E RT v - C_{v,E})
        auto basis_PkT_E_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2k_E);
        auto basis_PkmoT_E_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykmo, quad_2k_E);
        auto basis_PkmoE_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykmo, quad_2k_E);

        Eigen::MatrixXd mass_PkmoE = compute_gram_matrix(basis_PkmoE_E_quad, quad_2k_E);
        Eigen::MatrixXd gram_PkT_PkmoE = compute_gram_matrix(basis_PkT_E_quad, basis_PkmoE_E_quad, quad_2k_E);
        Eigen::MatrixXd gram_PkmoT_PkmoE = compute_gram_matrix(basis_PkmoT_E_quad, basis_PkmoE_E_quad, quad_2k_E);

        // Compute the projections on Pk-1(E) of the left and right "rotor" operators
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_mass_PkmoE(mass_PkmoE);
        Eigen::MatrixXd pikmoE_left_RT = lu_mass_PkmoE.solve(gram_PkT_PkmoE.transpose() * leftOp[offset_RT]);
        Eigen::MatrixXd pikmoE_right_RT = lu_mass_PkmoE.solve(gram_PkT_PkmoE.transpose() * rightOp[offset_RT]);
        
        L2P += w_hE3 * (pikmoE_left_RT - leftOp[offset_RE]).transpose() * mass_PkmoE * (pikmoE_right_RT - rightOp[offset_RE]);
      } // if
    } // for iE

    L2P *= penalty_factor;
  
    // Consistent (cell) term
    L2P += leftOp[offset_PT].transpose() * w_mass_Pk2_T * rightOp[offset_PT];
   
    return L2P;
  }

} // end of namespace HArDCore2D
#endif
