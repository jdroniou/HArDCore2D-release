#include <functional>
#include <xrotrot.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <operatorfillers.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XRotRot::XRotRot(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(
                   ddr_core.mesh(),
                   1,
                   PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree()) + PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree() - 1),
                   PolynomialSpaceDimension<Cell>::Roly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Cell>::RolyCompl(ddr_core.degree())        
                   ),
    m_ddr_core(ddr_core),
    m_use_threads(use_threads),    
    m_output(output),
    m_rot_dofs(
               ddr_core.mesh(),
               1,
               PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree() - 1),
               PolynomialSpaceDimension<Cell>::Poly(ddr_core.degree())
               ),
    m_cell_operators(ddr_core.mesh().n_cells())
{
  output << "[XRotRot] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XRotRot] Parallel execution" << std::endl;
  } else {
    m_output << "[XRotRot] Sequential execution" << std::endl;
  }

  // Construct cell rotors and potentials
  std::function<void(size_t, size_t)> construct_all_cell_rotors_potentials
    = [this](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        m_cell_operators[iT].reset( new LocalOperators(_compute_cell_rotor_potential(iT)) );
      } // for iT
    };

  m_output << "[XRotRot] Constructing cell rotors and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_rotors_potentials, use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XRotRot::interpolate(const FunctionType & v,
                                     const RotType & rot_v,
                                     const int deg_quad) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());
  
  // Degree of quadrature rules
  size_t dqr = (deg_quad >= 0 ? deg_quad : 2 * degree() + 3);

  // Interpolate at vertices
  std::function<void(size_t, size_t)> interpolate_vertices
    = [this, &vh, v, rot_v](size_t start, size_t end)->void
    {
      for (size_t iV = start; iV < end; iV++) {
        vh(iV) = rot_v(mesh().vertex(iV)->coords());
      } // for iV
    };
  parallel_for(mesh().n_vertices(), interpolate_vertices, m_use_threads);
  
  // Interpolate at edges
  std::function<void(size_t, size_t)> interpolate_edges
    = [this, &vh, v, rot_v, &dqr](size_t start, size_t end)->void
    {
      for (size_t iE = start; iE < end; iE++) {
        const Edge & E = *mesh().edge(iE);

        Eigen::Vector2d tE = E.tangent();
        auto v_dot_tE = [&tE, v](const Eigen::Vector2d & x)->double {
          return v(x).dot(tE);
        };

        Eigen::Index offset_E = globalOffset(E);
        
        QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr);
        auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_dqr_E);
        vh.segment(offset_E, edgeBases(iE).Polyk->dimension())
          = l2_projection(v_dot_tE, *edgeBases(iE).Polyk, quad_dqr_E, basis_Pk_E_quad);

        if (degree() > 0) {
          offset_E += PolynomialSpaceDimension<Edge>::Poly(m_ddr_core.degree());
          auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_dqr_E);
          vh.segment(offset_E, edgeBases(iE).Polykmo->dimension())
            = l2_projection(rot_v, *edgeBases(iE).Polykmo, quad_dqr_E, basis_Pkmo_E_quad);
        }
      } // for iE
    };
  parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);

  if (degree() > 0 ) {
    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &vh, v, &dqr](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          const Cell & T = *mesh().cell(iT);

          QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr );

          Eigen::Index offset_T = globalOffset(T);
          auto basis_Rkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Rolykmo, quad_dqr_T);
          vh.segment(offset_T, PolynomialSpaceDimension<Cell>::Roly(degree() - 1))
            = l2_projection(v, *cellBases(iT).Rolykmo, quad_dqr_T, basis_Rkmo_T_quad);

          offset_T += PolynomialSpaceDimension<Cell>::Roly(degree() - 1);
          auto basis_Rck_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).RolyComplk, quad_dqr_T);
          vh.segment(offset_T, PolynomialSpaceDimension<Cell>::RolyCompl(degree()))
            = l2_projection(v, *cellBases(iT).RolyComplk, quad_dqr_T, basis_Rck_T_quad);
        } // for iT
      };
    parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  } // if degree() > 0
  
  return vh;
}

//------------------------------------------------------------------------------
// Rot and potential reconstruction
//------------------------------------------------------------------------------

XRotRot::LocalOperators XRotRot::_compute_cell_rotor_potential(size_t iT)
{
  const Cell & T = *mesh().cell(iT);
  
  //------------------------------------------------------------------------------
  // Rot
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix for the element component of the rot

  MonomialCellIntegralsType int_mono_2kpo_T = IntegrateCellMonomials(T, 2 * degree() + 1);
  Eigen::MatrixXd MCT = GramMatrix(T, *cellBases(iT).Polyk, int_mono_2kpo_T);
    
  //------------------------------------------------------------------------------
  // Right-hand side matrix for the element component of the rot

  Eigen::MatrixXd BCT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk->dimension(), dimensionCell(iT));

  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    BCT.block(0, localOffset(T, E), cellBases(iT).Polyk->dimension(), edgeBases(E).Polyk->dimension())
      -= T.edge_orientation(iE) * compute_gram_matrix(
                                                      evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2k_E),
                                                      evaluate_quad<Function>::compute(*edgeBases(E).Polyk, quad_2k_E),
                                                      quad_2k_E
                                                      );
  } // for iE

  if (degree() > 0) {
    CurlBasis<DDRCore::PolyBasisCellType> rot_Pk_T(*cellBases(iT).Polyk);
    BCT.block(0, localOffset(T), cellBases(iT).Polyk->dimension(), cellBases(iT).Rolykmo->dimension())
      += GramMatrix(T, rot_Pk_T, *cellBases(iT).Rolykmo, int_mono_2kpo_T);
  } // if degree() > 0

  // Assemble the rotor
  Eigen::MatrixXd CT = Eigen::MatrixXd::Zero(m_rot_dofs.dimensionCell(iT), dimensionCell(iT));
  CT.block(m_rot_dofs.localOffset(T), 0, PolynomialSpaceDimension<Cell>::Poly(degree()), dimensionCell(iT))
    = MCT.ldlt().solve(BCT);

  if (m_ddr_core.degree() > 0 ) {
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
      const Edge & E = *T.edge(iE);
      size_t row_offset_E = m_rot_dofs.localOffset(T, E);
      size_t col_offset_E = localOffset(T, E) + PolynomialSpaceDimension<Edge>::Poly(m_ddr_core.degree());;
      size_t dim_Pkmo_E = PolynomialSpaceDimension<Edge>::Poly(m_ddr_core.degree() - 1);
      CT.block(row_offset_E, col_offset_E, dim_Pkmo_E, dim_Pkmo_E)
        = Eigen::MatrixXd::Identity(dim_Pkmo_E, dim_Pkmo_E);
    } // for iE
  } // if

  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);
    CT(m_rot_dofs.localOffset(T, V), localOffset(T, V)) = 1.;
  } // for iV

  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  auto basis_Pkpo0_T = ShiftedBasis<DDRCore::PolyBasisCellType>(*cellBases(iT).Polykpo, 1);

  size_t dim_Pk2_T = cellBases(iT).Polyk2->dimension();
  
  Eigen::MatrixXd MPT = Eigen::MatrixXd::Zero(dim_Pk2_T, dim_Pk2_T);
  Eigen::MatrixXd BPT = Eigen::MatrixXd::Zero(dim_Pk2_T, dimensionCell(iT));

  CurlBasis<ShiftedBasis<DDRCore::PolyBasisCellType> > rot_Pkpo0_T(basis_Pkpo0_T);
  MPT.topLeftCorner(rot_Pkpo0_T.dimension(), dim_Pk2_T)
    = GramMatrix(T, rot_Pkpo0_T, *cellBases(iT).Polyk2, int_mono_2kpo_T);

  if (degree() > 0) {
    MPT.bottomLeftCorner(cellBases(iT).RolyComplk->dimension(), dim_Pk2_T)
      = GramMatrix(T, *cellBases(iT).RolyComplk, *cellBases(iT).Polyk2, int_mono_2kpo_T);
    BPT.bottomRightCorner(cellBases(iT).RolyComplk->dimension(), cellBases(iT).RolyComplk->dimension())
      += GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2kpo_T);    
  } // if degree() > 0

  
  BPT.topLeftCorner(basis_Pkpo0_T.dimension(), dimensionCell(iT))
    += GramMatrix(T, basis_Pkpo0_T, *cellBases(iT).Polyk, int_mono_2kpo_T)
    * CT.block(m_rot_dofs.localOffset(T), 0, PolynomialSpaceDimension<Cell>::Poly(degree()), dimensionCell(iT));

  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * degree() + 1);
    BPT.block(0, localOffset(T, E), basis_Pkpo0_T.dimension(), edgeBases(E).Polyk->dimension())
      += T.edge_orientation(iE) * compute_gram_matrix(
                                                      evaluate_quad<Function>::compute(basis_Pkpo0_T, quad_2kpo_E),
                                                      evaluate_quad<Function>::compute(*edgeBases(E).Polyk, quad_2kpo_E),
                                                      quad_2kpo_E
                                                      );    
  } // for iE  

  return LocalOperators(CT, MPT.partialPivLu().solve(BPT));
}

//------------------------------------------------------------------------------
//        Functions to compute matrices for local L2 products
//------------------------------------------------------------------------------

Eigen::MatrixXd XRotRot::computeL2Product(
                                          size_t iT,
                                          const double & penalty_factor,
                                          const IntegralWeight & weight
                                          ) const
{
  auto potential_filler = [this](size_t iT) { return XRotRotDetail::_fill_potential_operators(iT, this); };
  
  return computeL2ProductWithOperatorFillers(
                                             iT,
                                             penalty_factor,
                                             weight,
                                             potential_filler,
                                             potential_filler
                                             );
}

//------------------------------------------------------------------------------

Eigen::MatrixXd XRotRot::computeGradientPotentialL2Product(
                                                           size_t iT,
                                                           const XGrad * x_grad,
                                                           const double & penalty_factor,
                                                           const IntegralWeight & weight
                                                           ) const
{
  auto gradient_filler = [x_grad](size_t iT) { return XRotRotDetail::_fill_gradient_operators(iT, x_grad); };
  auto potential_filler = [this](size_t iT) { return XRotRotDetail::_fill_potential_operators(iT, this); };
 
  return computeL2ProductWithOperatorFillers(
                                             iT,
                                             penalty_factor,
                                             weight,
                                             gradient_filler,
                                             potential_filler
                                             );
}

//------------------------------------------------------------------------------

Eigen::MatrixXd XRotRot::computeGradientL2Product(
                                                  size_t iT,
                                                  const XGrad * x_grad,
                                                  const double & penalty_factor,
                                                  const IntegralWeight & weight
                                                  ) const
{
  auto gradient_filler = [x_grad](size_t iT) { return XRotRotDetail::_fill_gradient_operators(iT, x_grad); };
  
  return computeL2ProductWithOperatorFillers(
                                             iT,
                                             penalty_factor,
                                             weight,
                                             gradient_filler,
                                             gradient_filler
                                             );
}

//------------------------------------------------------------------------------
// Compute local and global norms
//------------------------------------------------------------------------------

double XRotRot::_compute_squared_l2_norm(
                                         size_t iT,
                                         const Eigen::VectorXd & vT
                                         ) const
{
  return vT.transpose() * computeL2Product(iT) * vT;
}

//------------------------------------------------------------------------------

double XRotRot::_compute_squared_gradient_l2_norm(
                                                  size_t iT,
                                                  const XGrad * x_grad,
                                                  const Eigen::VectorXd & vT
                                                  ) const
{
  return vT.transpose() * computeGradientL2Product(iT, x_grad) * vT;
}

//------------------------------------------------------------------------------

double XRotRot::computeL2Norm(const Eigen::VectorXd & v) const
{
  Eigen::VectorXd l2 = Eigen::VectorXd::Zero(m_mesh.n_cells());
  
  std::function<void(size_t, size_t)> compute_squared_l2_norm_all
    = [this, &v, &l2](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        l2(iT) = _compute_squared_l2_norm(iT, restrictCell(iT, v));
      } // for iT
    };

  parallel_for(m_ddr_core.mesh().n_cells(), compute_squared_l2_norm_all, m_use_threads);

  return std::sqrt(l2.sum());
}

//------------------------------------------------------------------------------

double XRotRot::computeGradientL2Norm(
                                      const Eigen::VectorXd & v,
                                      const XGrad * x_grad
                                      ) const
{
  Eigen::VectorXd grad_l2 = Eigen::VectorXd::Zero(m_mesh.n_cells());
  
  std::function<void(size_t, size_t)> compute_squared_grad_l2_norm_all
    = [this, &v, &grad_l2, x_grad](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        grad_l2(iT)
          = _compute_squared_gradient_l2_norm(iT, x_grad, x_grad->restrictCell(iT, v));
      } // for iT
    };

  parallel_for(m_ddr_core.mesh().n_cells(), compute_squared_grad_l2_norm_all, m_use_threads);

  return std::sqrt(grad_l2.sum());
}
