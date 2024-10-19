#include <xrot.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XRot::XRot(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(ddr_core.mesh(),
                   1,
                   PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree() - 1),
                   PolynomialSpaceDimension<Cell>::Poly(ddr_core.degree())
                   ),
    m_ddr_core(ddr_core),
    m_use_threads(use_threads),
    m_output(output),
    m_edge_potentials(ddr_core.mesh().n_edges()),
    m_cell_operators(ddr_core.mesh().n_cells())
{
  m_output << "[XRot] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XRot] Parallel execution" << std::endl;
  } else {
    m_output << "[XRot] Sequential execution" << std::endl;
  }
  
  // Construct potentials
  std::function<void(size_t, size_t)> construct_all_edge_potentials
    = [this](size_t start, size_t end)->void
    {
      for (size_t iE = start; iE < end; iE++) {
        m_edge_potentials[iE].reset( new Eigen::MatrixXd(_compute_edge_potential(iE)) );
      } // for iE
    };

  m_output << "[XRot] Constructing edge potentials" << std::endl;
  parallel_for(mesh().n_edges(), construct_all_edge_potentials, use_threads);

  // Construct cell rotors
  std::function<void(size_t, size_t)> construct_all_cell_operators
    = [this](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        m_cell_operators[iT].reset( new LocalOperators(_compute_cell_operators(iT)) );
      } // for iT
    };

  m_output << "[XRot] Constructing cell rotors" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_operators, use_threads);
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XRot::interpolate(const FunctionType & q, const int deg_quad) const
{
  Eigen::VectorXd qh = Eigen::VectorXd::Zero(dimension());

  // Degree of quadrature rules
  size_t dqr = (deg_quad >= 0 ? deg_quad : 2 * degree() + 3);
  
  // Interpolate at vertices
  std::function<void(size_t, size_t)> interpolate_vertices
    = [this, &qh, q](size_t start, size_t end)->void
    {
      for (size_t iV = start; iV < end; iV++) {
        qh(iV) = q(mesh().vertex(iV)->coords());
      } // for iV
    };
  parallel_for(mesh().n_vertices(), interpolate_vertices, m_use_threads);

  if (degree() > 0) {    
    // Interpolate at edges
    std::function<void(size_t, size_t)> interpolate_edges
      = [this, &qh, q, &dqr](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          const Edge & E = *mesh().edge(iE);
          QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr);
          auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_dqr_E);
          qh.segment(globalOffset(E), PolynomialSpaceDimension<Edge>::Poly(degree() - 1)) 
            = l2_projection(q, *edgeBases(iE).Polykmo, quad_dqr_E, basis_Pkmo_E_quad);
        } // for iE
      };
    parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);    
  } // if degree()

  // Interpolate at cells
  std::function<void(size_t, size_t)> interpolate_cells
    = [this, &qh, q, &dqr](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        const Cell & T = *mesh().cell(iT);
        QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr);
        auto basis_Pk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_dqr_T);
        qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree())) 
          = l2_projection(q, *cellBases(iT).Polyk, quad_dqr_T, basis_Pk_T_quad);
      } // for iT
    };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);

  return qh;
}

//------------------------------------------------------------------------------
// Edge potential reconstructions
//------------------------------------------------------------------------------

Eigen::MatrixXd XRot::_compute_edge_potential(size_t iE)
{
  const Edge & E = *mesh().edge(iE);
  
  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  MonomialEdgeIntegralsType int_mono_2kpo_E = IntegrateEdgeMonomials(E, 2*degree()+1);
  auto MGE = GramMatrix(E, *edgeBases(iE).Polyk, int_mono_2kpo_E);

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BGE
    = Eigen::MatrixXd::Zero(edgeBases(iE).Polyk->dimension(), dimensionEdge(iE));
  for (size_t i = 0; i < edgeBases(iE).Polyk->dimension(); i++) {
    BGE(i, 0) = -edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(0)->coords());
    BGE(i, 1) = edgeBases(iE).Polyk->function(i, mesh().edge(iE)->vertex(1)->coords());
  } // for i

  QuadratureRule quad_2kmo_E = generate_quadrature_rule(E, 2 * (degree() - 1));
  
  if (degree() > 0) {    
    GradientBasis<DDRCore::PolyBasisEdgeType> grad_Pk_E(*edgeBases(iE).Polyk);
    BGE.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1)) 
      = -GramMatrix(E, grad_Pk_E, *edgeBases(iE).Polykmo, int_mono_2kpo_E);
  }
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BPE
    = Eigen::MatrixXd::Zero(PolynomialSpaceDimension<Edge>::Poly(degree()) + 1, dimensionEdge(iE));

  // Enforce the gradient of the potential reconstruction
  BPE.topRows(PolynomialSpaceDimension<Edge>::Poly(degree())) = BGE;

  // Enforce the average value of the potential reconstruction
  if (degree() == 0) {
    // We set the average equal to the mean of vertex values
    BPE.bottomRows(1)(0, 0) = 0.5 * E.measure();
    BPE.bottomRows(1)(0, 1) = 0.5 * E.measure();
  } else {
    QuadratureRule quad_kmo_E = generate_quadrature_rule(E, degree() - 1);
    auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_kmo_E);
    
    // We set the average value of the potential equal to the average of the edge unknown
    for (size_t i = 0; i < PolynomialSpaceDimension<Edge>::Poly(degree() - 1); i++) {
      for (size_t iqn = 0; iqn < quad_kmo_E.size(); iqn++) {
        BPE.bottomRows(1)(0, 2 + i) += quad_kmo_E[iqn].w * basis_Pkmo_E_quad[i][iqn];
      } // for iqn
    } // for i
  }
  
  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  Eigen::MatrixXd MPE
    = Eigen::MatrixXd::Zero(PolynomialSpaceDimension<Edge>::Poly(degree()) + 1, PolynomialSpaceDimension<Edge>::Poly(degree() + 1));

  GradientBasis<DDRCore::PolyBasisEdgeType> grad_Pkpo_E(*edgeBases(iE).Polykpo);                                    
  MPE.topRows(PolynomialSpaceDimension<Edge>::Poly(degree()))
    = GramMatrix(E, *edgeBases(iE).Polyk, grad_Pkpo_E, int_mono_2kpo_E);

  MonomialScalarBasisEdge basis_P0_E(E, 0);
  MPE.bottomRows(1) = GramMatrix(E, basis_P0_E, *edgeBases(iE).Polykpo, int_mono_2kpo_E);
  
  return MPE.partialPivLu().solve(BPE);
}

//------------------------------------------------------------------------------

XRot::LocalOperators XRot::_compute_cell_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  //------------------------------------------------------------------------------
  // Rotor
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  MonomialCellIntegralsType int_mono_2k_T = IntegrateCellMonomials(T, 2 * degree());
  auto MRT = GramMatrix(T, *cellBases(iT).Polyk2, int_mono_2k_T);
  
  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BRT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk2->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    QuadratureRule quad_2kp2_E = generate_quadrature_rule(E, 2 * (degree() + 1));
    auto basis_Pk2_tE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2kp2_E), E.tangent());
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kp2_E);
    Eigen::MatrixXd PE = extendOperator(T, E, edgePotential(E));
    BRT += T.edge_orientation(iE) * compute_gram_matrix(basis_Pk2_tE_E_quad, basis_Pkpo_E_quad, quad_2kp2_E) * PE;
  } // for iE

  // Cell contribution
  if (degree() > 0) {
    RotorBasis<DDRCore::Poly2BasisCellType> rot_Pk2_T(*cellBases(iT).Polyk2);
    QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
    BRT.block(0, localOffset(T), cellBases(iT).Polyk2->dimension(), cellBases(iT).Polyk->dimension())
      += compute_gram_matrix(
                             evaluate_quad<Function>::compute(rot_Pk2_T, quad_2k_T),
                             evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2k_T),
                             quad_2k_T
                             );
  } // if degree() > 0

  Eigen::MatrixXd RT = MRT.ldlt().solve(BRT);

  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  RotorBasis<DDRCore::GolyComplBasisCellType> rot_Gckp2_T(*cellBases(iT).GolyComplkp2);
  auto quad_2kp2_T = generate_quadrature_rule(T, 2 * degree() + 2);
  Eigen::MatrixXd MPT
    = compute_gram_matrix(
                          evaluate_quad<Rotor>::compute(*cellBases(iT).GolyComplkp2, quad_2kp2_T),
                          evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kp2_T),
                          quad_2kp2_T
                          );
  
  //------------------------------------------------------------------------------
  // Right-hand side matrix

  // Cell contribution
  auto quad_2kp3_T = generate_quadrature_rule(T, 2 * degree() + 3);
  Eigen::MatrixXd BPT
    = compute_gram_matrix(
                          evaluate_quad<Function>::compute(*cellBases(iT).GolyComplkp2, quad_2kp3_T),
                          evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2kp3_T),
                          quad_2kp3_T
                          ) * RT;

  // Boundary contribution
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    QuadratureRule quad_2kp3_E = generate_quadrature_rule(E, 2 * degree() + 3);
    auto basis_Gckp2_tE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).GolyComplkp2, quad_2kp3_E), E.tangent());
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kp3_E);
    Eigen::MatrixXd PE = extendOperator(T, E, edgePotential(E));
    BPT -= T.edge_orientation(iE) * compute_gram_matrix(basis_Gckp2_tE_E_quad, basis_Pkpo_E_quad, quad_2kp3_E) * PE;
  } // for iE

  return LocalOperators(RT, BRT, MPT.partialPivLu().solve(BPT));
}

//------------------------------------------------------------------------------
// Compute rotor L2-product

Eigen::MatrixXd XRot::computeRotorL2Product(
                                            const size_t iT,
                                            const double & penalty_factor_cell,
                                            const double & penalty_factor_edge
                                            ) const
{
  const Cell & T = *mesh().cell(iT);

  double hT = T.diam();
  
  Eigen::MatrixXd PT = cellOperators(iT).potential;

  Eigen::MatrixXd RL2P
    = cellOperators(iT).rotor.transpose() * cellOperators(iT).rotor_rhs;

  // Penalisation of the element unknown, which is one degree higher with respect
  // to XGrad
  { 
    MonomialCellIntegralsType int_mono_2kpo_T = IntegrateCellMonomials(T, 2 * degree() + 1);
    Eigen::MatrixXd MTk = GramMatrix(T, *cellBases(iT).Polyk, int_mono_2kpo_T);
    Eigen::MatrixXd gram_Pk_Pkpo = GramMatrix(T, *cellBases(iT).Polyk, *cellBases(iT).Polykpo, int_mono_2kpo_T);
    size_t dim_Pk = cellBases(iT).Polyk->dimension();

    // Projection on Pk(T) of the element potential
    Eigen::MatrixXd diff_T = MTk.ldlt().solve(gram_Pk_Pkpo * PT);
    // Subtract the element component    
    diff_T.block(0, localOffset(T), dim_Pk, dim_Pk) -= Eigen::MatrixXd::Identity(dim_Pk, dim_Pk);
   
    RL2P += ( penalty_factor_cell / std::pow(hT, 2) ) * diff_T.transpose() * MTk * diff_T;
  }

  // Edge penalty terms
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);

    QuadratureRule quad_2kpt_E = generate_quadrature_rule(E, 2 * (degree() + 1));

    double hE = E.measure();

    auto basis_Pkpo_T_quad_E = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpt_E);
    auto basis_Pkpo_E_quad_E = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kpt_E);
    
    Eigen::MatrixXd gram_PkpoT_PkpoE = compute_gram_matrix(basis_Pkpo_T_quad_E, basis_Pkpo_E_quad_E, quad_2kpt_E);

    Eigen::MatrixXd PE = extendOperator(T, E, edgePotential(E));

    RL2P += ( penalty_factor_edge / hE ) * (
                                            PT.transpose() * compute_gram_matrix(basis_Pkpo_T_quad_E, quad_2kpt_E) * PT
                                            - PT.transpose() * gram_PkpoT_PkpoE * PE
                                            - PE.transpose() * gram_PkpoT_PkpoE.transpose() * PT
                                            + PE.transpose() * compute_gram_matrix(basis_Pkpo_E_quad_E, quad_2kpt_E) * PE
                                            );

  } // for iE

  return RL2P;
}

//------------------------------------------------------------------------------
// Compute rotor L2-norm

double XRot::computeRotorL2Norm(const Eigen::VectorXd & v) const
{
  Eigen::VectorXd squared_rotor_l2_norm = Eigen::VectorXd::Zero(m_mesh.n_cells());

  std::function<void(size_t, size_t)> compute_squared_rotor_l2_norm_all
    = [this, &v, &squared_rotor_l2_norm](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        Eigen::VectorXd vT = restrictCell(iT, v);
        squared_rotor_l2_norm(iT) = vT.transpose() * computeRotorL2Product(iT) * vT;
      } // for iT
    };

  parallel_for(m_ddr_core.mesh().n_cells(), compute_squared_rotor_l2_norm_all, m_use_threads);

  return std::sqrt(squared_rotor_l2_norm.sum());
}
