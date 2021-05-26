
#include <xgrad.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XGrad::XGrad(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : DDRSpace(ddr_core.mesh(),
	     1,
	     PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree() - 1),
	     PolynomialSpaceDimension<Cell>::Poly(ddr_core.degree() - 1)
	     ),
    m_ddr_core(ddr_core),
    m_use_threads(use_threads),
    m_output(output),
    m_edge_operators(ddr_core.mesh().n_edges()),
    m_cell_operators(ddr_core.mesh().n_cells())
{
  m_output << "[XGrad] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XGrad] Parallel execution" << std::endl;
  } else {
    m_output << "[XGrad] Sequential execution" << std::endl;
  }
  
  // Construct edge gradients and potentials
  std::function<void(size_t, size_t)> construct_all_edge_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          m_edge_operators[iE].reset( new LocalOperators(_compute_edge_gradient_potential(iE)) );
        } // for iE
      };

  m_output << "[XGrad] Constructing edge gradients and potentials" << std::endl;
  parallel_for(mesh().n_edges(), construct_all_edge_gradients_potentials, use_threads);

  // Construct cell gradients and potentials
  std::function<void(size_t, size_t)> construct_all_cell_gradients_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_gradient_potential(iT)) );
        } // for iT
      };

  m_output << "[XGrad] Constructing cell gradients and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_gradients_potentials, use_threads);
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XGrad::interpolate(const FunctionType & q, const int deg_quad) const
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
    
    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &qh, q, &dqr](size_t start, size_t end)->void
        {
          for (size_t iT = start; iT < end; iT++) {
            const Cell & T = *mesh().cell(iT);
            QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr);
            auto basis_Pkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykmo, quad_dqr_T);
            qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree() - 1)) 
              = l2_projection(q, *cellBases(iT).Polykmo, quad_dqr_T, basis_Pkmo_T_quad);
          } // for iT
        };
    parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  } // if degree() > 0 

  return qh;
}

//------------------------------------------------------------------------------
// Gradient and potential reconstructions
//------------------------------------------------------------------------------

XGrad::LocalOperators XGrad::_compute_edge_gradient_potential(size_t iE)
{
  const Edge & E = *mesh().edge(iE);
  
  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
  auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_2k_E);
  auto MGE = compute_gram_matrix(basis_Pk_E_quad, basis_Pk_E_quad, quad_2k_E, "sym");

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
    auto grad_Pk_tE_E_quad = scalar_product(
                                            evaluate_quad<Gradient>::compute(*edgeBases(iE).Polyk, quad_2kmo_E),
                                            E.tangent()
                                            );
    auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_2kmo_E);
    BGE.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1))
      = -compute_gram_matrix(grad_Pk_tE_E_quad, basis_Pkmo_E_quad, quad_2kmo_E);
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
  auto grad_Pkpo_tE_E_quad = scalar_product(
					    evaluate_quad<Gradient>::compute(*edgeBases(iE).Polykpo, quad_2k_E),
					    E.tangent()
					    );
  MPE.topRows(PolynomialSpaceDimension<Edge>::Poly(degree()))
    = compute_gram_matrix(basis_Pk_E_quad, grad_Pkpo_tE_E_quad, quad_2k_E);
  
  QuadratureRule quad_kpo_E = generate_quadrature_rule(E, degree() + 1);  
  auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykpo, quad_kpo_E);  
  for (size_t i = 0; i < PolynomialSpaceDimension<Edge>::Poly(degree() + 1); i++) {
    for (size_t iqn = 0; iqn < quad_kpo_E.size(); iqn++) {
      MPE.bottomRows(1)(0, i) += quad_kpo_E[iqn].w * basis_Pkpo_E_quad[i][iqn];
    } // for iqn
  } // for i
  
  return LocalOperators(MGE.ldlt().solve(BGE), MPE.partialPivLu().solve(BPE));
}

//------------------------------------------------------------------------------

XGrad::LocalOperators XGrad::_compute_cell_gradient_potential(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
  auto basis_Pk2_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2k_T);
  auto MGT = compute_gram_matrix(basis_Pk2_T_quad, basis_Pk2_T_quad, quad_2k_T, "sym");

  //------------------------------------------------------------------------------
  // Right-hand side matrix
  
  Eigen::MatrixXd BGT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk2->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * (degree() + 1));
    auto basis_Pk2_nTE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2kpo_E), T.edge_normal(iE));
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kpo_E);
    Eigen::MatrixXd BGT_E
      = compute_gram_matrix(basis_Pk2_nTE_E_quad, basis_Pkpo_E_quad, quad_2kpo_E) * edgeOperators(E).potential;

    // Assemble local contribution
    BGT.col(localOffset(T, *E.vertex(0))) += BGT_E.col(0);
    BGT.col(localOffset(T, *E.vertex(1))) += BGT_E.col(1);
    if (degree() > 0) {
      BGT.block(0, localOffset(T, E), cellBases(iT).Polyk2->dimension(), PolynomialSpaceDimension<Edge>::Poly(degree() - 1))
        = BGT_E.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1));
    } // if degree() > 0
  } // for iE

  // Cell contribution
  if (degree() > 0) {
    QuadratureRule quad_2kmo_T = generate_quadrature_rule(T, 2 * (degree() - 1));
    auto div_Pk2_T_quad = evaluate_quad<Divergence>::compute(*cellBases(iT).Polyk2, quad_2kmo_T);
    auto basis_Pkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykmo, quad_2kmo_T);

    BGT.rightCols(PolynomialSpaceDimension<Cell>::Poly(degree() - 1))
      -= compute_gram_matrix(div_Pk2_T_quad, basis_Pkmo_T_quad, quad_2kmo_T);
  } // if degree() > 0

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix
  
  QuadratureRule quad_2kpo_T = generate_quadrature_rule(T, 2 * (degree() + 1));
  auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_T);
  auto div_Rckp2_T_quad = evaluate_quad<Divergence>::compute(*cellBases(iT).RolyComplkp2, quad_2kpo_T);
  Eigen::MatrixXd MPT = compute_gram_matrix(div_Rckp2_T_quad, basis_Pkpo_T_quad, quad_2kpo_T);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  // Cell contribution
  Eigen::MatrixXd BPT
    = -compute_gram_matrix(
			   evaluate_quad<Function>::compute(*cellBases(iT).RolyComplkp2, quad_2kpo_T),
			   evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2kpo_T),
			   quad_2kpo_T
			   ) * GT;
  
  // Boundary contribution
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    QuadratureRule quad_2kp2_E = generate_quadrature_rule(E, 2 * (degree() + 2));
    auto basis_Rckp2_nTE_E_quad
      = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).RolyComplkp2, quad_2kp2_E), T.edge_normal(iE));
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykpo, quad_2kp2_E);
    Eigen::MatrixXd BPT_E
      = compute_gram_matrix(basis_Rckp2_nTE_E_quad, basis_Pkpo_E_quad, quad_2kp2_E) * edgeOperators(E).potential;

    // Assemble local contribution
    BPT.col(localOffset(T, *E.vertex(0))) += BPT_E.col(0);
    BPT.col(localOffset(T, *E.vertex(1))) += BPT_E.col(1);
    if (degree() > 0) {
      BPT.block(0, localOffset(T, E), cellBases(iT).RolyComplkp2->dimension(), PolynomialSpaceDimension<Edge>::Poly(degree() - 1))
        += BPT_E.rightCols(PolynomialSpaceDimension<Edge>::Poly(degree() - 1));
    } // if degree() > 0
  } // for iE
  
  return LocalOperators(GT, MPT.partialPivLu().solve(BPT));
}


//-----------------------------------------------------------------------------
// local L2 inner product
//-----------------------------------------------------------------------------
Eigen::MatrixXd XGrad::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pkpo_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pkpo_T;
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pkpo_T.rows()==1){
      // We have to compute the mass matrix
      QuadratureRule quad_2kpo_T = generate_quadrature_rule(T, 2 * (degree()+1));
      w_mass_Pkpo_T = weight.value(T, T.center_mass()) * compute_gram_matrix(evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_T), quad_2kpo_T);
    }else{
      w_mass_Pkpo_T = weight.value(T, T.center_mass()) * mass_Pkpo_T;
    }
  }else{
    // weight is not constant, we create a weighted mass matrix
    QuadratureRule quad_2kpo_pw_T = generate_quadrature_rule(T, 2 * (degree() + 1) + weight.deg(T));
    auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_pw_T);
    std::function<double(const Eigen::Vector2d &)> weight_T 
              = [&T, &weight](const Eigen::Vector2d &x)->double {
                  return weight.value(T, x);
                };
    w_mass_Pkpo_T = compute_weighted_gram_matrix(weight_T, basis_Pkpo_T_quad, basis_Pkpo_T_quad, quad_2kpo_pw_T, "sym");
  }

  // Compute matrix of L2 product  
  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));

  // We need the potential in the cell
  Eigen::MatrixXd Potential_T = cellOperators(iT).potential;

  // Edge penalty terms
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
        
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * (degree()+1) );
    
    // weight and scaling hE
    double max_weight_quad_E = weight.value(T, quad_2kpo_E[0].vector());
    // If the weight is not constant, we want to take the largest along the edge
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2kpo_E.size(); iqn++) {
        max_weight_quad_E = std::max(max_weight_quad_E, weight.value(T, quad_2kpo_E[iqn].vector()));
      } // for
    }
    double w_hE = max_weight_quad_E * E.measure();

    // The penalty term int_E (PT q - q_E) * (PT r - r_E) is computed by developping.
    auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, quad_2kpo_E);
    auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polykpo, quad_2kpo_E);
    Eigen::MatrixXd gram_PkpoT_PkpoE = compute_gram_matrix(basis_Pkpo_T_quad, basis_Pkpo_E_quad, quad_2kpo_E);
    
    Eigen::MatrixXd Potential_E = extendOperator(T, E, edgeOperators(E).potential);

    // Contribution of edge E
    L2P += w_hE * ( Potential_T.transpose() * compute_gram_matrix(basis_Pkpo_T_quad, quad_2kpo_E) * Potential_T
                   - Potential_T.transpose() * gram_PkpoT_PkpoE * Potential_E
                   - Potential_E.transpose() * gram_PkpoT_PkpoE.transpose() * Potential_T
                   + Potential_E.transpose() * compute_gram_matrix(basis_Pkpo_E_quad, quad_2kpo_E) * Potential_E );
  } // for iE

  L2P *= penalty_factor;

  // Cell term
  L2P += Potential_T.transpose() * w_mass_Pkpo_T * Potential_T;

  return L2P;

}


//-----------------------------------------------------------------------------
// Evaluate the potential at a point
//-----------------------------------------------------------------------------
double XGrad::evaluatePotential(const size_t iT, const Eigen::VectorXd & vT, const VectorRd & x) const {
  // Ancestor of basis of P^{k+1}(T) and values at x
  MonomialScalarBasisCell monomial_ancestor = cellBases(iT).Polykpo->ancestor();
  Eigen::VectorXd anc_values = Eigen::VectorXd::Zero(monomial_ancestor.dimension());
  for (size_t i=0; i<monomial_ancestor.dimension(); i++){
    anc_values(i) = monomial_ancestor.function(i, x);  
  }

  // return value
  return anc_values.transpose() * (cellBases(iT).Polykpo->matrix()).transpose() * cellOperators(iT).potential * vT;

}


