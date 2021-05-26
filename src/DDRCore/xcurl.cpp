#include <xcurl.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XCurl::XCurl(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : DDRSpace(
	     ddr_core.mesh(),
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree()),
	     PolynomialSpaceDimension<Cell>::Roly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Cell>::RolyCompl(ddr_core.degree())	     
	     ),
    m_ddr_core(ddr_core),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_operators(ddr_core.mesh().n_cells())
{
  output << "[XCurl] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XCurl] Parallel execution" << std::endl;
  } else {
    m_output << "[XCurl] Sequential execution" << std::endl;
  }

  // Construct cell curls and potentials
  std::function<void(size_t, size_t)> construct_all_cell_curls_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_curl_potential(iT)) );
        } // for iT
      };

  m_output << "[XCurl] Constructing cell curls and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_curls_potentials, use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XCurl::interpolate(const FunctionType & v, const int deg_quad) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());
  
  // Degree of quadrature rules
  size_t dqr = (deg_quad >= 0 ? deg_quad : 2 * degree() + 3);
  
  // Interpolate at edges
  std::function<void(size_t, size_t)> interpolate_edges
    = [this, &vh, v, &dqr](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        const Edge & E = *mesh().edge(iE);

	        Eigen::Vector2d tE = E.tangent();
	        auto v_dot_tE = [&tE, v](const Eigen::Vector2d & x)->double {
			          return v(x).dot(tE);
			        };

	        QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr );
	        auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_dqr_E);
	        vh.segment(globalOffset(E), edgeBases(iE).Polyk->dimension())
	          = l2_projection(v_dot_tE, *edgeBases(iE).Polyk, quad_dqr_E, basis_Pk_E_quad);
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
// Curl and potential reconstruction
//------------------------------------------------------------------------------

XCurl::LocalOperators XCurl::_compute_cell_curl_potential(size_t iT)
{
  const Cell & T = *mesh().cell(iT);
  
  //------------------------------------------------------------------------------
  // Curl
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());

  auto basis_Pk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2k_T);
  Eigen::MatrixXd MCT = compute_gram_matrix(basis_Pk_T_quad, quad_2k_T);
  
  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BCT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk->dimension(), dimensionCell(iT));

  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    BCT.block(0, localOffset(T, E), cellBases(iT).Polyk->dimension(), edgeBases(E.global_index()).Polyk->dimension())
      -= T.edge_orientation(iE) * compute_gram_matrix(
						      evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2k_E),
						      evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2k_E),
						      quad_2k_E
						      );
  } // for iE

  if (degree() > 0) {
    QuadratureRule quad_2kmo_T = generate_quadrature_rule(T, 2 * (degree() - 1));

    BCT.block(0, localOffset(T), cellBases(iT).Polyk->dimension(), cellBases(iT).Rolykmo->dimension())
      += compute_gram_matrix(
			     evaluate_quad<Curl>::compute(*cellBases(iT).Polyk, quad_2kmo_T),
			     evaluate_quad<Function>::compute(*cellBases(iT).Rolykmo, quad_2kmo_T),
			     quad_2kmo_T
			     );
  } // if degree() > 0
 
  Eigen::MatrixXd CT = MCT.ldlt().solve(BCT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  auto basis_Pkpo0_T = ShiftedBasis<typename DDRCore::PolyBasisCellType>(*cellBases(iT).Polykpo, 1);

  Eigen::MatrixXd MPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk2->dimension(), cellBases(iT).Polyk2->dimension());
  Eigen::MatrixXd BPT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polyk2->dimension(), dimensionCell(iT));

  auto basis_Pk2_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2k_T);
  MPT.topLeftCorner(basis_Pkpo0_T.dimension(), cellBases(iT).Polyk2->dimension())
    = compute_gram_matrix(
			  evaluate_quad<Curl>::compute(basis_Pkpo0_T, quad_2k_T),
			  basis_Pk2_T_quad, quad_2k_T
			  );

  if (degree() > 0) {
    auto basis_Rck_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).RolyComplk, quad_2k_T);
    MPT.bottomLeftCorner(cellBases(iT).RolyComplk->dimension(), cellBases(iT).Polyk2->dimension())
      = compute_gram_matrix(basis_Rck_T_quad, basis_Pk2_T_quad, quad_2k_T);
    BPT.bottomRightCorner(cellBases(iT).RolyComplk->dimension(), cellBases(iT).RolyComplk->dimension())
      += compute_gram_matrix(basis_Rck_T_quad, quad_2k_T);    
  } // if degree() > 0
 
  auto quad_2kpo_T = generate_quadrature_rule(T, 2 * degree() + 1);
  BPT.topLeftCorner(basis_Pkpo0_T.dimension(), dimensionCell(iT))
    += compute_gram_matrix(
			   evaluate_quad<Function>::compute(basis_Pkpo0_T, quad_2kpo_T),
			   evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_2kpo_T),
			   quad_2kpo_T
			   ) * CT;
 
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * degree() + 1);
    BPT.block(0, localOffset(T, E), basis_Pkpo0_T.dimension(), edgeBases(E.global_index()).Polyk->dimension())
      += T.edge_orientation(iE) * compute_gram_matrix(
						      evaluate_quad<Function>::compute(basis_Pkpo0_T, quad_2kpo_E),
						      evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2kpo_E),
						      quad_2kpo_E
						      );    
  } // for iE

  return LocalOperators(CT, MPT.partialPivLu().solve(BPT));
}


//------------------------------------------------------------------------------
//        Functions to compute matrices for local L2 products on Xcurl
//------------------------------------------------------------------------------

Eigen::MatrixXd XCurl::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk2_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pk2_T;
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pk2_T.rows()==1){
      // We have to compute the mass matrix
      QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
      w_mass_Pk2_T = weight.value(T, T.center_mass()) * compute_gram_matrix(evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2k_T), quad_2k_T);
    }else{
      w_mass_Pk2_T = weight.value(T, T.center_mass()) * mass_Pk2_T;
    }
  }else{
    // weight is not constant, we create a weighted mass matrix
    QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
    auto basis_Pk2_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2kpw_T);
    std::function<double(const Eigen::Vector2d &)> weight_T 
              = [&T, &weight](const Eigen::Vector2d &x)->double {
                  return weight.value(T, x);
                };
    w_mass_Pk2_T = compute_weighted_gram_matrix(weight_T, basis_Pk2_T_quad, basis_Pk2_T_quad, quad_2kpw_T, "sym");
  }

  
  // The leftOp and rightOp will come from the potentials
  std::vector<Eigen::MatrixXd> potentialOp(T.n_edges()+1);
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge & E = *T.edge(iE);
    potentialOp[iE] = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E),dimensionEdge(E)));
  }
  potentialOp[T.n_edges()] = m_cell_operators[iT]->potential;

  return computeL2Product_with_Ops(iT, potentialOp, potentialOp, penalty_factor, w_mass_Pk2_T, weight);

}

Eigen::MatrixXd XCurl::computeL2ProductGradient(
                                        const size_t iT,
                                        const XGrad & x_grad,
                                        const std::string & side,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk2_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT);

  // create the weighted mass matrix, with simple product if weight is constant
  Eigen::MatrixXd w_mass_Pk2_T;
  if (weight.deg(T)==0){
    // constant weight
    if (mass_Pk2_T.rows()==1){
      // We have to compute the mass matrix
      QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * degree());
      w_mass_Pk2_T = weight.value(T, T.center_mass()) * compute_gram_matrix(evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2k_T), quad_2k_T);
    }else{
      w_mass_Pk2_T = weight.value(T, T.center_mass()) * mass_Pk2_T;
    }
  }else{
    // weight is not constant, we create a weighted mass matrix
    QuadratureRule quad_2kpw_T = generate_quadrature_rule(T, 2 * degree() + weight.deg(T));
    auto basis_Pk2_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2kpw_T);
    std::function<double(const Eigen::Vector2d &)> weight_T 
              = [&T, &weight](const Eigen::Vector2d &x)->double {
                  return weight.value(T, x);
                };
    w_mass_Pk2_T = compute_weighted_gram_matrix(weight_T, basis_Pk2_T_quad, basis_Pk2_T_quad, quad_2kpw_T, "sym");
  }


  // list of full gradients
  std::vector<Eigen::MatrixXd> gradientOp(T.n_edges()+1);
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge & E = *T.edge(iE);
    gradientOp[iE] = x_grad.extendOperator(T, E, x_grad.edgeOperators(E).gradient);
  }
  gradientOp[T.n_edges()] = x_grad.cellOperators(iT).gradient;
  
  // If we apply the gradient on one side only we'll need the potentials
  if (side != "both"){
    // list of potentials
    std::vector<Eigen::MatrixXd> potentialOp(T.n_edges()+1);
    for (size_t iE = 0; iE < T.n_edges(); iE++){
      const Edge & E = *T.edge(iE);
      potentialOp[iE] = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E),dimensionEdge(E)));
    }
    potentialOp[T.n_edges()] = m_cell_operators[iT]->potential;
  
    // Depending on side of gradient
    if (side == "left"){
      return computeL2Product_with_Ops(iT, gradientOp, potentialOp, penalty_factor, w_mass_Pk2_T, weight);
    }else{
      return computeL2Product_with_Ops(iT, potentialOp, gradientOp, penalty_factor, w_mass_Pk2_T, weight);
    }
    
  }

  // Default: gradient on both sides
  return computeL2Product_with_Ops(iT, gradientOp, gradientOp, penalty_factor, w_mass_Pk2_T, weight);

}


Eigen::MatrixXd XCurl::computeL2Product_with_Ops(
                                        const size_t iT,
                                        const std::vector<Eigen::MatrixXd> & leftOp,
                                        const std::vector<Eigen::MatrixXd> & rightOp,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & w_mass_Pk2_T,
                                        const IntegralWeight & weight
                                        ) const
{
  const Cell & T = *mesh().cell(iT); 

  // leftOp and rightOp must list the operators acting on the DOFs, and which we want to
  // use for the L2 product. Specifically, each one lists operators (matrices) returning
  // values in edges space P^k(E), faces space P^k(F)^2 (tangent) and element space P^k(T)^3.
  // For the standard Xcurl L2 product, these will respectively be identity (for each edge),
  // gamma_tF (for each F) and PT. 
  // To compute the Xcurl L2 product applied (left or right) to the discrete gradient,
  // leftOp or rightOp must list the edge, face and element (full) gradient operators.
  // All these operators must have the same domain, so possibly being extended appropriately
  // using extendOperator from ddrspace.

  Eigen::MatrixXd L2P = Eigen::MatrixXd::Zero(leftOp[0].cols(), rightOp[0].cols());
  
  size_t offset_T = T.n_edges();

  // Edge penalty terms
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    VectorRd tE = E.tangent();
        
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    
    // weight and scaling hE
    double max_weight_quad_E = weight.value(T, quad_2k_E[0].vector());
    // If the weight is not constant, we want to take the largest along the edge
    if (weight.deg(T)>0){
      for (size_t iqn = 1; iqn < quad_2k_E.size(); iqn++) {
        max_weight_quad_E = std::max(max_weight_quad_E, weight.value(T, quad_2k_E[iqn].vector()));
      } // for
    }
    double w_hE = max_weight_quad_E * E.measure();

    // The penalty term int_E (PT w . tE - w_E) * (PT v . tE - v_E) is computed by developping.
    auto basis_Pk2_T_dot_tE_quad = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, quad_2k_E), tE);
    auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(E.global_index()).Polyk, quad_2k_E);
    Eigen::MatrixXd gram_Pk2T_dot_tE_PkE = compute_gram_matrix(basis_Pk2_T_dot_tE_quad, basis_Pk_E_quad, quad_2k_E);
    
    // Contribution of edge E
    L2P += w_hE * ( leftOp[offset_T].transpose() * compute_gram_matrix(basis_Pk2_T_dot_tE_quad, quad_2k_E) * rightOp[offset_T]
                   - leftOp[offset_T].transpose() * gram_Pk2T_dot_tE_PkE * rightOp[iE]
                   - leftOp[iE].transpose() * gram_Pk2T_dot_tE_PkE.transpose() * rightOp[offset_T]
                   + leftOp[iE].transpose() * compute_gram_matrix(basis_Pk_E_quad, quad_2k_E) * rightOp[iE]);

  } // for iE

  L2P *= penalty_factor;
  
  // Consistent (cell) term
  L2P += leftOp[offset_T].transpose() * w_mass_Pk2_T * rightOp[offset_T];
  
   
  return L2P;
}


