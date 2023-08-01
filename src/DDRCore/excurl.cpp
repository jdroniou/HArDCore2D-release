#include "excurl.hpp"
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

EXCurl::EXCurl(const DDRCore & ddr_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(
	     ddr_core.mesh(),
	     0,
	     2*PolynomialSpaceDimension<Edge>::Poly(ddr_core.degree()),
	     PolynomialSpaceDimension<Cell>::Roly(ddr_core.degree() - 1) + PolynomialSpaceDimension<Cell>::RolyCompl(ddr_core.degree())	     
	     ),
    m_ddr_core(ddr_core),
    m_xcurl(XCurl(ddr_core, use_threads, output)),
    m_Polyk2x2(ddr_core.mesh().n_cells()),
    m_Polykpo2(ddr_core.mesh().n_cells()),
    m_reduction(ddr_core.mesh().n_cells()),
    m_use_threads(use_threads),
    m_output(output),
    m_hho_operators(ddr_core.mesh().n_cells())
{
  // Construct reduction matrices
  for (size_t iT = 0; iT < mesh().n_cells(); iT++){
    const Cell & T = *mesh().cell(iT);
    m_reduction[iT] = Eigen::MatrixXd::Zero(m_xcurl.dimensionCell(iT), dimensionCell(iT));
    for (size_t iE = 0; iE < T.n_edges(); iE++){
      size_t dofsE = edgeBases(iE).Polyk->dimension();
      m_reduction[iT].block(iE*dofsE, (2*iE+1)*dofsE, dofsE, dofsE) = 
            Eigen::MatrixXd::Identity(dofsE, dofsE);
    }
    m_reduction[iT].bottomRightCorner(numLocalDofsCell(), numLocalDofsCell()) = 
            Eigen::MatrixXd::Identity(numLocalDofsCell(), numLocalDofsCell());
  }
  
  // Construct Matricial Polyk basis and hho operators
  std::function<void(size_t, size_t)> construct_all_hho_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_Polyk2x2[iT].reset( new Polyk2x2Type(*(m_ddr_core.cellBases(iT).Polyk)) );
          m_Polykpo2[iT].reset( new Polykpo2Type(*(m_ddr_core.cellBases(iT).Polykpo)) );
          m_hho_operators[iT].reset( new hhoLocalOperators(_compute_hho_operators(iT)) );
        } // for iT
      };

  m_output << "[EXCurl] Constructing hho operators" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_hho_operators, use_threads);  

}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd EXCurl::interpolate(const FunctionType & v, const int deg_quad) const
{
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(dimension());
  
  // Degree of quadrature rules
  size_t dqr = (deg_quad >= 0 ? deg_quad : 2 * degree() + 3);
  
  // Interpolant for tangential components on edges, and cell components
  Eigen::VectorXd vh_xcurl = m_xcurl.interpolate(v, deg_quad);
  
  // Interpolate normal component to the edges, and plug in tangential components
  std::function<void(size_t, size_t)> interpolate_normal_edges
    = [this, &vh, &vh_xcurl, v, &dqr](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        const Edge & E = *mesh().edge(iE);

	        Eigen::Vector2d nE = E.normal();
	        auto v_dot_nE = [&nE, v](const Eigen::Vector2d & x)->double {
			          return v(x).dot(nE);
			        };

	        QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr);
	        auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_dqr_E);
	        vh.segment(globalOffset(E), edgeBases(iE).Polyk->dimension())
	          = l2_projection(v_dot_nE, *edgeBases(iE).Polyk, quad_dqr_E, basis_Pk_E_quad);
	        vh.segment(globalOffset(E) + edgeBases(iE).Polyk->dimension(), edgeBases(iE).Polyk->dimension())
	          = vh_xcurl.segment(m_xcurl.globalOffset(E), edgeBases(iE).Polyk->dimension());
	      } // for iE
      };
  parallel_for(mesh().n_edges(), interpolate_normal_edges, m_use_threads);
  
  // Interpolate cell component
  for (size_t iT = 0; iT < mesh().n_cells(); iT++){
    const Cell & T = *mesh().cell(iT);
    vh.segment(globalOffset(T), numLocalDofsCell() ) 
          = vh_xcurl.segment(m_xcurl.globalOffset(T), m_xcurl.numLocalDofsCell() );
  }
  
  return vh;
}

//------------------------------------------------------------------------------
// HHO operators
//------------------------------------------------------------------------------

EXCurl::hhoLocalOperators EXCurl::_compute_hho_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);
  auto qr_2k_T = generate_quadrature_rule(T, 2*degree());
  boost::multi_array<double, 2> basis_PkT_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, qr_2k_T);
  
  //
  // GRADIENT
  //
  // Left-hand side
  Eigen::MatrixXd MGT = compute_gram_matrix(Polyk2x2(iT), basis_PkT_quad, qr_2k_T);

  // Right-hand side, starting with volumetric contribution
  Eigen::MatrixXd BGT = 
        - compute_gram_matrix(evaluate_quad<Divergence>::compute(Polyk2x2(iT), qr_2k_T), basis_PkT_quad, qr_2k_T) * ddrPotential(T);

  // edge contributions
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge & E = *T.edge(iE);
    VectorRd nE = E.normal();
    VectorRd tE = E.tangent();
    auto qr_2k_E = generate_quadrature_rule(E, 2*degree());
    size_t dimPkE = edgeBases(E).Polyk->dimension();
    
    // Basis of Polyk2x2 and Polyk(E) at quadrature nodes
    boost::multi_array<MatrixRd, 2> basis_Polyk2x2_quad = evaluate_quad<Function>::compute(Polyk2x2(iT), qr_2k_E);
    boost::multi_array<double, 2> basis_PolykE_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polyk, qr_2k_E);

    // normal component: compute nE.(<basis>nE) by <basis>:(nE x nE)
    MatrixRd nE_tens_nE = nE*(nE.transpose());
    boost::multi_array<double, 2> basis_Polyk2x2_nEnE_quad = scalar_product(basis_Polyk2x2_quad, nE_tens_nE);
    BGT.block(0, localOffset(T, E),  Polyk2x2(iT).dimension(), dimPkE) 
        += T.edge_orientation(iE) * compute_gram_matrix(basis_Polyk2x2_nEnE_quad, basis_PolykE_quad, qr_2k_E);
  
    // tangential component: compute tE.(<basis>nE) by <basis>:(tE x nE)
    MatrixRd tE_tens_nE = tE*(nE.transpose());
    boost::multi_array<double, 2> basis_Polyk2x2_tEnE_quad = scalar_product(basis_Polyk2x2_quad, tE_tens_nE);
    BGT.block(0, localOffset(T, E) + dimPkE, Polyk2x2(iT).dimension(), dimPkE) 
        += T.edge_orientation(iE) * compute_gram_matrix(basis_Polyk2x2_tEnE_quad, basis_PolykE_quad, qr_2k_E);
  }

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);
  Eigen::MatrixXd GsT = Polyk2x2(iT).symmetriseOperator() * GT;
  
  //
  // POTENTIAL
  //
  // left-hand side, starting with stiffness
  boost::multi_array<MatrixRd, 2> grad_Pkpo2_T_quad = evaluate_quad<Gradient>::compute(Polykpo2(iT), qr_2k_T);
  boost::multi_array<MatrixRd, 2> gradS_Pkpo2_T_quad = transform_values_quad<MatrixRd>(grad_Pkpo2_T_quad, symmetrise_matrix);
  Eigen::MatrixXd MPT = compute_gram_matrix(gradS_Pkpo2_T_quad, qr_2k_T);
  // closure equations (only average is scaled): integral of potential (only over cell for k>=1 here, we'll do the integral
  // over the boundary when k=0 below), and integral of skew-symmetric gradient (computed at a higher degree than 
  // required because useful later).
  // The scaling parameter for the closure equation depends if it's an average in the cell or on the boundary
  double closure_average_scaling = 1.;
  auto qr_2kpo_T = generate_quadrature_rule(T, 2*degree()+1);
  boost::multi_array<VectorRd, 2> basis_Polykpo2_quad = evaluate_quad<Function>::compute(Polykpo2(iT), qr_2kpo_T);
  if (degree()>0){
    // closure for average over cell
    closure_average_scaling = 1./ (std::pow(T.diam(),2) * T.measure());
    MPT +=  closure_average_scaling * compute_closure_matrix(basis_Polykpo2_quad, qr_2kpo_T);
  }
  boost::multi_array<MatrixRd, 2> gradSS_Pkpo2_T_quad = transform_values_quad<MatrixRd>(grad_Pkpo2_T_quad, skew_symmetrise_matrix);
  MPT += compute_closure_matrix(gradSS_Pkpo2_T_quad, qr_2k_T);
  // If k=0, the closure equation involves integrals of basis of P^{k+1}(T)^2 over the boundary of T
  std::vector<VectorRd> int_partialT_Pkpo2(Polykpo2(iT).dimension(), VectorRd::Zero());
  if (degree()==0){
    closure_average_scaling = 1. / T.measure();
    for (Edge *E : T.get_edges()){
      auto qr_1_E = generate_quadrature_rule(*E, 1);
      for (size_t i=0; i<Polykpo2(iT).dimension(); i++){
        for (size_t iqn=0; iqn<qr_1_E.size(); iqn++){
          int_partialT_Pkpo2[i] += qr_1_E[iqn].w * Polykpo2(iT).function(i, qr_1_E[iqn].vector());
        }
      }
    }

    // Closure equation
    for (size_t i=0; i<Polykpo2(iT).dimension(); i++){
      for (size_t j=0; j<Polykpo2(iT).dimension(); j++){
        MPT(i, j) += closure_average_scaling * scalar_product(int_partialT_Pkpo2[i], int_partialT_Pkpo2[j]);
      }
    }
  }
  
  // right-hand side, starting with projection of GsT
  boost::multi_array<MatrixRd, 2> gradS_Polyk2_quad = 
      transform_values_quad<MatrixRd>(evaluate_quad<Gradient>::compute(*cellBases(iT).Polyk2, qr_2k_T), symmetrise_matrix);
  Eigen::MatrixXd BPT = compute_gram_matrix(gradS_Pkpo2_T_quad, evaluate_quad<Function>::compute(Polyk2x2(iT), qr_2k_T), qr_2k_T) * GsT;
  // Closure: integral of cell value (we compute with a higher degree than required because the values will be useful later)
  boost::multi_array<VectorRd, 2> basis_Polyk2_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, qr_2kpo_T);
  if (degree()>0){
    BPT += closure_average_scaling * compute_closure_matrix(basis_Polykpo2_quad, basis_Polyk2_quad, qr_2kpo_T) * ddrPotential(iT);
  }
  // Closure for symmetric gradient: contribution of edges
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge &E = *T.edge(iE);
    auto qr_k_E = generate_quadrature_rule(E, degree());
    boost::multi_array<double, 2> basis_Polyk_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polyk, qr_k_E);
    size_t dim_PkE = edgeBases(E).Polyk->dimension();

    // In the term vF x nTF - nTF x vF, only the tangential components of vF count so we compute 0.5*(tE x nTE - nTE x tE)
    MatrixRd SS_tExnTE = skew_symmetrise_matrix(E.tangent() * T.edge_normal(iE).transpose());

    // Take the scalar product of 0.5*(tE x nTE - nTE x tE) with nabla_ss w on quadrature nodes
    boost::multi_array<double, 2> SS_tExnTE_gradSS_Pkpo2_T_quad = 
          transform_values_quad<double>(gradSS_Pkpo2_T_quad, [&SS_tExnTE](const MatrixRd &x)->double { return scalar_product(x,SS_tExnTE);});

    // Compute contribution of edge to left-hand side (closure of nabla_ss pT)
    BPT.block(0, localOffset(T, E) + dim_PkE, Polykpo2(iT).dimension(), dim_PkE) 
        += compute_closure_matrix(SS_tExnTE_gradSS_Pkpo2_T_quad, basis_Polyk_E_quad, qr_2k_T, qr_k_E);
        
    // If k=0, closure equation involving average over the boundary
    if (degree()==0){
      // Integrals of basis functions on E is trivial because it is constant
      double int_basis_P0E = E.measure() * (*edgeBases(E).Polyk).function(0, E.center_mass());
      for (size_t i=0; i<Polykpo2(iT).dimension(); i++){
        // normal component
        BPT(i, localOffset(T,E)) += closure_average_scaling * scalar_product(int_partialT_Pkpo2[i], int_basis_P0E * E.normal());
        // tangential component
        BPT(i, localOffset(T,E)+1) += closure_average_scaling * scalar_product(int_partialT_Pkpo2[i], int_basis_P0E * E.tangent());
      }
    }
  }
  Eigen::MatrixXd PT = MPT.ldlt().solve(BPT);
  
  //
  // STABILISATION
  //
  // Operators delta_T and delta_TE. For delta_T we need to compute Pcurl Icurl (purely Xcurl, not EXcurl)
  // from P^{k+1} to P^k. We therefore compute Icurl:P^{k+1}->XcurlT
  Eigen::MatrixXd IcurlT = Eigen::MatrixXd::Zero(m_xcurl.dimensionCell(iT), Polykpo2(iT).dimension() );
  size_t dim_Pkpo2 = Polykpo2(iT).dimension();
  if (degree()>=1){
    size_t offset_T = m_xcurl.localOffset(T);
    size_t dim_Rolykmo = cellBases(iT).Rolykmo->dimension();
    size_t dim_RolyComplk = cellBases(iT).RolyComplk->dimension();
    // Icurl: Components on R^{k-1}(T)
    MonomialCellIntegralsType int_mono_2kpo_T = IntegrateCellMonomials(T, 2*degree()+1);
    Eigen::LDLT<Eigen::MatrixXd> cholesky_mass_Rolykmo(GramMatrix(T, *cellBases(iT).Rolykmo, int_mono_2kpo_T));
    IcurlT.block(offset_T, 0, dim_Rolykmo, dim_Pkpo2)
      = cholesky_mass_Rolykmo.solve(GramMatrix(T, *cellBases(iT).Rolykmo, Polykpo2(iT), int_mono_2kpo_T));
    // Icurl: Components R^{c,k}(T)
    Eigen::LDLT<Eigen::MatrixXd> cholesky_mass_RolyComplk(GramMatrix(T, *cellBases(iT).RolyComplk, int_mono_2kpo_T));
    IcurlT.block(offset_T + dim_Rolykmo, 0, dim_RolyComplk, dim_Pkpo2)
      = cholesky_mass_RolyComplk.solve(GramMatrix(T, *cellBases(iT).RolyComplk, Polykpo2(iT), int_mono_2kpo_T));

  }
  // Edge contributions to Icurl and difference operators
  // delta_TE is made of two components: tangential and normal to the edge
  std::vector<Eigen::MatrixXd> deltaTE_tE(T.n_edges());
  std::vector<Eigen::MatrixXd> deltaTE_nE(T.n_edges());
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge &E = *T.edge(iE);
    auto qr_2kpo_E = generate_quadrature_rule(E, 2*degree()+1);
    const VectorRd tE = E.tangent();
    const VectorRd nE = E.normal();
    const size_t dim_PkE = edgeBases(E).Polyk->dimension();
    // Icurl: contribution edge E
    boost::multi_array<double, 2> basis_PkE_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polyk, qr_2kpo_E);
    boost::multi_array<VectorRd, 2> basis_Polykpo2_quad_E = evaluate_quad<Function>::compute(Polykpo2(iT), qr_2kpo_E);
    boost::multi_array<double, 2> basis_Polykpo2_tE_quad_E 
            = transform_values_quad<double>(basis_Polykpo2_quad_E, [&tE](const VectorRd &x)->double { return x.dot(tE);});
    Eigen::LDLT<Eigen::MatrixXd> cholesky_mass_PolykE(compute_gram_matrix(basis_PkE_quad, qr_2kpo_E)); 
    Eigen::MatrixXd MkE_inv_MkEkpoT_tE = cholesky_mass_PolykE.solve(compute_gram_matrix(basis_PkE_quad, basis_Polykpo2_tE_quad_E, qr_2kpo_E));
    IcurlT.block(m_xcurl.localOffset(T, E), 0, dim_PkE, dim_Pkpo2) = MkE_inv_MkEkpoT_tE; 

    // delta_TE: made of tangential and normal components
    boost::multi_array<double, 2> basis_Polykpo2_nE_quad_E 
            = transform_values_quad<double>(basis_Polykpo2_quad_E, [&nE](const VectorRd &x)->double { return x.dot(nE);});    
    Eigen::MatrixXd MkE_inv_MkEkpoT_nE = cholesky_mass_PolykE.solve(compute_gram_matrix(basis_PkE_quad, basis_Polykpo2_nE_quad_E, qr_2kpo_E));
    deltaTE_nE[iE] = MkE_inv_MkEkpoT_nE * PT;
    deltaTE_nE[iE].block(0, localOffset(T, E), dim_PkE, dim_PkE) -= Eigen::MatrixXd::Identity(dim_PkE, dim_PkE);
    
    deltaTE_tE[iE] = MkE_inv_MkEkpoT_tE * PT;
    deltaTE_tE[iE].block(0, localOffset(T, E) + dim_PkE, dim_PkE, dim_PkE) -= Eigen::MatrixXd::Identity(dim_PkE, dim_PkE);
  }
  
  Eigen::MatrixXd deltaT = m_xcurl.cellOperators(iT).potential * IcurlT * PT - ddrPotential(iT);
  
  // Create stabilisation matrix
  Eigen::MatrixXd ST = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  for (size_t iE = 0; iE < T.n_edges(); iE++){
    const Edge &E = *T.edge(iE);
    auto qr_2k_E = generate_quadrature_rule(E, 2*degree());
    const VectorRd tE = E.tangent();
    const VectorRd nE = E.normal();
   
    boost::multi_array<double, 2> basis_PkE_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polyk, qr_2k_E);
    Eigen::MatrixXd MEE = compute_gram_matrix(basis_PkE_quad, qr_2k_E);
    boost::multi_array<VectorRd, 2> basis_Polyk2_quad_E = evaluate_quad<Function>::compute(*cellBases(iT).Polyk2, qr_2k_E);
    boost::multi_array<double, 2> basis_Polyk2_nE_quad_E 
        = transform_values_quad<double>(basis_Polyk2_quad_E, [&nE](const VectorRd &x)->double { return x.dot(nE);});
    boost::multi_array<double, 2> basis_Polyk2_tE_quad_E 
        = transform_values_quad<double>(basis_Polyk2_quad_E, [&tE](const VectorRd &x)->double { return x.dot(tE);});
    
    Eigen::MatrixXd MET_nE = compute_gram_matrix(basis_PkE_quad, basis_Polyk2_nE_quad_E, qr_2k_E);
    Eigen::MatrixXd MET_tE = compute_gram_matrix(basis_PkE_quad, basis_Polyk2_tE_quad_E, qr_2k_E);

    Eigen::MatrixXd deltaTE_m_deltaT_nE = deltaTE_nE[iE] - MEE.ldlt().solve(MET_nE)*deltaT;
    ST +=  deltaTE_m_deltaT_nE.transpose() * MEE * deltaTE_m_deltaT_nE;

    Eigen::MatrixXd deltaTE_m_deltaT_tE = deltaTE_tE[iE] - MEE.ldlt().solve(MET_tE)*deltaT;
    ST +=  deltaTE_m_deltaT_tE.transpose() * MEE * deltaTE_m_deltaT_tE;

  }
  ST /= T.diam();
  
  return hhoLocalOperators(GT, GsT, PT, ST);
  
}
  
//------------------------------------------------------------------------------
//        Functions to compute matrices for local L2 products on EXcurl
//------------------------------------------------------------------------------

Eigen::MatrixXd EXCurl::computeL2Product(
                                        const size_t iT,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk2_T,
                                        const IntegralWeight & weight
                                        ) const
{
  return m_reduction[iT].transpose() * m_xcurl.computeL2Product(iT, penalty_factor, mass_Pk2_T, weight) * m_reduction[iT];
}

Eigen::MatrixXd EXCurl::computeL2ProductGradient(
                                        const size_t iT,
                                        const XGrad & x_grad,
                                        const std::string & side,
                                        const double & penalty_factor,
                                        const Eigen::MatrixXd & mass_Pk2_T,
                                        const IntegralWeight & weight
                                        ) const
{
  if (side == "left"){
    return m_xcurl.computeL2ProductGradient(iT, x_grad, side, penalty_factor, mass_Pk2_T, weight) * m_reduction[iT];
  }else if (side == "right"){
    return m_reduction[iT].transpose() * m_xcurl.computeL2ProductGradient(iT, x_grad, side, penalty_factor, mass_Pk2_T, weight);
  }else{
    return m_xcurl.computeL2ProductGradient(iT, x_grad, side, penalty_factor, mass_Pk2_T, weight);
  }
}



