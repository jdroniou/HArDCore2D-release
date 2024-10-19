#include <cassert>

#include <vhhospace.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------

VHHOSpace::VHHOSpace(const Mesh & mesh, size_t K, const CellSelection & BoundaryStab, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(mesh,
	     0,
	     dimspace * PolynomialSpaceDimension<Edge>::Poly(K),
	     dimspace * PolynomialSpaceDimension<Cell>::Poly(K)
	     ),
	  m_mesh(mesh),
    m_K(K),
    m_boundary_stab(BoundaryStab),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_edge_bases(mesh.n_edges()),
    m_operators(mesh.n_cells())
{
  m_output << "[VHHOSpace] Initializing" << std::endl;
  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[VHHOSpace] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct edge bases
  std::function<void(size_t, size_t)> construct_all_edge_bases
    = [this](size_t start, size_t end)->void
      {
          for (size_t iE = start; iE < end; iE++) {
	        this->m_edge_bases[iE].reset( new EdgeBases(_construct_edge_bases(iE)) );
          } // for iE
      };
  
  m_output << "[VHHOSpace] Constructing edge bases" << std::endl;
  parallel_for(mesh.n_edges(), construct_all_edge_bases, use_threads);

  // Construct gradients, potentials and stabilisation
  std::function<void(size_t, size_t)> construct_all_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_operators[iT].reset( new LocalOperators(_compute_operators(iT)) );
        } // for iT
      };

  m_output << "[VHHOSpace] Constructing operators" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_operators, use_threads);
}

//------------------------------------------------------------------------------
// Polynomial bases
//------------------------------------------------------------------------------

VHHOSpace::CellBases VHHOSpace::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;
  
  MonomialCellIntegralsType int_monoT_2kp2 = IntegrateCellMonomials(T, 2*(m_K+1));
  
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T), Pk(T), and vector versions
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, GramMatrix(T, basis_Pkpo_T, int_monoT_2kp2))) );
  bases_T.Polykpod.reset( new PolydBasisCellType(*bases_T.Polykpo) );

  MonomialScalarBasisCell basis_Pk_T(T, m_K);
  bases_T.Polyk.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pk_T, GramMatrix(T, basis_Pk_T, int_monoT_2kp2))) );
  bases_T.Polykd.reset( new PolydBasisCellType(*bases_T.Polyk) );

  // Check that we got the dimensions right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polykpod->dimension() == dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polyk->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K) );
  assert( bases_T.Polykd->dimension() == dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Pk(T)^{dxd}
  //------------------------------------------------------------------------------

  bases_T.Polykdxd.reset( new PolydxdBasisCellType(*bases_T.Polyk) );
  assert( bases_T.Polykdxd->dimension() == dimspace * dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K) );

  return bases_T;
}

//------------------------------------------------------------------------------

VHHOSpace::EdgeBases VHHOSpace::_construct_edge_bases(size_t iE)
{
  const Edge & E = *m_mesh.edge(iE);
  
  EdgeBases bases_E;

  MonomialEdgeIntegralsType int_monoF_2k = IntegrateEdgeMonomials(E, 2*m_K);

  //------------------------------------------------------------------------------
  // Basis for Pk(E)
  //------------------------------------------------------------------------------
  MonomialScalarBasisEdge basis_Pk_E(E, m_K);
  bases_E.Polyk.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pk_E, GramMatrix(E, basis_Pk_E, int_monoF_2k))) );

  //------------------------------------------------------------------------------
  // Basis for Pk(E)^d
  //------------------------------------------------------------------------------
  // Old version: tensorized on the canonical basis of R^3 (requires to change the class PolydBasisFaceType)
  // bases_F.Polykd.reset( new PolydBasisFaceType(*bases_F.Polyk) );
  
  // For internal edges, the basis of Pk(E)^d$ is just a tensorized one, that we deal as a (trivial) Family<Tensorized> to match the expected class.
  // For boundary edges, the basis of Pk(E)^d is designed such that, if r=dim Pk(E), the first r polynomials are tangent to E and the last r polynomials are orthogonal to E
  if (!E.is_boundary()){
    TensorizedVectorFamily<PolyBasisEdgeType, dimspace> basis_tens_PolykdE(*bases_E.Polyk);
    bases_E.Polykd.reset(
              new PolydBasisEdgeType(basis_tens_PolykdE,
                                     Eigen::MatrixXd::Identity(basis_tens_PolykdE.dimension(), basis_tens_PolykdE.dimension()) )
                        );                                     
  }else{
    bases_E.Polykd.reset(
        new PolydBasisEdgeType(
            GenericTensorization<PolyBasisEdgeType, dimspace>(*bases_E.Polyk, std::vector<Eigen::VectorXd> {E.tangent(), E.normal()})
            )
          );
  }
  
  // Check that we got the dimensions right
  assert( bases_E.Polyk->dimension() == PolynomialSpaceDimension<Edge>::Poly(m_K) );
  assert( bases_E.Polykd->dimension() == dimspace * PolynomialSpaceDimension<Edge>::Poly(m_K) );

  return bases_E;
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd VHHOSpace::interpolate(const FunctionType & q, const int doe_cell, const int doe_edge) const
{
  Eigen::VectorXd qh = Eigen::VectorXd::Zero(dimension());

  // Degrees of quadrature rules
  size_t dqr_cell = (doe_cell >= 0 ? doe_cell : 2 * degree() + 3);
  size_t dqr_edge = (doe_edge >= 0 ? doe_edge : 2 * degree() + 3);
      
  // Interpolate at edges
  std::function<void(size_t, size_t)> interpolate_edges
    = [this, &qh, q, &dqr_edge](size_t start, size_t end)->void
      {
        for (size_t iE = start; iE < end; iE++) {
          const Edge & E = *mesh().edge(iE);
          QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr_edge);
          auto basis_Pkd_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykd, quad_dqr_E);
          qh.segment(globalOffset(E), numLocalDofsEdge())
            = l2_projection(q, *edgeBases(iE).Polykd, quad_dqr_E, basis_Pkd_E_quad, GramMatrix(E, *edgeBases(iE).Polykd));
        } // for iE
      };
  parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);

  // Interpolate at cells
  std::function<void(size_t, size_t)> interpolate_cells
    = [this, &qh, q, &dqr_cell](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          const Cell & T = *mesh().cell(iT);
          QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);
          auto basis_Pkd_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, quad_dqr_T);
          qh.segment(globalOffset(T), numLocalDofsCell()) 
            = l2_projection(q, *cellBases(iT).Polykd, quad_dqr_T, basis_Pkd_T_quad, GramMatrix(T, *cellBases(iT).Polykd));
        } // for iT
      };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);

  return qh;
}

Eigen::VectorXd VHHOSpace::local_interpolate(size_t iT, const FunctionType & q, const int doe_cell, const int doe_edge) const
{
    const Cell & T = *mesh().cell(iT);
    Eigen::VectorXd qT = Eigen::VectorXd::Zero(dimensionCell(T));

    // Global interpolate
    Eigen::VectorXd qh = interpolate(q, doe_cell, doe_edge);

    // Interpolate at edges
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
        const Edge & E = *T.edge(iE);
        qT.segment(localOffset(T,E), numLocalDofsEdge()) = qh.segment(globalOffset(E), numLocalDofsEdge());
    }

    // Interpolate at cells
    qT.segment(localOffset(T),numLocalDofsCell()) = qh.segment(globalOffset(T), numLocalDofsCell());
    
    return qT;
}


//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

VHHOSpace::LocalOperators VHHOSpace::_compute_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  // Dimension
  size_t dim_Pkpod = dimspace * PolynomialSpaceDimension<Cell>::Poly(degree()+1);

  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  // Compute all integrals of monomial powers to degree 2k+3 and the mass matrix
  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
  Eigen::MatrixXd MGT = GramMatrix(T, *cellBases(iT).Polykdxd, int_mono_2kp2);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BGT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polykdxd->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    DecomposePoly dec(E, TensorizedVectorFamily<MonomialScalarBasisEdge, dimspace>(MonomialScalarBasisEdge(E, degree())));
    VectorRd nTE = T.edge_normal(iE);
    auto PkdxdT_nTE_nodes = transform_values_quad<VectorRd>(evaluate_quad<Function>::compute(*cellBases(iT).Polykdxd, dec.get_nodes()), [&nTE](const MatrixRd &M)->VectorRd { return M*nTE;});
    auto PkdxdT_nTE_family_PkdE = dec.family(PkdxdT_nTE_nodes);
    // PE extracts the edge unknowns corresponding to E when read among all the element & edge unknowns of T
    Eigen::MatrixXd PE = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E), dimensionEdge(E)));
    MonomialEdgeIntegralsType int_mono_2k_E = IntegrateEdgeMonomials(E, 2*degree());
    BGT += GramMatrix(E, PkdxdT_nTE_family_PkdE, *edgeBases(E).Polykd, int_mono_2k_E) * PE;
    // Replaces the block above, without using DecomposePoly
    /*
      QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree() );
      VectorRd nTE = T.edge_normal(iE);
      auto basis_Pkdxd_nTE_E_quad = transform_values_quad<VectorRd>(
				         evaluate_quad<Function>::compute(*cellBases(iT).Polykdxd, quad_2k_E),
				         [&nTE](const MatrixRd &M)->VectorRd { return M*nTE;}
				         );
      auto basis_Pkd_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykd, quad_2k_E);
      Eigen::MatrixXd PE = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E), dimensionEdge(E)));
      BGT += compute_gram_matrix(basis_Pkdxd_nTE_E_quad, basis_Pkd_E_quad, quad_2k_E) * PE;
    */
    
  } // for iE

  // Cell contribution
  DivergenceBasis<VHHOSpace::PolydxdBasisCellType> div_Pkdxd_basis(*cellBases(iT).Polykdxd);
  BGT.rightCols(numLocalDofsCell()) -= GramMatrix(T, div_Pkdxd_basis, *cellBases(iT).Polykd, int_mono_2kp2);

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);

  //------------------------------------------------------------------------------
  // HHO Potential
  //------------------------------------------------------------------------------

  // We write \nabla pT = projection on \nabla P^{k+1} of GT, and add the closure relation:
  //    (nabla pT v, nabla w)_T + lambda_T(p_T v,1)_T(w,1)_T = (GT v,\nabla w)_T + lambda_T(v_T,1)_T(w,1)_T

  // LHS
  GradientBasis<VHHOSpace::PolydBasisCellType> GradPolykpod = GradientBasis(*cellBases(iT).Polykpod);
  Eigen::MatrixXd StiffT = GramMatrix(T, GradPolykpod, int_mono_2kp2);

  // RHS, starting with volumetric term
  Eigen::MatrixXd RHS_PT = GramMatrix(T, GradPolykpod, *cellBases(iT).Polykdxd, int_mono_2kp2) * GT;

  // Closure matrix (\int phi_i).(\int phi_j) computed as \sum_k \int phi_i.e_k \int phi_j.e_k, with e_k basis of R^d represented by a polynomial basis for P^0(T)^d
  auto basis_P0d = TensorizedVectorFamily<MonomialScalarBasisCell, dimspace>(MonomialScalarBasisCell(T, 0));
  Eigen::MatrixXd Integrate_Pkpod_P0d = GramMatrix(T, *cellBases(iT).Polykpod, basis_P0d, int_mono_2kp2);
  Eigen::MatrixXd Integrate_Pkd_P0d = GramMatrix(T, *cellBases(iT).Polykd, basis_P0d, int_mono_2kp2);

  Eigen::MatrixXd Closure = Integrate_Pkpod_P0d * Integrate_Pkpod_P0d.transpose();
  double scalT = StiffT.trace() / Closure.trace();

  RHS_PT.block(
              0, localOffset(T),
              dim_Pkpod, numLocalDofsCell()
              ) += scalT * (Integrate_Pkpod_P0d * Integrate_Pkd_P0d.transpose());

  Eigen::MatrixXd PT = ((StiffT + scalT * Closure).ldlt()).solve(RHS_PT);

  //------------------------------------------------------------------------------
  // Stabilisations: HHO and coming from the divergence.
  //    They are both based on bilinear forms on the local space, and operators Id - I_T P_T
  //    (with P_T the HHO or divergence potential)
  //  The bilinear forms for HHO built here are (depending on choice_bilinear):
  //        1) 'components_norm': ||v_T||_T^2 + h_T\sum_E ||v_E||_E^2 with scaling depending on regT
  //        2) 'original_hho': original hho one, h_T \sum_E ||v_E-v_T||_E^2
  //  For the divergence stabilisation, it is always 'components_norm'
  //------------------------------------------------------------------------------
  std::string choice_bilinear = "components_norm";

  Eigen::MatrixXd BilinearForm = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  Eigen::MatrixXd Id_minus_ITPT = Eigen::MatrixXd::Identity(dimensionCell(iT), dimensionCell(iT));

  // Element contributions
  Eigen::MatrixXd GramTkd = GramMatrix(T, *cellBases(iT).Polykd, int_mono_2kp2);
  Eigen::MatrixXd GramTkd_kpod = GramMatrix(T, *cellBases(iT).Polykd, *cellBases(iT).Polykpod, int_mono_2kp2);
  double regT = std::pow(T.diam(), dimspace)/T.measure() * T.n_edges();
  // Element contribution to bilinear form only for the components norm
  if (choice_bilinear == "components_norm"){
    BilinearForm.bottomRightCorner(numLocalDofsCell(), numLocalDofsCell()) +=  regT * GramTkd;
  }
  Id_minus_ITPT.bottomRows(numLocalDofsCell()) -= GramTkd.ldlt().solve(GramTkd_kpod * PT);

  // Edge contributions
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);

    // For boundary faces, we include the contribution to the inner product only if m_boundary_stab says so (we do not complete
    // the construction of Id_minus_ITPT for these faces because that part will not be used by BilinearForm; so Id_minus_ITPT
    // is not empty on these faces but incorrect, only contain u_F)
    if ( !E.is_boundary() || m_boundary_stab(T) ){

      DecomposePoly dec(E, TensorizedVectorFamily<MonomialScalarBasisEdge, dimspace>(MonomialScalarBasisEdge(E, degree()+1)));
      auto PkpodT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykpod, dec.get_nodes());
      auto PkpodT_family_PkdE = dec.family(PkpodT_nodes);
      auto PkdT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, dec.get_nodes());
      auto PkdT_family_PkdE = dec.family(PkdT_nodes);

      MonomialEdgeIntegralsType int_mono_2kpo_E = IntegrateEdgeMonomials(E, 2*degree()+1);
      Eigen::MatrixXd GramEkd = GramMatrix(E, *edgeBases(E).Polykd, int_mono_2kpo_E);
      Eigen::MatrixXd GramEkd_Tkpod = GramMatrix(E, *edgeBases(E).Polykd, PkpodT_family_PkdE, int_mono_2kpo_E);
      Eigen::MatrixXd GramEkd_Tkd = GramMatrix(E, *edgeBases(E).Polykd, PkdT_family_PkdE, int_mono_2kpo_E);

      // Replaces the block above, without using DecomposePoly
      /*
        QuadratureRule quad_2kpo_F = generate_quadrature_rule(F, 2 * degree() + 1);
        auto PkpodT_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykpod, quad_2kpo_F);
        auto PkdT_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, quad_2kpo_F);
        auto PkdF_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykd, quad_2kpo_F);
        Eigen::MatrixXd GramFkd_Tkpod = compute_gram_matrix(PkdF_quad, PkpodT_quad, quad_2kpo_F);
        Eigen::MatrixXd GramFkd_Tkd = compute_gram_matrix(PkdF_quad, PkdT_quad, quad_2kpo_F);
        Eigen::MatrixXd GramFkd = GramMatrix(F, *faceBases(F).Polykd);
      */

    Eigen::LDLT<Eigen::MatrixXd> GramEkd_inv(GramEkd);
    GramEkd_inv.compute(GramEkd);

    // PE extracts the edge unknowns corresponding to E when read among all the element & edge unknowns of T
    Eigen::MatrixXd PE = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E), dimensionEdge(E)));

    // For the original hho stabilisation, PE must correspond to v_E - v_T (the latter being projected on P^k(E)^d)
    if (choice_bilinear == "original_hho"){
      PE.rightCols(numLocalDofsCell()) -= GramEkd_inv.solve(GramEkd_Tkd);
    }
      
    BilinearForm += T.diam() * PE.transpose() * GramEkd * PE;
    Id_minus_ITPT.middleRows(localOffset(T, E), numLocalDofsEdge()) -= GramEkd_inv.solve(GramEkd_Tkpod * PT);

    }
  } // for iE
  // Creation of the stabilisations
  Eigen::MatrixXd ST = std::pow(T.diam(), -2) * Id_minus_ITPT.transpose() * BilinearForm * Id_minus_ITPT;

  return LocalOperators (GT, PT, ST);
}



//------------------------------------------------------------------------------
// Norms
//------------------------------------------------------------------------------

std::vector<std::pair<double,double>> VHHOSpace::computeNorms( const std::vector<Eigen::VectorXd> & list_dofs ) const
{
  size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_L2_sqnorms(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));
  std::vector<Eigen::VectorXd> local_H1_sqnorms(nb_vectors, Eigen::VectorXd::Zero(mesh().n_cells()));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &list_dofs, &local_L2_sqnorms, &local_H1_sqnorms, &nb_vectors](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *mesh().cell(iT);

        // Mass matrices
	      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*degree());
        Eigen::MatrixXd mass_PkdT = GramMatrix(T, *cellBases(iT).Polykd, int_mono_2k);
        Eigen::MatrixXd mass_GradPkdT = GramMatrix(T, GradientBasis<VHHOSpace::PolydBasisCellType>(*cellBases(iT).Polykd), int_mono_2k);

        // Local vectors
        std::vector<Eigen::VectorXd> uT(nb_vectors, Eigen::VectorXd::Zero(dimensionCell(iT)));
        for (size_t i=0; i < nb_vectors; i++){
          uT[i] = restrict(T, list_dofs[i]);
        }

        // Cell contributions
        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT_cell = uT[i].tail(numLocalDofsCell());

          // L2 norm
          local_L2_sqnorms[i](iT) += uT_cell.transpose() * mass_PkdT * uT_cell;
          // H1 norm
          local_H1_sqnorms[i](iT) += uT_cell.transpose() * mass_GradPkdT * uT_cell;
        }
    
        // Face contributions
        for (size_t iE=0; iE<T.n_edges(); iE++){
          Edge & E = *T.edge(iE);

          MonomialEdgeIntegralsType int_mono_2kE = IntegrateEdgeMonomials(E, 2*degree());
          Eigen::MatrixXd mass_PkdE = GramMatrix(E, *edgeBases(E).Polykd, int_mono_2kE);

          DecomposePoly dec(E, TensorizedVectorFamily<MonomialScalarBasisEdge, dimspace>(MonomialScalarBasisEdge(E, degree())));
          auto PkdT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, dec.get_nodes());
          auto PkdT_family_PkdE = dec.family(PkdT_nodes);
          Eigen::MatrixXd mass_PkdE_PkdT = GramMatrix(E, *edgeBases(E).Polykd, PkdT_family_PkdE, int_mono_2kE);

          // Replaces the block above, without using DecomposePoly
          /*
            QuadratureRule quad_2k_F = generate_quadrature_rule(F, 2 * degree());
            auto PkdT_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykd, quad_2k_F);
            auto PkdF_quad = evaluate_quad<Function>::compute(*faceBases(F).Polykd, quad_2k_F);
            Eigen::MatrixXd mass_PkdF_PkdT = compute_gram_matrix(PkdF_quad, PkdT_quad, quad_2k_F);          
          */
          
          for (size_t i=0; i<nb_vectors; i++){
            // Decompose u_E-u_T on Polyk basis of E
            Eigen::VectorXd uE_minus_uT =
                    uT[i].segment(iE*numLocalDofsEdge(), numLocalDofsEdge())
                    -
                    mass_PkdE.ldlt().solve(mass_PkdE_PkdT * uT[i].tail(numLocalDofsCell()));
            // Contribution of gradients, without any weight (no permeability)
            double sqnorm_uE_minus_uT = uE_minus_uT.transpose() * mass_PkdE * uE_minus_uT;
            local_H1_sqnorms[i](iT) += sqnorm_uE_minus_uT/T.diam();
          }
        } // for iE
      }
    };
  parallel_for(mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  // Vector of outputs
  std::vector<std::pair<double,double>> list_norms(nb_vectors);
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i].first = std::sqrt(std::abs(local_L2_sqnorms[i].sum()));
    list_norms[i].second = std::sqrt(std::abs(local_H1_sqnorms[i].sum()));
  }
  
  return list_norms;
}

//------------------------------------------------------------------------------
// Vertex values 
//------------------------------------------------------------------------------

std::vector<Eigen::VectorXd> VHHOSpace::computeVertexValues(const Eigen::VectorXd & u) const
{
  std::vector<Eigen::VectorXd> values(dimspace, Eigen::VectorXd::Zero(m_mesh.n_vertices()));
  
  // Value at each vertex obtained averaging the values from all the cells around
  for (Vertex * V : mesh().get_vertices()) {
    size_t iV = V->global_index();

    VectorRd uV = VectorRd::Zero();
 
    for (Cell * T : V->get_cells()) {
      size_t iT = T->global_index();
      Eigen::VectorXd pTuT = operators(iT).potential * restrict(*T, u);
      
      for (size_t i=0; i < cellBases(iT).Polykpod->dimension(); i++) {
        uV += pTuT(i) * cellBases(iT).Polykpod->function(i, V->coords());
      }
    }
    
    uV /= V->n_cells();
    for (size_t i = 0; i < dimspace; i++)
      values[i](iV) = uV(i);
  }
  
  return values;
}                  


