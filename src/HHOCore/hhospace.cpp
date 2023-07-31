#include <cassert>

#include <hhospace.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------

HHOSpace::HHOSpace(const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(mesh,
	     0,
	     PolynomialSpaceDimension<Edge>::Poly(K),
	     PolynomialSpaceDimension<Cell>::Poly(K)
	     ),
	  m_mesh(mesh),
    m_K(K),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_edge_bases(mesh.n_edges()),
    m_operators(mesh.n_cells())
{
  m_output << "[HHOSpace] Initializing" << std::endl;
  
  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[HHOSpace] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct edge bases
  std::function<void(size_t, size_t)> construct_all_edge_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        this->m_edge_bases[iE].reset( new EdgeBases(_construct_edge_bases(iE)) );
	      } // for iE
      };
  
  m_output << "[HHOSpace] Constructing edge bases" << std::endl;
  parallel_for(mesh.n_edges(), construct_all_edge_bases, use_threads);

  // Construct gradients, potentials and stabilisation
  std::function<void(size_t, size_t)> construct_all_operators
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_operators[iT].reset( new LocalOperators(_compute_operators(iT)) );
        } // for iT
      };

  m_output << "[HHOSpace] Constructing operators" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_operators, use_threads);

}

//------------------------------------------------------------------------------
// Polynomial bases
//------------------------------------------------------------------------------

HHOSpace::CellBases HHOSpace::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;
  
  MonomialCellIntegralsType int_monoT_2kp2 = IntegrateCellMonomials(T, 2*(m_K+1));
  
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T), Pk(T)
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, GramMatrix(T, basis_Pkpo_T, int_monoT_2kp2))) );  

  MonomialScalarBasisCell basis_Pk_T(T, m_K);
  bases_T.Polyk.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pk_T, GramMatrix(T, basis_Pk_T, int_monoT_2kp2))) );  

  // Check that we got the dimensions right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polyk->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Pk(T)^d
  //------------------------------------------------------------------------------

  bases_T.Polykd.reset( new PolydBasisCellType(*bases_T.Polyk) );
  assert( bases_T.Polykd->dimension() == dimspace * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  return bases_T;
}

//------------------------------------------------------------------------------

HHOSpace::EdgeBases HHOSpace::_construct_edge_bases(size_t iE)
{
  const Edge & E = *m_mesh.edge(iE);
  
  EdgeBases bases_E;

  MonomialEdgeIntegralsType int_monoE_2k = IntegrateEdgeMonomials(E, 2*m_K);

  //------------------------------------------------------------------------------
  // Basis for Pk(E)
  //------------------------------------------------------------------------------
  MonomialScalarBasisEdge basis_Pk_E(E, m_K);
  bases_E.Polyk.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pk_E, GramMatrix(E, basis_Pk_E, int_monoE_2k))) );
  
  // Check that we got the dimensions right
  assert( bases_E.Polyk->dimension() == PolynomialSpaceDimension<Edge>::Poly(m_K) );

  return bases_E;
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd HHOSpace::interpolate(const FunctionType & q, const int doe_cell, const int doe_edge) const
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
          auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polyk, quad_dqr_E);
          qh.segment(globalOffset(E), PolynomialSpaceDimension<Edge>::Poly(degree())) 
            = l2_projection(q, *edgeBases(iE).Polyk, quad_dqr_E, basis_Pk_E_quad, GramMatrix(E, *edgeBases(iE).Polyk));
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
          auto basis_Pk_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, quad_dqr_T);
          qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree())) 
            = l2_projection(q, *cellBases(iT).Polyk, quad_dqr_T, basis_Pk_T_quad, GramMatrix(T, *cellBases(iT).Polyk));
        } // for iT
      };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);


  return qh;
}

//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

HHOSpace::LocalOperators HHOSpace::_compute_operators(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  // Dimension
  size_t dim_Pkpo = PolynomialSpaceDimension<Cell>::Poly(degree()+1);

  //------------------------------------------------------------------------------
  // Gradient
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  // Compute all integrals of monomial powers to degree 2k+3 and the mass matrix
  MonomialCellIntegralsType int_mono_2kp2 = IntegrateCellMonomials(T, 2*degree()+2);
  Eigen::MatrixXd MGT = GramMatrix(T, *cellBases(iT).Polykd, int_mono_2kp2);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BGT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polykd->dimension(), dimensionCell(iT));

  // Boundary contribution
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    DecomposePoly dec(E, MonomialScalarBasisEdge(E, degree()));
    auto PkdT_dot_nTE_nodes = scalar_product(evaluate_quad<Function>::compute(*cellBases(iT).Polykd, dec.get_nodes()), T.edge_normal(iE));
    auto PkdT_dot_nTE_family_PkE = dec.family(PkdT_dot_nTE_nodes);
    // PE extracts the edge unknowns corresponding to E when read among all the element & edge unknowns of T
    Eigen::MatrixXd PE = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E), dimensionEdge(E)));
    MonomialEdgeIntegralsType int_mono_2k_E = IntegrateEdgeMonomials(E, 2*degree());
    BGT += GramMatrix(E, PkdT_dot_nTE_family_PkE, *edgeBases(E).Polyk, int_mono_2k_E) * PE;

    // Following commented block could replace the block above, without using DecomposePoly (more expensive, but sometimes better rounding errors)
    /*
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree() );
    auto basis_Pkd_nTE_E_quad = scalar_product(
				         evaluate_quad<Function>::compute(*cellBases(iT).Polykd, quad_2k_E),
				         T.edge_normal(iE)
				         );
    auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polyk, quad_2k_E);
    Eigen::MatrixXd PE = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E), dimensionEdge(E)));
    BGT += compute_gram_matrix(basis_Pkd_nTE_E_quad, basis_Pk_E_quad, quad_2k_E) * PE;
    */
        
  } // for iE

  // Cell contribution
  DivergenceBasis<HHOSpace::PolydBasisCellType> div_Pkd_basis(*cellBases(iT).Polykd);
  BGT.rightCols(numLocalDofsCell()) -= GramMatrix(T, div_Pkd_basis, *cellBases(iT).Polyk, int_mono_2kp2);

  Eigen::MatrixXd GT = MGT.ldlt().solve(BGT);
  
  //------------------------------------------------------------------------------
  // Potential
  //------------------------------------------------------------------------------

  // We write \nabla pT = projection on \nabla P^{k+1} of GT, and add the closure relation:
  //    (nabla pT v, nabla w)_T + lambda_T(p_T v,1)_T(w,1)_T = (GT v,\nabla w)_T + lambda_T(v_T,1)_T(w,1)_T

  // LHS
  GradientBasis<HHOSpace::PolyBasisCellType> GradPolykpo = GradientBasis(*cellBases(iT).Polykpo);
  Eigen::MatrixXd StiffT = GramMatrix(T, GradPolykpo, int_mono_2kp2);

  // RHS, starting with volumetric term
  Eigen::MatrixXd RHS_PT = GramMatrix(T, GradPolykpo, *cellBases(iT).Polykd, int_mono_2kp2) * GT;
  
  // Vector LT of (phi_j,1)_T for phi_j up to degree K+1, and LT^t*LT, for the closure relation
  Eigen::VectorXd LT = GramMatrix(T, *cellBases(iT).Polykpo, RestrictedBasis(*cellBases(iT).Polykpo, 1), int_mono_2kp2);
  Eigen::MatrixXd LTtLT = LT * (LT.transpose());
  double scalT = StiffT.trace() / std::pow(LT.norm(), 2);
 
  // Add closure relation and compute PT
  RHS_PT.block(
              0, localOffset(T), 
              dim_Pkpo, numLocalDofsCell() 
              ) += 
      scalT * LTtLT.topLeftCorner( dim_Pkpo, numLocalDofsCell() );
  Eigen::MatrixXd PT = ((StiffT + scalT*LTtLT).ldlt()).solve(RHS_PT);

  //------------------------------------------------------------------------------
  // Stabilisation.
  //    It is based on bilinear forms on the local space, and the operators Id - I_T P_T
  //  The bilinear forms built here are (depending on choice_bilinear):
  //        1) 'components_norm': ||v_T||_T^2 + h_T\sum_E ||v_E||_E^2 with scaling depending on regT
  //        2) 'original_hho': original hho one, h_T \sum_E ||v_E-v_T||_E^2
  //------------------------------------------------------------------------------
  std::string choice_bilinear = "original_hho";

  Eigen::MatrixXd BilinearForm = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  Eigen::MatrixXd Id_minus_ITPT = Eigen::MatrixXd::Identity(dimensionCell(iT), dimensionCell(iT));

  // Element contributions
  Eigen::MatrixXd GramTk = GramMatrix(T, *cellBases(iT).Polyk, int_mono_2kp2);
  Eigen::MatrixXd GramTk_kpo = GramMatrix(T, *cellBases(iT).Polyk, *cellBases(iT).Polykpo, int_mono_2kp2);

  // Element contribution to bilinear form only for the components norm
  if (choice_bilinear == "components_norm"){
    double regT = std::pow(T.diam(), dimspace)/T.measure() * T.n_edges();
    BilinearForm.bottomRightCorner(numLocalDofsCell(), numLocalDofsCell()) +=  regT * GramTk;
  }
  Id_minus_ITPT.bottomRows(numLocalDofsCell()) -= GramTk.ldlt().solve(GramTk_kpo * PT);

  // Edge contributions
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    DecomposePoly dec(E, MonomialScalarBasisEdge(E, degree()+1));
    auto PkpoT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polykpo, dec.get_nodes());
    auto PkpoT_family_PkE = dec.family(PkpoT_nodes);

    MonomialEdgeIntegralsType int_mono_2kpo_E = IntegrateEdgeMonomials(E, 2*degree()+1);
    Eigen::MatrixXd GramEk = GramMatrix(E, *edgeBases(E).Polyk, int_mono_2kpo_E);
    Eigen::MatrixXd GramEk_Tkpo = GramMatrix(E, *edgeBases(E).Polyk, PkpoT_family_PkE, int_mono_2kpo_E);
    Eigen::LDLT<Eigen::MatrixXd> GramEk_inv(GramEk);
    GramEk_inv.compute(GramEk);
    
    // PE extracts the edge unknowns corresponding to E when read among all the element & edge unknowns of T
    Eigen::MatrixXd PE = extendOperator(T, E, Eigen::MatrixXd::Identity(dimensionEdge(E), dimensionEdge(E)));

    // For the original hho stabilisation, PE must correspond to v_E - v_T (the latter being projected on P^k(E)^d)
    if (choice_bilinear == "original_hho"){
      auto PkT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, dec.get_nodes());
      auto PkT_family_PkE = dec.family(PkT_nodes);
      Eigen::MatrixXd GramEk_Tk = GramMatrix(E, *edgeBases(E).Polyk, PkT_family_PkE, int_mono_2kpo_E);

      PE.rightCols(numLocalDofsCell()) -= GramEk_inv.solve(GramEk_Tk);
    }
    
    BilinearForm += T.diam() * PE.transpose() * GramEk * PE;
    Id_minus_ITPT.middleRows(localOffset(T, E), numLocalDofsEdge()) -= GramEk_inv.solve(GramEk_Tkpo * PT);

  } // for iE

  // Creation of the stabilisation
  Eigen::MatrixXd ST = std::pow(T.diam(), -2) * Id_minus_ITPT.transpose() * BilinearForm * Id_minus_ITPT;
  

  return LocalOperators(GT, PT, ST);
}




//------------------------------------------------------------------------------
// Norms
//------------------------------------------------------------------------------

std::vector<std::pair<double,double>> HHOSpace::computeNorms( const std::vector<Eigen::VectorXd> & list_dofs ) const
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
        Eigen::MatrixXd mass_PkT = GramMatrix(T, *cellBases(iT).Polyk, int_mono_2k);
        Eigen::MatrixXd mass_GradPkT = GramMatrix(T, GradientBasis<HHOSpace::PolyBasisCellType>(*cellBases(iT).Polyk), int_mono_2k);

        // Cell contributions
        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT_cell = restrict(T, list_dofs[i]).tail(numLocalDofsCell());

          // L2 norm
          local_L2_sqnorms[i](iT) += uT_cell.transpose() * mass_PkT * uT_cell;
          // H1 norm
          local_H1_sqnorms[i](iT) += uT_cell.transpose() * mass_GradPkT * uT_cell;
        }
    
        // Edge contributions
        for (size_t iE=0; iE<T.n_edges(); iE++){
          Edge & E = *T.edge(iE);

          MonomialEdgeIntegralsType int_mono_2kE = IntegrateEdgeMonomials(E, 2*degree());
          Eigen::MatrixXd mass_PkE = GramMatrix(E, *edgeBases(E).Polyk, int_mono_2kE);
          
          DecomposePoly dec(E, MonomialScalarBasisEdge(E, degree()));
          auto PkT_nodes = evaluate_quad<Function>::compute(*cellBases(iT).Polyk, dec.get_nodes());
          auto PkT_family_PkE = dec.family(PkT_nodes);
          Eigen::MatrixXd mass_PkE_PkT = GramMatrix(E, *edgeBases(E).Polyk, PkT_family_PkE, int_mono_2kE);

          for (size_t i=0; i<nb_vectors; i++){
            Eigen::VectorXd uT = restrict(T, list_dofs[i]);

            // Decompose u_E-u_T on Polyk basis of E
            Eigen::VectorXd uE_minus_uT = 
                    uT.segment(iE*numLocalDofsEdge(), numLocalDofsEdge())
                    -
                    mass_PkE.ldlt().solve(mass_PkE_PkT * uT.tail(numLocalDofsCell()));
            // Contribution of gradients, without any weight (no permeability)
            double sqnorm_uE_minus_uT = uE_minus_uT.transpose() * mass_PkE * uE_minus_uT;
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
Eigen::VectorXd HHOSpace::computeVertexValues(const Eigen::VectorXd & u) const
{
  Eigen::VectorXd values = Eigen::VectorXd::Zero(m_mesh.n_vertices());
  
  // Value at each vertex obtained averaging the values from all the cells around
  for (Vertex * V : mesh().get_vertices()){
    size_t iV = V->global_index();
 
    for (Cell * T : V->get_cells()){
      size_t iT = T->global_index();
      Eigen::VectorXd pTuT = operators(iT).potential * restrict(*T, u);
      
      for (size_t i=0; i < cellBases(iT).Polykpo->dimension(); i++){
        values[iV] += pTuT(i) * cellBases(iT).Polykpo->function(i, V->coords());
      }
    }
    
    values[iV] /= V->n_cells();

  }
  
  return values;
}


