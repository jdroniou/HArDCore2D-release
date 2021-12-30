#include <xdivdiv.hpp>
#include <basis.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_edge.hpp>

#include <functional>

using namespace HArDCore2D;


//------------------------------------------------------------------------------
// Constructor of SymmetricMatrixBasisVertex
//------------------------------------------------------------------------------

XDivDiv::SymmetricMatrixBasisVertex::SymmetricMatrixBasisVertex()
{
  m_basis[0] <<
    1, 0,
    0, 0;
  m_basis[1] <<
    0, 1,
    1, 0;
  m_basis[2] <<
    0, 0,
    0, 1;
}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

XDivDiv::XDivDiv(const PlatesCore & plates_core, bool use_threads, std::ostream & output)
  : GlobalDOFSpace(
                   plates_core.mesh(),
                   3,
                   PolynomialSpaceDimension<Edge>::Poly(plates_core.degree() - 3) + PolynomialSpaceDimension<Edge>::Poly(plates_core.degree() - 2),
                   PolynomialSpaceDimension<Cell>::HolyCompl(plates_core.degree() - 1) + PolynomialSpaceDimension<Cell>::Holy(plates_core.degree() - 4)
                   ),
    m_plates_core(plates_core),
    m_use_threads(use_threads),
    m_output(output),
    m_cell_operators(plates_core.mesh().n_cells()),
    m_edge_potentials(plates_core.mesh().n_edges())  
{
  output << "[XDivDiv] Initializing" << std::endl;
  if (use_threads) {
    m_output << "[XDivDiv] Parallel execution" << std::endl;
  } else {
    m_output << "[XDivDiv] Sequential execution" << std::endl;
  }

  // Construct edge potentials
  std::function<void(size_t, size_t)> construct_all_edge_potentials
    = [this](size_t start, size_t end)->void
    {
      for (size_t iE = start; iE < end; iE++) {
        m_edge_potentials[iE] = _compute_edge_potential(iE);
      } // for iE
    };
  m_output << "[XDivDiv] Constructing edge potentials" << std::endl;
  parallel_for(mesh().n_edges(), construct_all_edge_potentials, use_threads);

  // Construct cell divdiv and potentials
  std::function<void(size_t, size_t)> construct_all_cell_divdiv_potentials
    = [this](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          m_cell_operators[iT].reset( new LocalOperators(_compute_cell_divdiv_potential(iT)) );
        } // for iT
      };
  m_output << "[XDivDiv] Constructing cell divdiv operators and potentials" << std::endl;
  parallel_for(mesh().n_cells(), construct_all_cell_divdiv_potentials, use_threads);  
}

//------------------------------------------------------------------------------
// Interpolator
//------------------------------------------------------------------------------

Eigen::VectorXd XDivDiv::interpolate(
                                     const FunctionType & tau,
                                     const EdgeFunctionType & Dtau,
                                     const int deg_quad
                                     ) const
{
  Eigen::VectorXd tauh = Eigen::VectorXd::Zero(dimension());

  // Degree of quadrature rules
  size_t dqr = (deg_quad >= 0 ? deg_quad : 2 * degree() + 3);

  // Interpolate at vertices
  std::function<void(size_t, size_t)> interpolate_vertices
    = [this, &tauh, tau](size_t start, size_t end)->void
    {
      for (size_t iV = start; iV < end; iV++) {
        const Vertex & V = *mesh().vertex(iV);
        
        Eigen::Matrix2d tau_V = tau(V.coords());

        tauh.segment(globalOffset(V), 3) << tau_V(0,0), tau_V(0,1), tau_V(1,1);
      } // for iV
    };
  parallel_for(mesh().n_vertices(), interpolate_vertices, m_use_threads);

  // Interpolate at edges
  std::function<void(size_t, size_t)> interpolate_edges
    = [this, &tauh, tau, Dtau, &dqr](size_t start, size_t end)->void
    {
      for (size_t iE = start; iE < end; iE++) {
        const Edge & E = *mesh().edge(iE);
        
        QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr);
        auto basis_Pkm3_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykm3, quad_dqr_E);
        auto basis_Pkm2_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykm2, quad_dqr_E);
      
        Eigen::Vector2d nE = E.normal();

        auto tau_nn_E = [&nE, tau](const Eigen::Vector2d & x)->double {
          return (tau(x) * nE).dot(nE);
        };
                
        auto DtauE = [&Dtau, &E](const Eigen::Vector2d & x)->double {
          return Dtau(x, E);
        };

        tauh.segment(globalOffset(E), edgeBases(iE).Polykm3->dimension())
          = l2_projection(tau_nn_E, *edgeBases(iE).Polykm3, quad_dqr_E, basis_Pkm3_quad);

        tauh.segment(globalOffset(E) + edgeBases(iE).Polykm3->dimension(), edgeBases(iE).Polykm2->dimension())          
          = l2_projection(DtauE, *edgeBases(iE).Polykm2, quad_dqr_E, basis_Pkm2_quad);
      } // for iE
    };
  parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);

  // Interpolate at cells
  std::function<void(size_t, size_t)> interpolate_cells
    = [this, &tauh, tau, &dqr](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        const Cell & T = *mesh().cell(iT);

        QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr);

        auto basis_Hckm1_quad = evaluate_quad<Function>::compute(*cellBases(iT).HolyComplkm1, quad_dqr_T);
        tauh.segment(globalOffset(T), cellBases(iT).HolyComplkm1->dimension())
          = l2_projection(tau, *cellBases(iT).HolyComplkm1, quad_dqr_T, basis_Hckm1_quad);
                     
        if (degree() > 3) {
          auto basis_Hkm4_quad = evaluate_quad<Function>::compute(*cellBases(iT).Holykm4, quad_dqr_T);
          tauh.segment(globalOffset(T) + cellBases(iT).HolyComplkm1->dimension(), cellBases(iT).Holykm4->dimension())
            = l2_projection(tau, *cellBases(iT).Holykm4, quad_dqr_T, basis_Hkm4_quad);
        } // if degree > 3
      } // for iT
    };
  parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  
  return tauh;
}

//------------------------------------------------------------------------------
//    Compute operators
//------------------------------------------------------------------------------

Eigen::MatrixXd XDivDiv::_compute_edge_potential(size_t iE)
{
  const Edge & E = *mesh().edge(iE);

  const Vertex & V1 = *E.vertex(0);
  const Vertex & V2 = *E.vertex(1);

  Vector2d nE = E.normal();
  
  Eigen::MatrixXd MPE
    = Eigen::MatrixXd::Zero(edgeBases(iE).Polykm1->dimension(), edgeBases(iE).Polykm1->dimension());
  Eigen::MatrixXd BPE
    = Eigen::MatrixXd::Zero(edgeBases(iE).Polykm1->dimension(), dimensionEdge(iE));

  // Enforce vertex values
  for (size_t j = 0; j < edgeBases(iE).Polykm1->dimension(); j++) {
    MPE(0, j) = edgeBases(iE).Polykm1->function(j, V1.coords());
    MPE(1, j) = edgeBases(iE).Polykm1->function(j, V2.coords());
  } // for j
  
  for (size_t j = 0; j < m_symm_vertex.size(); j++) {
    BPE(0, localOffset(E, V1) + j) = (m_symm_vertex.function(j)*nE).dot(nE);
    BPE(1, localOffset(E, V2) + j) = (m_symm_vertex.function(j)*nE).dot(nE);
  } // for j

  // Enforce projection
  MPE.block(2, 0, edgeBases(iE).Polykm3->dimension(), edgeBases(iE).Polykm1->dimension())
    = GramMatrix(E, *edgeBases(iE).Polykm3, *edgeBases(iE).Polykm1);

  BPE.block(2, localOffset(E), edgeBases(iE).Polykm3->dimension(), edgeBases(iE).Polykm3->dimension())
    = Eigen::MatrixXd::Identity(edgeBases(iE).Polykm3->dimension(), edgeBases(iE).Polykm3->dimension());

  return MPE.partialPivLu().solve(BPE);
}

//------------------------------------------------------------------------------

XDivDiv::LocalOperators XDivDiv::_compute_cell_divdiv_potential(size_t iT)
{
  const Cell & T = *mesh().cell(iT);

  //------------------------------------------------------------------------------
  // Div-div
  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  MonomialCellIntegralsType int_monoT_2km2 = IntegrateCellMonomials(T, 2 * (degree() - 1));
  Eigen::MatrixXd MDDT = GramMatrix(T, *cellBases(iT).Polykm2, int_monoT_2km2);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  Eigen::MatrixXd BDDT
    = Eigen::MatrixXd::Zero(cellBases(iT).Polykm2->dimension(), dimensionCell(iT));

  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);

    Vector2d tE = E.tangent();
    Vector2d nE = E.normal();
    double omegaTE = T.edge_orientation(iE);

    // Vertex terms
    for (size_t iV = 0; iV < 2; iV++){
      const Vertex & V = *E.vertex(iV);
      for (size_t i = 0; i < cellBases(iT).Polykm2->dimension(); i++) {
        for (size_t j = 0; j < m_symm_vertex.size(); j++) {
          BDDT(i, localOffset(T, V) + j)
            -= omegaTE * E.vertex_orientation(iV) * (m_symm_vertex.function(j) * nE).dot(tE) * m_plates_core.cellBases(iT).Polykm2->function(i, V.coords());
        } // for j
      } // for i
    } // for iV
    
    // Edge terms
    QuadratureRule quad_2km4_E = generate_quadrature_rule(E, 2 * (degree() - 2));
    
    auto basis_grad_Pkm2_nE_quad
      = scalar_product(evaluate_quad<Gradient>::compute(*cellBases(iT).Polykm2, quad_2km4_E), nE);

    BDDT.block(0, localOffset(T, E), cellBases(iT).Polykm2->dimension(), edgeBases(E).Polykm3->dimension())      
      += -omegaTE * compute_gram_matrix(
                                        basis_grad_Pkm2_nE_quad,
                                        evaluate_quad<Function>::compute(*edgeBases(E).Polykm3, quad_2km4_E),
                                        quad_2km4_E
                                        );

    BDDT.block(0, localOffset(T, E) + edgeBases(E).Polykm3->dimension(), cellBases(iT).Polykm2->dimension(), edgeBases(E).Polykm2->dimension())
      += omegaTE * compute_gram_matrix(
                                       evaluate_quad<Function>::compute(*cellBases(iT).Polykm2, quad_2km4_E),
                                       evaluate_quad<Function>::compute(*edgeBases(E).Polykm2, quad_2km4_E),
                                       quad_2km4_E
                                       );
  } // for iE

  // Element term
  if (degree() > 3) {
    QuadratureRule quad_2km8_T = generate_quadrature_rule(T, 2 * (degree() - 4));

    BDDT.block(0, localOffset(T) + cellBases(iT).HolyComplkm1->dimension(), cellBases(iT).Polykm2->dimension(), cellBases(iT).Holykm4->dimension())
      += compute_gram_matrix(
                             evaluate_quad<Hessian>::compute(*cellBases(iT).Polykm2, quad_2km8_T),
                             evaluate_quad<Function>::compute(*cellBases(iT).Holykm4, quad_2km8_T),
                             quad_2km8_T
                             );

  } // if degree() > 3

  Eigen::MatrixXd DDT = MDDT.ldlt().solve(BDDT);

  //------------------------------------------------------------------------------
  // Potentials
  //------------------------------------------------------------------------------  

  //------------------------------------------------------------------------------
  // Left-hand side matrix

  Eigen::MatrixXd MPT = Eigen::MatrixXd::Zero(cellBases(iT).PolySymkm1->dimension(), cellBases(iT).PolySymkm1->dimension());
  QuadratureRule quad_2km2_T = generate_quadrature_rule(T, 2 * (degree() - 1));
  auto basis_PSymkm1_quad = evaluate_quad<Function>::compute(*cellBases(iT).PolySymkm1, quad_2km2_T);
  auto basis_HComplkm1_quad = evaluate_quad<Function>::compute(*cellBases(iT).HolyComplkm1, quad_2km2_T);

  ShiftedBasis<Family<MonomialScalarBasisCell>> basis_Pkp1_minus_P1(*cellBases(iT).Polykp1, PolynomialSpaceDimension<Cell>::Poly(1));
  auto basis_Hkm1_quad = evaluate_quad<Hessian>::compute(basis_Pkp1_minus_P1, quad_2km2_T);
  
  MPT.topRows(basis_Pkp1_minus_P1.dimension())
    = compute_gram_matrix(basis_Hkm1_quad, basis_PSymkm1_quad, quad_2km2_T);
  
  MPT.bottomRows(cellBases(iT).HolyComplkm1->dimension()) 
    = compute_gram_matrix(basis_HComplkm1_quad, basis_PSymkm1_quad, quad_2km2_T);

  //------------------------------------------------------------------------------
  // Right-hand side matrix

  // DDT contribution
  Eigen::MatrixXd BPT = Eigen::MatrixXd::Zero(cellBases(iT).PolySymkm1->dimension(), dimensionCell(iT));
  BPT.topRows(basis_Pkp1_minus_P1.dimension()) 
    += GramMatrix(T, basis_Pkp1_minus_P1, *cellBases(iT).Polykm2) * DDT;

  // HolyCompl contribution
  BPT.block(basis_Pkp1_minus_P1.dimension(), localOffset(T), cellBases(iT).HolyComplkm1->dimension(), cellBases(iT).HolyComplkm1->dimension())
    += compute_gram_matrix(basis_HComplkm1_quad, basis_HComplkm1_quad, quad_2km2_T);
  
  // Edge and vertex contributions
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);

    Vector2d tE = E.tangent();
    Vector2d nE = E.normal();
    double omegaTE = T.edge_orientation(iE);

    // Vertex terms
    for (size_t iV = 0; iV < 2; iV++){
      const Vertex & V = *E.vertex(iV);
       for (size_t i = 0; i < basis_Pkp1_minus_P1.dimension(); i++) {
        for (size_t j = 0; j < m_symm_vertex.size(); j++) {
          BPT(i, localOffset(T, V) + j)
            += omegaTE * E.vertex_orientation(iV) * (m_symm_vertex.function(j) * nE).dot(tE) * basis_Pkp1_minus_P1.function(i, V.coords());
        } // for j
      } // for i
    } // for iV
    
    // Edge terms
    QuadratureRule quad_2k_E = generate_quadrature_rule(E, 2 * degree());
    
    auto basis_grad_Pkp1_minus_P1_nE_quad
      = scalar_product(evaluate_quad<Gradient>::compute(basis_Pkp1_minus_P1, quad_2k_E), nE);

    // Edge potential contribution
    BPT.topRows(basis_Pkp1_minus_P1.dimension())     
      += omegaTE * compute_gram_matrix(
                                        basis_grad_Pkp1_minus_P1_nE_quad,
                                        evaluate_quad<Function>::compute(*edgeBases(E).Polykm1, quad_2k_E),
                                        quad_2k_E
                                        ) * extendOperator(T, E, m_edge_potentials[E.global_index()]);

    // Dtau_E contribution
    BPT.block(0, localOffset(T, E) + edgeBases(E).Polykm3->dimension(), basis_Pkp1_minus_P1.dimension(), edgeBases(E).Polykm2->dimension())
      += -omegaTE * compute_gram_matrix(
                                       evaluate_quad<Function>::compute(basis_Pkp1_minus_P1, quad_2k_E),
                                       evaluate_quad<Function>::compute(*edgeBases(E).Polykm2, quad_2k_E),
                                       quad_2k_E
                                       );
  } // for iE

 
    return LocalOperators(DDT, BDDT, MPT.partialPivLu().solve(BPT));
}

//------------------------------------------------------------------------------
//        L2 product
//------------------------------------------------------------------------------

Eigen::MatrixXd XDivDiv::computeL2Product(
                                          const size_t iT,
                                          const ConstitutiveLawType & A,
                                          const double & penalty_factor
                                          ) const
{
  const Cell & T = *mesh().cell(iT); 
  
  Eigen::MatrixXd L2T = Eigen::MatrixXd::Zero(dimensionCell(iT), dimensionCell(iT));
  
  double hT = T.diam();
  const Eigen::MatrixXd & pot_T = cellOperators(iT).potential;

  // Stabilisation: edges  
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);

    VectorRd tE = E.tangent();
    VectorRd nE = E.normal();
    
    QuadratureRule quad_2km2_E = generate_quadrature_rule(E, 2*(degree()-1));

    // Basis P^{k-1}(T;Sym) at the quadrature nodes, transformed by ()*nE.nE
    auto basis_PSymkm1_nnE_quad 
      = transform_values_quad<double>(
                                      evaluate_quad<Function>::compute(*cellBases(iT).PolySymkm1, quad_2km2_E), 
                                      [&nE](const MatrixRd & M)->double { return (M*nE).dot(nE);}
                                      );
    auto basis_Pkm1_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykm1, quad_2km2_E);
    
    //
    // First contribution of edge: hT * int (PT*nE.nE - PE)(PT*nE.nE - PE), developed
    //
    Eigen::MatrixXd M1TTnnE = compute_gram_matrix(basis_PSymkm1_nnE_quad, quad_2km2_E);
    Eigen::MatrixXd M1TnnE_E = compute_gram_matrix(basis_PSymkm1_nnE_quad, basis_Pkm1_E_quad, quad_2km2_E);
    MonomialEdgeIntegralsType int_mono_2km2_E = IntegrateEdgeMonomials(E, 2*degree()-2);
    Eigen::MatrixXd M1EE = GramMatrix(E, *edgeBases(E).Polykm1, int_mono_2km2_E);
    
    Eigen::MatrixXd pot_E = extendOperator(T, E, m_edge_potentials[E.global_index()]);
    Eigen::MatrixXd cross_PTnnE_PE = pot_T.transpose() * M1TnnE_E * pot_E;
    L2T += penalty_factor * hT * (pot_T.transpose() * M1TTnnE * pot_T 
                                  - cross_PTnnE_PE - cross_PTnnE_PE.transpose() 
                                  + pot_E.transpose() * M1EE * pot_E);
  
    //
    // Second contribution of edge from term d_tE(PT nE.tE) + div PT .nE - DvE
    //
    // Values of P^{k-1}(T;Sym) at quad nodes
    boost::multi_array<MatrixGradient<2>, 2> basis_PSymkm1_quadE
      = evaluate_quad<Gradient>::compute(*cellBases(iT).PolySymkm1, quad_2km2_E);
    // Transformation of these values to obtain those of d_tE(PT nE.tE) + (div PT).nE
    auto trace_E = [&tE, &nE](const MatrixGradient<2> & G)->double
                              { double dtEnEtE = ( (G.dx*tE.x() + G.dy*tE.y() )*nE).dot(tE);
                                double divnE = G.dx(0,0)*nE.x() + G.dy(0,1)*nE.x() + G.dx(1,0)*nE.y() + G.dy(1,1)*nE.y();
                                return dtEnEtE + divnE;
                              };
    boost::multi_array<double, 2> basis_dtnEtE_plus_divnE
      = transform_values_quad<double>(basis_PSymkm1_quadE, trace_E);
    
    //  Matrix to return the unknown DvE from all cell unknowns
    Eigen::MatrixXd unknown_DvE = Eigen::MatrixXd::Zero(edgeBases(E).Polykm2->dimension(), dimensionCell(iT));
    unknown_DvE.block(0, localOffset(T, E)+edgeBases(E).Polykm3->dimension(), edgeBases(E).Polykm2->dimension(), edgeBases(E).Polykm2->dimension())
      = Eigen::MatrixXd::Identity(edgeBases(E).Polykm2->dimension(), edgeBases(E).Polykm2->dimension());

    // Contribution by developing products of (dtE(PT nE).tE+divPT.nE) - DE
    boost::multi_array<double, 2> basis_Pkm2_E_quad = evaluate_quad<Function>::compute(*edgeBases(E).Polykm2, quad_2km2_E);
    Eigen::MatrixXd M2TTtrace = compute_gram_matrix(basis_dtnEtE_plus_divnE, quad_2km2_E);
    Eigen::MatrixXd M2Ttrace_E = compute_gram_matrix(basis_dtnEtE_plus_divnE, basis_Pkm2_E_quad, quad_2km2_E);
    Eigen::MatrixXd M2EE = M1EE.topLeftCorner(edgeBases(E).Polykm2->dimension(), edgeBases(E).Polykm2->dimension());
    Eigen::MatrixXd cross_PTtrace_DvE = pot_T.transpose() * M2Ttrace_E * unknown_DvE;

    L2T += penalty_factor * std::pow(hT, 3) * ( pot_T.transpose() * M2TTtrace * pot_T
                                                - cross_PTtrace_DvE - cross_PTtrace_DvE.transpose()
                                                + unknown_DvE.transpose() * M2EE * unknown_DvE
                                              ); 
    
  }  

  // Stabilisation: vertices
  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);

    // Matrix to return the unknowns at V from all cell unknowns
    Eigen::MatrixXd unknown_V = Eigen::MatrixXd::Zero(m_symm_vertex.size(), dimensionCell(iT));
    unknown_V.block(0, localOffset(T, V), m_symm_vertex.size(), m_symm_vertex.size())
      = Eigen::MatrixXd::Identity(m_symm_vertex.size(), m_symm_vertex.size());
    
    // We create a quadrature rule of one node (V) to use evaluate_quad, which is much more efficient than manual evaluations
    QuadratureRule quad_V(1, QuadratureNode(V.coords().x(), V.coords().y(), 1.));
    boost::multi_array<MatrixRd, 2> basis_PSymkm1_quad_V 
        = evaluate_quad<Function>::compute(*cellBases(iT).PolySymkm1, quad_V);
    Eigen::MatrixXd MTT_V = compute_gram_matrix(basis_PSymkm1_quad_V, quad_V);

    Eigen::MatrixXd MTV_V = Eigen::MatrixXd::Zero(cellBases(iT).PolySymkm1->dimension(), m_symm_vertex.size());
    for (size_t i = 0; i < cellBases(iT).PolySymkm1->dimension(); i++) {
      for (size_t j = 0; j < m_symm_vertex.size(); j++) {
        MTV_V(i, j) = scalar_product(basis_PSymkm1_quad_V[i][0], m_symm_vertex.function(j));
      } // for j
    } // for i
    
    Eigen::MatrixXd MVV_V = Eigen::MatrixXd(m_symm_vertex.size(), m_symm_vertex.size());
    for (size_t i = 0; i < m_symm_vertex.size(); i++) {
      for (size_t j = 0; j < m_symm_vertex.size(); j++) {
        MVV_V(i, j) = scalar_product(m_symm_vertex.function(i), m_symm_vertex.function(j));
      } // for j
    } // for i

    Eigen::MatrixXd cross_PTtauV = pot_T.transpose() * MTV_V * unknown_V;
    L2T += penalty_factor * pow(hT, 2) * (
                                          pot_T.transpose() * MTT_V * pot_T
                                          - cross_PTtauV
                                          - cross_PTtauV.transpose()
                                          + unknown_V.transpose() * MVV_V * unknown_V
                                          );
  } // for iV
  
  // Consistent contribution
  //    compute basis of P^{k-1}(T;Sym) at quadrature nodes, and transform it with A.
  QuadratureRule quad_2km2_T = generate_quadrature_rule(T, 2*(degree()-1));
  auto basis_PSymkm1_quad = evaluate_quad<Function>::compute(*cellBases(iT).PolySymkm1, quad_2km2_T);
  auto basis_A_PSymkm1_quad = transform_values_quad<Eigen::Matrix2d>(basis_PSymkm1_quad, A);
  
  L2T += pot_T.transpose() 
            * compute_gram_matrix(basis_PSymkm1_quad, basis_A_PSymkm1_quad, quad_2km2_T) 
                * pot_T;

  return L2T;

}

