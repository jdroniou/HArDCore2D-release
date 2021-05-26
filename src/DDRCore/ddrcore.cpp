#include <cassert>

#include <ddrcore.hpp>
#include <parallel_for.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------

DDRCore::DDRCore(const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
  : m_mesh(mesh),
    m_K(K),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_edge_bases(mesh.n_edges())
{
  m_output << "[DDRCore] Initializing" << std::endl;
  
  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[DDRCore] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct edge bases
  std::function<void(size_t, size_t)> construct_all_edge_bases   
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        this->m_edge_bases[iE].reset( new EdgeBases(_construct_edge_bases(iE)) );
	      } // for iF
      };
  
  m_output << "[DDRCore] Constructing edge bases" << std::endl;
  parallel_for(mesh.n_edges(), construct_all_edge_bases, use_threads);
}

//------------------------------------------------------------------------------

DDRCore::CellBases DDRCore::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;
  
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T)
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  QuadratureRule quad_2kpo_T = generate_quadrature_rule(T, 2 * (m_K + 1));
  boost::multi_array<double, 2> on_basis_Pkpo_T_quad = evaluate_quad<Function>::compute(basis_Pkpo_T, quad_2kpo_T);
  // Orthonormalize and store
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, quad_2kpo_T, on_basis_Pkpo_T_quad)) );   
  // Check that we got the dimension right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );

  //------------------------------------------------------------------------------
  // Basis for Pk(T), Pk-1(T) and Pk(T)^2
  //------------------------------------------------------------------------------

  // Given that the basis for Pk+1(T) is hierarchical, bases for Pk(T) and
  // Pk-1(T) can be obtained by restricting the former
  bases_T.Polyk.reset( new RestrictedBasis<PolyBasisCellType>(*bases_T.Polykpo, PolynomialSpaceDimension<Cell>::Poly(m_K)) );  
  bases_T.Polyk2.reset( new Poly2BasisCellType(*bases_T.Polyk) );
  if (PolynomialSpaceDimension<Cell>::Poly(m_K - 1) > 0) {
    bases_T.Polykmo.reset( new RestrictedBasis<PolyBasisCellType>(*bases_T.Polykpo, PolynomialSpaceDimension<Cell>::Poly(m_K - 1)) );
  }
  // Check dimension Pk(T)^2
  assert( bases_T.Polyk2->dimension() == 2 * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Rk-1(T)
  //------------------------------------------------------------------------------

  // Quadrature useful for various spaces to follow (degree might be too large in certain cases, but that is
  // not a major additional cost in 2D)
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_K);

  if (PolynomialSpaceDimension<Cell>::Roly(m_K - 1) > 0) {
    // Non-orthonormalised basis of Rk-1(T). 
    MonomialScalarBasisCell basis_Pk_T(T, m_K);
    ShiftedBasis<MonomialScalarBasisCell> basis_Pk0_T(basis_Pk_T,1);
    CurlBasis<ShiftedBasis<MonomialScalarBasisCell>> basis_Rkmo_T(basis_Pk0_T);
    // Orthonormalise, store and check dimension
    auto basis_Rkmo_T_quad = evaluate_quad<Function>::compute(basis_Rkmo_T, quad_2k_T);
    bases_T.Rolykmo.reset( new RolyBasisCellType(l2_orthonormalize(basis_Rkmo_T, quad_2k_T, basis_Rkmo_T_quad)) );
    assert( bases_T.Rolykmo->dimension() == PolynomialSpaceDimension<Cell>::Roly(m_K - 1) );
  }
  
  //------------------------------------------------------------------------------
  // Basis for Rck(T)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Cell>::RolyCompl(m_K) > 0) {
    // Non-orthonormalised
    RolyComplBasisCell basis_Rck_T(T, m_K);
    auto basis_Rck_T_quad = evaluate_quad<Function>::compute(basis_Rck_T, quad_2k_T);
    // Orthonormalise, store and check dimension
    bases_T.RolyComplk.reset( new RolyComplBasisCellType(l2_orthonormalize(basis_Rck_T, quad_2k_T, basis_Rck_T_quad)) );
    assert ( bases_T.RolyComplk->dimension() == PolynomialSpaceDimension<Cell>::RolyCompl(m_K) );
  }

  //------------------------------------------------------------------------------
  // Basis for Rck+2(T)
  //------------------------------------------------------------------------------

  // Non-orthonormalised
  RolyComplBasisCell basis_Rckp2_T(T, m_K+2);
  QuadratureRule quad_2kp2_T = generate_quadrature_rule(T, 2 * (m_K+2) );
  auto basis_Rckp2_T_quad = evaluate_quad<Function>::compute(basis_Rckp2_T, quad_2kp2_T);
  // Orthonormalise, store and check dimension
  bases_T.RolyComplkp2.reset( new RolyComplBasisCellType(l2_orthonormalize(basis_Rckp2_T, quad_2kp2_T, basis_Rckp2_T_quad)) );
  assert ( bases_T.RolyComplkp2->dimension() == PolynomialSpaceDimension<Cell>::RolyCompl(m_K+2) );

  
  return bases_T;
}


//------------------------------------------------------------------------------

DDRCore::EdgeBases DDRCore::_construct_edge_bases(size_t iE)
{
  const Edge & E = *m_mesh.edge(iE);

  EdgeBases bases_E;

  // Basis for Pk+1(E)
  MonomialScalarBasisEdge basis_Pkpo_E(E, m_K + 1);
  QuadratureRule quad_2kpo_E = generate_quadrature_rule(E, 2 * (m_K + 1));
  auto basis_Pkpo_E_quad = evaluate_quad<Function>::compute(basis_Pkpo_E, quad_2kpo_E);
  bases_E.Polykpo.reset( new PolyEdgeBasisType(l2_orthonormalize(basis_Pkpo_E, quad_2kpo_E, basis_Pkpo_E_quad)) );

  // Basis for Pk(E)
  bases_E.Polyk.reset( new RestrictedBasis<PolyEdgeBasisType>(*bases_E.Polykpo, PolynomialSpaceDimension<Edge>::Poly(m_K)) );
  
  // Basis for Pk-1(E)
  if (PolynomialSpaceDimension<Edge>::Poly(m_K - 1) > 0) {
    // Given that the basis for Pk+1(E) is hierarchical, a basis for Pk-1(E)
    // can be obtained by restricting the former
    bases_E.Polykmo.reset( new RestrictedBasis<PolyEdgeBasisType>(*bases_E.Polykpo, PolynomialSpaceDimension<Edge>::Poly(m_K - 1)) );
  }

  return bases_E;
}
