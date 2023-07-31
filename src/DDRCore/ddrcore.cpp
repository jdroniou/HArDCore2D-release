#include <cassert>

#include <ddrcore.hpp>
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_edge.hpp>

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
  
  MonomialCellIntegralsType int_monoT_2kp4 = IntegrateCellMonomials(T, 2*(m_K+2));
    
  //------------------------------------------------------------------------------
  // Basis for Pk+1(T), Pk(T), Pk-1(T) and Pk(T)^2
  //    We do not use 'RestrictedBasis' because the construction is quite 
  //    fast already with GramMatrix, and it simplifies the management
  //    of these bases afterwards
  //------------------------------------------------------------------------------
  
  MonomialScalarBasisCell basis_Pkpo_T(T, m_K + 1);
  bases_T.Polykpo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkpo_T, GramMatrix(T, basis_Pkpo_T, int_monoT_2kp4))) );  

  MonomialScalarBasisCell basis_Pk_T(T, m_K);
  bases_T.Polyk.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pk_T, GramMatrix(T, basis_Pk_T, int_monoT_2kp4))) );  

  // Check that we got the dimensions right
  assert( bases_T.Polykpo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K + 1) );
  assert( bases_T.Polyk->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K) );

  if (PolynomialSpaceDimension<Cell>::Poly(m_K - 1) > 0) {
    MonomialScalarBasisCell basis_Pkmo_T(T, m_K-1);
    bases_T.Polykmo.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkmo_T, GramMatrix(T, basis_Pkmo_T, int_monoT_2kp4))) );  
    assert( bases_T.Polykmo->dimension() == PolynomialSpaceDimension<Cell>::Poly(m_K - 1) );
  }
  
  //------------------------------------------------------------------------------
  // Basis Pk(T)^2
  //------------------------------------------------------------------------------

  bases_T.Polyk2.reset( new Poly2BasisCellType(*bases_T.Polyk) );

  // Check dimension
  assert( bases_T.Polyk2->dimension() == 2 * PolynomialSpaceDimension<Cell>::Poly(m_K) );
  
  //------------------------------------------------------------------------------
  // Basis for Rk-1(T)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Cell>::Roly(m_K - 1) > 0) {
    // Non-orthonormalised basis of Rk-1(T). 
    MonomialScalarBasisCell basis_Pk_T(T, m_K);
    ShiftedBasis<MonomialScalarBasisCell> basis_Pk0_T(basis_Pk_T,1);
    CurlBasis<ShiftedBasis<MonomialScalarBasisCell>> basis_Rkmo_T(basis_Pk0_T);
    // Orthonormalise, store and check dimension
    bases_T.Rolykmo.reset( new RolyBasisCellType(l2_orthonormalize(basis_Rkmo_T, GramMatrix(T, basis_Rkmo_T, int_monoT_2kp4))) );
    assert( bases_T.Rolykmo->dimension() == PolynomialSpaceDimension<Cell>::Roly(m_K - 1) );
  }
  
  //------------------------------------------------------------------------------
  // Basis for Rck(T)
  //------------------------------------------------------------------------------

  if (PolynomialSpaceDimension<Cell>::RolyCompl(m_K) > 0) {
    // Non-orthonormalised
    RolyComplBasisCell basis_Rck_T(T, m_K);
    // Orthonormalise, store and check dimension
    bases_T.RolyComplk.reset( new RolyComplBasisCellType(l2_orthonormalize(basis_Rck_T, GramMatrix(T, basis_Rck_T, int_monoT_2kp4))) );
    assert ( bases_T.RolyComplk->dimension() == PolynomialSpaceDimension<Cell>::RolyCompl(m_K) );
  }

  //------------------------------------------------------------------------------
  // Basis for Rck+2(T)
  //------------------------------------------------------------------------------

  // Non-orthonormalised
  RolyComplBasisCell basis_Rckp2_T(T, m_K+2);
  // Orthonormalise, store and check dimension
  bases_T.RolyComplkp2.reset( new RolyComplBasisCellType(l2_orthonormalize(basis_Rckp2_T, GramMatrix(T, basis_Rckp2_T, int_monoT_2kp4))) );
  assert ( bases_T.RolyComplkp2->dimension() == PolynomialSpaceDimension<Cell>::RolyCompl(m_K+2) );

  //------------------------------------------------------------------------------
  // Basis for Gck+2(T)
  //------------------------------------------------------------------------------

  // Non-orthonormalised
  GolyComplBasisCell basis_Gckp2_T(T, m_K+2);
  // Orthonormalise, store and check dimension
  bases_T.GolyComplkp2.reset( new GolyComplBasisCellType(l2_orthonormalize(basis_Gckp2_T, GramMatrix(T, basis_Gckp2_T, int_monoT_2kp4))) );
  assert ( bases_T.GolyComplkp2->dimension() == PolynomialSpaceDimension<Cell>::GolyCompl(m_K+2) );
  
  return bases_T;
}


//------------------------------------------------------------------------------

DDRCore::EdgeBases DDRCore::_construct_edge_bases(size_t iE)
{
  const Edge & E = *m_mesh.edge(iE);

  EdgeBases bases_E;

  MonomialEdgeIntegralsType int_monoE_2kp2 = IntegrateEdgeMonomials(E, 2*m_K+2);

  // Basis for Pk+1(E)
  MonomialScalarBasisEdge basis_Pkpo_E(E, m_K + 1);
  bases_E.Polykpo.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pkpo_E, GramMatrix(E, basis_Pkpo_E, int_monoE_2kp2))) );  

  // Basis for Pk(E)
  MonomialScalarBasisEdge basis_Pk_E(E, m_K);
  bases_E.Polyk.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pk_E, GramMatrix(E, basis_Pk_E, int_monoE_2kp2))) );  

  // Basis for Pk-1(E)
  if (PolynomialSpaceDimension<Edge>::Poly(m_K - 1) > 0) {
    MonomialScalarBasisEdge basis_Pkmo_E(E, m_K-1);
    bases_E.Polykmo.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pkmo_E, GramMatrix(E, basis_Pkmo_E, int_monoE_2kp2))) );  
  }

  return bases_E;
}
