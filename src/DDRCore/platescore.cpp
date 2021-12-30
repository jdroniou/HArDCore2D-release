#include <cassert>

#include <platescore.hpp>
#include <parallel_for.hpp>
#include <quadraturerule.hpp>
#include <GMpoly_cell.hpp>
#include <GMpoly_edge.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------

PlatesCore::PlatesCore(const Mesh & mesh, size_t K, bool use_threads, std::ostream & output)
  : m_mesh(mesh),
    m_K(K),
    m_output(output),
    m_cell_bases(mesh.n_cells()),
    m_edge_bases(mesh.n_edges())
{
  if (K < 3) {
    std::cerr << "Plates complex requires k >= 3" << std::endl;
    exit(1);
  }
  
  m_output << "[PlatesCore] Initializing" << std::endl;
  
  // Construct element bases
  std::function<void(size_t, size_t)> construct_all_cell_bases
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iT = start; iT < end; iT++) {
	        this->m_cell_bases[iT].reset( new CellBases(this->_construct_cell_bases(iT)) );
	      } // for iT
      };

  m_output << "[PlatesCore] Constructing element bases" << std::endl;
  parallel_for(mesh.n_cells(), construct_all_cell_bases, use_threads);
  
  // Construct edge bases
  std::function<void(size_t, size_t)> construct_all_edge_bases   
    = [this](size_t start, size_t end)->void
      {
	      for (size_t iE = start; iE < end; iE++) {
	        this->m_edge_bases[iE].reset( new EdgeBases(_construct_edge_bases(iE)) );
	      } // for iF
      };
  
  m_output << "[PlatesCore] Constructing edge bases" << std::endl;
  parallel_for(mesh.n_edges(), construct_all_edge_bases, use_threads);
}

//------------------------------------------------------------------------------

PlatesCore::CellBases PlatesCore::_construct_cell_bases(size_t iT)
{
  const Cell & T = *m_mesh.cell(iT);

  CellBases bases_T;

  MonomialCellIntegralsType int_monoT_2kp2 = IntegrateCellMonomials(T, 2 * (m_K + 1));

  // Basis for Pk+1 and Pk-2
  MonomialScalarBasisCell basis_Pkp1_T(T, m_K + 1);
  bases_T.Polykp1.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkp1_T, GramMatrix(T, basis_Pkp1_T, int_monoT_2kp2))) );
  MonomialScalarBasisCell basis_Pkm2_T(T, m_K - 2);
  bases_T.Polykm2.reset( new PolyBasisCellType(l2_orthonormalize(basis_Pkm2_T, GramMatrix(T, basis_Pkm2_T, int_monoT_2kp2))) );

  // Basis for (Pk-1)^2
  MonomialScalarBasisCell basis_Pkm1_T(T, m_K - 1);
  PolyBasisCellType on_basis_Pkm1_T = l2_orthonormalize(basis_Pkm1_T, GramMatrix(T, basis_Pkm1_T, int_monoT_2kp2));
  bases_T.Poly2km1.reset( new Poly2BasisCellType(on_basis_Pkm1_T) );

  // Basis for Pk-1(Sym): obtained as a family of the MatrixFamily of all 2x2 matrices
  MatrixFamily<PolyBasisCellType, 2> basis_Pkm1_2x2_T(on_basis_Pkm1_T);
  bases_T.PolySymkm1.reset( new PolySymBasisCellType(basis_Pkm1_2x2_T, basis_Pkm1_2x2_T.symmetricBasis()) );  

  // Basis for Hk-4
  if (m_K >= 4) {
    QuadratureRule qr_2km8 = generate_quadrature_rule(T, 2 * (m_K - 4));
    HessianBasis<ShiftedBasis<MonomialScalarBasisCell> > basis_Hkm4_T(ShiftedBasis<MonomialScalarBasisCell>(MonomialScalarBasisCell(T, m_K - 2), PolynomialSpaceDimension<Cell>::Poly(1)));
    auto basis_Hkm4_T_quad = evaluate_quad<Function>::compute(basis_Hkm4_T, qr_2km8);
    bases_T.Holykm4.reset( new HolyBasisCellType(l2_orthonormalize(basis_Hkm4_T, qr_2km8, basis_Hkm4_T_quad)) );
  }
  
  // Basis for Hck-1
  QuadratureRule qr_2km2 = generate_quadrature_rule(T, 2 * (m_K - 1));  
  HolyComplBasisCell basis_Hckm1_T(T, m_K - 1);
  auto basis_Hckm1_T_quad = evaluate_quad<Function>::compute(basis_Hckm1_T, qr_2km2);
  bases_T.HolyComplkm1.reset( new HolyComplBasisCellType(l2_orthonormalize(basis_Hckm1_T, qr_2km2, basis_Hckm1_T_quad)) );
    
  return bases_T;
}

//------------------------------------------------------------------------------

PlatesCore::EdgeBases PlatesCore::_construct_edge_bases(size_t iE)
{
  const Edge & E = *m_mesh.edge(iE);

  EdgeBases bases_E;

  MonomialEdgeIntegralsType int_monoE_2km2 = IntegrateEdgeMonomials(E, 2 * (m_K - 1));

  MonomialScalarBasisEdge basis_Pkm1_E(E, m_K -1 );
  bases_E.Polykm1.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pkm1_E, GramMatrix(E, basis_Pkm1_E, int_monoE_2km2))) );
  
  MonomialScalarBasisEdge basis_Pkm2_E(E, m_K - 2);
  bases_E.Polykm2.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pkm2_E, GramMatrix(E, basis_Pkm2_E, int_monoE_2km2))) );

  MonomialScalarBasisEdge basis_Pkm3_E(E, m_K - 3);
  bases_E.Polykm3.reset( new PolyBasisEdgeType(l2_orthonormalize(basis_Pkm3_E, GramMatrix(E, basis_Pkm3_E, int_monoE_2km2))) );

  return bases_E;
}
