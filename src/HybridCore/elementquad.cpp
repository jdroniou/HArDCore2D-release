// Class to create and store values of cell and edge basis functions on quadrature points
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "elementquad.hpp"
#include <string>
#include <mesh.hpp>
#include <quadraturerule.hpp>

using namespace HArDCore2D;

// Creation class

ElementQuad::ElementQuad(const HybridCore& hcore, const size_t iT, const size_t doeT, const size_t doeF )
  : m_hcore(hcore),
    m_iT(iT),
    m_doeT(doeT),
    m_doeF(doeF),
    m_T(*hcore.get_mesh()->cell(iT)) {
  // Create cell quadratures and values of basis functions and gradients
  m_quadT = generate_quadrature_rule(m_T, m_doeT);
  const auto& basisT = m_hcore.CellBasis(m_iT);

  m_phiT_quadT.resize( boost::extents[basisT.dimension()][m_quadT.size()] );
  m_phiT_quadT = evaluate_quad<Function>::compute(basisT, m_quadT);

  m_dphiT_quadT.resize( boost::extents[basisT.dimension()][m_quadT.size()] );
  m_dphiT_quadT = evaluate_quad<Gradient>::compute(basisT, m_quadT);

  // Create edge quadrature and values of basis/gradient at nodes
  size_t nedges = m_T.n_edges();

  m_quadF.resize(nedges);
  m_phiT_quadF.resize(nedges);
  m_phiF_quadF.resize(nedges);
  m_dphiT_quadF.resize(nedges);
  for (size_t ilF = 0; ilF < nedges ; ilF++){
    Edge* edge = m_T.edge(ilF);
    const size_t iF = edge->global_index();
    const auto& basisF = m_hcore.EdgeBasis(iF);
    m_quadF[ilF] = generate_quadrature_rule(*edge, m_doeF);
    size_t nqF = m_quadF[ilF].size();
    m_phiT_quadF[ilF].resize( boost::extents[basisT.dimension()][nqF] ); 
    m_phiT_quadF[ilF] = evaluate_quad<Function>::compute(basisT, m_quadF[ilF]);

    m_phiF_quadF[ilF].resize( boost::extents[basisF.dimension()][nqF] );
    m_phiF_quadF[ilF] = evaluate_quad<Function>::compute(basisF, m_quadF[ilF]);

    m_dphiT_quadF[ilF].resize( boost::extents[basisT.dimension()][nqF] );
    m_dphiT_quadF[ilF] = evaluate_quad<Gradient>::compute(basisT, m_quadF[ilF]);
  }

}

//------------------------------------------------------------------------------
// Class methods
//------------------------------------------------------------------------------

boost::multi_array<VectorRd, 2> ElementQuad::get_vec_phiT_quadT(size_t degree) const
{
  const size_t n_scalar_basis = DimPoly<Cell>(degree);
  assert(n_scalar_basis <= m_phiT_quadT.shape()[0]);

  const size_t dim = m_hcore.get_mesh()->dim(); 
  const size_t nqT = m_quadT.size();
  boost::multi_array<VectorRd, 2> m_vec_phiT_quadT;
  m_vec_phiT_quadT.resize( boost::extents[dim * n_scalar_basis][nqT] );

  for (size_t k = 0; k < dim; k++) {
    VectorRd ek = VectorRd::Zero();
    ek(k) = 1.;
    for (size_t i = 0; i < n_scalar_basis; i++) {
      for (size_t iqn = 0; iqn < nqT; iqn++){
        m_vec_phiT_quadT[k * n_scalar_basis + i][iqn] = m_phiT_quadT[i][iqn] * ek;
      } // for iqn  
    } // for i
  } // for k

  return m_vec_phiT_quadT;
}


boost::multi_array<VectorRd, 2> ElementQuad::get_vec_phiT_quadF(size_t ilF, size_t degree) const
{
  const size_t n_scalar_basis = DimPoly<Cell>(degree);
  assert(n_scalar_basis <= m_phiT_quadF[ilF].shape()[0]);

  const size_t dim = m_hcore.get_mesh()->dim();
  const size_t nqF = m_quadF[ilF].size();
  boost::multi_array<VectorRd, 2> m_vec_phiT_quadF;
  m_vec_phiT_quadF.resize( boost::extents[dim * n_scalar_basis][nqF] );

  for (size_t k = 0; k < dim; k++) {
    VectorRd ek = VectorRd::Zero();
    ek(k) = 1.;
    for (size_t i = 0; i < n_scalar_basis; i++) {
      for (size_t iqn = 0; iqn < nqF; iqn++){
        m_vec_phiT_quadF[k * n_scalar_basis + i][iqn] = m_phiT_quadF[ilF][i][iqn] * ek;
      } // for iqn
    } // for i
  } // for k

  return m_vec_phiT_quadF;
}


boost::multi_array<VectorRd, 2> ElementQuad::get_vec_phiF_quadF(size_t ilF, size_t degree) const
{
  const size_t n_scalar_basis = DimPoly<Edge>(degree);
  assert(n_scalar_basis <= m_phiF_quadF[ilF].shape()[0]);

  const size_t dim = m_hcore.get_mesh()->dim();
  const size_t nqF = m_quadF[ilF].size();
  boost::multi_array<VectorRd, 2> m_vec_phiF_quadF;
  m_vec_phiF_quadF.resize( boost::extents[dim * n_scalar_basis][nqF] ); 

  for (size_t k = 0; k < dim; k++) {
    VectorRd ek = VectorRd::Zero();
    ek(k) = 1.;
    for (size_t i = 0; i < n_scalar_basis; i++) {
      for (size_t iqn = 0; iqn < nqF; iqn++){
        m_vec_phiF_quadF[k * n_scalar_basis + i][iqn] = m_phiF_quadF[ilF][i][iqn] * ek;
      } // for iqn
    } // for i
  } // for k

  return m_vec_phiF_quadF;
}



