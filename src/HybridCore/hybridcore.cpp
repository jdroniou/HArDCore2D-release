// Data structures and methods to implement hybrid schemes in 2D, with polynomial unknowns
// in the cells and on the edges, such as Hybrid High-order (HHO) schemes.
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "hybridcore.hpp"
#include "mesh.hpp"
#include <quad2d.hpp>
#include <quad1d.hpp>
#include <Eigen/Dense>
#include <iostream>

using namespace HArDCore2D;

//----------------------------------------------------------------------------
//                UVector
//----------------------------------------------------------------------------
UVector::UVector(const Eigen::VectorXd & values, const Mesh & mesh, const int cell_deg, const size_t edge_deg)
  : m_values(values),
    m_mesh(mesh), 
    m_cell_deg(cell_deg),
    m_edge_deg(edge_deg)
  {
    // Do nothing
  }

Eigen::VectorXd UVector::restr(size_t iT) const{
  Cell* cell = m_mesh.cell(iT);
  size_t nedgesT = cell->n_edges();

  Eigen::VectorXd XTF = Eigen::VectorXd::Zero(n_cell_dofs() + nedgesT * n_edge_dofs());

  XTF.head(n_cell_dofs()) = m_values.segment(iT * n_cell_dofs(), n_cell_dofs());
  for (size_t ilE = 0; ilE < nedgesT; ilE++){
    size_t offset_E = n_total_cell_dofs() + cell->edge(ilE)->global_index() * n_edge_dofs();
    XTF.segment(n_cell_dofs() + ilE * n_edge_dofs(), n_edge_dofs()) = m_values.segment(offset_E, n_edge_dofs());
  }

  return XTF;

}


//----------------------------------------------------------------------------
//                HybridCore
//----------------------------------------------------------------------------

// Creation class

HybridCore::HybridCore(const Mesh* mesh_ptr, const int cell_deg, const size_t edge_deg, const bool use_threads, std::ostream & output, const bool ortho)
  : m_mesh(mesh_ptr),
    m_cell_deg(cell_deg),
    m_cell_deg_pos(std::max(cell_deg,0)),
    m_edge_deg(edge_deg),
    m_use_threads(use_threads),
    m_output(output),
    m_ortho(ortho),
    m_cell_basis(mesh_ptr->n_cells()),
    m_edge_basis(mesh_ptr->n_edges())
    {
  m_output << "[HybridCore] Construction" << (m_use_threads ? " (multi-threading)" : "") << (m_ortho ? " (orthonormalised basis functions)" : "") << "\n";
  // Create cell bases
  std::function<void(size_t, size_t)> construct_all_cell_basis
    =[&] (size_t start, size_t end)->void 
    {  
      for (size_t iT = start; iT < end; iT++){
        m_cell_basis[iT].reset( new PolyCellBasisType(_construct_cell_basis(iT)) );
      }
    };
  parallel_for(mesh_ptr->n_cells(), construct_all_cell_basis, m_use_threads);

  // Create edge bases
  std::function<void(size_t, size_t)> construct_all_edge_basis
    =[&] (size_t start, size_t end)->void 
    {  
      for (size_t iF = start; iF < end; iF++){
        m_edge_basis[iF].reset( new PolyEdgeBasisType(_construct_edge_basis(iF)) );
      }
    };
  parallel_for(mesh_ptr->n_edges(), construct_all_edge_basis, m_use_threads);

}

// -------------------------------------------------------------
// ------- Construction of cell and edge basis functions
// -------------------------------------------------------------


HybridCore::PolyCellBasisType HybridCore::_construct_cell_basis(size_t iT)
{
  const Cell & T = *m_mesh->cell(iT);

  // Basis of monomials
  MonomialScalarBasisCell basis_mono_T(T, m_cell_deg_pos);

  // We construct the non-orthonormalised basis (same as above but as a "Family")
  Eigen::MatrixXd B = Eigen::MatrixXd::Identity(basis_mono_T.dimension(), basis_mono_T.dimension());
  PolyCellBasisType CellBasis(basis_mono_T, B);
  
  // Orthonormalisation
  if (m_ortho){
    QuadratureRule quadT = generate_quadrature_rule(T, 2 * m_cell_deg_pos);
    auto basis_mono_T_quadT = evaluate_quad<Function>::compute(basis_mono_T, quadT);
    CellBasis = l2_orthonormalize(basis_mono_T, quadT, basis_mono_T_quadT);
  }

  return CellBasis;
}

HybridCore::PolyEdgeBasisType HybridCore::_construct_edge_basis(size_t iE)
{
  const Edge & E = *m_mesh->edge(iE);

  // Basis of monomials
  MonomialScalarBasisEdge basis_mono_E(E, m_edge_deg);

  // We construct the non-orthonormalised basis (same as above but as a "Family")
  Eigen::MatrixXd B = Eigen::MatrixXd::Identity(basis_mono_E.dimension(), basis_mono_E.dimension());
  PolyEdgeBasisType EdgeBasis(basis_mono_E, B);

  // Orthonormalisation
  if (m_ortho){
    QuadratureRule quadE = generate_quadrature_rule(E, 2 * m_edge_deg);
    auto basis_mono_E_quadE = evaluate_quad<Function>::compute(basis_mono_E, quadE);
    EdgeBasis = l2_orthonormalize(basis_mono_E, quadE, basis_mono_E_quadE);
  }

  return EdgeBasis;
}


// -----------------------------------------------------------------
// ------- Weights to eliminate cell unknows if L=-1
// -----------------------------------------------------------------

Eigen::VectorXd HybridCore::compute_weights(size_t iT) const {
  Cell* iCell = m_mesh->cell(iT);
  size_t nlocal_edges = iCell->n_edges();
  Eigen::VectorXd barycoefT = Eigen::VectorXd::Zero(nlocal_edges);

  // Rule degree 0: all coefficients identical
  //    barycoefT = (1/double(nlocal_edges)) * Eigen::VectorXd::Ones(nlocal_edges);

  // Rule degree 1: coefficient is |F|d_{TF}/(d|T|)
  for (size_t ilF=0; ilF < nlocal_edges; ilF++){
    VectorRd normalTF = iCell->edge_normal(ilF);
    VectorRd xFxT = iCell->edge(ilF)->center_mass() - iCell->center_mass();

    double dTF = xFxT.dot(normalTF); 

    barycoefT(ilF) = iCell->edge(ilF)->measure() * dTF / (m_mesh->dim() * iCell->measure());

  }

  return barycoefT;
}

// ----------------------------------------------------------
// ---------- Norms of discrete unknowns
// ----------------------------------------------------------

double HybridCore::L2norm(const UVector &Xh) const {
  size_t n_cell_dofs = Xh.n_cell_dofs();
  Eigen::ArrayXd cell_norm = Eigen::ArrayXd::Zero(m_mesh->n_cells());

  std::function<void(size_t, size_t)> compute_local_norms
      = [&](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          // L2 norm computed using the mass matrix
          // Compute cell quadrature nodes and values of cell basis functions at these nodes
          Cell* cell = m_mesh->cell(iT);
          QuadratureRule quadT = generate_quadrature_rule(*cell, 2*m_cell_deg);
          boost::multi_array<double, 2> phiT_quadT = evaluate_quad<Function>::compute(CellBasis(iT), quadT);

          Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, n_cell_dofs, n_cell_dofs, "sym");
          Eigen::VectorXd XT = Xh.asVectorXd().segment(iT*n_cell_dofs, n_cell_dofs);

          cell_norm[iT] = XT.dot(MTT*XT);
        }
      };

  parallel_for(m_mesh->n_cells(), compute_local_norms, m_use_threads);

  return sqrt( cell_norm.sum() );
}

double HybridCore::H1norm(const UVector &Xh) const {

  Eigen::ArrayXd cell_norm = Eigen::ArrayXd::Zero(m_mesh->n_cells());
  size_t n_cell_dofs = Xh.n_cell_dofs();
  size_t n_edge_dofs = Xh.n_edge_dofs();
  int cell_deg = Xh.get_cell_deg();
  int edge_deg = Xh.get_edge_deg();

  std::function<void(size_t, size_t)> compute_local_norms
      = [&](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          Cell* cell = m_mesh->cell(iT);
          size_t nedgesT = cell->n_edges();
          size_t ndofs_T = n_cell_dofs + nedgesT * n_edge_dofs;

          // Local matrix of the bilinear form corresponding to the cell contribution in the discrete H1 norm
          Eigen::MatrixXd H1aT = Eigen::MatrixXd::Zero(ndofs_T, ndofs_T);

          // Compute cell quadrature nodes and values of gradients of cell basis functions
          QuadratureRule quadT = generate_quadrature_rule(*cell, 2*(cell_deg-1));
          boost::multi_array<VectorRd, 2> dphiT_quadT = evaluate_quad<Gradient>::compute(CellBasis(iT), quadT);

          // CELL CONTRIBUTION
          //
          // Stiffness matrix
          H1aT.topLeftCorner(n_cell_dofs, n_cell_dofs) = 
              compute_gram_matrix(dphiT_quadT, dphiT_quadT, quadT, n_cell_dofs, n_cell_dofs, "sym");


          // EDGES CONTRIBUTION               
          for (size_t ilF=0; ilF < nedgesT ; ilF++){
            Edge* edge = cell->edge(ilF);
            size_t iF = edge->global_index();
            size_t offset_F = n_cell_dofs + ilF * n_edge_dofs;
            // Different scalings for boundary terms
            //                double dTF = m_mesh->edge(iF)->diam();
            double dTF = cell->measure() / m_mesh->edge(iF)->measure();
            // Edge quadrature nodes and values of cell and edge basis functions at these nodes
            size_t doeF = std::max(2*edge_deg, 2*cell_deg);
            QuadratureRule quadF = generate_quadrature_rule(*edge, doeF);
            boost::multi_array<double, 2> phiT_quadF = evaluate_quad<Function>::compute(CellBasis(iT), quadF);
            boost::multi_array<double, 2> phiF_quadF = evaluate_quad<Function>::compute(EdgeBasis(iF), quadF);

            // Edge, edge-cell and edge-edge Gram matrix on F
            Eigen::MatrixXd MFF = compute_gram_matrix(phiF_quadF, phiF_quadF, quadF, n_edge_dofs, n_edge_dofs, "sym");
            Eigen::MatrixXd MFT = compute_gram_matrix(phiF_quadF, phiT_quadF, quadF, n_edge_dofs, n_cell_dofs, "nonsym");
            Eigen::MatrixXd MTT_on_F = compute_gram_matrix(phiT_quadF, phiT_quadF, quadF, n_cell_dofs, n_cell_dofs, "sym");

            // Contribution of the edge to the local bilinear form
            H1aT.block(offset_F, offset_F, n_edge_dofs, n_edge_dofs) += MFF / dTF;
            H1aT.block(offset_F, 0, n_edge_dofs, n_cell_dofs) -= MFT / dTF;
            H1aT.block(0, offset_F, n_cell_dofs, n_edge_dofs) -= MFT.transpose() / dTF;
            H1aT.topLeftCorner(n_cell_dofs, n_cell_dofs) += MTT_on_F / dTF;
          }

          Eigen::VectorXd XTF = Xh.restr(iT);
          cell_norm[iT] = XTF.transpose() * H1aT * XTF;
        }
       };

  parallel_for(m_mesh->n_cells(), compute_local_norms, m_use_threads);

  return sqrt( cell_norm.sum() );
}


// -----------------------------------------------------------
// --------- Evaluate discrete functions in cell or edges
// -----------------------------------------------------------

double HybridCore::evaluate_in_cell(const UVector Xh, size_t iT, VectorRd x) const {
  double value = 0.0;
  const auto& T_basis = CellBasis(iT);
  for (size_t i = 0; i < Xh.n_cell_dofs(); i++) {
    const size_t index = i + iT * Xh.n_cell_dofs();
    value += Xh(index) * T_basis.function(i, x);

  }
  return value;
}

double HybridCore::evaluate_in_edge(const UVector Xh, size_t iF, VectorRd x) const {
  double value = 0.0;
  const auto& F_basis = EdgeBasis(iF);
  for (size_t i = 0; i < Xh.n_edge_dofs(); i++) {
    const size_t index = Xh.n_total_cell_dofs() + iF * Xh.n_edge_dofs() + i;
    value += Xh(index) * F_basis.function(i, x);

  }
  return value;
}


// ---------------------------------------------------------------
// ------- Compute vertex values from a hybrid function
// ---------------------------------------------------------------

Eigen::VectorXd HybridCore::VertexValues(const UVector Xh, const std::string from_dofs) {
  Eigen::VectorXd function_vertex = Eigen::VectorXd::Zero(m_mesh->n_vertices());

  if (from_dofs == "cell"){
    // compute from the cell polynomials, averaging all values coming from the cells around each vertex
    for (size_t iV = 0; iV < m_mesh->n_vertices(); iV++){
      auto xV = m_mesh->vertex(iV)->coords();
      auto cList = m_mesh->vertex(iV)->get_cells();
      auto weight = cList.size();
      for (size_t ilT = 0; ilT < weight; ilT++){
        auto temp = this->evaluate_in_cell(Xh, cList[ilT]->global_index(), xV);
        function_vertex(iV) += temp;
      }
      function_vertex(iV) = function_vertex(iV)/(weight);
    }

  }else{
    // compute from the edge polynomials, averaging all values coming from the edges around each vertex
    for (size_t iV = 0; iV < m_mesh->n_vertices(); iV++){
      auto xV = m_mesh->vertex(iV)->coords();
      auto eList = m_mesh->vertex(iV)->get_edges();
      double weight = eList.size();
      for (size_t ilF = 0; ilF < weight; ilF++){
        auto temp = this->evaluate_in_edge(Xh, eList[ilF]->global_index(), xV);
        function_vertex(iV) += temp;
      }
      function_vertex(iV) = function_vertex(iV)/(weight);
    }
  }

  return function_vertex;
}


