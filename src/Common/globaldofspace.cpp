#include "globaldofspace.hpp"

using namespace HArDCore2D;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------

GlobalDOFSpace::GlobalDOFSpace(
                   const Mesh & mesh,
                   size_t n_local_vertex_dofs,
                   size_t n_local_edge_dofs,
                   size_t n_local_cell_dofs		   
                   )
  : LocalDOFSpace(mesh, n_local_vertex_dofs, n_local_edge_dofs, n_local_cell_dofs),
    m_labelDOF(dimension(),0)
{
  // Do nothing
}

//------------------------------------------------------------------------------
// Restrictions
//------------------------------------------------------------------------------

Eigen::VectorXd GlobalDOFSpace::restrictEdge(size_t iE, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vE = Eigen::VectorXd::Zero(dimensionEdge(iE));
  const Edge & E = *m_mesh.edge(iE);

  if (m_n_local_vertex_dofs > 0) {
    vE.head(m_n_local_vertex_dofs)
      = vh.segment(globalOffset(*E.vertex(0)), m_n_local_vertex_dofs);
    vE.segment(m_n_local_vertex_dofs, m_n_local_vertex_dofs)
      = vh.segment(globalOffset(*E.vertex(1)), m_n_local_vertex_dofs);
  }

  if (m_n_local_edge_dofs > 0) {
    vE.tail(m_n_local_edge_dofs) = vh.segment(globalOffset(E), m_n_local_edge_dofs);
  }
  
  return vE;
}

//------------------------------------------------------------------------------

Eigen::VectorXd GlobalDOFSpace::restrictCell(size_t iT, const Eigen::VectorXd & vh) const
{
  Eigen::VectorXd vT = Eigen::VectorXd::Zero(dimensionCell(iT));
  const Cell & T = *m_mesh.cell(iT);

  if (m_n_local_vertex_dofs > 0) {
    for (size_t iV = 0; iV < T.n_vertices(); iV++) {
      const Vertex & V = *T.vertex(iV);
      vT.segment(localOffset(T, V), m_n_local_vertex_dofs)
        = vh.segment(globalOffset(V), m_n_local_vertex_dofs);
    } // for iV
  }

  if (m_n_local_edge_dofs > 0) {
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
      const Edge & E = *T.edge(iE);      
      vT.segment(localOffset(T, E), m_n_local_edge_dofs)
        = vh.segment(globalOffset(E), m_n_local_edge_dofs);
    } // for iE
  }

  vT.tail(m_n_local_cell_dofs)
    = vh.segment(globalOffset(T), m_n_local_cell_dofs);
  
  return vT;
}

//------------------------------------------------------------------------------

Eigen::MatrixXd GlobalDOFSpace::extendOperator(const Cell & T, const Edge & E, const Eigen::MatrixXd & opE) const
{
  Eigen::MatrixXd opT = Eigen::MatrixXd::Zero(opE.rows(), dimensionCell(T));

  // Vertex DOFs
  if (m_n_local_vertex_dofs > 0) {
    for (size_t iV = 0; iV < 2; iV++) { 
      const Vertex & V = *E.vertex(iV);     
      opT.block(0, localOffset(T, V), opE.rows(), m_n_local_vertex_dofs)
        = opE.block(0, localOffset(E, V), opE.rows(), m_n_local_vertex_dofs);
    } // for iV
  }
  
  // Edge DOFs
  if (m_n_local_edge_dofs > 0) {
      opT.block(0, localOffset(T, E), opE.rows(), m_n_local_edge_dofs)
        = opE.block(0, localOffset(E), opE.rows(), m_n_local_edge_dofs);
  }
  
  return opT;
}
//------------------------------------------------------------------------------

std::vector<size_t> GlobalDOFSpace::globalDOFIndices(const Cell & T) const
{
  std::vector<size_t> I(dimensionCell(T));

  size_t dof_index = 0;
  
  // Vertex DOFs
  if (m_n_local_vertex_dofs) {
    for (size_t iV = 0; iV < T.n_vertices(); iV++) {
      size_t offset_V = globalOffset(*T.vertex(iV));
      for (size_t i = 0; i < m_n_local_vertex_dofs; i++, dof_index++) {
        I[dof_index] = offset_V + i;
      } // for i
    } // for iV
  }

  // Edge DOFs
  if (m_n_local_edge_dofs) {
    for (size_t iE = 0; iE < T.n_edges(); iE++) {
      size_t offset_E = globalOffset(*T.edge(iE));
      for (size_t i = 0; i < m_n_local_edge_dofs; i++, dof_index++) {
        I[dof_index] = offset_E + i;
      } // for i
    } // for iE
  }

  // Cell DOFs
  if (m_n_local_cell_dofs) {
    size_t offset_T = globalOffset(T);
    for (size_t i = 0; i < m_n_local_cell_dofs; i++, dof_index++) {
      I[dof_index] = offset_T + i;
    } // for i
  }

  assert( dimensionCell(T) == dof_index );
  
  return I;
}


