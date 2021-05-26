#include "dofspace.hpp"

using namespace HArDCore2D;

DOFSpace::DOFSpace(
		               const Mesh & mesh,
		               size_t n_local_vertex_dofs,
                   size_t n_local_edge_dofs,
                   size_t n_local_cell_dofs		   
                   )
  : m_mesh(mesh),
    m_n_local_vertex_dofs(n_local_vertex_dofs),
    m_n_local_edge_dofs(n_local_edge_dofs),
    m_n_local_cell_dofs(n_local_cell_dofs)
{
  // Do nothing
}
