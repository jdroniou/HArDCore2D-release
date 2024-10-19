#ifndef LOCALDOFSPACE_HPP
#define LOCALDOFSPACE_HPP

#include <mesh.hpp>
//#include <vertex.hpp>
//#include <edge.hpp>
//#include <cell.hpp>

namespace HArDCore2D {

  /*!
   *	\addtogroup Common
   * @{
   */
  
  /// Base class for local DOF spaces: functions to access local DOFs (organised from the smallest dimension to the largest) associated with each geometric entity.
  /** In a cell T, for example, the DOFs are presented in the order: DOFs of vertices of T, DOFs of edges of T, DOFs of T. */
  class LocalDOFSpace {
  public:
    /// Constructor
    LocalDOFSpace(
                  const Mesh & mesh,
                  size_t n_local_vertex_dofs,
                  size_t n_local_edge_dofs,
                  size_t n_local_cell_dofs
                  );

    //------------------------------------------------------------------------------
    // Accessors
    //------------------------------------------------------------------------------
    
    /// Returns the mesh
    const Mesh & mesh() const
    {
      return m_mesh;
    }

    /// Returns the number of local vertex DOFs
    inline size_t numLocalDofsVertex() const
    {
      return m_n_local_vertex_dofs;
    }
    
    /// Returns the number of local edge DOFs
    inline size_t numLocalDofsEdge() const
    {
      return m_n_local_edge_dofs;
    }
        
    /// Returns the number of local cell DOFs
    inline size_t numLocalDofsCell() const
    {
      return m_n_local_cell_dofs;
    }
    
    //------------------------------------------------------------------------------
    // Dimensions
    //------------------------------------------------------------------------------
    
    /// Returns the dimension of the global space (all DOFs for all geometric entities)
    inline size_t dimension() const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs
        + m_mesh.n_edges() * m_n_local_edge_dofs
        + m_mesh.n_cells() * m_n_local_cell_dofs;
    }

    /// Returns the total number of vertex DOFs
    inline size_t numTotalDofsVertices() const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs;
    }

    /// Returns the total number of vertex DOFs
    inline size_t numTotalDofsEdges() const
    {
      return m_mesh.n_edges() * m_n_local_edge_dofs;
    }

    /// Returns the total number of vertex DOFs
    inline size_t numTotalDofsCells() const
    {
      return m_mesh.n_cells() * m_n_local_cell_dofs;
    }
    
    /// Returns the dimension of the local space on the vertex V
    inline size_t dimensionVertex(const Vertex & V) const
    {
      return m_n_local_vertex_dofs;
    }

    /// Returns the dimension of the local space on the vertex of index iV
    inline size_t dimensionVertex(size_t iV) const
    {
      return dimensionVertex(*m_mesh.vertex(iV));
    }

    /// Returns the dimension of the local space on the edge E (including vertices)
    inline size_t dimensionEdge(const Edge & E) const
    {
      return 2 * m_n_local_vertex_dofs
        + m_n_local_edge_dofs;
    }

    /// Returns the dimension of the local space on the edge of index iE (including vertices)
    inline size_t dimensionEdge(size_t iE) const
    {
      return dimensionEdge(*m_mesh.edge(iE));
    }

    /// Returns the dimension of the local space on the cell T (including faces, edges and vertices)
    inline size_t dimensionCell(const Cell & T) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
        + T.n_edges() * m_n_local_edge_dofs
        + m_n_local_cell_dofs;
    }

    /// Returns the dimension of the local space on the boundary of T
    inline size_t dimensionCellBoundary(const Cell & T) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
        + T.n_edges() * m_n_local_edge_dofs;
    }
    
    /// Returns the dimension of the local space on the cell of index iT (including faces, edges and vertices)
    inline size_t dimensionCell(size_t iT) const
    {
      return dimensionCell(*m_mesh.cell(iT));
    }

    /// Returns the dimension of the local space on the boundary of T
    inline size_t dimensionCellBoundary(size_t iT) const
    {
      return dimensionCellBoundary(*m_mesh.cell(iT));
    }
    
    //------------------------------------------------------------------------------
    // Local offsets
    //------------------------------------------------------------------------------

    /// Returns the local offset of the vertex V with respect to the edge E
    inline size_t localOffset(const Edge & E, const Vertex & V) const
    {
      return E.index_vertex(&V) * m_n_local_vertex_dofs;
    }
    
    /// Returns the local offset of the unknowns attached to the edge E
    inline size_t localOffset(const Edge & E) const
    {
      return 2 * m_n_local_vertex_dofs;
    }
    
    /// Returns the local offset of the vertex V with respect to the cell T
    inline size_t localOffset(const Cell & T, const Vertex & V) const
    {
      return T.index_vertex(&V) * m_n_local_vertex_dofs;
    }

    /// Returns the local offset of the edge E with respect to the cell T
    inline size_t localOffset(const Cell & T, const Edge & E) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
        + T.index_edge(&E) * m_n_local_edge_dofs;
    }

    /// Returns the local offset of the unknowns attached to the element T
    inline size_t localOffset(const Cell & T) const
    {
      return T.n_vertices() * m_n_local_vertex_dofs
        + T.n_edges() * m_n_local_edge_dofs;
    }

  protected:
    const Mesh & m_mesh;
    size_t m_n_local_vertex_dofs;
    size_t m_n_local_edge_dofs;
    size_t m_n_local_cell_dofs;
  };
  
  //@}

} // namespace HArDCore2D
#endif
