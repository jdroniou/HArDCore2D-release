#ifndef DDRSPACE_HPP
#define DDRSPACE_HPP

#include <dofspace.hpp>

namespace HArDCore2D
{
  /*!
   * \addtogroup DDRCore
   * @{
   */

  /// Base class for DDR spaces. Provides functions to manipulate global DOFs (the local version being provided by DOFSpace).
  /** The DOFs are organised by increasing geometric entities dimensions: DOFs of vertices, DOFs of edges, DOFs of cells. */
  
  class DDRSpace : public DOFSpace {
  public:
    /// Constructor
    DDRSpace(
             const Mesh & mesh,
             size_t n_local_vertex_dofs,
             size_t n_local_edge_dofs,
             size_t n_local_cell_dofs
             );

    //------------------------------------------------------------------------------
    // Global offsets
    //------------------------------------------------------------------------------
    
    /// Return the global offset for the unknowns on the vertex V
    inline size_t globalOffset(const Vertex & V) const
    {
      return V.global_index() * m_n_local_vertex_dofs;

    }

    /// Return the global offset for the unknowns on the edge E
    inline size_t globalOffset(const Edge & E) const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs
        + E.global_index() * m_n_local_edge_dofs;
    }

    /// Return the global offset for the unknowns on the cell T
    inline size_t globalOffset(const Cell & T) const
    {
      return m_mesh.n_vertices() * m_n_local_vertex_dofs
        + m_mesh.n_edges() * m_n_local_edge_dofs
        + T.global_index() * m_n_local_cell_dofs;
    }
    
    //------------------------------------------------------------------------------
    // Restrictions
    //------------------------------------------------------------------------------

    /// Restrict to the edge (including its vertices) of index iE
    Eigen::VectorXd restrictEdge(size_t iE, const Eigen::VectorXd & vh) const;

    /// Restrict to the cell (including vertices and edges) of index iT
    Eigen::VectorXd restrictCell(size_t iT, const Eigen::VectorXd & vh) const;

    /// Restrict to an edge
    inline Eigen::VectorXd restrict(const Edge & E, const Eigen::VectorXd vh) const
    {
      return restrictEdge(E.global_index(), vh);
    }

    /// Restrict to a cell
    inline Eigen::VectorXd restrict(const Cell & T, const Eigen::VectorXd vh) const
    {
      return restrictCell(T.global_index(), vh);
    }
    
    //------------------------------------------------------------------------------
    // Extensions
    //------------------------------------------------------------------------------

    /// Extend an edge operator to a cell
    Eigen::MatrixXd extendOperator(const Cell & T, const Edge & E, const Eigen::MatrixXd & opE) const;

    //------------------------------------------------------------------------------
    // Global DOF indices for an element T
    //------------------------------------------------------------------------------

    std::vector<size_t> globalDOFIndices(const Cell & T) const;

    //------------------------------------------------------------------------------
    // Functions to handle labels
    //------------------------------------------------------------------------------
    /// Set a label to the DOF number i (default label is -1)
    inline void setLabelDOF(const size_t i, const int label)
    {
      m_labelDOF[i] = label;
    }
    
    /// Get label of DOF number i
    inline int getLabelDOF(const size_t i)
    {
      return m_labelDOF[i];
    }
    
  private:
    std::vector<int> m_labelDOF;    /// Vector to store labels for the DOFs  
    
  };

} // namespace HArDCore2D

#endif
