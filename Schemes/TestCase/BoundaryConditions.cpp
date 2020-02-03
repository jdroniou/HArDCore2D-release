// Class to provide description of boundary conditions
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "BoundaryConditions.hpp"
#include "vertex.hpp"
#include "basis.hpp" // for the VectorRd type

using namespace HArDCore2D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class
BoundaryConditions::BoundaryConditions(const std::string bc_id, Mesh& mesh)
  : m_bc_id(bc_id),
    m_mesh(mesh),
    m_n_dir_edges(0)
  {
    // Compute number of Dirichlet edges
    for (Edge* edge : m_mesh.get_b_edges()){
      if (type(*edge) == "dir"){
        m_n_dir_edges++;
      }
    }
  }



// Returns the type of BC of an edge
const std::string BoundaryConditions::type(const Edge& edge) const
  {
    if (edge.is_boundary()){
      if (m_bc_id == "D"){
        return "dir";
      }else if (m_bc_id == "N"){
        return "neu";
      }else if (m_bc_id == "M0"){
        // Dirichlet at x=0, Neumann everywhere else.
        VectorRd v0 = edge.vertex(0)->coords();
        VectorRd v1 = edge.vertex(1)->coords();
        double eps = 1e-8;
        if ( (std::abs(v0.x())<eps) && (std::abs(v1.x())<eps)){
          return "dir";
        }else{
          return "neu";
        }
      }
    }
    return "int";
  }

// Reorder edges
void BoundaryConditions::reorder_edges()
  {
    // Create vector with all non-Dirichlet edges first, and all Dirichlet edges at the end
    std::vector<size_t> new_to_old(m_mesh.n_edges(), 0);
    // Index for non-Dirichlet edges start at 0 and increases, index for Dirichlet edges start at n_edges()-1
    // and decreases
    size_t idx_nondir = 0;
    size_t idx_dir = m_mesh.n_edges()-1;
    for (Edge* edge : m_mesh.get_edges()){
      if (idx_nondir > idx_dir){
        std::cout << "Error during creation vector to renumber edges: " << idx_nondir << ", " << idx_dir << "\n";
        exit(1);
      }
      if (type(*edge) == "dir"){
        new_to_old[idx_dir] = edge->global_index();
        idx_dir--;
      }else{
        new_to_old[idx_nondir] = edge->global_index();
        idx_nondir++;
      }
    }
    // Check: idx_dir and idx_nondir must just have crossed
    if (idx_nondir != idx_dir + 1){
      std::cout << "Error in creating vector to renumber edges: " << idx_nondir << "/" << m_mesh.n_edges() - m_n_dir_edges << " || " << idx_dir << "/" << m_mesh.n_edges() - m_n_dir_edges << "\n";
      exit(1);
    }

    // Reordering
    m_mesh.renum('E', new_to_old);

 }

