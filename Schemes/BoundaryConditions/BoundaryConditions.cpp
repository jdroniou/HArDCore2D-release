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
    m_n_dir_edges(0),
    m_n_dir_vertices(0)
  {
    // Compute number of Dirichlet edges
    for (Edge* edge : m_mesh.get_b_edges()){
      if (type(*edge) == "dir"){
        m_n_dir_edges++;
      }
    }
    // Compute number of Dirichlet vertices
    for (Vertex* vertex : m_mesh.get_b_vertices()){
      if (type(*vertex) == "dir"){
        m_n_dir_vertices++;
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

// Returns the type of BC of a vertex
const std::string BoundaryConditions::type(const Vertex& vertex) const
  {
    if (vertex.is_boundary()){
      if (m_bc_id == "D"){
        return "dir";
      }else if (m_bc_id == "N"){
        return "neu";
      }else if (m_bc_id == "M0"){
        // Dirichlet at x=0, Neumann everywhere else.
        VectorRd v0 = vertex.coords();
        double eps = 1e-8;
        if ( (std::abs(v0.x())<eps) ){
          return "dir";
        }else{
          return "neu";
        }
      }
    }
    return "int";
  }

// Reorder edges
void BoundaryConditions::reorder_edges(const std::string pos)
  {
    // Create a vector of Dirichlet boundary edges, and a vector of other edges
    std::vector<size_t> dir_edges(m_mesh.n_edges(), 0);
    std::vector<size_t> nondir_edges(m_mesh.n_edges(), 0);
    size_t dir_idx = 0;
    size_t nondir_idx = 0;
    for (Edge* edge : m_mesh.get_edges()){
      if (type(*edge) == "dir"){
        dir_edges[dir_idx] = edge->global_index();
        dir_idx++;
      }else{
        nondir_edges[nondir_idx] = edge->global_index();
        nondir_idx++;
      }
    }
    // check
    if (dir_idx + nondir_idx != m_mesh.n_edges()){
     std::cout << "Error during renumbering edges: " << dir_idx << ", " << nondir_idx << ", " << m_mesh.n_edges() << "\n";
     exit(1);
    }
    
    // Depending on "pos" we put the Dirichlet edges at the end or the start
    std::vector<size_t> new_to_old(m_mesh.n_edges(), 0);
    if (pos=="end"){
      for (size_t i=0; i < nondir_idx; i++){
        new_to_old[i] = nondir_edges[i];
      }     
      for (size_t i=nondir_idx; i < m_mesh.n_edges(); i++){
        new_to_old[i] = dir_edges[i-nondir_idx];
      }
    }else{
      for (size_t i=0; i < dir_idx; i++){
        new_to_old[i] = dir_edges[i];
      }     
      for (size_t i=dir_idx; i < m_mesh.n_edges(); i++){
        new_to_old[i] = nondir_edges[i-dir_idx];
      }
    }

    // Reordering
    m_mesh.renum('E', new_to_old);


//////    // Create vector with all non-Dirichlet edges first, and all Dirichlet edges at the end
//////    std::vector<size_t> new_to_old(m_mesh.n_edges(), 0);
//////    // Index for non-Dirichlet edges start at 0 and increases, index for Dirichlet edges start at n_edges()-1
//////    // and decreases
//////    size_t idx_nondir = 0;
//////    size_t idx_dir = m_mesh.n_edges()-1;
//////    for (Edge* edge : m_mesh.get_edges()){
//////      if (idx_nondir > idx_dir){
//////        std::cout << "Error during creation vector to renumber edges: " << idx_nondir << ", " << idx_dir << "\n";
//////        exit(1);
//////      }
//////      if (type(*edge) == "dir"){
//////        new_to_old[idx_dir] = edge->global_index();
//////        idx_dir--;
//////      }else{
//////        new_to_old[idx_nondir] = edge->global_index();
//////        idx_nondir++;
//////      }
//////    }
//////    // Check: idx_dir and idx_nondir must just have crossed
//////    if (idx_nondir != idx_dir + 1){
//////      std::cout << "Error in creating vector to renumber edges: " << idx_nondir << "/" << m_mesh.n_edges() - m_n_dir_edges << " || " << idx_dir << "/" << m_mesh.n_edges() - m_n_dir_edges << "\n";
//////      exit(1);
//////    }

//////    // Reordering
//////    m_mesh.renum('E', new_to_old);

 }

// Reorder vertices
void BoundaryConditions::reorder_vertices(const std::string pos)
  {
    // Create a vector of Dirichlet boundary vertices, and a vector of other vertices
    std::vector<size_t> dir_vertices(m_mesh.n_vertices(), 0);
    std::vector<size_t> nondir_vertices(m_mesh.n_vertices(), 0);
    size_t dir_idx = 0;
    size_t nondir_idx = 0;
    for (Vertex* vertex : m_mesh.get_vertices()){
      if (type(*vertex) == "dir"){
        dir_vertices[dir_idx] = vertex->global_index();
        dir_idx++;
      }else{
        nondir_vertices[nondir_idx] = vertex->global_index();
        nondir_idx++;
      }
    }
    // check
    if (dir_idx + nondir_idx != m_mesh.n_vertices()){
     std::cout << "Error during renumbering vertices: " << dir_idx << ", " << nondir_idx << ", " << m_mesh.n_vertices() << "\n";
     exit(1);
    }
    
    // Depending on "pos" we put the Dirichlet vertices at the end or the start
    std::vector<size_t> new_to_old(m_mesh.n_vertices(), 0);
    if (pos=="end"){
      for (size_t i=0; i < nondir_idx; i++){
        new_to_old[i] = nondir_vertices[i];
      }     
      for (size_t i=nondir_idx; i < m_mesh.n_vertices(); i++){
        new_to_old[i] = dir_vertices[i-nondir_idx];
      }
    }else{
      for (size_t i=0; i < dir_idx; i++){
        new_to_old[i] = dir_vertices[i];
      }     
      for (size_t i=dir_idx; i < m_mesh.n_vertices(); i++){
        new_to_old[i] = nondir_vertices[i-dir_idx];
      }
    }

    // Reordering
    m_mesh.renum('V', new_to_old);

 }

