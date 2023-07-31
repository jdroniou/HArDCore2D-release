// Class to provide description of boundary conditions
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
 *
 *      This library was developed around HHO methods, although some parts of it have a more
 * general purpose. If you use this code or part of it in a scientific publication, 
 * please mention the following book as a reference for the underlying principles
 * of HHO schemes:
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
 * D. A. Di Pietro and J. Droniou. 2019, 516p. 
 * url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 */


#ifndef _BOUNDARY_CONDITIONS_HPP
#define _BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <string>
#include <mesh.hpp>

using namespace HArDCore2D;

/*!
 * @defgroup BoundaryConditions
 *      @brief Classes and functions to handle boundary conditions and associated numbering of unknowns
 */

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/*!
 *      \addtogroup BoundaryConditions
 * @{
 */

/// The BoundaryConditions class provides definition of boundary conditions
class BoundaryConditions {

public:
  /// Initialise data
  BoundaryConditions(
                     const std::string bc_id,  ///< The identifier for the boundary condition (D, N or Mx)
                     Mesh& mesh          ///< reference to the mesh
                     );

  /// Test the boundary condition of an edge
  const std::string type(
                         const Edge& edge    ///< Edge on which to check the boundary condition
                         ) const ;
  /**< @returns "dir" or "neu" depending if the edge is a Dirichlet or Neumann boundary condition, as determined by m_bc_id below. For an internal edge, returns "int".
     m_bc_id = "D": all edges are Dirichlet
     m_bc_id = "N": all edges are Neumann
     m_bc_id = "Mx" (x=0,1,...): various types of mixed boundary conditions, some edges will be Dirichlet and some will be Neumman.
  **/

  /// Test the boundary condition of a vertex
  const std::string type(
                         const Vertex& vertex    ///< Vertex on which to check the boundary condition
                         ) const ;

  /// Returns the number of Dirichlet edges
  inline const size_t n_dir_edges() const {
    return m_n_dir_edges;
  };

  /// Returns the number of Dirichlet vertices
  inline const size_t n_dir_vertices() const {
    return m_n_dir_vertices;
  };

  /// Returns the complete name of the boundary condition
  inline const std::string name() const {
    if (m_bc_id == "D"){
      return "Dirichlet";
    }else if (m_bc_id == "N"){
      return "Neumann";
    }else if (m_bc_id == "M0"){
      return "Mixed #0";
    }
    std::cout << "Unknown boundary conditions: " << m_bc_id << "\n";
    exit(1);
  };

  /// Re-order edges of the mesh to put the Dirichlet edges at the position "pos" (=end or start)
  void reorder_edges(const std::string pos = "end");

  /// Re-order vertices of the mesh to put the Dirichlet vertices at the position "pos" (=end or start)
  void reorder_vertices(const std::string pos = "end");

private:
  // Parameters: id of boundary condition, reference to mesh
  const std::string m_bc_id;
  Mesh& m_mesh;

  // Number of Dirichlet edges and vertices
  size_t m_n_dir_edges;
  size_t m_n_dir_vertices;

};


//@}

#endif //_BOUNDARY_CONDITION_HPP
