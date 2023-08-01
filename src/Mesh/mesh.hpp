// Class to describe a mesh.
//    Members: cells, vertices, edges...
//    Methods: h_max, add cells and edges...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
* This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
*  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
*  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
*  url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/



#ifndef MESH_HPP
#define MESH_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "cell.hpp"
namespace HArDCore2D {  // forward declaration
  class Edge;
  class Vertex;
}

/*!  
*  @defgroup Mesh 
* @brief Classes to construct and describe a 2D mesh
*/

namespace HArDCore2D {
/**
* Class which represents a 2D mesh. Contains cells, vertices, edges.
*/

using Eigen::Vector2d;

/*!
*  @addtogroup Mesh
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The Mesh class provides description of a mesh

class Mesh {
public:
    /**
    * default constructor for an empty mesh
    */
    Mesh();
    ~Mesh();

    inline void set_name(std::string name);  ///< set the name of the mesh
    inline std::string get_name();  ///< getter for the edge name
    inline size_t n_cells() const;     ///< number of cells in the mesh
    inline size_t n_edges() const;     ///< number of edges in the mesh
    inline size_t n_vertices() const;  ///< number of vertices in the mesh
    inline double h_max() const;       ///< max of diameter of cells
    inline size_t dim() const;         ///< dimension of the mesh (2)
    size_t n_b_cells() const;           ///< number of boundary cells
    size_t n_b_edges() const;           ///< number of boundary edges
    size_t n_b_vertices() const;        ///< number of boundary vertices
    size_t n_i_cells() const;           ///< number of interior cells
    size_t n_i_edges() const;           ///< number of interior edges
    size_t n_i_vertices() const;        ///< number of interior vertices

    inline std::vector<Cell*> get_cells() const;  ///< lists the cells in the mesh.
    inline std::vector<Edge*> get_edges() const;  ///< lists the edges in the mesh.
    inline std::vector<Vertex*> get_vertices() const;  ///< lists the vertices in the mesh.
    Cell* cell(size_t iC) const;  ///< get a constant pointer to a cell using its global index
    Edge* edge(size_t iE) const;   ///< get a constant pointer to an edge using its global index
    Vertex* vertex(size_t iV) const;   ///< get a constant pointer to a vertex using its global index

    inline std::vector<Cell*> get_b_cells() const;  ///< lists the boundary cells in the mesh.
    inline std::vector<Edge*> get_b_edges() const;  ///< lists the boundary edges in the mesh.
    inline std::vector<Vertex*> get_b_vertices() const;  ///< lists the boundary vertices in the mesh.
    Cell* b_cell(size_t iC) const;  ///< get a constant pointer to the iC-th boundary cell
    Edge* b_edge(size_t iE) const;   ///< get a constant pointer to the iE-th boundary edge
    Vertex* b_vertex(size_t iV) const;   ///< get a constant pointer to the iV-th boundary vertex

    inline std::vector<Cell*> get_i_cells() const;  ///< lists the interior cells in the mesh.
    inline std::vector<Edge*> get_i_edges() const;  ///< lists the interior edges in the mesh.
    inline std::vector<Vertex*> get_i_vertices() const;  ///< lists the interior vertices in the mesh.
    Cell* i_cell(size_t iC) const;  ///< get a constant pointer to the iC-th interior cell
    Edge* i_edge(size_t iE) const;   ///< get a constant pointer to the iE-th interior edge
    Vertex* i_vertex(size_t iV) const;   ///< get a constant pointer to the iV-th interior vertex

    size_t find_cell(const Vector2d & x);   ///< returns the index of a cell that contains x (this function is a bit expensive)

    inline bool add_cell(Cell* cell);  ///<  adds a cell to the mesh
    inline bool add_vertex(Vertex* vertex);  ///<  adds a vertex to the mesh
    Edge* add_edge(std::vector<size_t> vertex_ids, Cell* cell);  ///< add an edge to the mesh

    inline bool add_b_cell(Cell* cell);  ///<  adds a boundary cell to the mesh
    inline bool add_b_edge(Edge* edge);  ///<  adds a boundary edge to the mesh
    inline bool add_b_vertex(Vertex* vertex);  ///<  adds a boundary vertex to the mesh

    inline bool add_i_cell(Cell* cell);  ///<  adds an interior cell to the mesh
    inline bool add_i_edge(Edge* edge);  ///<  adds an interior edge to the mesh
    inline bool add_i_vertex(Vertex* vertex);  ///<  adds an interior vertex to the mesh

    inline size_t next_edge_idx();  ///< gets the next global edge index

    std::vector<double> regularity(); ///< returns regularity factors
   
    /// Re-index the cells, edges or vertices
    void renum(                      
      const char B,                         ///< T for cells, E for edges, V for vertices
      const std::vector<size_t> new_to_old   ///< Vector of new indices to old ones (new_to_old[i]=j: the index formerly j will become i)
      );

    
private:
    std::string _mesh_name;
    size_t _next_edge_idx = 0;

    // primary data: list of cells, edges, vertices...
    std::vector<Cell*> _cells;
    std::vector<Edge*> _edges;
    std::vector<Vertex*> _vertices;
    std::vector<Cell*> _b_cells;
    std::vector<Edge*> _b_edges;
    std::vector<Vertex*> _b_vertices;
    std::vector<Cell*> _i_cells;
    std::vector<Edge*> _i_edges;
    std::vector<Vertex*> _i_vertices;
    double _h_max;
    const size_t _dim = 2;
  
};



// ----------------------------------------------------------------------------
//                            Implementations
// ----------------------------------------------------------------------------

size_t Mesh::n_cells() const { return _cells.size(); }
size_t Mesh::n_edges() const { return _edges.size(); }
size_t Mesh::n_vertices() const { return _vertices.size(); }
double Mesh::h_max() const { return _h_max; }
size_t Mesh::dim() const { return _dim; }
size_t Mesh::next_edge_idx() {
    _next_edge_idx++;
    return _next_edge_idx - 1;
}
void Mesh::set_name(std::string name) {
    _mesh_name = name;
    return;
}
std::string Mesh::get_name() { return _mesh_name; }

bool Mesh::add_cell(Cell* cell) {
    // Add the cell to the mesh and update mesh size
    _cells.push_back(cell);
    _h_max = std::max(_h_max, cell->diam());

    return true;
}

bool Mesh::add_b_cell(Cell* cell) {
    _b_cells.push_back(cell);

    return true;
}

bool Mesh::add_i_cell(Cell* cell) {
    _i_cells.push_back(cell);

    return true;
}

bool Mesh::add_b_edge(Edge* edge) {
    _b_edges.push_back(edge);

    return true;
}

bool Mesh::add_i_edge(Edge* edge) {
    _i_edges.push_back(edge);

    return true;
}

bool Mesh::add_vertex(Vertex* vertex) {
    _vertices.push_back(vertex);

    return true;
}

bool Mesh::add_b_vertex(Vertex* vertex) {
    _b_vertices.push_back(vertex);

    return true;
}

bool Mesh::add_i_vertex(Vertex* vertex) {
    _i_vertices.push_back(vertex);

    return true;
}

std::vector<Cell*> Mesh::get_cells() const { return _cells; }
std::vector<Edge*> Mesh::get_edges() const { return _edges; }
std::vector<Vertex*> Mesh::get_vertices() const { return _vertices; }
std::vector<Cell*> Mesh::get_b_cells() const { return _b_cells; }
std::vector<Edge*> Mesh::get_b_edges() const { return _b_edges; }
std::vector<Vertex*> Mesh::get_b_vertices() const { return _b_vertices; }
std::vector<Cell*> Mesh::get_i_cells() const { return _i_cells; }
std::vector<Edge*> Mesh::get_i_edges() const { return _i_edges; }
std::vector<Vertex*> Mesh::get_i_vertices() const { return _i_vertices; }
/*@}*/
}

#endif /* MESH_HPP */

