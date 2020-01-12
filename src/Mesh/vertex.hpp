// Class to define a vertex in 2D
//    Members: cells, edges, connected vertices...
//    Methods: index, coordinates
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef VERTEX_HPP
#define VERTEX_HPP
#include <vector>
#include <Eigen/Dense>

namespace HArDCore2D {  // forward declaration
  class Mesh;
  class Cell;
  class Edge;
}

namespace HArDCore2D {

/*!
*  @addtogroup Mesh
* @{
*/

using Eigen::Vector2d;

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The Vertex class provides description of a vertex
class Vertex {
    /**
    * A class representing a vertex of a cell for a 2D mesh. Contains a
    * pointer to the mesh, the vertex coordinates, and cells, edges and vertices linked to that vertex.
    */
public:
    /**
    * Default constructor
    *
    * @param iV global vertex number
    * @param coord  coordinates of the vertex
    * @param mesh pointer to the mesh
    */
    Vertex(size_t iV, Vector2d coords, Mesh *mesh);
    ~Vertex(); // destructor, nothing special

    inline size_t global_index() const;  ///< returns the edges global index
    inline Vector2d coords() const;    ///< returns the coordinates of the vertex

    inline size_t n_cells() const;    ///< returns the number of cells that contain the vertex
    inline size_t n_edges() const;    ///< returns the number of edges that contain the vertex
    inline size_t n_vertices() const;    ///< returns the number of vertices connected by an edge to the vertex
    inline std::vector<Cell *> get_cells() const;    ///< returns the list of cells containing the vertex      
    inline std::vector<Edge *> get_edges() const;    ///< returns the list of edges containint the vertex      
    inline std::vector<Vertex *> get_vertices() const;  ///< returns the list of vertices linked to the vertex
    Cell *cell(size_t i) const;    ///< returns i-th cell containing the vertex
    Edge *edge(size_t i) const;    ///< returns i-th edge containing the vertex
    Vertex *vertex(size_t i) const;   ///< returns i-th vertex linked to the vertex

    inline bool is_boundary() const;   ///< returns true if vertex lies on the boundary

    void add_cell(Cell *cell);   ///< Add a new cell to the list of cells containing the vertex
    void add_edge(Edge *edge);   ///< Add a new edge to the list of edges containing the vertex
    void add_vertex(Vertex *vertex);   ///< Add a new vertex to the list of vertices connected by an edge to the vertex

    void set_boundary(bool val); ///< Set the _boundary value of the vertex to val
    void set_global_index(size_t idx);  ///< Set the global index of the vertex to idx. Used to re-index the vertices, should essentially only be used inside Mesh::renum

private:
    size_t _iV;
    Vector2d _coords;
    Mesh *_mesh;
    std::vector<Cell *> _cells;
    std::vector<Edge *> _edges;
    std::vector<Vertex *> _vertices;
    bool _boundary;
  
};

// ----------------------------------------------------------------------------
//                            Implementations
// ----------------------------------------------------------------------------

size_t Vertex::global_index() const { return _iV; }
Vector2d Vertex::coords() const { return _coords;}

size_t Vertex::n_cells() const { return _cells.size();}
size_t Vertex::n_edges() const { return _edges.size();}
size_t Vertex::n_vertices() const { return _vertices.size();}
std::vector<Cell *> Vertex::get_cells() const { return _cells;}
std::vector<Edge *> Vertex::get_edges() const { return _edges;}
std::vector<Vertex *> Vertex::get_vertices() const { return _vertices;}
bool Vertex::is_boundary() const { return _boundary; }

/*@}*/
}
#endif /* VERTEX_HPP */


