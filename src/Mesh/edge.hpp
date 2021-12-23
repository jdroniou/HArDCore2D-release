// Class to define an edge
//    Members: cells, vertices...
//    Methods: index, diameter, measure, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef EDGE_HPP
#define EDGE_HPP
#include <vector>
#include <Eigen/Dense>
namespace HArDCore2D {  // forward declaration
  class Mesh;
  class Cell;
  class Vertex;
}

namespace HArDCore2D {
using Eigen::Vector2d;
/*!
*  @addtogroup Mesh
* @{
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The Edge class provides description of an edge
class Edge {
    /**
    * A class representing an edge of a 2D mesh. Contains a
    * pointer to the mesh, as well as to the vertices and cells connected to that edge
    */
public:
    /**
    * Default constructor
    *
    * @param iE global edge number
    * @param vertices list of vertices making up edge
    * @param mesh pointer to the mesh
    * @param cell pointer to the cell the edge belongs to
    */
    Edge(size_t iE, std::vector<size_t> vertices, Mesh *mesh, Cell *cell);
    ~Edge(); //  destructor, nothing special

    inline size_t global_index() const;  ///< returns the edge global index
    inline size_t n_cells() const;      ///< returns the number of cells neighbouring the edge

    std::vector<Cell *> get_cells() const; ///< list of cells that are neighbours of the edge
    std::vector<Vertex *> get_vertices() const; ///< list of vertices of the edge
    Cell *cell(size_t i) const;  ///< returns pointer to the i-th cell neighbour of the edge
    Vertex *vertex(size_t i) const;  ///< returns a pointer to the i-th vertex of the edge
    size_t index_vertex(const Vertex* V) const; ///< reciprocal of vertex(i): returns the local index of vertex V in the edge

    double measure() const;   ///< length of the edge
    double diam() const;   ///< length of the edge
    inline Vector2d center_mass() const;  ///< get the midpoint of the edge
    inline Vector2d tangent() const;   ///< get a normalised tangent to the edge, from the first vertex to the second
    inline Vector2d normal() const;   ///< get a unit normal to the edge, rotation by +pi/2 of tangent()
    inline bool is_boundary() const;  ///< getter to see if edge is boundary edge
    int vertex_orientation(size_t i) const; ///< returns the relative orientation of the i-th vertex with respect to the edge (that is, +1 if the tangent to the edge points toward the vertex, -1 otherwise).


    void add_cell(Cell *cell);      ///< Add a new cell to the edge
    void set_global_index(size_t idx);  ///< Set the global index of the edge to idx. Used to re-index the edges, should essentially only be used inside Mesh::renum

private:
    size_t _iE;
    std::vector<size_t> _vertex_ids;
    Mesh *_mesh;
    bool _boundary;
    std::vector<Cell *> _cells;
    std::vector<Vertex *> _vertices;
    Vector2d _mp;
    Vector2d _line;
  
};


// ----------------------------------------------------------------------------
//                            Implementations
// ----------------------------------------------------------------------------

size_t Edge::global_index() const { return _iE; }
size_t Edge::n_cells() const { return _cells.size();}
bool Edge::is_boundary() const { return _boundary; }

inline Vector2d Edge::center_mass() const { return _mp; }
inline Vector2d Edge::tangent() const { return _line.normalized(); }
inline Vector2d Edge::normal() const { return Vector2d(-tangent().y(), tangent().x()); }


/*@}*/
}
#endif /* EDGE_HPP */
