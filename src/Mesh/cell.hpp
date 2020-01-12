// Class to define a cell
//    Members: vertices, edges, neighbouring cells...
//    Methods: index, diameter, area, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef CELL_HPP
#define CELL_HPP
#include <vector>
#include <Eigen/Dense>
namespace HArDCore2D {  // forward declaration
  class Mesh;
  class Edge;
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

/// The Cell class provides description of a cell
class Cell {
public:
    /**
    * Constructor from global indexes and a mesh
    * @param iG is the global index of the Cell (basically its just the location in the mesh cells vector
    * @param vertex_ids vector containing the global indices of the cell vertices
    * @param mesh pointer to the mesh that the cell is contained in
    */
    Cell(size_t iC, std::vector<size_t> vertex_ids, Mesh *mesh);
    ~Cell();

    inline size_t global_index() const;    ///< cell index
    inline size_t n_edges() const;  ///<  returns number of edges of the cell
    inline size_t n_vertices() const;  ///< returns number of vertices of the cell

    std::vector<Edge *> get_edges() const;    ///< returns the list of edges of the cell
    std::vector<Vertex *> get_vertices() const;    ///< returns the list of vertices of the cell
    std::vector<Cell *> get_neighbours() const;     ///< returns the list of neighbours of the cell
    Edge *edge(size_t iL) const;  ///< returns the iL-th edge of the cell
    Vertex *vertex(size_t iL) const;   ///< returns the iL-th edge of the cell
    Cell *neighbour(size_t iL) const; ///< returns the iL-th neighbour of the cell

    inline bool is_boundary() const; ///< returns true if cell touches the boundary

    inline double measure() const; ///< returns area of cell
    inline double diam() const;  ///< returns diameter of cell
    Vector2d edge_normal(size_t i);  ///< returns the outer normal to the i-th edge
    inline Vector2d center_mass() const;  ///< returns the center of mass of the cell
    bool calc_cell_geometry_factors();  ///< calculate cell diam, area etc

    bool add_neighbour(Cell *neigh);  ///< add a cell to the neighbour
    bool is_neighbour(const Cell *rhs) const;  ///< true if shares edge

    Vector2d ari_coords() const;  ///< return average of edge midpoints

    size_t shared_edge(size_t i);  ///< Returns the global index to the cell
                             /// that shares the
                             /// edge defined by the local coordinates with
                             /// this cell if there is no edge it will
                             /// return a null pointer, it would be best to
                             /// check if the cell is a boundary before
    /// calling this function. Returns -1 if there is no cell

    void set_boundary(bool val); ///< Set the _boundary value of the cell to val
    void set_global_index(size_t idx);  ///< Set the global index of the cell to idx. Used to re-index the cells, should essentially only be used inside Mesh::renum


private:
    size_t _iC;    ///< cell global index
    const std::vector<size_t> _vertex_ids;  ///< the global ids of the cell vertices
    Mesh *_mesh;                     ///< pointer to the owner mesh
    std::vector<Edge *> _edges;      ///< list of cell edges
    std::vector<Vertex *> _vertices;      ///< a list of cell vertices
    std::vector<Cell *> _neighbours;  ///< list of cell neighbours

    bool _boundary;                    ///< flag is cell boundary?

    Vector2d _center_mass;                ///< center of mass of the cell
    double _cell_diam;                 ///< diameter of the cell
    double _cell_area;                 ///< area of the cell

};


// ----------------------------------------------------------------------------
//                            Implementations
// ----------------------------------------------------------------------------

Vector2d Cell::center_mass() const { return _center_mass; }
bool Cell::is_boundary() const { return _boundary; }
double Cell::diam() const { return _cell_diam; }
double Cell::measure() const { return _cell_area; }
size_t Cell::global_index() const { return _iC; }
size_t Cell::n_edges() const { return _edges.size(); }
size_t Cell::n_vertices() const { return _vertices.size(); }

/*@}*/
}
#endif /* CELL_HPP */
