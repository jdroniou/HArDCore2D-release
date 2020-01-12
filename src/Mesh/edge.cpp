// Class to define an edge
//    Members: cells, vertices...
//    Methods: index, diameter, measure, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "edge.hpp"
#include "cell.hpp"
#include "vertex.hpp"
#include "mesh.hpp"
#include <iostream>

using namespace HArDCore2D;
Edge::Edge(size_t iE, std::vector<size_t> vertices, Mesh *mesh, Cell *cell)
    : _iE(iE), 
      _vertex_ids(vertices), 
      _mesh(mesh), 
      _boundary(true),
      _cells(0) {
        _cells.push_back(cell);
        _vertices.push_back(_mesh->vertex(vertices[0]));
        _vertices.push_back(_mesh->vertex(vertices[1]));
        Vector2d p1 = _mesh->vertex(vertices[0])->coords();
        Vector2d p2 = _mesh->vertex(vertices[1])->coords();

        if (_cells[0] == NULL) {
          std::cout << "cell nil" << std::endl;
        }
        _line = p2 - p1;
        _mp = (p1 + p2)/2;
}

Edge::~Edge() {}

Cell *Edge::cell(size_t i) const {
    if (i < _cells.size()) {
        return _cells[i];
    } else {
        throw "No cell at edge local index";
    }
}

Vertex *Edge::vertex(size_t i) const {
    if (i < _vertices.size()) {
        return _vertices[i];
    } else {
        throw "No vertex at edge local index";
    }
}

/// Measure and diameter are the same for an edge, but for consistency with the way schemes
/// are usually defined, we provide two functions. For translation of each scheme's code to 3D (with
/// edges becoming faces), it is better to use the proper one depending if the measure or diameter is expected
double Edge::measure() const { return _line.norm(); }
/// Measure and diameter are the same for an edge, but for consistency with the way schemes
/// are usually defined, we provide two functions. For translation of each scheme's code to 3D (with
/// edges becoming faces), it is better to use the proper one depending if the measure or diameter is expected
double Edge::diam() const { return _line.norm(); }
std::vector<Cell *> Edge::get_cells() const { return _cells; }
std::vector<Vertex *> Edge::get_vertices() const { return _vertices; }

void Edge::add_cell(Cell *cell) {
    if (_cells.size() < 2) {
        _cells.push_back(cell);
        // if the edge has 2 cells its not a boundary
        _boundary = false;
    } else {
        throw "An edge cannot have three cells!";
    }
}

void Edge::set_global_index(size_t idx) {
    _iE = idx;
  }


