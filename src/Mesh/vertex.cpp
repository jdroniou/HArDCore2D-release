// Class to define a vertex in 2D
//    Members: cells, edges, connected vertices...
//    Methods: index, coordinates
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "vertex.hpp"
#include "mesh.hpp"
#include "cell.hpp"
#include "edge.hpp"
#include <iostream>

using namespace HArDCore2D;
Vertex::Vertex(size_t iV, Vector2d coords)
    : _iV(iV), 
      _coords(coords), 
      _cells(0),
      _edges(0),
      _vertices(0),
      _boundary(false) {
        // Do nothing
        }

Vertex::~Vertex() {}

Cell *Vertex::cell(size_t i) const {
    if (i < _cells.size()) {
        return _cells[i];
    } else {
        throw "No cell at vertex local index";
    }
}

Edge *Vertex::edge(size_t i) const {
    if (i < _edges.size()) {
        return _edges[i];
    } else {
        throw "No edge at vertex local index";
    }
}

Vertex *Vertex::vertex(size_t i) const {
    if (i < _vertices.size()) {
        return _vertices[i];
    } else {
        throw "No vertex at vertex local index";
    }
}

void Vertex::add_cell(Cell *cell) {
    _cells.push_back(cell);
}

void Vertex::add_edge(Edge *edge) {
    _edges.push_back(edge);
}

void Vertex::add_vertex(Vertex *vertex) {
    _vertices.push_back(vertex);
}

void Vertex::set_boundary(bool val) {
    _boundary = val;
  }

void Vertex::set_global_index(size_t idx) {
    _iV = idx;
  }

