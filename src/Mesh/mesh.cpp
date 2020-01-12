// Class to describe a mesh.
//    Members: cells, vertices, edges...
//    Methods: h_max, add cells and edges...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "mesh.hpp"
#include "cell.hpp"
#include "edge.hpp"
#include "vertex.hpp"
#include <iostream>
#include <set>
#include <Eigen/Dense>  //Vector2d

using namespace HArDCore2D;
using Eigen::Vector2d;
Mesh::Mesh() : _mesh_name("mesh-2d"),
          _next_edge_idx(0),
          _cells(0),
          _edges(0),
          _vertices(0),
          _b_cells(0),
          _b_edges(0),
          _b_vertices(0),
          _i_cells(0),
          _i_edges(0),
          _i_vertices(0),
          _h_max(0)
          {
            // do nothing
          }

Mesh::~Mesh() {
      for (auto& cell : _cells){
        delete cell;
      }
      for (auto& edge : _edges){
        delete edge;
      }
      for (auto& vertex : _vertices){
        delete vertex;
      }
    }

Cell* Mesh::cell(size_t iC) const {
    if (iC < n_cells()) {
        return _cells[iC];
    } else {
        throw "Trying to access cell at global index which does not exist";
    }
}

Edge* Mesh::edge(size_t iE) const {
    if (iE < n_edges()) {
        return _edges[iE];
    } else {
        throw "Trying to access edge at global index which does not exist";
    }
}

Vertex* Mesh::vertex(size_t iV) const {
    if (iV < n_vertices()) {
        return _vertices[iV];
    } else {
        throw "Trying to access vertex at global index which does not exist";
    }
}

Cell* Mesh::b_cell(size_t iC) const {
    if (iC < n_b_cells()) {
        return _b_cells[iC];
    } else {
        throw "Trying to access boundary cell at global index which does not exist";
    }
}

Edge* Mesh::b_edge(size_t iE) const {
    if (iE < n_b_edges()) {
        return _b_edges[iE];
    } else {
        throw "Trying to access boundary edge at global index which does not exist";
    }
}

Vertex* Mesh::b_vertex(size_t iV) const {
    if (iV < n_b_vertices()) {
        return _b_vertices[iV];
    } else {
        throw "Trying to access boundary vertex at global index which does not exist";
    }
}

Cell* Mesh::i_cell(size_t iC) const {
    if (iC < n_i_cells()) {
        return _i_cells[iC];
    } else {
        throw "Trying to access interior cell at global index which does not exist";
    }
}

Edge* Mesh::i_edge(size_t iE) const {
    if (iE < n_i_edges()) {
        return _i_edges[iE];
    } else {
        throw "Trying to access interior edge at global index which does not exist";
    }
}

Vertex* Mesh::i_vertex(size_t iV) const {
    if (iV < n_i_vertices()) {
        return _i_vertices[iV];
    } else {
        throw "Trying to access interior vertex at global index which does not exist";
    }
}

Edge* Mesh::add_edge(std::vector<size_t> vertex_ids, Cell* cell) {
    // vertex_ids are the ids of two consecutive vertices in counter-clockwise order in cell
    // These vertices define an edge of cell

    // We look if the edge already exists, by looking at all the vertices already linked to vertex_ids[0]
    auto vlist = this->vertex(vertex_ids[0])->get_vertices();
    for (size_t j = 0; j < vlist.size(); j++) {
      if (vlist[j]->global_index() == vertex_ids[1]) {
        // The edge exists, we grab it together with the cell it already belongs to
        auto elist = this->vertex(vertex_ids[0])->get_edges();
        Edge* e = elist[j];
        auto cell_of_e = e->cell(0);

        // We add cell as a neighbour of the edge, and the two cells as neighbours of each other
        e->add_cell(cell);
        cell_of_e->add_neighbour(cell);
        cell->add_neighbour(cell_of_e);

        return e;
      }
    }

    // The edge does not exist: get the new global number and create the edge
    size_t iE = next_edge_idx();  
    Edge* edge = new Edge(iE, vertex_ids, this, cell);
    _edges.push_back(edge);

    // add each vertex to each other's list, and corresponding edge
    Vertex* vertex1 = vertex(vertex_ids[0]);
    Vertex* vertex2 = vertex(vertex_ids[1]);
    vertex1->add_edge(edge);
    vertex2->add_edge(edge);
    vertex1->add_vertex(vertex2);
    vertex2->add_vertex(vertex1);

    return edge;
}


size_t Mesh::n_b_cells() const {
    return _b_cells.size();
}
size_t Mesh::n_b_edges() const {
    return _b_edges.size();
}
size_t Mesh::n_b_vertices() const {
    return _b_vertices.size();
}

size_t Mesh::n_i_cells() const {
    return _i_cells.size();
}
size_t Mesh::n_i_edges() const {
    return _i_edges.size();
}
size_t Mesh::n_i_vertices() const {
    return _i_vertices.size();
}

std::vector<double> Mesh::regularity(){
    /// Regularity factor = 
    ///   1st component: maximum of
    ///      * diameter of cell / (measure of cell)^{1/2}
    ///      * diameter of cell / diameter of edge  [for each edge of the cell]
    ///
    ///   2nd component: evaluation of ratio "diam of cell / radius ball inscribed in cell"

    std::vector<double> value(2,0.0);
    for (auto& icell : get_cells()){
      double hT = icell->diam();
      Eigen::Vector2d xT = icell->center_mass();

      value[0] = std::max(value[0], hT / pow(icell->measure(), 1/this->dim()));
 
      double rhoT = hT;
      std::vector<Edge *> edges = icell->get_edges();
      for (size_t i=0; i < edges.size(); i++){
        Edge* iedge = edges[i];
        double hF = iedge->diam();
        Eigen::Vector2d xF = iedge->center_mass();
        Eigen::Vector2d nTF = icell->edge_normal(i);

        value[0] = std::max(value[0], hT / hF);

        rhoT = std::min(rhoT, std::abs( (xT-xF).dot(nTF) ) );
      }
      value[1] = std::max(value[1], hT/rhoT);

    }

    return value;

}

void Mesh::renum(const char B, const std::vector<size_t> new_to_old){
  
  switch (B) {
    case 'C': {
      std::vector<Cell*> old_index = _cells;
      for (size_t i=0; i < n_cells(); i++){
        old_index[new_to_old[i]]->set_global_index(i);
        _cells[i] = old_index[new_to_old[i]];
        }
      break;
      }

    case 'E': {
      std::vector<Edge*> old_index = _edges;
      for (size_t i=0; i < n_edges(); i++){
        old_index[new_to_old[i]]->set_global_index(i);
        _edges[i] = old_index[new_to_old[i]];
        }
      break;
      }

    case 'V': {
      std::vector<Vertex*> old_index = _vertices;
      for (size_t i=0; i < n_vertices(); i++){
        old_index[new_to_old[i]]->set_global_index(i);
        _vertices[i] = old_index[new_to_old[i]];
        }
      break;
      }

  }

}

