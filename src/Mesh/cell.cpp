// Class to define a cell
//    Members: vertices, edges, neighbouring cells...
//    Methods: index, diameter, area, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "cell.hpp"
#include "mesh.hpp"
#include "edge.hpp"
#include "vertex.hpp"
#include <cmath>
#include <math.h>
#include <iostream>

using namespace HArDCore2D;

Cell::Cell(size_t iC, std::vector<size_t> vertex_ids, Mesh *mesh)
    : _iC(iC),
      _vertex_ids(vertex_ids),
      _mesh(mesh),
       _edges(0),
       _vertices(0),
      _neighbours(0),
      _boundary(false) {
        // need to initialise the cells edge - while we're at it why not calculate
        // the edge midpoints
        std::vector<size_t> tmp = {1,1};
        size_t nbV = vertex_ids.size();
        for (size_t i = 0; i < nbV; i++) {
          // grab vertex i and i+1 (modulo nb of vertices in cell)
          size_t inext = (i + 1) % nbV;
          // Vertices i,j appear counter-clockwise in cell
          tmp[0] = vertex_ids[i];
          tmp[1] = vertex_ids[inext];

          Edge *edge = _mesh->add_edge(tmp, this);
          _edges.push_back(edge);

          // Add vertex to cell
          Vertex *vertex = _mesh->vertex(vertex_ids[i]);
          _vertices.push_back(vertex);

        }  
        calc_cell_geometry_factors();

      }

Cell::~Cell() {}

std::vector<Edge *> Cell::get_edges() const { return _edges; }

std::vector<Vertex *> Cell::get_vertices() const { return _vertices; }

std::vector<Cell *> Cell::get_neighbours() const { return _neighbours; }

Edge *Cell::edge(size_t i) const {
    if (i < _edges.size()) {
        return _edges[i];
    } else {
      std::cerr << "[Cell] No edge at local index" << std::endl;
      exit(1);
    }
}

Vertex *Cell::vertex(size_t i) const {
    if (i < _vertices.size()) {
        return _vertices[i];
    } else {
      std::cerr << "[Cell] No vertex at local index" << std::endl;
      exit(1);
    }
}

Cell *Cell::neighbour(size_t i) const {
    if (i < _neighbours.size()) {
        return _neighbours[i];
    } else {
      std::cerr << "[Cell] No neighbour at local index" << std::endl;
      exit(1);
    }
}

size_t Cell::index_edge(const Edge* E) const {
  size_t i = 0;
  size_t nedg = n_edges();
  while(i < nedg && edge(i) != E){
    i++;
  }
  if (i >= nedg || edge(i) != E){
    std::cerr << "[Cell] Edge does not belong to cell" << std::endl;
    exit(1);
  }
  return i;
}

size_t Cell::index_vertex(const Vertex* V) const {
  size_t i = 0;
  size_t nvert = n_vertices();
  while(i < nvert && vertex(i) != V){
    i++;
  }
  if (i >= nvert || vertex(i) != V){
    std::cerr << "[Cell] Vertex does not belong to cell" << std::endl;
    exit(1);
  }
  return i;
}


Vector2d Cell::edge_normal(size_t i) const {
    //because we need to know the order of the indices in the cell we should just recreate the edge indicies here!
    size_t k = i;
    size_t j = i+1;
    if (j>=n_vertices()){
        j = 0;
        }
    Vector2d v1 = vertex(k)->coords();
    Vector2d v2 = vertex(j)->coords();

    Vector2d normal = Vector2d((v2 - v1).y(), -(v2 - v1).x());
    return normal.normalized();
}

int Cell::edge_orientation(size_t i) const { 
    return ( (_edges[i]->normal()).dot(edge_normal(i)) > 0 ? 1 : -1);
}

Vector2d Cell::ari_coords() const {
    double x = 0.0;
    double y = 0.0;
    // calculate the arithmetic coordinates using the calculated edge midpoints.
    // Currently this is just computed when needed but it could be saved or even
    // precomputed when the midpoints are calculated
    for (auto& e : _edges) {
        auto mp = e->center_mass();
        x += mp.x();
        y += mp.y();
    }

    x /= n_edges();
    y /= n_edges();
    Vector2d ari(x, y);
    return ari;
}

bool Cell::is_neighbour(const Cell *rhs) const {
    // check if the cell rhs is a neighbour, by looking at the number of vertices
    // it has in common with the current cell
    size_t n = 0;
    if (rhs->global_index() == global_index()) {
        return false;
    }
    for (auto& v : _vertices) {
        for (auto& v_rhs : rhs->get_vertices()) {
            if (v == v_rhs) {
                n++;  // nb of vertices in common
            }
        }
    }
    if (n >= 2) {
        return true;
    }
    return false;
}

bool Cell::add_neighbour(Cell *neigh) {
    _neighbours.push_back(neigh);
    return true;
}

bool Cell::calc_cell_geometry_factors() {
    std::vector<Vertex *> vlist = get_vertices();
    size_t nFV = vlist.size();
    _cell_diam = 0.0;
    for (size_t iVl = 0; iVl < nFV; iVl++) {
        Vector2d p1 = vlist[iVl]->coords();
        for (size_t jVl = iVl + 1; jVl < nFV; jVl++) {
            Vector2d p2 = vlist[jVl]->coords();
            _cell_diam = std::max(_cell_diam, (p1 - p2).norm());
        }
    }
    Vector2d xc = Vector2d::Zero();
    double de(0.), ar(0.);

    if (nFV == 3) {  // cell is a triangle
        Vector2d v0 = vlist[1]->coords() - vlist[0]->coords();
        Vector2d v1 = vlist[2]->coords() - vlist[0]->coords();
        de = v0(0) * v1(1) - v0(1) * v1(0);
        ar = sqrt(de * de);
        // be careful that ar is twice the triangle area!
        for (auto& v : vlist) {
          xc += v->coords() / 3.;
        }
        _center_mass = xc;

    } else {  // F is a planar polygon
        // get sub-triangle vertices
        Vector2d ari_center = ari_coords();
        for (size_t ilV = 0; ilV < nFV; ++ilV) {
            size_t ilVnext = (ilV + 1) % nFV;
            Vector2d v1 = vlist[ilV]->coords();
            Vector2d v2 = vlist[ilVnext]->coords();
            Vector2d x = v1 - ari_center;
            Vector2d y = v2 - ari_center;

            double tmp_de = x(0) * y(1) - x(1) * y(0);
            double tmp_ar = sqrt(tmp_de * tmp_de);
            ar += tmp_ar;
            // be careful that: tmp_ar is twice the triangle area!
            // (the coefficient should be sub_area_triangle/6)
            xc += (ari_center + v1 + v2) * tmp_ar / 3;
        }
        _center_mass = xc / ar;
    }
    // --------------------------
    _cell_area = ar / 2.;

    return true;  
}

void Cell::set_boundary(bool val) {
    _boundary = val;
  }

void Cell::set_global_index(size_t idx) {
    _iC = idx;
  }

size_t Cell::shared_edge(size_t i) {
    Edge *e = edge(i);
    size_t neighbour_cell = _mesh->n_cells() + 1;
    if (e->is_boundary()) {
        neighbour_cell = _mesh->n_cells();
    }
    //get the index for both cells that belong to this edge
    size_t iC1 = e->cell(0)->global_index();
    size_t iC2 = e->cell(1)->global_index();
    //check which one isn't the current cell 
    if (iC1 != global_index()) {
        neighbour_cell = iC1;
    }
    if (iC2 != global_index()) {
        neighbour_cell = iC2;
    }

    return neighbour_cell;
}





