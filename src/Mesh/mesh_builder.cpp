// Class to build the mesh data after having read the mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "mesh_builder.hpp"
#include <iostream>
#include <deque>

using namespace HArDCore2D;
MeshBuilder::MeshBuilder() {}
MeshBuilder::MeshBuilder(const std::string mesh_file) : _mesh_file(mesh_file) {}
std::unique_ptr<Mesh> MeshBuilder::build_the_mesh(
    std::vector<std::vector<double> > vertices,
    std::vector<std::vector<size_t> > cells) {
    if (vertices.size() > 0 && cells.size() > 0) {
        std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>();  // make a pointer to the mesh so that it outlives the builder
        std::cout << "Mesh: ";

        // Create vertices
        size_t iG = 0;
        for (auto& v : vertices) {
            Vertex* vertex = new Vertex(iG, Vector2d(v[0], v[1]), mesh.get());
            mesh->add_vertex(vertex);
            iG++;
        }

        // Create cells
        double total_area = 0.0;
        iG = 0;
        for (auto& c : cells) {
            // first entry is the number of nodes so skip it
            std::vector<size_t> vertex_ids;
            for (size_t i = 1; i < c.size(); i++) {
                vertex_ids.push_back(c[i]);
            }
            Cell* cell = new Cell(iG, vertex_ids, mesh.get());
            total_area += cell->measure();
            mesh->add_cell(cell);
            iG++;

            // add the cell to all its vertices
            for (auto& vid : vertex_ids){
              mesh->vertex(vid)->add_cell(cell);
            }
        }


        // build boundary
        build_boundary(mesh.get());

        std::cout << "added " << mesh->n_cells() << " cells; Total area= " << total_area << std::endl;
        return mesh;
    } else {
        if (vertices.size() <= 0) {
            throw "Can't build mesh vertices is empty. Check the input file";
        }
        if (cells.size() <= 0) {
            throw "Can't build mesh cells is empty. Check the input file";
        }
    }
    return NULL;
}

std::unique_ptr<Mesh> MeshBuilder::build_the_mesh()
{
    MeshReaderTyp2 mesh_ptr(_mesh_file);

    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<size_t>> cells;
    std::vector<std::vector<double>> centers;

    mesh_ptr.read_mesh(vertices, cells, centers);

    return build_the_mesh(vertices, cells);
}

void MeshBuilder::build_boundary(Mesh* mesh) {
    // When the mesh is built, boundary edges are already identified (see Edge::add_cell)
    // Here we fill in the _boundary variables of the cells and vertices, and the lists of boundary
    // edges, cells and vertices
    for (auto& cell : mesh->get_cells()) {
      for (auto& edge : cell->get_edges()) {
        if (edge->is_boundary()) {
          // The cell has a boundary edge, so it is a boundary cell
          cell->set_boundary(true);
          mesh->add_b_cell(cell);
          // We also add the edge to the boundary edges
          mesh->add_b_edge(edge);
        }
      }
      // If we have a boundary cell, we explore its vertices and those connected to 
      // boundary edges are boundary vertices
      if (cell->is_boundary()) {
        for (auto& vertex : cell->get_vertices())  {
          for (auto& edge : vertex->get_edges()){
            if (edge->is_boundary()){
              vertex->set_boundary(true);
              mesh->add_b_vertex(vertex);
            }
          }
        }

      }
    }

    // Pass to fill in interior elements
    for (auto& cell : mesh->get_cells()){
      if ( !(cell->is_boundary()) ){
        mesh->add_i_cell(cell);
      }  
    }
    for (auto& edge : mesh->get_edges()){
      if ( !(edge->is_boundary()) ){
        mesh->add_i_edge(edge);
      }  
    }
    for (auto& vertex : mesh->get_vertices()){
      if ( !(vertex->is_boundary()) ){
        mesh->add_i_vertex(vertex);
      }  
    }
}


