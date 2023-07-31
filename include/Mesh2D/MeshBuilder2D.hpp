// WARNING - this code assumes that cell vertices are listed anti clockwise and that cells do not intersect themselves
// If these conditions are not met, behaviour is unknown

#ifndef MESHBUILDER2D_HPP
#define MESHBUILDER2D_HPP

#include <memory> //std::unique_ptr
#include <MeshReaderTyp2.hpp>
#include <Polytope2D.hpp>
#include <Mesh2D.hpp>

namespace Mesh2D
{
    inline double shortest_distance(VectorRd p, VectorRd A, VectorRd B) // computes min distance from a point p to a line segment AB
    {
        VectorRd line_vec = B - A;
        VectorRd pnt_vec = p - A;

        double t = line_vec.dot(pnt_vec) / line_vec.squaredNorm();

        if (t < 0.0)
            t = 0.0;
        if (t > 1.0)
            t = 1.0;

        VectorRd nearest = t * line_vec;

        return (nearest - pnt_vec).norm();
    }

    inline Simplices<2> simplexify(std::vector<VectorRd> vertices)
    {
        assert(vertices.size() >= 3);

        Simplices<2> simplices;

        std::size_t pos = 0; // Start at the first vertex (because why not)

        double tolerance = 0.5;
        // We start with an initially high tolerance to generate as isotropic simplices as possible.
        // Each time we complete a full loop of the remaining vertices without generating a new simplex, we reduce tolerance
        // so to allow for less isotropic simplices. Continue until fully simplexified.

        bool simplex_made = true;
        std::size_t its_without_new_simplex = 0;
        while (vertices.size() > 3)
        {
            if (!simplex_made)
            {
                ++its_without_new_simplex;
                if (its_without_new_simplex > vertices.size())
                {
                    tolerance /= 2.0; // reduce the tolerance if we have tried every remaining vertex but still haven't created a simplex
                    its_without_new_simplex = 0;
                }
            }
            else
            {
                simplex_made = false;
                its_without_new_simplex = 0;
            }

            std::size_t pos_next = (pos + 1) % vertices.size();
            std::size_t pos_last = (pos_next + 1) % vertices.size();

            VectorRd v1 = vertices[pos];
            VectorRd v2 = vertices[pos_next];
            VectorRd v3 = vertices[pos_last];

            // If we know cell is convex and contains no hanging nodes, then all of the following checks are pointless
            // Perhaps in a future iteration we should add an 'is_convex' flag and if true skip the checks
            // However, it would require prior knowledge of convexity (i.e. an is_convex line in the data file)

            // Candidate simplex is v1,v2,v3. Need to check validity

            if (signed_area(v1, v2, v3) <= 0) // If negative area, trying to simplexify at non-convex vertex - simplex will be outside cell.
            {
                pos = pos_next; // skip to next vertex and continue
                continue;
            }

            // Next, we make sure we are not generating a more 'squished' triangle then necessary
            // In the limiting case, this also ensure the simplex isn't trivial (a line)

            double new_edge_length = (v1 - v3).norm();

            if ((new_edge_length > (v1 - v2).norm()) && (new_edge_length > (v2 - v3).norm()))
            {
                if (shortest_distance(v2, v1, v3) < tolerance * new_edge_length)
                {
                    // If new edge is longest edge, and second vert is very near it, triangle is bad
                    // If new edge is not longest edge, then triangle is best as can be (even if it is bad)
                    pos = pos_next; // skip to next vertex and continue
                    continue;
                }
            }

            // Finally,  we check if any of the remaining verts are in or on the new triangle
            // or if any of the remaining verts are unneccessarily close to the new triangle

            bool vert_in_triangle = false;
            for (std::size_t i = 0; i < vertices.size(); ++i)
            {
                if ((i == pos) || (i == pos_next) || (i == pos_last))
                {
                    // skip all vertices of the triangle
                    continue;
                }
                VectorRd vert = vertices[i];

                if (signed_area(vert, v3, v1) >= 0.0)
                {
                    // vert is to the right of (or in line with) the line v1v3 (new_edge)
                    // might be in or on triangle - check
                    // If not, we continue. We don't care about closeness as it must already be closer to a pre existing edge
                    if ((signed_area(vert, v2, v3) > 0.0) && (signed_area(vert, v1, v2) > 0.0))
                    {
                        vert_in_triangle = true;
                        break;
                    }
                    else
                    {
                        continue;
                    }
                }
                else
                {
                    // vert is to the left of the line v1v3 (new_edge)
                    // If vert is too close, triangle is bad
                    if (shortest_distance(vert, v1, v3) < tolerance * new_edge_length)
                    {
                        vert_in_triangle = true; // not actually in triangle, but too close for our liking
                        break;
                    }
                }
            }

            if (vert_in_triangle)
            {
                pos = pos_next; // skip to next vertex and continue
                continue;
            }

            simplices.push_back(Simplex<2>({v1, v2, v3}));
            vertices.erase(vertices.begin() + pos_next); // erase vertex located at pos_next

            // Want to move to vertex v3 (used to be located at pos_last). However, vertex at pos_next was deleted, so v3 is now located at 
            // pos_next unless pos_last < pos_next (pos_last = 0). In this case pos_next = vertices.size() - 1 so we set pos = pos_last (= 0)

            pos = pos_next % (vertices.size() - 1);
            simplex_made = true;

            assert(tolerance > 1E-16); // If tolerance gets too small, code has failed - exit.
        }

        simplices.push_back(Simplex<2>({vertices[0], vertices[1], vertices[2]})); // final simplex is the remaining vertices

        return simplices;
    }

    /*!
     *  @addtogroup Mesh
     * @{
     */

    // ----------------------------------------------------------------------------
    //                            Class definition
    // ----------------------------------------------------------------------------

    /// The MeshBuilder class provides build tools to create a full mesh with all connectivities
    class MeshBuilder
    {
    public:
        /**
         * Constructor for MeshBuilder.
         */
        MeshBuilder() {}

        /**
         * Overloaded constructor for MeshBuilder so read_mesh can be called from build_the_mesh().
         */
        MeshBuilder(const std::string mesh_file) : _mesh_file(mesh_file) {}

        /**
         *  Build mesh
         */
        std::unique_ptr<Mesh> build_the_mesh()
        {
            MeshReaderTyp2 mesh_reader(_mesh_file);

            std::vector<std::array<double, 2>> vertices;
            std::vector<std::vector<size_t>> cells;

            mesh_reader.read_mesh(vertices, cells);

            if (vertices.size() > 0 && cells.size() > 0)
            {
                std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>(); // make a pointer to the mesh so that it outlives the builder
                std::cout << "[MeshBuilder] ";

                // Create vertices
                for (auto &v : vertices)
                {
                    VectorRd vert(v[0], v[1]);
                    Vertex *vertex = new Vertex(mesh->n_vertices(), vert);
                    mesh->add_vertex(vertex);
                }

                // Create cells
                double total_area = 0.0;
                for (auto &c : cells)
                {
                    std::vector<std::size_t> vertex_ids;
                    std::vector<VectorRd> vertex_coords;
                    for (size_t i = 0; i < c.size(); i++)
                    {
                        vertex_ids.push_back(c[i] - 1); // data file starts at 1, so subtract 1
                        VectorRd coord(vertices[c[i] - 1][0], vertices[c[i] - 1][1]);
                        vertex_coords.push_back(coord);
                    }

                    Simplices<2> simplices = simplexify(vertex_coords);

                    Cell *cell = new Cell(mesh->n_cells(), simplices);
                    total_area += cell->measure();
                    mesh->add_cell(cell);

                    for (std::size_t i = 0; i < vertex_ids.size(); ++i)
                    {
                        std::size_t i_next = (i + 1) % vertex_ids.size();

                        bool edge_exists = false;

                        std::vector<Vertex *> vlist = mesh->vertex(vertex_ids[i])->get_vertices();
                        for (size_t j = 0; j < vlist.size(); j++)
                        {
                            if (vlist[j]->global_index() == vertex_ids[i_next]) // The edge exists
                            {
                                std::vector<Edge *> elist = mesh->vertex(vertex_ids[i])->get_edges();
                                Edge *edge = elist[j];

                                cell->add_edge(edge);
                                edge->add_cell(cell);
                                edge_exists = true;
                                break;
                            }
                        }

                        if (!edge_exists)
                        {
                            // edge does not exists - so create it

                            Edge *edge = new Edge(mesh->n_edges(), Simplex<1>({vertex_coords[i], vertex_coords[i_next]}));

                            mesh->add_edge(edge);
                            cell->add_edge(edge);
                            edge->add_cell(cell);

                            Vertex *vertex1 = mesh->vertex(vertex_ids[i]);
                            Vertex *vertex2 = mesh->vertex(vertex_ids[i_next]);

                            edge->add_vertex(vertex1);
                            edge->add_vertex(vertex2);
                            vertex1->add_edge(edge);
                            vertex2->add_edge(edge);
                            vertex1->add_vertex(vertex2);
                            vertex2->add_vertex(vertex1);
                        }
                        mesh->vertex(vertex_ids[i])->add_cell(cell);
                        cell->add_vertex(mesh->vertex(vertex_ids[i]));
                    }

                    cell->construct_face_normals();
                }

                // build boundary
                build_boundary(mesh.get());
                std::cout << "added " << mesh->n_cells() << " cells; Total area = " << total_area << "\n";

                return mesh;
            }
            else
            {
                throw "Cannot build mesh. Check input file";
            }
            return NULL;
        }

    private:
        void build_boundary(Mesh *mesh) ///< identifies boundary cells and vertices, and compute lists of boundary cells, edges and vertices
        {
            for (auto &face : mesh->get_faces())
            {
                if (face->n_cells() == 1) // If face has only one cell, it is a boundary face
                {
                    face->set_boundary(true);                 // set face boundary to true
                    face->get_cells()[0]->set_boundary(true); // set cell boundary to true

                    for (auto &vert : face->get_vertices())
                    {
                        vert->set_boundary(true); // set vertex boundary to true
                    }

                    mesh->add_b_face(face); // add face to list of boundary faces
                }
                else
                {
                    // internal face
                    mesh->add_i_face(face); // add face to list of internal faces
                }
            }

            for (auto &cell : mesh->get_cells())
            {
                if (cell->is_boundary())
                {
                    mesh->add_b_cell(cell);
                }
                else
                {
                    mesh->add_i_cell(cell);
                }
            }

            for (auto &vert : mesh->get_vertices())
            {
                if (vert->is_boundary())
                {
                    mesh->add_b_vertex(vert);
                }
                else
                {
                    mesh->add_i_vertex(vert);
                }
            }
        }
        const std::string _mesh_file;
    };

    /*@}*/
}
#endif /* MESHBUILDER2D_HPP */
