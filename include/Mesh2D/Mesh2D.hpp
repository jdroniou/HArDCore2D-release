#include "Polytope2D.hpp"

#ifndef _MESH2D_HPP
#define _MESH2D_HPP

namespace Mesh2D
{
    inline double signed_area(VectorRd A, VectorRd B, VectorRd C) // returns signed area of trangle ABC. Positive if anticlock, negative otherwise
    {
        return (A(0) * (B(1) - C(1)) + B(0) * (C(1) - A(1)) + C(0) * (A(1) - B(1))) / 2.0;
    }


     /*!
     * @defgroup Mesh
     * @brief Classes to construct and describe a 2D mesh
     */


    /*!
     *  @addtogroup Mesh
     * @{
     */

    class Mesh
    {
    public:
        Mesh() {}
        ~Mesh()
        {
            for (auto &vertex : _vertices)
            {
                delete vertex;
            }
            for (auto &edge : _edges)
            {
                delete edge;
            }
            for (auto &cell : _cells)
            {
                delete cell;
            }
        }

        void set_name(std::string name) { _mesh_name = name; } ///< set the name of the mesh
        inline std::string get_name() { return _mesh_name; }   ///< getter for the mesh name

        double h_max() const ///< max diameter of cells
        {
            double val = 0.0;
            for (auto &cell : _cells)
            {
                val = std::max(val, cell->diam());
            }
            return val;
        }

        inline std::size_t dim() const { return DIMENSION; } ///< dimension of the mesh

        inline std::size_t n_vertices() const { return _vertices.size(); } ///< number of vertices in the mesh.
        inline std::size_t n_edges() const { return _edges.size(); }       ///< number of edges in the mesh.
        inline std::size_t n_faces() const { return _edges.size(); }       ///< number of faces in the mesh.
        inline std::size_t n_cells() const { return _cells.size(); }       ///< number of cells in the mesh.

        inline std::size_t n_b_vertices() const { return _b_vertices.size(); } ///< number of boundary vertices in the mesh.
        inline std::size_t n_b_edges() const { return _b_edges.size(); }       ///< number of boundary edges in the mesh.
        inline std::size_t n_b_faces() const { return _b_edges.size(); }       ///< number of boundary faces in the mesh.
        inline std::size_t n_b_cells() const { return _b_cells.size(); }       ///< number of boundary cells in the mesh.

        inline std::size_t n_i_vertices() const { return _i_vertices.size(); } ///< number of internal vertices in the mesh.
        inline std::size_t n_i_edges() const { return _i_edges.size(); }       ///< number of internal edges in the mesh.
        inline std::size_t n_i_faces() const { return _i_edges.size(); }       ///< number of internal faces in the mesh.
        inline std::size_t n_i_cells() const { return _i_cells.size(); }       ///< number of internal cells in the mesh.

        inline std::vector<Vertex *> get_vertices() const { return _vertices; } ///< lists the vertices in the mesh.
        inline std::vector<Edge *> get_edges() const { return _edges; }         ///< lists the edges in the mesh.
        inline std::vector<Face *> get_faces() const { return _edges; }         ///< lists the faces in the mesh.
        inline std::vector<Cell *> get_cells() const { return _cells; }         ///< lists the cells in the mesh.

        inline std::vector<Vertex *> get_b_vertices() const { return _b_vertices; } ///< lists the boundary vertices in the mesh.
        inline std::vector<Edge *> get_b_edges() const { return _b_edges; }         ///< lists the boundary edges in the mesh.
        inline std::vector<Face *> get_b_faces() const { return _b_edges; }         ///< lists the boundary faces in the mesh.
        inline std::vector<Cell *> get_b_cells() const { return _b_cells; }         ///< lists the boundary cells in the mesh.

        inline std::vector<Vertex *> get_i_vertices() const { return _i_vertices; } ///< lists the internal vertices in the mesh.
        inline std::vector<Edge *> get_i_edges() const { return _i_edges; }         ///< lists the internal edges in the mesh.
        inline std::vector<Face *> get_i_faces() const { return _i_edges; }         ///< lists the internal faces in the mesh.
        inline std::vector<Cell *> get_i_cells() const { return _i_cells; }         ///< lists the internal cells in the mesh.

        void add_vertex(Vertex *vertex) ///<  adds a vertex to the mesh
        {
            assert(std::find(_vertices.begin(), _vertices.end(), vertex) == _vertices.end());
            _vertices.push_back(vertex);
        }

        void add_edge(Edge *edge) ///<  adds a edge to the mesh
        {
            assert(std::find(_edges.begin(), _edges.end(), edge) == _edges.end());
            _edges.push_back(edge);
        }

        void add_face(Face *face) ///<  adds a face to the mesh
        {
            assert(std::find(_edges.begin(), _edges.end(), face) == _edges.end());
            _edges.push_back(face);
        }

        void add_cell(Cell *cell) ///<  adds a cell to the mesh
        {
            assert(std::find(_cells.begin(), _cells.end(), cell) == _cells.end());
            _cells.push_back(cell);
        }

        void add_b_vertex(Vertex *vertex) ///<  adds a boundary vertex to the mesh
        {
            assert(std::find(_b_vertices.begin(), _b_vertices.end(), vertex) == _b_vertices.end());
            _b_vertices.push_back(vertex);
        }

        void add_b_edge(Edge *edge) ///<  adds a boundary edge to the mesh
        {
            assert(std::find(_b_edges.begin(), _b_edges.end(), edge) == _b_edges.end());
            _b_edges.push_back(edge);
        }

        void add_b_face(Face *face) ///<  adds a boundary face to the mesh
        {
            assert(std::find(_b_edges.begin(), _b_edges.end(), face) == _b_edges.end());
            _b_edges.push_back(face);
        }

        void add_b_cell(Cell *cell) ///<  adds a boundary cell to the mesh
        {
            assert(std::find(_b_cells.begin(), _b_cells.end(), cell) == _b_cells.end());
            _b_cells.push_back(cell);
        }

        void add_i_vertex(Vertex *vertex) ///<  adds an internal vertex to the mesh
        {
            assert(std::find(_i_vertices.begin(), _i_vertices.end(), vertex) == _i_vertices.end());
            _i_vertices.push_back(vertex);
        }

        void add_i_edge(Edge *edge) ///<  adds an internal edge to the mesh
        {
            assert(std::find(_i_edges.begin(), _i_edges.end(), edge) == _i_edges.end());
            _i_edges.push_back(edge);
        }

        void add_i_face(Face *face) ///<  adds an internal face to the mesh
        {
            assert(std::find(_i_edges.begin(), _i_edges.end(), face) == _i_edges.end());
            _i_edges.push_back(face);
        }

        void add_i_cell(Cell *cell) ///<  adds an internal cell to the mesh
        {
            assert(std::find(_i_cells.begin(), _i_cells.end(), cell) == _i_cells.end());
            _i_cells.push_back(cell);
        }

        // Note that all these assume that a MeshObject's index is equal to its position in the mesh!!
        inline Vertex *vertex(std::size_t index) const
        {
            assert(index < _vertices.size());
            return _vertices[index];
        } ///<  get a constant pointer to a vertex using its global index
        inline Edge *edge(std::size_t index) const
        {
            assert(index < _edges.size());
            return _edges[index];
        } ///<  get a constant pointer to a edge using its global index
        inline Face *face(std::size_t index) const
        {
            assert(index < _edges.size());
            return _edges[index];
        } ///<  get a constant pointer to a face using its global index
        inline Cell *cell(std::size_t index) const
        {
            assert(index < _cells.size());
            return _cells[index];
        } ///<  get a constant pointer to a cell using its global index

        inline Vertex *b_vertex(std::size_t index) const
        {
            assert(index < _b_vertices.size());
            return _b_vertices[index];
        } ///<  get a constant pointer to a boundary vertex using an index
        inline Edge *b_edge(std::size_t index) const
        {
            assert(index < _b_edges.size());
            return _b_edges[index];
        } ///<  get a constant pointer to boundary a edge using an index
        inline Face *b_face(std::size_t index) const
        {
            assert(index < _b_edges.size());
            return _b_edges[index];
        } ///<  get a constant pointer to boundary a face using an index
        inline Cell *b_cell(std::size_t index) const
        {
            assert(index < _b_cells.size());
            return _b_cells[index];
        } ///<  get a constant pointer to boundary a cell using an index

        inline Vertex *i_vertex(std::size_t index) const
        {
            assert(index < _i_vertices.size());
            return _i_vertices[index];
        } ///<  get a constant pointer to an internal vertex using an index
        inline Edge *i_edge(std::size_t index) const
        {
            assert(index < _i_edges.size());
            return _i_edges[index];
        } ///<  get a constant pointer to an internal edge using an index
        inline Face *i_face(std::size_t index) const
        {
            assert(index < _i_edges.size());
            return _i_edges[index];
        } ///<  get a constant pointer to an internal face using an index
        inline Cell *i_cell(std::size_t index) const
        {
            assert(index < _i_cells.size());
            return _i_cells[index];
        } ///<  get a constant pointer to an internal cell using an index

        std::vector<double> regularity()
        {
            /// Regularity factor =
            ///   1st component: maximum of
            ///      * diameter of cell / (measure of cell)^{1/dim}
            ///      * diameter of cell / diameter of face  [for each face of the cell]
            ///
            ///   2nd component: evaluation of max of ratio "diam of cell / radius ball inscribed in cell"

            std::vector<std::vector<double>> reg_cell(n_cells(), {0.0, 0.0});
            std::size_t count = 0;
            for (auto &T : _cells)
            {
                double hT = T->diam();
                VectorRd xT = T->center_mass();

                reg_cell[count][0] = hT / pow(T->measure(), 1.0 / DIMENSION);

                double rhoT = hT;
                std::vector<Face *> faces = T->get_faces();
                for (auto &F : faces)
                {
                    double hF = F->diam();
                    VectorRd xF = F->center_mass();
                    VectorRd nTF = F->normal(); // sign does not matter

                    reg_cell[count][0] = std::max(reg_cell[count][0], hT / hF);

                    rhoT = std::min(rhoT, std::abs((xT - xF).dot(nTF))); // If xT is not in T, is this really a good measure?
                }
                reg_cell[count][1] = hT / rhoT;
                ++count; // could just use iterators
            }

            std::vector<double> value(2, 0.0);
            for (size_t iT = 0; iT < n_cells(); iT++)
            {
                value[0] = std::max(value[0], reg_cell[iT][0]);
                value[1] = std::max(value[1], reg_cell[iT][1]);
            }

            return value;
        }

        void renum(const char B, const std::vector<size_t> new_to_old)
        {

            switch (B)
            {
            case 'C':
            {
                std::vector<Cell *> old_index = _cells;
                for (size_t i = 0; i < _cells.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _cells[i] = old_index[new_to_old[i]];
                }
                break;
            }

            case 'F':
            {
                std::vector<Face *> old_index = _edges;
                for (size_t i = 0; i < _edges.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _edges[i] = old_index[new_to_old[i]];
                }
                break;
            }

            case 'E':
            {
                std::vector<Edge *> old_index = _edges;
                for (size_t i = 0; i < _edges.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _edges[i] = old_index[new_to_old[i]];
                }
                break;
            }

            case 'V':
            {
                std::vector<Vertex *> old_index = _vertices;
                for (size_t i = 0; i < _vertices.size(); i++)
                {
                    old_index[new_to_old[i]]->set_global_index(i);
                    _vertices[i] = old_index[new_to_old[i]];
                }
                break;
            }
            }
        }

        size_t find_cell(const VectorRd x) const ///<  returns the index of the cell containing the point x
        {

            // Locate neighbouring cells
            std::vector<size_t> neighbouring_cells;
            for (Cell *T : get_cells())
            {
                if ((T->center_mass() - x).norm() <= T->diam())
                {
                    neighbouring_cells.push_back(T->global_index());
                }
            }

            // In neighbouring cells, find one such that x is contained in one of the subtriangles
            size_t i = 0;
            bool found = false;
            size_t iT = 0;
            while (i < neighbouring_cells.size() && !found)
            {
                iT = neighbouring_cells[i];
                Cell *T = cell(iT);
    
                for(auto& simplex : T->get_simplices())
                {
                    double area1 = signed_area(x, simplex[0], simplex[1]);
                    double area2 = signed_area(x, simplex[1], simplex[2]);
                    double area3 = signed_area(x, simplex[2], simplex[0]);

                    found = (area1 >= 0.0 && area2 >= 0.0 && area3 >= 0.0) || (area1 <= 0.0 && area2 <= 0.0 && area3 <= 0.0); // if all non neg or non pos must be in or on triangle
                    if(found) break;
                }
                ++i;
            }
            assert(found);

            return iT;
        }

        void plot_simplices(std::ofstream *out)
        {
            for (auto &cell : _cells)
            {
                cell->plot_simplices(out);
            }
        }

        void plot_mesh(std::ofstream *out)
        {
            for (auto &edge : _edges)
            {
                edge->plot(out);
            }
        }

    private:
        std::string _mesh_name;

        std::vector<Vertex *> _vertices;
        std::vector<Edge *> _edges;
        // std::vector<Face*> _faces;
        std::vector<Cell *> _cells;

        std::vector<Vertex *> _b_vertices;
        std::vector<Edge *> _b_edges;
        // std::vector<Face*> _b_faces;
        std::vector<Cell *> _b_cells;

        std::vector<Vertex *> _i_vertices;
        std::vector<Edge *> _i_edges;
        // std::vector<Face*> _i_faces;
        std::vector<Cell *> _i_cells;
    };
} // namespace MeshND

#endif
