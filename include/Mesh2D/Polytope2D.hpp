// standard libraries
#include <vector>    // std::vector
#include <algorithm> // std::find
#include <fstream> // std::ofstream
 
// external linking to Eigen is required
#include <Eigen/Dense> // Eigen::Matrix, Eigen::MatrixXd, Eigen::FullPivLU

// path to Math specified so that external linking is NOT required
#include <../Math/math.hpp> // Math::factorial, Math::sgn

#ifndef _POLYTOPE2D_HPP
#define _POLYTOPE2D_HPP

namespace Mesh2D
{
    inline constexpr int DIMENSION = 2;

    using VectorRd = Eigen::Matrix<double, DIMENSION, 1>;
    using VectorZd = Eigen::Matrix<int, DIMENSION, 1>;

    template <size_t object_dim>
    using Simplex = std::array<VectorRd, object_dim + 1>;

    template <size_t object_dim>
    using Simplices = std::vector<Simplex<object_dim>>;

    ///< Method to find the center mass of an arbitrary simplex in arbitrary space
    template <size_t object_dim>
    VectorRd simplex_center_mass(Simplex<object_dim> simplex);

    ///< Method to find the Lebesgue measure of an arbitrary simplex in arbitrary space
    template <size_t object_dim>
    double simplex_measure(Simplex<object_dim> simplex);


    /*!
     * \addtogroup Mesh
     * @brief Class providing a template for all objects (Vertex, Edge, Face, Cell) that appear in a Mesh
     * @{
     */

    // ----------------------------------------------------------------------------
    //                         Polytope class definition
    // ----------------------------------------------------------------------------

    /** Polytope is a templated class describing polytopes in two dimensional space.
        It takes in one template parameter - object_dim - the dimension of the polytope.
        object_dim must be less than or equal to two
     **/

    template <size_t object_dim>
    class Polytope
    {
    public:
        ///< Constructor for a Polytope defined by its simplices
        Polytope(size_t index, Simplices<object_dim> simplices);

        ///< Constructor for a Simplex
        Polytope(size_t index, Simplex<object_dim> simplex);

        ///< Constructor for a Polytope defined by a single coordinate (Vertex)
        Polytope(size_t index, VectorRd vertex);

        ///< Null constructor
        Polytope();

        ///< Destructor
        ~Polytope();

        // inline size_t dim() const { return object_dim; }

        inline size_t global_index() const { return _index; }                     ///< Return the global index of the Polytope
        inline double diam() const { return _diameter; }                          ///< Return the diameter of the Polytope
        inline VectorRd center_mass() const { return _center_mass; }              ///< Return the center mass of the Polytope
        inline double measure() const { return _measure; }                        ///< Return the Lebesgue measure of the Polytope
        inline Simplices<object_dim> get_simplices() const { return _simplices; } ///< Return the simplices making up the Polytope
        inline void set_global_index(const size_t idx) { _index = idx; }          ///< Set the global index
        inline bool is_boundary() const { return _is_boundary; }                  ///< Return true if Polytope is a boundary object, false otherwise
        inline void set_boundary(bool val) { _is_boundary = val; }                ///< Set the boundary value of the Polytope

        inline std::vector<Polytope<0> *> get_vertices() const { return _vertices; }       ///< Return the vertices of the Polytope
        inline std::vector<Polytope<1> *> get_edges() const { return _edges; }             ///< Return the edges of the Polytope
        inline std::vector<Polytope<DIMENSION - 1> *> get_faces() const { return _edges; } ///< Return the faces of the Polytope
        inline std::vector<Polytope<DIMENSION> *> get_cells() const { return _cells; }     ///< Return the cells of the Polytope

        inline size_t n_vertices() const { return _vertices.size(); } ///< Return the number of vertices of the Polytope
        inline size_t n_edges() const { return _edges.size(); }       ///< Return the number of edges of the Polytope
        inline size_t n_faces() const { return _edges.size(); }       ///< Return the number of faces of the Polytope
        inline size_t n_cells() const { return _cells.size(); }       ///< Return the number of cells of the Polytope

        Polytope<0> *vertex(const size_t i) const;           ///< Return the i-th vertex of the Polytope
        Polytope<1> *edge(const size_t i) const;             ///< Return the i-th edge of the Polytope
        Polytope<DIMENSION - 1> *face(const size_t i) const; ///< Return the i-th face of the Polytope
        Polytope<DIMENSION> *cell(const size_t i) const;     ///< Return the i-th cell of the Polytope

        void add_vertex(Polytope<0> *vertex);         ///< Add a vertex to the Polytope
        void add_edge(Polytope<1> *edge);             ///< Add an edge to the Polytope
        void add_face(Polytope<DIMENSION - 1> *face); ///< Add a face to the Polytope
        void add_cell(Polytope<DIMENSION> *cell);     ///< Add a cell to the Polytope

        int index_vertex(const Polytope<0> *vertex) const;         ///< Returns the local index of a vertex
        int index_edge(const Polytope<1> *edge) const;             ///< Returns the local index of an edge
        int index_face(const Polytope<DIMENSION - 1> *face) const; ///< Returns the local index of a face
        int index_cell(const Polytope<DIMENSION> *cell) const;     ///< Returns the local index of a cell

        VectorRd coords() const; ///< Return the coordinates of a Vertex

        VectorRd face_normal(const size_t face_index) const; ///< Return the outer normal of a Cell towards the Face located at face_index
        VectorRd edge_normal(const size_t edge_index) const; ///< Return the edge normal of a 2D object

        int face_orientation(const size_t face_index) const; ///< Return the orientation of a Face
        int edge_orientation(const size_t edge_index) const; ///< Return the orientation of a Edge
        int vertex_orientation(const size_t vertex_index) const; ///< Return the orientation of a Vertex

        VectorRd normal() const;  ///< Return the normal of a Face
        VectorRd tangent() const; ///< Return the tangent of a Edge

        void construct_face_normals(); ///< Set the directions of the face normals of a cell

        void plot_simplices(std::ofstream *out) const; ///< Plot the simplices to out
        void plot(std::ofstream *out) const; ///< Plot the polytope to out

    private:
        size_t _index;
        VectorRd _center_mass;
        double _measure;
        double _diameter;
        bool _is_boundary;
        Simplices<object_dim> _simplices;
        VectorRd _normal;                  // uninitialised unless object_dim == DIMENSION - 1 (face)
        VectorRd _tangent;                  // uninitialised unless object_dim == 1 (edge)
        std::vector<int> _face_directions; // empty unless object_dim == DIMENSION (cell)

        std::vector<Polytope<0> *> _vertices;
        std::vector<Polytope<1> *> _edges;
        // std::vector<Polytope<DIMENSION - 1>*> _faces;
        std::vector<Polytope<DIMENSION> *> _cells;
    };

    ///< A Vertex is a Polytope with object_dim = 0
    using Vertex = Polytope<0>;

    ///< An Edge is a Polytope with object_dim = 1
    using Edge = Polytope<1>;

    ///< A Face is a Polytope with object_dim = DIMENSION - 1
    using Face = Polytope<DIMENSION - 1>;

    ///< A Cell is a Polytope with object_dim = DIMENSION
    using Cell = Polytope<DIMENSION>;

    //@}

    // -----------------------------------------------------------------
    // ----------- Implementation of templated functions ---------------
    // -----------------------------------------------------------------

    template <size_t object_dim>
    VectorRd simplex_center_mass(Simplex<object_dim> simplex)
    {
        VectorRd center_mass = VectorRd::Zero();

        for (auto &coord : simplex)
        {
            for (size_t i = 0; i < DIMENSION; ++i)
            {
                center_mass(i) += coord(i);
            }
        }
        return center_mass / (object_dim + 1);
    }

    template <size_t object_dim>
    double simplex_measure(Simplex<object_dim> simplex)
    {
        Eigen::MatrixXd CM_mat = Eigen::MatrixXd::Zero(object_dim + 2, object_dim + 2);
        for (size_t i = 0; i < object_dim + 1; ++i)
        {
            VectorRd i_coords = simplex[i];
            for (size_t j = i + 1; j < object_dim + 1; ++j)
            {
                VectorRd j_coords = simplex[j];
                double norm = (j_coords - i_coords).norm();
                CM_mat(i, j) = norm * norm;
                CM_mat(j, i) = CM_mat(i, j);
            }
            CM_mat(i, object_dim + 1) = 1;
            CM_mat(object_dim + 1, i) = 1;
        }

        // Calculate Cayley-Menger determinant
        double det = CM_mat.determinant();

        double scaling = std::pow(Math::factorial(object_dim), 2) * std::pow(2, object_dim);

        return std::sqrt(std::abs(det / scaling));
    }

    template <size_t object_dim>
    Polytope<object_dim>::Polytope(size_t index, Simplices<object_dim> simplices)
        : _index(index), _simplices(simplices)
    {
        assert(object_dim <= DIMENSION);

        if (object_dim == 0 || object_dim == 1)
        {
            assert(simplices.size() == 1);
        }

        _is_boundary = false; // by default, is_boundary is set to false. If polytope is on boundary, set is_boundary to true upon construction of the mesh.

        _measure = 0.0;
        _center_mass = VectorRd::Zero();
        std::vector<VectorRd> vertex_coords;
        for (auto &simplex : simplices)
        {
            // assert(simplex.size() == _dim + 1);
            double measure_of_simplex = simplex_measure<object_dim>(simplex);
            _measure += measure_of_simplex;
            _center_mass += measure_of_simplex * simplex_center_mass<object_dim>(simplex);
            for (auto &coord : simplex)
            {
                if (std::find(vertex_coords.begin(), vertex_coords.end(), coord) == vertex_coords.end())
                {
                    vertex_coords.push_back(coord); // if coord not already in vertex_coords, add it
                }
            }
        }
        _center_mass /= _measure;

        _diameter = 0.0;
        for (auto it = vertex_coords.begin(); it != vertex_coords.end(); ++it)
        {
            for (auto jt = it; jt != vertex_coords.end(); ++jt)
            {
                _diameter = std::max(_diameter, (*it - *jt).norm()); 
            }
        }

        if (object_dim == DIMENSION - 1) // find the tangent and normal
        {
            _tangent = (this->get_simplices()[0][1] - this->get_simplices()[0][0]).normalized();
            _normal = VectorRd(-_tangent(1), _tangent(0));
        }
    }

    template <size_t object_dim>
    Polytope<object_dim>::Polytope(size_t index, Simplex<object_dim> simplex) // simplex
        : Polytope(index, Simplices<object_dim>(1, simplex))
    {
    }

    template <size_t object_dim>
    Polytope<object_dim>::Polytope(size_t index, VectorRd vertex) // vertex
        : Polytope(index, Simplex<object_dim>({vertex}))
    {
        assert(object_dim == 0); // must be a zero dimensional object to use vertex constructor
    }

    template <size_t object_dim>
    Polytope<object_dim>::Polytope() {}

    template <size_t object_dim>
    Polytope<object_dim>::~Polytope() {}

    template <size_t object_dim>
    void Polytope<object_dim>::add_vertex(Polytope<0> *vertex)
    {
        assert(std::find(_vertices.begin(), _vertices.end(), vertex) == _vertices.end()); // ensure vertex does not already exist in _vertices
        _vertices.push_back(vertex);
    }

    template <size_t object_dim>
    void Polytope<object_dim>::add_edge(Polytope<1> *edge)
    {
        assert(std::find(_edges.begin(), _edges.end(), edge) == _edges.end());
        _edges.push_back(edge);
    }

    template <size_t object_dim>
    void Polytope<object_dim>::add_face(Polytope<DIMENSION - 1> *face)
    {
        assert(std::find(_edges.begin(), _edges.end(), face) == _edges.end());
        _edges.push_back(face);
    }

    template <size_t object_dim>
    void Polytope<object_dim>::add_cell(Polytope<DIMENSION> *cell)
    {
        assert(std::find(_cells.begin(), _cells.end(), cell) == _cells.end());
        _cells.push_back(cell);
    }

    template <size_t object_dim>
    Polytope<0> *Polytope<object_dim>::vertex(const size_t i) const
    {
        assert(i < _vertices.size());
        return _vertices[i];
    }

    template <size_t object_dim>
    Polytope<1> *Polytope<object_dim>::edge(const size_t i) const
    {
        assert(i < _edges.size());
        return _edges[i];
    }

    template <size_t object_dim>
    Polytope<DIMENSION - 1> *Polytope<object_dim>::face(const size_t i) const
    {
        assert(i < _edges.size());
        return _edges[i];
    }

    template <size_t object_dim>
    Polytope<DIMENSION> *Polytope<object_dim>::cell(const size_t i) const
    {
        assert(i < _cells.size());
        return _cells[i];
    }

    template <size_t object_dim>
    void Polytope<object_dim>::construct_face_normals() // not very efficient - probably room for improvement
    {
        assert(object_dim == DIMENSION);
        for (size_t iF = 0; iF < _edges.size(); ++iF)
        {
            VectorRd normal = _edges[iF]->normal();
            Simplex<object_dim - 1> face_simplex = _edges[iF]->get_simplices()[0];
            VectorRd center = VectorRd::Zero();
            double count;
            for (auto &cell_simplex : this->get_simplices())
            {
                count = 0;
                for (size_t i = 0; (i < cell_simplex.size()) && count < 2; ++i)
                {
                    if (std::find(face_simplex.begin(), face_simplex.end(), cell_simplex[i]) == face_simplex.end()) // requires numerical precision
                    {
                        ++count;
                    }
                }
                if (count == 1) // only don't share one coordinate
                {
                    center = simplex_center_mass<object_dim>(cell_simplex);
                    break;
                }
            }
            assert(count == 1);
            //    _face_directions.push_back(Math::sgn((_faces[iF]->center_mass() - _center_mass).dot(normal))); // star shaped wrt center mass
            _face_directions.push_back(Math::sgn((_edges[iF]->center_mass() - center).dot(normal)));
            assert(_face_directions[iF] != 0);
        }
    }

    template <size_t object_dim>
    VectorRd Polytope<object_dim>::face_normal(const size_t face_index) const
    {
        assert(object_dim == DIMENSION);
        assert(face_index < _edges.size());
        return _face_directions[face_index] * _edges[face_index]->normal();
    }

    template <size_t object_dim>
    VectorRd Polytope<object_dim>::edge_normal(const size_t edge_index) const
    {
        return this->face_normal(edge_index);
    }

    template <size_t object_dim>
    VectorRd Polytope<object_dim>::coords() const // only for vertices
    {
        assert(object_dim == 0);            // can only return coordinate of a vertex
        return this->get_simplices()[0][0]; // only has one simplex, with one coordinate
    }

    template <size_t object_dim>
    VectorRd Polytope<object_dim>::normal() const
    {
        assert(object_dim == DIMENSION - 1);
        return _normal;
    }

    template <size_t object_dim>
    VectorRd Polytope<object_dim>::tangent() const
    {
        assert(object_dim == 1);
        return _tangent;
        // return (this->get_vertices()[1]->coords() - this->get_vertices()[0]->coords()).normalized();
    }

    template <size_t object_dim>
    int Polytope<object_dim>::index_vertex(const Polytope<0> *vertex) const
    {
        auto itr = std::find(_vertices.begin(), _vertices.end(), vertex);
        if (itr != _vertices.end())
        {
            return itr - _vertices.begin();
        }
        else
        {
            throw "Vertex not found";
        }
    }

    template <size_t object_dim>
    int Polytope<object_dim>::index_edge(const Polytope<1> *edge) const
    {
        auto itr = std::find(_edges.begin(), _edges.end(), edge);
        if (itr != _edges.end())
        {
            return itr - _edges.begin();
        }
        else
        {
            throw "Edge not found";
        }
    }

    template <size_t object_dim>
    int Polytope<object_dim>::index_face(const Polytope<DIMENSION - 1> *face) const
    {
        auto itr = std::find(_edges.begin(), _edges.end(), face);
        if (itr != _edges.end())
        {
            return itr - _edges.begin();
        }
        else
        {
            throw "Face not found";
        }
    }

    template <size_t object_dim>
    int Polytope<object_dim>::index_cell(const Polytope<DIMENSION> *cell) const
    {
        auto itr = std::find(_cells.begin(), _cells.end(), cell);
        if (itr != _cells.end())
        {
            return itr - _cells.begin();
        }
        else
        {
            throw "Cell not found";
        }
    }

    template <size_t object_dim>
    int Polytope<object_dim>::face_orientation(const size_t face_index) const
    {
        assert(object_dim == DIMENSION);
        assert(face_index < _edges.size());
        return _face_directions[face_index];
    }

    template <size_t object_dim>
    int Polytope<object_dim>::edge_orientation(const size_t edge_index) const
    {
        assert(object_dim == 2);
        assert(edge_index < _edges.size());
        return _face_directions[edge_index];
    }

    template <size_t object_dim>
    int Polytope<object_dim>::vertex_orientation(const size_t vertex_index) const
    {
        assert(object_dim == 1);
        assert(vertex_index < _vertices.size());
        return ( (_vertices[vertex_index]->coords()-_center_mass).dot(this->tangent()) > 0 ? 1 : -1);
    }      
    
    template <size_t object_dim>
    void Polytope<object_dim>::plot_simplices(std::ofstream *out) const 
    {
        assert(object_dim > 0);
        for (auto &simplex : _simplices)
        {
            for (size_t i = 0; i < object_dim + 1; ++i)
            {
                size_t i_next = (i + 1) % (object_dim + 1);
                *out << simplex[i](0) << " " << simplex[i](1) << std::endl;
                *out << simplex[i_next](0) << " " << simplex[i_next](1) << std::endl;
                *out << std::endl;
            }
        }
    }

    template <size_t object_dim>
    void Polytope<object_dim>::plot(std::ofstream *out) const 
    {
        assert(object_dim > 0);
        for (size_t i = 0; i < _vertices.size(); ++i)
        {
            size_t i_next = (i + 1) % _vertices.size();
            *out << _vertices[i]->coords()(0) << " " << _vertices[i]->coords()(1) << std::endl;
            *out << _vertices[i_next]->coords()(0) << " " << _vertices[i_next]->coords()(1) << std::endl;
            *out << std::endl;
        }
    }

} // namespace Mesh2D

#endif
