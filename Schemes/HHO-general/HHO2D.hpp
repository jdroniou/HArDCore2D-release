// Helper methods and assemble and solve routines for implementing 2D Hybrid High Order schemes.

#include <basis.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <parallel_for.hpp>

#include <mesh_builder.hpp>
#include <hybridcore.hpp>
#include <elementquad.hpp>

#include <TestCase/TestCase.hpp>
#include <BoundaryConditions/BoundaryConditions.hpp>

#include <boost/timer/timer.hpp>
#include "vtu_writer.hpp"

#ifndef _HHO2D_HPP
#define _HHO2D_HPP

/*!	
 * @defgroup HHO2D 
 * @brief Class providing helper methods and assemble and solve routines for general 2D HHO schemes
 */

/** The HHO2D class contains methods applicable to a general 2D HHO model. It contains assemble and
  * solve routines usable by any scheme. The global operator and load vector are passed to the model
  * through setters and are then accessible by assemble and solve. The class also contains helper
  * methods for writing HHO schemes. 
 **/

/*!
 * \addtogroup HHO2D
 * @{
 */

// ----------------------------------------------------------------------------
//                            HHO2D class definition
// ----------------------------------------------------------------------------


typedef std::function<Eigen::MatrixXd(Cell *, ElementQuad &)> MatrixFType; ///< type for the global operator as a function of a Cell and an ElementQuad
typedef std::function<Eigen::VectorXd(Cell *, ElementQuad &)> VectorFType; ///< type for the load vector as a function of a Cell and an ElementQuad

class HHO2D
{
public:
    /**@brief Class constructor: initialises the model by providing a HybridCore object, and the 
      * exact solution and boundary conditions of the model.
     **/
    HHO2D(
        HybridCore &,                  ///< A reference to the HybridCore object containing mesh data
        const size_t,                  ///< Cell polynomial degree
        const size_t,                  ///< Edge polynomial degree
        const bool use_threads = true, ///< Optional argument to indicate if threads should be used
        size_t doeT = 0,               ///< Optional argument to set cell quadrature degree. Default is L + K + 1
        size_t doeF = 0                ///< Optional argument to set edge quadrature degree. Default is 2 * K + 1
    );

    /// A general assemble routine that calculates the statically condensed matrices required by solve
    void assemble();

    /// Solves the statically condensed system
    UVector solve();

    /// Solves the system when the model is ill posed (not yet running)
    UVector neumann_solve();

    /// Returns the energy norm of a given UVector
    double energy_norm(const UVector);

    /// Set the global operator
    void set_global_operator(const MatrixFType &);

    /// Set the load vector
    void set_load_vector(const VectorFType &);

    /// Plot the numerical and exact solutions
    void plot(
        const std::string,    ///< Plot file
        const UVector &,      ///< Numerical solution
        const FType<double> & ///< Exact solution
    );

    /// Returns the standard load vector (f, v_T)_T with no Neumann boundary conditions
    VectorFType standard_load_vector(
        const CellFType<double> & ///< Source term
    );

    /// Returns the standard load vector (f, v_T)_T
    VectorFType standard_load_vector(
        const CellFType<double> &,   ///< Source term
        const CellFType<VectorRd> &, ///< Gradient to be dotted with the normal vector on each Neumann edge
        const BoundaryConditions &   ///< A reference to the BC to determine Neumann edges
    );

    /// Returns the standard load vector (f, v_T)_T
    VectorFType standard_load_vector(
        const CellFType<double> &, ///< Source term
        const FType<double> &,     ///< Function on the Neumann BC edges
        const BoundaryConditions & ///< A reference to the BC to determine Neumann edges
    );

    /// Set the Dirichlet boundary conditions
    void set_dirichlet(
        const FType<double> &, ///< Function on the Dirichlet BC edges
        const size_t           ///< Number of Dirichlet edges
    );

    /// Set the Dirichlet boundary condition to zero
    void set_dirichlet(
        const size_t ///< Number of Dirichlet edges
    );

    /// Return the (statically condensed) matrix system
    inline Eigen::SparseMatrix<double> get_SysMat()
    {
        return SysMat;
    };

    /// CPU time to assemble the scheme
    inline double get_assembly_time() const
    {
        return double(assembly_time) * pow(10, -9);
    };

    /// CPU time to solve the scheme
    inline double get_solving_time() const
    {
        return double(solving_time) * pow(10, -9);
    };

    /// Residual after solving the scheme
    inline double get_solving_error() const
    {
        return solving_error;
    };

private:
    HybridCore &m_hho;
    const size_t m_L;
    const size_t m_K;
    const bool m_use_threads;
    size_t m_doeT;
    size_t m_doeF;

    MatrixFType global_operator;
    VectorFType load_vector;

    const Mesh *mesh_ptr = m_hho.get_mesh();

    const size_t n_cells = mesh_ptr->n_cells();
    const size_t n_edges = mesh_ptr->n_edges();

    const size_t n_local_cell_dofs = DimPoly<Cell>(m_L);
    const size_t n_local_edge_dofs = DimPoly<Edge>(m_K);

    const size_t n_total_cell_dofs = n_local_cell_dofs * n_cells;
    const size_t n_total_edge_dofs = n_local_edge_dofs * n_edges;
    const size_t n_total_dofs = n_total_cell_dofs + n_total_edge_dofs;

    double solving_error;
    size_t solving_time;
    size_t assembly_time;

    std::vector<Eigen::MatrixXd> AT;

    Eigen::VectorXd GlobRHS = Eigen::VectorXd::Zero(n_total_edge_dofs);
    Eigen::VectorXd ScRHS = Eigen::VectorXd::Zero(n_total_cell_dofs);
    Eigen::VectorXd UDir;

    Eigen::SparseMatrix<double> GlobMat;
    Eigen::SparseMatrix<double> SysMat;
    Eigen::SparseMatrix<double> ScBeMat;
};

#endif
