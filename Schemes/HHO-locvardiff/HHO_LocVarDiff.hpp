// Implementation of the HHO scheme in 2D for the diffusion equation
//
//   { -div(K \grad(u)) = f,       inside Omega
//   { K \grad(u) . nTF = g,       on GammaN
//   {                 u = g,       on GammaD
// 
//  At the moment, only pure Neumann or pure Dirichlet
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
 *
 *  This implementation of HHO was developped following the principles described in 
 * Appendix B of the book
 *
 * The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications. 
 *  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
 *  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
 *  url: https://hal.archives-ouvertes.fr/hal-02151813.
 *
 * If you use this code or part of it for a scientific publication, please cite the book
 * above as a reference for the implementation.
 *
 */

#ifndef _HHO_LOCVARDIFF_HPP
#define _HHO_LOCVARDIFF_HPP

#include <functional>
#include <utility>
#include <iostream>

#include <boost/timer/timer.hpp>

// Matrices and linear solvers
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include "Eigen/MA41.cpp"

#include <vertex.hpp>
#include <hybridcore.hpp>
#include <elementquad.hpp>
#include <parallel_for.hpp>
#include <TestCase/BoundaryConditions.hpp>
#include "TestCase/TestCase.hpp"

/*!
 * @defgroup HHO_LocVarDiff
 * @brief HHO scheme for diffusion equation -div(Diff grad u)=f, with Diff possibly varying in each cell
 */

namespace HArDCore2D {

  /*!
   *  @addtogroup HHO_LocVarDiff
   * @{
   */
  // ----------------------------------------------------------------------------
  //                            Class definition
  // ----------------------------------------------------------------------------

  /// The HHO_LocVarDiff class provides tools to implement the HHO method for the diffusion problem

  class HHO_LocVarDiff {

  public:

    ///@brief Constructor of the class
    HHO_LocVarDiff(
       HybridCore& hho,        ///< reference to the mesh
       size_t K,              ///< degree of polynomials on edges
       int L,                 ///< degree of polynomials in cells
       CellFType<MatrixRd> kappa,   ///< diffusion tensor
       size_t deg_kappa,              ///< polynomial degree of the diffusion tensor
       CellFType<double> source,  ///< source term
       BoundaryConditions BC,                     ///< type of boundary conditions
       FType<double> exact_solution,   ///< exact solution
       CellFType<VectorRd> grad_exact_solution,   ///< gradient of the exact solution
       std::string solver_type,    ///< type of solver to use for the global system (bicgstab at the moment)
       bool use_threads,    ///< optional argument determining if local parallelisation is to be used
       std::ostream & output = std::cout    ///< optional argument for output of messages
       );

    /// Assemble and solve the scheme
    void assemble(HybridCore& hho);
    UVector solve(HybridCore& hho);

    /// Discrete energy norm (associated to the diffusion operator) of an hybrid function
    double EnergyNorm(HybridCore& hho, const UVector Xh); 

    /// Return the (statically condensed) matrix system
    Eigen::SparseMatrix<double> get_SysMat(){
      return SysMat;
    }

    /// cpu time to assemble the scheme
    inline double get_assembly_time() const
      {
        return double(_assembly_time) * pow(10, -9);
      }; 
    /// cpu time to solve the scheme
    inline double get_solving_time() const
      {
        return double(_solving_time) * pow(10, -9);
      };  
    /// residual after solving the scheme
    inline double get_solving_error() const
      {
        return _solving_error;
      };
    /// various intermediate assembly times
    inline double get_itime(size_t idx) const      
      {
        return double(_itime[idx]) * pow(10, -9);
      };   

    /// Number of DOFs in each cell
    const inline size_t get_nlocal_cell_dofs() { return m_nlocal_cell_dofs; }
    /// Number of DOFs on each edge
    const inline size_t get_nlocal_edge_dofs() { return m_nlocal_edge_dofs; }
    /// Number of DOFs per cell for high-order (K+1) polynomials
    const inline size_t get_nhighorder_dofs() { return m_nhighorder_dofs; }
    /// Total number of cell DOFs over the entire mesh
    const inline size_t get_ntotal_cell_dofs() { return m_ntotal_cell_dofs; }
    /// Total number of edge DOFs over the entire mesh
    const inline size_t get_ntotal_edge_dofs() { return m_ntotal_edge_dofs; }
    /// Total number of edge DOFs for Dirichlet edges
    const inline size_t get_ndir_edge_dofs() { return m_ndir_edge_dofs; }
    /// Total number of degrees of freedom over the entire mesh
    const inline size_t get_ntotal_dofs() { return m_ntotal_dofs; }

  private:
    /// Compute the local diffusion operator in the cell iT
    Eigen::MatrixXd diffusion_operator(HybridCore& hho, const size_t iT, const ElementQuad& elquad) const;

    /// Compute the local load operator (the source term) in the cell iT
    Eigen::VectorXd load_operator(HybridCore& hho, const size_t iT, const ElementQuad &elquad) const;

    // Reference to the HybridCore structure
    HybridCore& m_hho;

    // Degrees on edges and in cells
    size_t m_K;
    int m_L;
    size_t m_Ldeg;

    // Data
    const CellFType<MatrixRd> kappa;
    size_t _deg_kappa;
    const CellFType<double> source;
    const BoundaryConditions m_BC;
    const FType<double> exact_solution;
    const CellFType<VectorRd> grad_exact_solution;
    const std::string solver_type;
    const bool m_use_threads;
    std::ostream & m_output;

    // DOFs
    const size_t m_nlocal_cell_dofs;
    const size_t m_nlocal_edge_dofs;
    const size_t m_nhighorder_dofs;
    const size_t m_ntotal_cell_dofs;
    const size_t m_ntotal_edge_dofs;
    const size_t m_ndir_edge_dofs;
    const size_t m_nnondir_edge_dofs;
    const size_t m_ntotal_dofs;

    // Local bilinear forms
    std::vector<Eigen::MatrixXd> aT;
    // Global matrix (without BC accounted for), and system matrix
    Eigen::SparseMatrix<double> GlobMat;
    Eigen::SparseMatrix<double> SysMat;
    // If static condensation (L>=0): matrix to recover cell unknowns
    // If barycentric elimination (L=-1): matrix to recover cell unknowns
    Eigen::SparseMatrix<double> ScBeMat;
    // Source terms for the system, and for recovering cell unknowns from static condensation
    Eigen::VectorXd GlobRHS;
    Eigen::VectorXd ScRHS;

    // Computation statistics
    size_t _assembly_time;
    size_t _solving_time;
    double _solving_error;
    mutable std::vector<size_t> _itime = std::vector<size_t>(10, 0);

  };

  HHO_LocVarDiff::HHO_LocVarDiff(HybridCore& hho, size_t K, int L, CellFType<MatrixRd> kappa, size_t deg_kappa, CellFType<double> source, BoundaryConditions BC, FType<double> exact_solution, CellFType<VectorRd> grad_exact_solution, std::string solver_type, bool use_threads, std::ostream & output)
    : m_hho(hho),
      m_K(K),
      m_L(L),
      m_Ldeg(std::max(L,0)),
      kappa(kappa),
      _deg_kappa(deg_kappa),
      source(source),
      m_BC(BC),
      exact_solution(exact_solution),
      grad_exact_solution(grad_exact_solution),
      solver_type(solver_type),
      m_use_threads(use_threads),
      m_output(output),
      m_nlocal_cell_dofs(DimPoly<Cell>(m_Ldeg)),
      m_nlocal_edge_dofs(DimPoly<Edge>(m_K)),
      m_nhighorder_dofs(DimPoly<Cell>(m_K+1)),
      m_ntotal_cell_dofs(m_nlocal_cell_dofs * m_hho.get_mesh()->n_cells()),
      m_ntotal_edge_dofs(m_nlocal_edge_dofs * m_hho.get_mesh()->n_edges()),
      m_ndir_edge_dofs(m_nlocal_edge_dofs * m_BC.n_dir_edges()),
      m_nnondir_edge_dofs(m_nlocal_edge_dofs * m_hho.get_mesh()->n_edges() - m_ndir_edge_dofs),
      m_ntotal_dofs(m_ntotal_cell_dofs + m_ntotal_edge_dofs),
      GlobRHS(Eigen::VectorXd::Zero(m_ntotal_edge_dofs)),
      ScRHS(Eigen::VectorXd::Zero(m_ntotal_cell_dofs))
 {
      GlobMat.resize(m_ntotal_edge_dofs, m_ntotal_edge_dofs);        
      ScBeMat.resize(m_ntotal_cell_dofs, m_ntotal_edge_dofs);
    // Do nothing
  }

  void HHO_LocVarDiff::assemble(HybridCore &hho) {

    boost::timer::cpu_timer timer;  // Time the matrix assembly
    timer.start();
    const Mesh* mesh = hho.get_mesh();

    //--------------- PREPARE SYSTEM ------------------------//

    // Global triplets for: system matrix, static condensaion/barycentric elimination
    std::vector<Eigen::Triplet<double>> triplets_GlobMat;
    std::vector<Eigen::Triplet<double>> triplets_ScBe;

    // Local bilinear form, triplets and source term (one per cell)
    aT.resize(mesh->n_cells());
    std::vector<std::vector<Eigen::Triplet<double>>> cell_triplets_GlobMat;
    std::vector<std::vector<Eigen::Triplet<double>>> cell_triplets_ScBe;
    std::vector<Eigen::VectorXd> cell_source(mesh->n_cells());
    cell_triplets_GlobMat.resize(mesh->n_cells());
    cell_triplets_ScBe.resize(mesh->n_cells());
    size_t size_triplets_GlobMat = 0;
    size_t size_triplets_ScBe = 0;

    //-------------- ASSEMBLE LOCAL CONTRIBUTIONS -------------//
    
    // Function to create local contribution between cell start and cell end-1
    std::function<void(size_t, size_t)> construct_all_local_contributions
      = [&](size_t start, size_t end)->void
      {
        for (size_t iT = start; iT < end; iT++) {
          Cell* iCell = mesh->cell(iT);

          // Total number of edge degrees of freedom local to this cell (adjacent edges to the cell)
          size_t nlocal_edges = iCell->n_edges();
          size_t edge_dofs = nlocal_edges * m_nlocal_edge_dofs;

          // Local bilinear form and source term
          size_t doeT = std::max( std::max(m_K,m_Ldeg) + m_K+1 , 2*m_K + _deg_kappa );
          size_t doeF = 2*m_K + 1;
          ElementQuad elquad(hho, iT, doeT, doeF);

          aT[iT] = diffusion_operator(hho, iT, elquad);
          Eigen::VectorXd bT = load_operator(hho, iT, elquad);

          // Local matrix and right-hand side on the edge unknowns
          Eigen::MatrixXd MatF = Eigen::MatrixXd::Zero(edge_dofs,edge_dofs);
          cell_source[iT] = Eigen::VectorXd::Zero(edge_dofs);

          if (m_L>=0) {
            // STATIC CONDENSATION OF ELEMENT UNKNOWNS

            // Perform static condensation
            Eigen::MatrixXd ATT = aT[iT].topLeftCorner(m_nlocal_cell_dofs, m_nlocal_cell_dofs);
            Eigen::MatrixXd ATF = aT[iT].topRightCorner(m_nlocal_cell_dofs, edge_dofs);
            Eigen::MatrixXd AFF = aT[iT].bottomRightCorner(edge_dofs, edge_dofs);

            Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
            invATT.compute(ATT);
                
            Eigen::MatrixXd invATT_ATF = invATT.solve(ATF);
            Eigen::VectorXd invATT_bTcell = invATT.solve(bT.head(m_nlocal_cell_dofs));
            MatF = AFF - ATF.transpose() * invATT_ATF;
                
            cell_source[iT] = bT.tail(edge_dofs) - ATF.transpose() * invATT_bTcell;
                      
            // Assemble local triplets for static condensation operator
            ScRHS.segment(iT * m_nlocal_cell_dofs, m_nlocal_cell_dofs) = invATT_bTcell;
            for (size_t i = 0; i < m_nlocal_cell_dofs; i++) {
              for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
                const size_t jF = iCell->edge(jlF)->global_index();
                for (size_t jk = 0; jk < m_nlocal_edge_dofs; jk++) {
                  const size_t jLocal = jlF * m_nlocal_edge_dofs + jk;
                  const size_t jGlobal = jF * m_nlocal_edge_dofs + jk;
                  cell_triplets_ScBe[iT].emplace_back(iT * m_nlocal_cell_dofs + i, jGlobal, invATT_ATF(i, jLocal));
                }
              }
            }
            size_triplets_ScBe += cell_triplets_ScBe[iT].size();

          } else {
            // BARYCENTRIC ELIMINATION OF ELEMENT UNKNOWNS
            // Create reduction matrix: 1+nlocal_edges * nlocal_edges matrix with the coefficients on the first row, and the identity below. When multiplied by the edge unknowns, return cell and edge unknowns
            // Note that the basis functions are constant, but not necessarily assumed to be one (which is not the case after orthonormalisation for example), which is why we have to adjust the first row.
            Eigen::MatrixXd red_matT = Eigen::MatrixXd::Zero(1+nlocal_edges,nlocal_edges);
            red_matT.row(0) = hho.compute_weights(iT);
            VectorRd xT = iCell->center_mass();
            double phiT_cst = hho.CellBasis(iT).function(0, xT);
            for (size_t ilF = 0; ilF < nlocal_edges; ilF++){
              VectorRd xF = iCell->edge(ilF)->center_mass();
              size_t iF = iCell->edge(ilF)->global_index();
              double phiF_cst = hho.EdgeBasis(iF).function(0, xF);
              red_matT(0,ilF) *= phiF_cst / phiT_cst;
            }
            red_matT.bottomRightCorner(nlocal_edges,nlocal_edges) = Eigen::MatrixXd::Identity(nlocal_edges,nlocal_edges);

            cell_source[iT] = red_matT.transpose() * bT;
            MatF = red_matT.transpose() * aT[iT] * red_matT;

            // Assemble local triplets for barycentric combination to recover cell unknown
            for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
              const size_t jF = iCell->edge(jlF)->global_index();
              const size_t jGlobal = jF * m_nlocal_edge_dofs;
              cell_triplets_ScBe[iT].emplace_back(iT, jGlobal, red_matT(0,jlF));
            }
            size_triplets_ScBe += cell_triplets_ScBe[iT].size();

          }

          // Assemble local triplets for scheme's matrix and source term
          for (size_t ilF = 0; ilF < nlocal_edges; ilF++) {
            const size_t iF = iCell->edge(ilF)->global_index();
            for (size_t ik = 0; ik < m_nlocal_edge_dofs; ik++) {
              const size_t iLocal = ilF * m_nlocal_edge_dofs + ik;
              const size_t iGlobal = iF * m_nlocal_edge_dofs + ik;
              for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
                const size_t jF = iCell->edge(jlF)->global_index();
                for (size_t jk = 0; jk < m_nlocal_edge_dofs; jk++) {
                  const size_t jLocal = jlF * m_nlocal_edge_dofs + jk;
                  const size_t jGlobal = jF * m_nlocal_edge_dofs + jk;
                  cell_triplets_GlobMat[iT].emplace_back(iGlobal, jGlobal, MatF(iLocal, jLocal));
                }
              }
            }
          }
          size_triplets_GlobMat += cell_triplets_GlobMat[iT].size();
      
        }
    };    // End function to construct local contributions

    // Running the local constructions in parallel
    parallel_for(mesh->n_cells(), construct_all_local_contributions, m_use_threads);

    // Assemble local contribution into global matrix
    triplets_ScBe.reserve(size_triplets_ScBe);
    triplets_GlobMat.reserve(size_triplets_GlobMat);
    for (size_t iT = 0; iT < mesh->n_cells(); iT++){
      for (size_t i = 0; i < cell_triplets_ScBe[iT].size(); i++){
        triplets_ScBe.push_back(cell_triplets_ScBe[iT][i]);
      }
      for (size_t i = 0; i < cell_triplets_GlobMat[iT].size(); i++){
        triplets_GlobMat.push_back(cell_triplets_GlobMat[iT][i]);
      }
      Cell& T = *mesh->cell(iT);      
      for (size_t ilF = 0; ilF < T.n_edges(); ilF++) {
        const size_t iF = T.edge(ilF)->global_index();
        for (size_t ik = 0; ik < m_nlocal_edge_dofs; ik++) {
          const size_t iLocal = ilF * m_nlocal_edge_dofs + ik;
          const size_t iGlobal = iF * m_nlocal_edge_dofs + ik;
          GlobRHS(iGlobal) += cell_source[iT](iLocal);
        }
      }
    }

    if (m_BC.name()=="Neumann"){
      // Neumann BC: remove a row in the matrix and fix the first degree of freedom
      triplets_GlobMat.erase(std::remove_if(std::begin(triplets_GlobMat), std::end(triplets_GlobMat),
              [](const auto& x) { return (x.row() == 0); }), std::end(triplets_GlobMat));
      triplets_GlobMat.emplace_back(0, 0, 1);
      GlobRHS(0) = 0;
    }

    // Assemble the global linear system (without BC), and matrix to recover statically-condensed cell dofs
    GlobMat.setFromTriplets(std::begin(triplets_GlobMat), std::end(triplets_GlobMat));
    ScBeMat.setFromTriplets(std::begin(triplets_ScBe), std::end(triplets_ScBe));

    // Record assembly time 
  //    _assembly_time = timer.elapsed().user + timer.elapsed().system;
    _assembly_time = timer.elapsed().wall;

  }

  UVector HHO_LocVarDiff::solve(HybridCore& hho) {
    const Mesh* mesh = hho.get_mesh();
    boost::timer::cpu_timer timer;  // Time the matrix assembly
    timer.start();

    //-------------- TREATMENT OF BOUNDARY CONDITIONS -------------//

    // If Dirichlet, the final system is only posed on the interior edge unknowns and we have to subtract from the source
    //    term the contribution of the boundary values
    // If Neumann, the final system is posed on all edge unknowns

    size_t n_unknowns = 0;
    size_t n_fixed_dofs = 0;
    Eigen::VectorXd B;
    Eigen::VectorXd UDir;

    if (m_BC.name() != "Neumann"){
      // Dirichlet boundary conditions
      n_unknowns = m_nnondir_edge_dofs;
      n_fixed_dofs = m_ndir_edge_dofs;
      SysMat = GlobMat.topLeftCorner(n_unknowns, n_unknowns);
    
      // Boundary value: UDir corresponds to the L2 projection of the exact solution on the polynomial spaces on the Dirichlet edges (last BC.n_dir_edges() edges)
      UDir = Eigen::VectorXd::Zero(n_fixed_dofs);
      size_t n_dir_edges = m_BC.n_dir_edges();
      size_t n_nondir_edges = mesh->n_edges() - n_dir_edges; 
      for (size_t idF = 0; idF < n_dir_edges; idF++){
        Edge* edge = mesh->edge(n_nondir_edges + idF);
        size_t iF = edge->global_index();
        QuadratureRule quadF = generate_quadrature_rule(*edge, 2*m_K+2);
        boost::multi_array<double, 2> phiF_quadF = evaluate_quad<Function>::compute(hho.EdgeBasis(iF), quadF);
        UDir.segment(idF * m_nlocal_edge_dofs, m_nlocal_edge_dofs) = l2_projection<HybridCore::PolyEdgeBasisType>(exact_solution, hho.EdgeBasis(iF), quadF, phiF_quadF);
      }

      B = GlobRHS.segment(0, n_unknowns) - GlobMat.topRightCorner(n_unknowns, n_fixed_dofs) * UDir;

    } else {
      // We will solve the complete system
      n_unknowns = m_ntotal_edge_dofs;
      SysMat = GlobMat;
      B = GlobRHS;
      UDir = Eigen::VectorXd::Zero(n_fixed_dofs);
    }

    //-------------- SOLVE CONDENSED SYSTEM -------------//

    Eigen::VectorXd xF = Eigen::VectorXd::Zero(n_unknowns);

    //  if (solver_type == "ma41") {
    //    Eigen::MA41<Eigen::SparseMatrix<double>, Eigen::VectorXd> solver;
    //    solver.analyzePattern(SysMat);
    //    solver.factorize(SysMat);
    //    xF = solver.solve(B);
    //  } else {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(SysMat);
    xF = solver.solve(B);
    std::cout << "  [solver] #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << std::endl;
    //  }
    _solving_error = (SysMat * xF - B).norm();
    // Recover the fixed boundary values, cell unknowns (from static condensation/barycentric elimination)
    Eigen::VectorXd Xh = Eigen::VectorXd::Zero(m_ntotal_dofs);
    Xh.tail(n_fixed_dofs) = UDir;
    Xh.segment(m_ntotal_cell_dofs, n_unknowns) = xF;
    if (m_L>=0) {
      Xh.head(m_ntotal_cell_dofs) = ScRHS - ScBeMat * Xh.tail(m_ntotal_edge_dofs);
    } else {
      Xh.head(m_ntotal_cell_dofs) = ScBeMat * Xh.tail(m_ntotal_edge_dofs);
    }

    // Only Neumann: translate to get the proper average
    if (m_BC.name()=="Neumann"){
      // Compute average to translate
      double average = 0.;
      double total_measure = 0;
      for (size_t iT = 0; iT < mesh->n_cells(); iT++){
        Cell& T = *mesh->cell(iT);
        total_measure += T.measure();
        HybridCore::PolyCellBasisType basisT = hho.CellBasis(iT);
        QuadratureRule quadT = generate_quadrature_rule(T, 2*m_K);
        boost::multi_array<double, 2> phiT_quadT = evaluate_quad<Function>::compute(basisT, quadT);
        for (size_t i = 0; i < basisT.dimension(); i++){
          for (size_t iqn = 0; iqn < quadT.size(); iqn++){
            average += quadT[iqn].w * Xh(iT * m_nlocal_cell_dofs + i) * basisT.function(i, quadT[iqn].vector());
          }
        }
      }
      double average_exact_sol = 0.0;
      for (auto& T : mesh->get_cells()){
        QuadratureRule quadT = generate_quadrature_rule(*T, 2 * m_Ldeg + 2);
        for (QuadratureNode& qT : quadT){
          average_exact_sol += qT.w * exact_solution(qT.vector());
        }
      }
      average_exact_sol /= total_measure;

      // Translate the cells and edges
      // We compute the interpolant of the constant function "average_exact_sol - average" 
      // and we translate Xh by that amount
      std::function<double(VectorRd)> AveDiff = [&average_exact_sol,&average](VectorRd x)->double
  { return average_exact_sol - average;};

      UVector Cst = hho.interpolate(AveDiff, m_L, m_K, 2*m_K+3);
      Xh += Cst.asVectorXd();
    }

    _solving_time = timer.elapsed().user + timer.elapsed().system;  // Record the final solving time

    return UVector(Xh, *hho.get_mesh(), m_L, m_K);
  }

  //******************************** 
  //    local diffusion matrix 
  //********************************

  Eigen::MatrixXd HHO_LocVarDiff::diffusion_operator(HybridCore &hho, const size_t iT, const ElementQuad &elquad) const 
  {

    boost::timer::cpu_timer timeint;

    const auto mesh = hho.get_mesh();
    const size_t dimPKcell = DimPoly<Cell>(m_K);
    const size_t dimPKcell_vec = mesh->dim() * dimPKcell;
    Cell* cell = mesh->cell(iT);
    const size_t nedgesT = cell->n_edges();

    // Total number of degrees of freedom local to this cell (cell and its adjacent edges)
    size_t local_dofs = m_nlocal_cell_dofs + nedgesT * m_nlocal_edge_dofs;

    //-------------------  Initialisatons: quadratures, mass matrices... --------------------//

    // Note: representing the gradient reconstruction GT supposes a basis for (P^k)^d. This basis can be
    //   built from the basis (phi_i)_i of P^k, by working component-by-component:
    //    (phi_1 e_1, ..., phi_N e_1, phi_1 e_2, ..., phi_N e_2, ...),
    //    where e_1=[1 0], e_2=[0 1]... is the canonical basis of R^d
    //  Hence, the first dimPKcell functions are on the first component of R^d, the next dimPKcell on the second
    //    component etc.
    //  We do not explicitly build this basis function, but we compute the values of these functions at the
    //    quadrature nodes, using the values computed for the scalar basis functions


    // QUADRATURES
    // Cell quadrature nodes, and values of cell basis functions (up to degree K+1) and gradients thereof.
    QuadratureRule quadT = elquad.get_quadT();
    boost::multi_array<double, 2> phiT_quadT = elquad.get_phiT_quadT();
    boost::multi_array<VectorRd, 2> dphiT_quadT = elquad.get_dphiT_quadT();

    // Vector basis functions (up to degree K) at the quadrature nodes. 
    boost::multi_array<VectorRd, 2> vec_phiT_quadT = elquad.get_vec_phiT_quadT(m_K);
    
    // Diffusion tensor at the quadrature nodes
    std::vector<MatrixRd> kappaT_quadT(quadT.size());
    std::transform(quadT.begin(), quadT.end(), kappaT_quadT.begin(),
       [this,&cell](QuadratureNode qr) -> Eigen::MatrixXd { return kappa(qr.vector(), cell); });

    // MASS MATRICES
    // Scalar cell mass matrix (phi_i,phi_j)_T up to degree max(L,K) * (K+1), and 
    // Vector cell mass matrix (Phi_i,Phi_j)_T up to degree K*K [this one is block diagonal]
    size_t maxdimPKL = std::max(m_nlocal_cell_dofs, dimPKcell);

    _itime[0] += timeint.elapsed().user + timeint.elapsed().system;
    timeint.start();

    Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, maxdimPKL, m_nhighorder_dofs, "sym");
    Eigen::MatrixXd VecMTT = Eigen::MatrixXd::Zero(dimPKcell_vec, dimPKcell_vec);
    for (size_t r=0; r < mesh->dim(); r++){
      VecMTT.block(r*dimPKcell, r*dimPKcell, dimPKcell, dimPKcell) = MTT.topLeftCorner(dimPKcell, dimPKcell);
    } 

    // Edge mass matrices:
    // MFF[ilF]: edge-edge mass on edge with local number ilF, up to degree K*K
    // MFT[ilF]: edge-cell mass on edge with local number ilF, up to degree K*(K+1)
    // MTT_on_F[ilF]: cell-cell mass on edge with local number ilF, up to degree K*L
    std::vector<Eigen::MatrixXd> MFF(nedgesT, Eigen::MatrixXd::Zero(m_nlocal_edge_dofs, m_nlocal_edge_dofs));
    std::vector<Eigen::MatrixXd> MFT(nedgesT, Eigen::MatrixXd::Zero(m_nlocal_edge_dofs, m_nhighorder_dofs));
    std::vector<Eigen::MatrixXd> MTT_on_F(nedgesT, Eigen::MatrixXd::Zero(dimPKcell, m_nlocal_cell_dofs));

    for (size_t ilF = 0; ilF < nedgesT; ilF++) {

      // Edge quadrature nodes and values of cell and edge basis functions (and gradients) at these nodes
      auto quadF = elquad.get_quadF(ilF);
      boost::multi_array<double, 2> phiT_quadF = elquad.get_phiT_quadF(ilF);
      boost::multi_array<double, 2> phiF_quadF = elquad.get_phiF_quadF(ilF);

      // Mass matrices
      MFF[ilF] = compute_gram_matrix(phiF_quadF, phiF_quadF, quadF, "sym");
      MFT[ilF] = compute_gram_matrix(phiF_quadF, phiT_quadF, quadF, "nonsym");
      MTT_on_F[ilF] = compute_gram_matrix(phiT_quadF, phiT_quadF, quadF, dimPKcell, m_nlocal_cell_dofs, "nonsym");
    }

    _itime[1] += timeint.elapsed().user + timeint.elapsed().system;
    timeint.start();

    // STIFNESS mass-matrices (and the like): (\nabla phi_i,\nabla phi_j)_T up to degree (K+1)*(K+1)
    //  and (\nabla phi_i, Phi_j)_T up to degree (K+1)*K
    Eigen::MatrixXd StiffT = compute_gram_matrix(dphiT_quadT, dphiT_quadT, quadT, "sym");  
    Eigen::MatrixXd MdphiT_PhiT = compute_gram_matrix(dphiT_quadT, vec_phiT_quadT, quadT, "nonsym");

    _itime[2] += timeint.elapsed().user + timeint.elapsed().system;
    timeint.start();

    //-------------------- Compute GT, matrix of full gradient reconstruction ---------//

    // Right-hand side, starting with volumetric term (Phi_i, \nabla phi_j)_T
    Eigen::MatrixXd RHS_GT = Eigen::MatrixXd::Zero(dimPKcell_vec, local_dofs);
    RHS_GT.topLeftCorner(dimPKcell_vec, m_nlocal_cell_dofs) = (MdphiT_PhiT.topLeftCorner(m_nlocal_cell_dofs, dimPKcell_vec)).transpose();

    // Boundary terms
    for (size_t r=0; r < mesh->dim(); r++){
      for (size_t ilF = 0; ilF < nedgesT; ilF++) {
        // Offset for edge unknowns
        const size_t offset_F = m_nlocal_cell_dofs + ilF * m_nlocal_edge_dofs;
        const auto& nTF = cell->edge_normal(ilF);

        // Contribution of cell unknowns, and then edge unknowns on F
        RHS_GT.block(r*dimPKcell, 0, dimPKcell, m_nlocal_cell_dofs) -= nTF(r) * MTT_on_F[ilF];
        RHS_GT.block(r*dimPKcell, offset_F, dimPKcell, m_nlocal_edge_dofs) += 
          nTF(r) * (MFT[ilF].topLeftCorner(m_nlocal_edge_dofs, dimPKcell)).transpose();
      }
    }

    // Compute GT
    Eigen::MatrixXd GT = (VecMTT.ldlt()).solve(RHS_GT);

    //------------- Consistent contribution (K GT, GT)_T to the local bilinear form ------//


    //  Weighted mass matrix  (K Phi_i, Phi_j)_T of the basis (Phi_i)_i of (P^k)^d
    // We start by creating K Phi_i at the quadrature nodes
    boost::multi_array<VectorRd, 2> kappa_vec_phiT_quadT;
    kappa_vec_phiT_quadT.resize( boost::extents[dimPKcell_vec][quadT.size()] );
    for (size_t i = 0; i < dimPKcell_vec; i++){
      for (size_t iqn = 0; iqn < quadT.size(); iqn++){
        kappa_vec_phiT_quadT[i][iqn] = kappaT_quadT[iqn] * vec_phiT_quadT[i][iqn];
      }
    }
    _itime[3] += timeint.elapsed().user + timeint.elapsed().system;
    timeint.start();

    Eigen::MatrixXd kappaVecMTT = compute_gram_matrix(kappa_vec_phiT_quadT, vec_phiT_quadT, quadT, vec_phiT_quadT.shape()[0], vec_phiT_quadT.shape()[0], "sym");

    _itime[4] += timeint.elapsed().user + timeint.elapsed().system;
    timeint.start();

    Eigen::MatrixXd ATF = GT.transpose() * kappaVecMTT * GT;


    //------------- Compute PT, matrix of potential reconstruction, using GT ------//

    // We write that \nabla pT = projection on \nabla P^{k+1} of GT, and add the closure relation:
    // (nabla pT v, nabla w)_T + lambda_T(p_T v,1)_T(w,1)_T = (GT v,\nabla w)_T + lambda_T(v_T,1)_T(w,1)_T

    // Right-hand side, starting with volumetric term
    Eigen::MatrixXd RHS_PT = MdphiT_PhiT * GT;

    // Vector LT of (phi_j,1)_T for phi_j up to degree K+1, and LT^t*LT, for the closure relation
    Eigen::VectorXd LT = (MTT.row(0)).transpose();
    Eigen::MatrixXd LTtLT = LT * (LT.transpose());
    double scalT = StiffT.trace() / std::pow(LT.norm(), 2);

    // Add closure relation and compute PT
    RHS_PT.topLeftCorner(m_nhighorder_dofs, m_nlocal_cell_dofs) += 
      scalT * LTtLT.topLeftCorner(m_nhighorder_dofs, m_nlocal_cell_dofs);
    Eigen::MatrixXd PT = ((StiffT + scalT*LTtLT).ldlt()).solve(RHS_PT);


    //-------------------- Compute stabilisation term sT ---------//

    Eigen::MatrixXd STF = Eigen::MatrixXd::Zero(local_dofs, local_dofs);

    // Cell residual delta_T^l = pi_T^l (rT uT) - u_T

    Eigen::MatrixXd MTT_LKp1 = MTT.topLeftCorner(m_nlocal_cell_dofs, m_nhighorder_dofs);
    Eigen::MatrixXd MTT_LL = MTT.topLeftCorner(m_nlocal_cell_dofs, m_nlocal_cell_dofs);
    Eigen::MatrixXd deltaTL = MTT_LL.ldlt().solve( MTT_LKp1 * PT );
    deltaTL.topLeftCorner(m_nlocal_cell_dofs, m_nlocal_cell_dofs) -= Eigen::MatrixXd::Identity(m_nlocal_cell_dofs, m_nlocal_cell_dofs);

    for (size_t ilF = 0; ilF < nedgesT; ilF++) {
      // Two options for stabilisation: diameter of edge, or ratio measure cell/measure edge
     double dTF = cell->edge(ilF)->diam();
    //   double dTF = cell->measure() / cell->edge(ilF)->measure();

      VectorRd xF = cell->edge(ilF)->center_mass();

    //   auto kappa_TF = kappa(xF, cell).trace();

      const VectorRd &nTF = cell->edge_normal(ilF);
      const double kappa_TF = (kappa(xF, cell) * nTF).dot(nTF);

      // Edge residual delta_TF^k = pi_F^k (rT uT) - u_F
      Eigen::MatrixXd MFFinv = MFF[ilF].inverse();
      Eigen::MatrixXd deltaTFK = MFFinv * MFT[ilF] * PT;
      deltaTFK.block(0, m_nlocal_cell_dofs + ilF * m_nlocal_edge_dofs, m_nlocal_edge_dofs, m_nlocal_edge_dofs) -=
        Eigen::MatrixXd::Identity(m_nlocal_edge_dofs, m_nlocal_edge_dofs);

      // Stabilisation term: here, we actually project deltaTL on P^k(F) so, for l=k+1, it actually corresponds to the stabilisation used in HDG methods (see Section 5.1.6 of HHO book)
      Eigen::MatrixXd deltaTFK_minus_deltaTL = deltaTFK - MFFinv * MFT[ilF].topLeftCorner(m_nlocal_edge_dofs, m_nlocal_cell_dofs) * deltaTL;

      STF += (kappa_TF / dTF) * deltaTFK_minus_deltaTL.transpose() * MFF[ilF] *  deltaTFK_minus_deltaTL;
    }

    _itime[5] += timeint.elapsed().user + timeint.elapsed().system;

    // Adjust local bilinear form with stabilisation term
    ATF += mesh->dim() * STF;

    return ATF;

  }


  //******************************** 
  //    local load term 
  //********************************

  Eigen::VectorXd HHO_LocVarDiff::load_operator(HybridCore &hho, const size_t iT, const ElementQuad &elquad) const {
    // Load for the cell DOFs (first indices) and edge DOFs (last indices)
    const auto mesh = hho.get_mesh();
    Cell* cell = mesh->cell(iT);
    size_t cell_edge_dofs = m_nlocal_cell_dofs + cell->n_edges()*m_nlocal_edge_dofs;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(cell_edge_dofs);

    // Quadrature nodes and values of cell basis functions at these nodes
    auto quadT = elquad.get_quadT();
    size_t nbq = quadT.size();
    boost::multi_array<double, 2> phiT_quadT = elquad.get_phiT_quadT();

    // Value of source times quadrature weights at the quadrature nodes
    Eigen::ArrayXd weight_source_quad = Eigen::ArrayXd::Zero(nbq);
    for (size_t iqn = 0; iqn < nbq; iqn++){
      weight_source_quad(iqn) = quadT[iqn].w * source(quadT[iqn].vector(), cell);
    }

    for (size_t i=0; i < m_nlocal_cell_dofs; i++){
      for (size_t iqn = 0; iqn < quadT.size(); iqn++){
        b(i) += weight_source_quad[iqn] * phiT_quadT[i][iqn];
      }
    }

    // Boundary values, if we have a boundary cell
    if (cell->is_boundary()){
      // Boundary values only on boundary Neumann edges
      for (size_t ilF = 0; ilF < cell->n_edges(); ilF++) {
        Edge* F = cell->edge(ilF);
        if (m_BC.type(*F)=="neu"){
          const size_t iF = F->global_index(); 
          // BC on boundary edges
          if (F->is_boundary()){
            // Offset for edge unknowns
            const size_t offset_F = m_nlocal_cell_dofs + ilF * m_nlocal_edge_dofs;
            // Normal to the edge and bases function
            const auto& nTF = cell->edge_normal(ilF);
            const auto& basisF = hho.EdgeBasis(iF);
            // for each DOF of the boundary edge
            for (size_t i = 0; i < m_nlocal_edge_dofs; i++){
              QuadratureRule quadF = generate_quadrature_rule(*F, 2*m_K+2);
              std::function<double(VectorRd)> Kgrad_n = [&](VectorRd p){
                return nTF.dot(kappa(p,cell) * grad_exact_solution(p,cell)) * basisF.function(i,p);
              };
              for (QuadratureNode& qF : quadF){
                b(offset_F + i) += qF.w * Kgrad_n(qF.vector());
              }
            }
          }
        }
      }
    }

    return b;
  }

  double HHO_LocVarDiff::EnergyNorm(HybridCore& hho, const UVector Xh) {
    const auto mesh = hho.get_mesh();
    double value = 0.0;

    for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
      Eigen::VectorXd XTF = Xh.restr(iT);
      value += XTF.transpose() * aT[iT] * XTF;
    }

    return sqrt(value);
  }


  //@}
} // end of namespace HArDCore2D

#endif //_HHO_LOCVARDIFF_HPP
