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
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
* If you use this code or part of it for a scientific publication, please cite the book
* above as a reference for the implementation.
*
*/


#ifndef _HHO_DIFFUSION_HPP
#define _HHO_DIFFUSION_HPP

#include <functional>
#include <utility>
#include <iostream>

#include <boost/timer/timer.hpp>

// Matrices and linear solvers
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include "Eigen/MA41.cpp"


#include <mesh.hpp>
#include <hybridcore.hpp>
#include <elementquad.hpp>

/*!
* @defgroup HHO_Diffusion
* @brief HHO scheme for diffusion equation -div(Diff grad u)=f, with Diff piecewise constant on the mesh
*/

namespace HArDCore2D {

/*!
*  @addtogroup HHO_Diffusion
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The HHO_Diffusion class provides tools to implement the HHO method for the diffusion problem

class HHO_Diffusion {

// Types
public:
  using solution_function_type = std::function<double(double,double)>;    ///< type for solution
  using source_function_type = std::function<double(double,double,Cell*)>;    ///< type for source
  using grad_function_type = std::function<Eigen::Vector2d(double,double,Cell*)>;    ///< type for gradient
  using tensor_function_type = std::function<Eigen::Matrix2d(double,double,Cell*)>;    ///< type for diffusion tensor

  ///@brief Constructor of the class
  HHO_Diffusion(
    tensor_function_type kappa,   ///< diffusion tensor
    source_function_type source,  ///< source term
    size_t BC,                     ///< type of boundary conditions (0 for Dirichlet, 1 for Neumann)
    solution_function_type exact_solution,   ///< exact solution
    grad_function_type grad_exact_solution,   ///< gradient of the exact solution
    std::string solver_type    ///< type of solver to use for the global system (bicgstab at the moment)
    );

  /// Assemble and solve the scheme
  Eigen::VectorXd solve(HybridCore& hho);

  /// Discrete energy norm (associated to the diffusion operator) of an hybrid function
  double EnergyNorm(HybridCore& hho, const Eigen::VectorXd Xh); 

  double get_assembly_time() const; ///< cpu time to assemble the scheme
  double get_solving_time() const;  ///< cpu time to solve the scheme
  double get_solving_error() const;  ///< residual after solving the scheme
  double get_itime(size_t idx) const;    ///< various intermediate assembly times

private:
  /// Compute the local diffusion operator in the cell iT
  Eigen::MatrixXd diffusion_operator(HybridCore& hho, const size_t iT, const ElementQuad &elquad) const;

  /// Compute the local load operator (the source term) in the cell iT
  Eigen::VectorXd load_operator(HybridCore& hho, const size_t iT, const ElementQuad &elquad) const;

  const tensor_function_type kappa;
  const source_function_type source;
  const size_t BC;
  const solution_function_type exact_solution;
  const grad_function_type grad_exact_solution;
  const std::string solver_type;

  // To store local bilinear forms
  std::vector<Eigen::MatrixXd> aT;

  // Computation statistics
  size_t _assembly_time;
  size_t _solving_time;
  double _solving_error;
  mutable std::vector<size_t> _itime = std::vector<size_t>(10, 0);

};

HHO_Diffusion::HHO_Diffusion(tensor_function_type kappa, source_function_type source, size_t BC, solution_function_type exact_solution, grad_function_type grad_exact_solution, std::string solver_type)
  : kappa(kappa),
    source(source),
    BC(BC),
    exact_solution(exact_solution),
    grad_exact_solution(grad_exact_solution),
    solver_type(solver_type) {
  // Do nothing
}

Eigen::VectorXd HHO_Diffusion::solve(HybridCore &hho) {

  boost::timer::cpu_timer timer;  // Time the matrix assembly
  boost::timer::cpu_timer timerint; 
  const auto mesh = hho.get_mesh_ptr();
  aT.resize(mesh->n_cells());

  //--------------- PREPARE SYSTEM ------------------------//

  // System matrix
  Eigen::SparseMatrix<double> GlobMat(hho.ntotal_edge_dofs(), hho.ntotal_edge_dofs());
  std::vector<Eigen::Triplet<double>> triplets_GlobMat;
  // If static condensation (L>=0): matrix to recover cell unknowns
  // If barycentric elimination (L=-1): matrix to recover cell unknowns
  Eigen::SparseMatrix<double> ScBeMat(hho.ntotal_cell_dofs(), hho.ntotal_edge_dofs());
  std::vector<Eigen::Triplet<double>> triplets_ScBe;

  // Source terms for the system, and for recovering cell unknowns from static condensation
  Eigen::VectorXd GlobRHS = Eigen::VectorXd::Zero(hho.ntotal_edge_dofs());
  Eigen::VectorXd ScRHS = Eigen::VectorXd::Zero(hho.ntotal_cell_dofs());

  // Global quadrature rule for the cells
  Eigen::VectorXd cell_quadrature = Eigen::VectorXd::Zero(hho.ntotal_cell_dofs());


  //-------------- ASSEMBLE LOCAL CONTRIBUTIONS -------------//
  auto total_measure = 0.0;

  for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
    Cell* iCell = mesh->cell(iT);
  
    total_measure += iCell->measure();

    // Total number of face degrees of freedom local to this cell (adjacent faces to the cell)
    size_t nlocal_edges = iCell->n_edges();
    size_t edge_dofs = nlocal_edges * hho.nlocal_edge_dofs();

    // Local bilinear form and source term
    size_t doeT = hho.Ldeg()+hho.K()+1;
    size_t doeF = 2*hho.K()+1;
    ElementQuad elquad(hho, iT, doeT, doeF);
    aT[iT] = diffusion_operator(hho, iT, elquad);
    Eigen::VectorXd bT = load_operator(hho, iT, elquad);

    // Local matrix and right-hand side on the face unknowns
    Eigen::MatrixXd MatF = Eigen::MatrixXd::Zero(edge_dofs,edge_dofs);
    Eigen::VectorXd bF = Eigen::VectorXd::Zero(edge_dofs);

    if (hho.L()>=0) {
      // STATIC CONDENSATION OF ELEMENT UNKNOWNS

      // Perform static condensation
      Eigen::MatrixXd ATT = aT[iT].topLeftCorner(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs());
      Eigen::MatrixXd ATF = aT[iT].topRightCorner(hho.nlocal_cell_dofs(), edge_dofs);
      Eigen::MatrixXd AFF = aT[iT].bottomRightCorner(edge_dofs, edge_dofs);

      Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
      invATT.compute(ATT);
      
      Eigen::MatrixXd invATT_ATF = invATT.solve(ATF);
      Eigen::VectorXd invATT_bTcell = invATT.solve(bT.head(hho.nlocal_cell_dofs()));
      MatF = AFF - ATF.transpose() * invATT_ATF;
      
      bF = bT.tail(edge_dofs) - ATF.transpose() * invATT_bTcell;
            
      // Assemble static condensation operator
      ScRHS.segment(iT * hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs()) = invATT_bTcell;
      for (size_t i = 0; i < hho.nlocal_cell_dofs(); i++) {
        for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
          const size_t jF = iCell->edge(jlF)->global_index();
          for (size_t jk = 0; jk < hho.nlocal_edge_dofs(); jk++) {
            const size_t jLocal = jlF * hho.nlocal_edge_dofs() + jk;
            const size_t jGlobal = jF * hho.nlocal_edge_dofs() + jk;
            triplets_ScBe.emplace_back(iT * hho.nlocal_cell_dofs() + i, jGlobal, invATT_ATF(i, jLocal));
          }
        }
      }
    } else {
      // BARYCENTRIC ELIMINATION OF ELEMENT UNKNOWNS
      // Create reduction matrix: 1+nlocal_edges * nlocal_edges matrix with the coefficients on the first row, and the identity below. When multiplied by the face unknowns, return cell and face unknowns
      // Note that the basis functions are constant, but not necessarily assumed to be one (which is not the case after orthonormalisation for example), which is why we have to adjust the first row.
      Eigen::MatrixXd red_matT = Eigen::MatrixXd::Zero(1+nlocal_edges,nlocal_edges);
      red_matT.row(0) = hho.compute_weights(iT);
      Eigen::Vector2d xT = iCell->center_mass();
      double phiT_cst = hho.cell_basis(iT,0)(xT.x(), xT.y());
      for (size_t ilE = 0; ilE < nlocal_edges; ilE++){
        Eigen::Vector2d xE = iCell->edge(ilE)->center_mass();
        size_t iE = iCell->edge(ilE)->global_index();
        double phiE_cst = hho.edge_basis(iE,0)(xE.x(), xE.y());
        red_matT(0,ilE) *= phiE_cst / phiT_cst;
      }
      red_matT.bottomRightCorner(nlocal_edges,nlocal_edges) = Eigen::MatrixXd::Identity(nlocal_edges,nlocal_edges);

      bF = red_matT.transpose() * bT;
      MatF = red_matT.transpose() * aT[iT] * red_matT;

      // To recover cell unknown
      for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
        const size_t jF = iCell->edge(jlF)->global_index();
        const size_t jGlobal = jF * hho.nlocal_edge_dofs();
        triplets_ScBe.emplace_back(iT, jGlobal, red_matT(0,jlF));
      }

    }

    // Assemble local contribution into global matrix
    for (size_t ilF = 0; ilF < nlocal_edges; ilF++) {
      const size_t iF = iCell->edge(ilF)->global_index();
      for (size_t ik = 0; ik < hho.nlocal_edge_dofs(); ik++) {
        const size_t iLocal = ilF * hho.nlocal_edge_dofs() + ik;
        const size_t iGlobal = iF * hho.nlocal_edge_dofs() + ik;
        for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
          const size_t jF = iCell->edge(jlF)->global_index();
          for (size_t jk = 0; jk < hho.nlocal_edge_dofs(); jk++) {
            const size_t jLocal = jlF * hho.nlocal_edge_dofs() + jk;
            const size_t jGlobal = jF * hho.nlocal_edge_dofs() + jk;
            triplets_GlobMat.emplace_back(iGlobal, jGlobal, MatF(iLocal, jLocal));
          }
        }
        GlobRHS(iGlobal) += bF(iLocal);
      }
    }
  
    // Record cell quadrature
    if (BC==1){
      for (size_t i = 0; i < hho.nlocal_cell_dofs(); i++) {
        const auto& phi_i = hho.cell_basis(iT, i);
        cell_quadrature(iT * hho.nlocal_cell_dofs() + i) = hho.integrate_over_cell(iT, phi_i);
      }
    }

  }

  if (BC==1){
    // Neumann BC: remove a row in the matrix and fix the first degree of freedom
    triplets_GlobMat.erase(std::remove_if(std::begin(triplets_GlobMat), std::end(triplets_GlobMat),
                                  [](const auto& x) { return (x.row() == 0); }), std::end(triplets_GlobMat));
    triplets_GlobMat.emplace_back(0, 0, 1);
    GlobRHS(0) = 0;
  }

  // Assemble the global linear system (without BC), and matrix to recover statically-condensed cell dofs
  GlobMat.setFromTriplets(std::begin(triplets_GlobMat), std::end(triplets_GlobMat));
  ScBeMat.setFromTriplets(std::begin(triplets_ScBe), std::end(triplets_ScBe));

  //-------------- TREATMENT OF BOUNDARY CONDITIONS -------------//

  // If Dirichlet, the final system is only posed on the interior edge unknowns and we have to subtract from the source
  //    term the contribution of the boundary values
  // If Neumann, the final system is posed on all edge unknowns

  size_t nb_unknowns = 0;
  size_t nb_fixed_dofs = 0;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd B;
  Eigen::VectorXd UDir;

  if (BC==0){
    // Dirichlet boundary conditions
    nb_unknowns = hho.ninternal_edge_dofs();
    nb_fixed_dofs = hho.nboundary_edge_dofs();
    A = GlobMat.topLeftCorner(nb_unknowns, nb_unknowns);
    
    // Boundary value: UDir corresponds to the L2 projection of the exact solution on the polynomial spaces on the edges
    UDir = Eigen::VectorXd::Zero(nb_fixed_dofs);
    std::vector<Edge *> b_edges = mesh->get_b_edges(); // List of boundary edges

    for (size_t ibF = 0; ibF < mesh->n_b_edges(); ibF++){
      Edge* edge = b_edges[ibF]; 
      size_t iF = edge->global_index();
      // Mass matrix and boundary values
      QuadratureRule quadF = generate_quadrature_rule(*edge, 2*hho.K());
      std::vector<Eigen::ArrayXd> phiF_quadF = hho.basis_quad("edge", iF, quadF, hho.K());
      Eigen::MatrixXd MFF = hho.gram_matrix(phiF_quadF, phiF_quadF, hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs(), quadF, true);

      // Compute (uexact, phi_i)_F for all edge basis functions phi_i
      Eigen::VectorXd RHS_UDirF = Eigen::VectorXd::Zero(hho.nlocal_edge_dofs());
      for (size_t j = 0; j < hho.nlocal_edge_dofs(); j++){
        const auto& phi_j = hho.edge_basis(iF, j);
        QuadratureRule quadF = generate_quadrature_rule(*edge, 2*hho.K()+2);
        for (QuadratureNode qF : quadF){
          RHS_UDirF(j) += qF.w * phi_j(qF.x, qF.y) * exact_solution(qF.x, qF.y);
        }

      }
      // Project exact solution
      UDir.segment(ibF * hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs()) = MFF.ldlt().solve(RHS_UDirF);

    }

    B = GlobRHS.segment(0, nb_unknowns) - GlobMat.topRightCorner(nb_unknowns, nb_fixed_dofs) * UDir;

  } else if (BC==1){
    // We will solve the complete system
    nb_unknowns = hho.ntotal_edge_dofs();
    A = GlobMat;
    B = GlobRHS;
    UDir = Eigen::VectorXd::Zero(nb_fixed_dofs);
  }

  _assembly_time = timer.elapsed().user + timer.elapsed().system;

  //-------------- SOLVE CONDENSED SYSTEM -------------//
  timer.start();

  Eigen::VectorXd xF = Eigen::VectorXd::Zero(nb_unknowns);

//  if (solver_type == "ma41") {
//    Eigen::MA41<Eigen::SparseMatrix<double>, Eigen::VectorXd> solver;
//    solver.analyzePattern(A);
//    solver.factorize(A);
//    xF = solver.solve(B);
//  } else {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    xF = solver.solve(B);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
//  }
  _solving_error = (A * xF - B).norm();
  // Recover the fixed boundary values, cell unknowns (from static condensation/barycentric elimination)
  Eigen::VectorXd Xh = Eigen::VectorXd::Zero(hho.ntotal_dofs());
  Xh.tail(nb_fixed_dofs) = UDir;
  Xh.segment(hho.ntotal_cell_dofs(), nb_unknowns) = xF;
  if (hho.L()>=0) {
    Xh.head(hho.ntotal_cell_dofs()) = ScRHS - ScBeMat * Xh.tail(hho.ntotal_edge_dofs());
  } else {
    Xh.head(hho.ntotal_cell_dofs()) = ScBeMat * Xh.tail(hho.ntotal_edge_dofs());
  }

  // Only Neumann: translate to get the proper average
  if (BC==1){
    // Compute average to translate
    auto average = cell_quadrature.dot(Xh.head(hho.ntotal_cell_dofs())) / total_measure;
    auto average_exact_sol = hho.integrate_over_domain([&](auto x, auto y) {
      return exact_solution(x,y);
      }) / total_measure;

    // Translate the cells and edges
    // We compute the interpolant of the constant function "average_exact_sol - average" 
    // and we translate Xh by that amount
    std::function<double(double, double)> AveDiff = [&average_exact_sol,&average](double x, double y)->double
        { return average_exact_sol - average;};
    Eigen::VectorXd Cst = hho.interpolate(AveDiff, 2*hho.K()+3);
    Xh += Cst;

  }

  _solving_time = timer.elapsed().user + timer.elapsed().system;  // Record the final solving time

  return Xh;
}

//******************************** 
//    local diffusion matrix 
//********************************

Eigen::MatrixXd HHO_Diffusion::diffusion_operator(HybridCore &hho, const size_t iT, const ElementQuad &elquad) const {

  boost::timer::cpu_timer timeint;

  const auto mesh = hho.get_mesh_ptr();
  Cell* cell = mesh->cell(iT);
  const size_t nedgesT = cell->n_edges();

  // Total number of degrees of freedom local to this cell (cell and its adjacent faces)
  size_t local_dofs = hho.nlocal_cell_dofs() + nedgesT * hho.nlocal_edge_dofs();

  //-------------------  Initialisatons ------------------------------//

timeint.start();
  
  // Diffusion in the cell.
  std::function<Eigen::Matrix2d(double,double)> kappaT = [&](double x, double y){
    // Constant in the cell
    Eigen::Vector2d cellcenter = cell->center_mass();
    return kappa(cellcenter.x(), cellcenter.y(), cell);
    // Variable in the cell - the scheme may not provide optimal convergence rate if the diffusion is actually variable in the cell!
    // If the diffusion is piecewise constant, choosing the previous version ensures that, for edge integrals, it's the value inside the cell that is computed (essential in case kappa is discontinuous across the edges)
//      return kappa(x, y, cell);
  };

  // Compute cell quadrature nodes and values of cell basis functions and gradients at these nodes
  QuadratureRule quadT = elquad.get_quadT();
  std::vector<Eigen::ArrayXd> phiT_quadT = elquad.get_phiT_quadT();
  std::vector<Eigen::ArrayXXd> dphiT_quadT = elquad.get_dphiT_quadT();
  // Diffusion tensor at the quadrature nodes
  std::vector<Eigen::Matrix2d> kappaT_quadT(quadT.size());
  std::transform(quadT.begin(), quadT.end(), kappaT_quadT.begin(),
      [&kappaT,&cell](QuadratureNode qr) -> Eigen::MatrixXd { return kappaT(qr.x, qr.y); });


_itime[0] += timeint.elapsed().user + timeint.elapsed().system;
timeint.start();

  // Preparing matrices to store face mass matrices (for face and cell basis functions), and face-cell mass matrices
  // (will be useful for the stabilisation term)
  std::vector<Eigen::MatrixXd> MFF(mesh->cell(iT)->n_edges(), Eigen::MatrixXd::Zero(hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs()));
  std::vector<Eigen::MatrixXd> MFT(mesh->cell(iT)->n_edges(), Eigen::MatrixXd::Zero(hho.nlocal_edge_dofs(), hho.nhighorder_dofs()));

_itime[1] += timeint.elapsed().user + timeint.elapsed().system;
timeint.start();


  //-------------------- Compute PT, matrix of potential reconstruction ---------//

  // Stiffness matrix (kappaT dphi_i,dphi_j)_T for phi_i, phi_j up to degree K+1
  Eigen::MatrixXd StiffT = hho.gram_matrix(dphiT_quadT, dphiT_quadT, hho.nhighorder_dofs(), hho.nhighorder_dofs(), quadT, true, kappaT_quadT);  

  // Mass matrix of (phi_i,phi_j)_T for phi_i up to degree L and phi_j up to degree K+1
  Eigen::MatrixXd MTT = hho.gram_matrix(phiT_quadT, phiT_quadT, hho.nlocal_cell_dofs(), hho.nhighorder_dofs(), quadT, true);

  // Row vector of (phi_j,1)_T for phi_j up to degree K+1, and LT^t*LT (used later to impose average condition on PT)
  Eigen::VectorXd LT = (MTT.row(0)).transpose();
  Eigen::MatrixXd LTtLT = LT * (LT.transpose());

  // Right-hand side: we start with volumetric terms (dphi_i,dphi_j)_T for phi_i up to degree K+1 and phi_j up to degree L
  Eigen::MatrixXd BP = Eigen::MatrixXd::Zero(hho.nhighorder_dofs(), local_dofs);
  BP.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs()) = StiffT.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs());

  // Boundary terms in BP
  for (size_t ilF = 0; ilF < mesh->cell(iT)->n_edges(); ilF++) {
    // Offset for face unknowns
    const size_t offset_F = hho.nlocal_cell_dofs() + ilF * hho.nlocal_edge_dofs();
    const auto& nTF = mesh->cell(iT)->edge_normal(ilF);

    // Compute face quadrature nodes and values of basis functions (and gradients) at these nodes
    QuadratureRule quadF = elquad.get_quadF(ilF);
    size_t nbqF = quadF.size();
    
    std::vector<Eigen::ArrayXd> phiT_quadF = elquad.get_phiT_quadF(ilF);
    std::vector<Eigen::ArrayXd> phiF_quadF = elquad.get_phiF_quadF(ilF);
    std::vector<Eigen::ArrayXXd> dphiT_quadF = elquad.get_dphiT_quadF(ilF);

    // Calculation of RHS BP
    //
    // kappaT*nTF on quadrature nodes
    std::vector<Eigen::Vector2d> quadF_kappaT_nTF(nbqF, Eigen::Vector2d::Zero());
    for (size_t iqn = 0; iqn < nbqF; iqn++){
      quadF_kappaT_nTF[iqn] = quadF[iqn].w * kappaT(quadF[iqn].x, quadF[iqn].y)*nTF;
    }
    for (size_t i = 1; i < hho.nhighorder_dofs(); i++) {
      // We do not need i=0 because it corresponds to dphi_i=0
      for (size_t iqn = 0; iqn < nbqF; iqn++){
        // Contribution of test function to BP(i, .)
        double contrib_i = (dphiT_quadF[i].col(iqn).matrix()).dot(quadF_kappaT_nTF[iqn]);
        // Integration
        for (size_t j = 0; j < hho.nlocal_cell_dofs(); j++) {
          // Cell unknown
          BP(i,j) -= contrib_i * phiT_quadF[j](iqn);
        }
        for (size_t j = 0; j < hho.nlocal_edge_dofs(); j++) {
          // Edge unknowns
          BP(i, offset_F + j) += contrib_i * phiF_quadF[j](iqn);
        }
      }
    }


    // Face mass matrix, and cell-face mass matrix, for stabilisation term below
    MFF[ilF] = hho.gram_matrix(phiF_quadF, phiF_quadF, hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs(), quadF, true);
    MFT[ilF] = hho.gram_matrix(phiF_quadF, phiT_quadF, hho.nlocal_edge_dofs(), hho.nhighorder_dofs(), quadF, false);
 
  }

  // Compute PT
  double scalT = StiffT.trace() / std::pow(LT.norm(), 2);
  BP.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs()) += 
          scalT * LTtLT.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs());
  Eigen::MatrixXd PT = ((StiffT+scalT*LTtLT).ldlt()).solve(BP);

  // Consistent component (K \nabla pT, \nabla pT)_T in local bilinear form
  Eigen::MatrixXd ATF = PT.transpose() * StiffT * PT;

  //-------------------- Compute stabilisation term sT ---------//

  Eigen::MatrixXd STF = Eigen::MatrixXd::Zero(local_dofs, local_dofs);

  // Cell residual delta_T^l = pi_T^l (rT uT) - u_T
  Eigen::MatrixXd MTTKinv = ( MTT.topLeftCorner(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs()) ).inverse();
  Eigen::MatrixXd deltaTL = MTTKinv * MTT * PT;
  deltaTL.topLeftCorner(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs()) -=
        Eigen::MatrixXd::Identity(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs());

  for (size_t ilF = 0; ilF < nedgesT; ilF++) {
    // Two options for stabilisation: diameter of edge, or ratio measure cell/measure edge
//    double dTF = cell->edge(ilF)->diam();
    double dTF = cell->measure() / cell->edge(ilF)->measure();

    Eigen::Vector2d xF = cell->edge(ilF)->center_mass();

    auto kappa_TF = kappaT(xF.x(), xF.y()).trace();

    // Face residual delta_TF^k = pi_F^k (rT uT) - u_F
    Eigen::MatrixXd MFFinv = MFF[ilF].inverse();
    Eigen::MatrixXd deltaTFK = MFFinv * MFT[ilF] * PT;
    deltaTFK.block(0, hho.nlocal_cell_dofs() + ilF * hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs()) -=
        Eigen::MatrixXd::Identity(hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs());

    // Stabilisation term
    Eigen::MatrixXd deltaTFK_minus_deltaTL = deltaTFK - MFFinv * MFT[ilF].topLeftCorner(hho.nlocal_edge_dofs(), hho.nlocal_cell_dofs()) * deltaTL;

    STF += (kappa_TF / dTF) * deltaTFK_minus_deltaTL.transpose() * MFF[ilF] *  deltaTFK_minus_deltaTL;

  }
_itime[2] += timeint.elapsed().user + timeint.elapsed().system;
timeint.start();

  // Adjust local bilinear form with stabilisation term
  ATF += STF;

  return ATF;

}


//******************************** 
//    local load term 
//********************************

Eigen::VectorXd HHO_Diffusion::load_operator(HybridCore &hho, const size_t iT, const ElementQuad &elquad) const {
  // Load for the cell DOFs (first indices) and face DOFs (last indices)
  const auto mesh = hho.get_mesh_ptr();
  Cell* cell = mesh->cell(iT);
  size_t cell_edge_dofs = hho.nlocal_cell_dofs() + cell->n_edges()*hho.nlocal_edge_dofs();
  Eigen::VectorXd b = Eigen::VectorXd::Zero(cell_edge_dofs);

  // Quadrature nodes and values of cell basis functions at these nodes
  QuadratureRule quadT = elquad.get_quadT();
  size_t nbq = quadT.size();
  std::vector<Eigen::ArrayXd> phiT_quadT = elquad.get_phiT_quadT();

  // Value of source times quadrature weights at the quadrature nodes
  Eigen::ArrayXd weight_source_quad = Eigen::ArrayXd::Zero(nbq);
  for (size_t iqn = 0; iqn < nbq; iqn++){
    weight_source_quad(iqn) = quadT[iqn].w * source(quadT[iqn].x, quadT[iqn].y, cell);
  }

  for (size_t i=0; i < hho.nlocal_cell_dofs(); i++){
    b(i) = (weight_source_quad * phiT_quadT[i]).sum();
  }
  // Boundary values, if we have a boundary cell
  if (cell->is_boundary()){
    if (BC==0){
      // Dirichlet BCs: no source terms on these edges
    } else if (BC==1) {
      // Neumann BCs
      for (size_t ilF = 0; ilF < cell->n_edges(); ilF++) {
        Edge* edge = cell->edge(ilF);
        const size_t iF = edge->global_index(); 
        // BC on boundary faces
        if (edge->is_boundary()){
          // Offset for face unknowns
          const size_t offset_F = hho.nlocal_cell_dofs() + ilF * hho.nlocal_edge_dofs();
          // Normal to the edge
          const auto& nTF = cell->edge_normal(ilF);
          // for each DOF of the boundary face
          for (size_t i = 0; i < hho.nlocal_edge_dofs(); i++){
            const auto& phi_i = hho.edge_basis(iF, i);
            QuadratureRule quadF = generate_quadrature_rule(*edge, 2*hho.K()+2);
            std::function<double(double,double)> Kgrad_n = [&](double x, double y){
                return nTF.dot(kappa(x,y,cell) * grad_exact_solution(x,y,cell)) * phi_i(x,y);
              };
            for (QuadratureNode qF : quadF){
              b(offset_F + i) += qF.w * Kgrad_n(qF.x,qF.y);
            }
          }
        }
      }
    }
  }
  
  return b;
}



//******************************** 
//    energy norm 
//********************************

double HHO_Diffusion::EnergyNorm(HybridCore& hho, const Eigen::VectorXd Xh) {
  const auto mesh = hho.get_mesh_ptr();
  double value = 0.0;

  for (size_t iT = 0; iT < mesh-> n_cells(); iT++) {
    Eigen::VectorXd XTF = hho.restr(Xh, iT);
    value += XTF.transpose() * aT[iT] * XTF;
  }

  return sqrt(value);

}



double HHO_Diffusion::get_assembly_time() const {
  return double(_assembly_time) * pow(10, -9);
}

double HHO_Diffusion::get_solving_time() const {
  return double(_solving_time) * pow(10, -9);
}

double HHO_Diffusion::get_itime(size_t idx) const {
  return double(_itime[idx]) * pow(10, -9);
}

double HHO_Diffusion::get_solving_error() const {
  return _solving_error;
}

//@}
} // end of namespace HArDCore2D

#endif //_HHO_DIFFUSION_HPP
