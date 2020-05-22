// Implementation of the NC Poly scheme in 2D for the Stefan and PME models
//
//   { u - div(K \grad(zeta(u))) = f,       inside Omega
//   { K \grad(zeta(u)) . nTF = g,       on GammaN
//   { 								u = g,			 on GammaD
// 
//	At the moment, only pure Neumann or pure Dirichlet
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef _LEPNC_STEFANPME_HPP
#define _LEPNC_STEFANPME_HPP

#include <functional>
#include <utility>
#include <iostream>

#include <boost/timer/timer.hpp>

// Matrices and linear solvers
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include "Eigen/MA41.cpp"

#include "mesh.hpp"
#include "lepnccore.hpp"
#include "quad2d.hpp"
#include "TestCase/TestCaseNonLinearity.hpp"
#include "TestCase/TestCaseStefanPME.hpp"
#include "TestCase/BoundaryConditions.hpp"

/*!
* @defgroup LEPNC
*/

namespace HArDCore2D {

/*!
*	@addtogroup LEPNC
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------
/* The LEPNC_StefanPME class provides an implementation of the LEPNC method for the stationnary Stefan and PME problems.

* If using this code in a scientific publication, please cite the reference for the LEPNC:
*
* Non-conforming finite elements on polytopal mesh,
* J. Droniou, R. Eymard, T. Gallouët and R. Herbin
* url:
*
*/
/// The vector Xh manipulated in the resolution has mixed components, corresponding either to the unknown u or to \f$\zeta(u)\f$, depending on the choice of weight of mass-lumping for the cell/edge unknowns. If no weight is put on the edges (resp. the cells), then the edge (resp. cell) unknowns represent \f$\zeta(u)\f$. Otherwise, they represent u.

class LEPNC_StefanPME { 

// Types
public:
  using solution_function_type = std::function<double(double,double)>;		///< type for solution
  using source_function_type = std::function<double(double,double,Cell*)>;		///< type for source
  using grad_function_type = std::function<VectorRd(double,double,Cell*)>;		///< type for gradient
  using tensor_function_type = std::function<Eigen::Matrix2d(double,double,Cell*)>;		///< type for diffusion tensor

	///@brief Constructor of the class
  LEPNC_StefanPME(
		LEPNCCore &nc,										///< instance of the core nonconforming class (with basis functions, etc.)
		tensor_function_type kappa, 	///< diffusion tensor
		size_t deg_kappa,							///< polynomial degree of the diffusion tensor
		source_function_type source,  ///< source term
		BoundaryConditions BC,					///< type of boundary conditions
		solution_function_type exact_solution, 	///< exact solution
		grad_function_type grad_exact_solution, 	///< gradient of the exact solution
		TestCaseNonLinearity::nonlinearity_function_type zeta,				///< function describing the nonlinearity
		double weight,				///< proportion of weight (between 0 and 1) of mass-lumping on the edges
		std::string solver_type,		///< type of solver to use for the global system (bicgstab at the moment)
    std::ostream & output = std::cout    ///< optional argument for output of messages
		);

  /// Assemble and solve the scheme
  Eigen::VectorXd solve();

  /// Compute non-linearity on vector (depends if weight=0, weight=1 or weight\in (0,1) )
  Eigen::VectorXd apply_nonlinearity(const Eigen::VectorXd& Y, const std::string type) const;

	/// Mass-lumped L2 norm of a function given by a vector
  double L2_MassLumped(const Eigen::VectorXd Xh) const; 

	/// Discrete energy norm (associated to the diffusion operator)
  double EnergyNorm(const Eigen::VectorXd Xh) const; 

  double get_assembly_time() const; ///< cpu time to assemble the scheme
  double get_solving_time() const;	///< cpu time to solve the scheme
  double get_solving_error() const;	///< residual after solving the scheme
  double get_itime(size_t idx) const;		///< various intermediate assembly times

private:
  /// Compute the local diffusion operator in the cell iT
  Eigen::MatrixXd diffusion_operator(const size_t iT) const;

  /// Compute the local load operator (the source term) in the cell iT. The second version computes the load using the mass-lumping (which, in particular, puts no load on basis functions associated with the edges)
  Eigen::VectorXd load_operator(const size_t iT) const;
  Eigen::VectorXd load_operator_ml(const size_t iT) const;

  /// Compute residual of non-linear scheme on vector Xh: source - scheme(Xh)
	/// The vectors DiffT, MassT and SourceT must have been calculated and stored before calling this function
	Eigen::VectorXd residual_scheme(const Eigen::VectorXd& Xh) const;

	const LEPNCCore& nc;
  const tensor_function_type kappa;
	size_t _deg_kappa;
  const source_function_type source;
	const BoundaryConditions m_BC;
  const solution_function_type exact_solution;
  const grad_function_type grad_exact_solution;
	const TestCaseNonLinearity::nonlinearity_function_type zeta;
	const double _weight;
	const std::string solver_type;
  std::ostream & m_output;

	// To store local diffusion, mass and source terms
	std::vector<Eigen::MatrixXd> DiffT;
	std::vector<Eigen::MatrixXd> MassT;
	std::vector<Eigen::VectorXd> SourceT;

  // Computation statistics
  size_t _assembly_time;
  size_t _solving_time;
  double _solving_error;
	mutable std::vector<size_t> _itime = std::vector<size_t>(10, 0);

};

LEPNC_StefanPME::LEPNC_StefanPME(LEPNCCore &nc, tensor_function_type kappa, size_t deg_kappa, source_function_type source, BoundaryConditions BC, solution_function_type exact_solution, grad_function_type grad_exact_solution, TestCaseNonLinearity::nonlinearity_function_type zeta, double weight, std::string solver_type, std::ostream & output)
  : nc(nc),
		kappa(kappa),
		_deg_kappa(deg_kappa),
    source(source),
		m_BC(BC),
		exact_solution(exact_solution),
		grad_exact_solution(grad_exact_solution),
		zeta(zeta),
		_weight(weight),
		solver_type(solver_type),
    m_output(output) {

		//---- Compute diffusion, mass matrix and source (do not change during Newton) ------//
		const Mesh* mesh = nc.get_mesh();
		DiffT.resize(mesh->n_cells());
		MassT.resize(mesh->n_cells());
		SourceT.resize(mesh->n_cells());
		for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
			size_t n_edgesT = mesh->cell(iT)->n_edges();
			size_t n_local_dofs = DimPoly<Cell>(1) + n_edgesT;
			double measT = mesh->cell(iT)->measure();
			DiffT[iT] = diffusion_operator(iT);
			MassT[iT] = Eigen::MatrixXd::Zero(n_local_dofs, n_local_dofs);
			MassT[iT].topLeftCorner(DimPoly<Cell>(1), DimPoly<Cell>(1)) = (1-_weight) * (measT / DimPoly<Cell>(1)) * Eigen::MatrixXd::Identity(DimPoly<Cell>(1), DimPoly<Cell>(1));
			MassT[iT].bottomRightCorner(n_edgesT, n_edgesT) = _weight * (measT / n_edgesT) * Eigen::MatrixXd::Identity(n_edgesT, n_edgesT);
			SourceT[iT] = load_operator_ml(iT);
		}
}

// ASSEMBLE AND SOLVE THE SCHEME
Eigen::VectorXd LEPNC_StefanPME::solve() {

	boost::timer::cpu_timer timer; 
	boost::timer::cpu_timer timerint; 
  const auto mesh = nc.get_mesh();
	size_t n_edges_dofs = mesh->n_edges();
	size_t n_cell_dofs = mesh->n_cells() * DimPoly<Cell>(1);

	//----- PARAMETERS FOR NEWTON ------------//
	constexpr double tol = 1e-6;
	constexpr size_t maxiter = 400;
  // relaxation parameter
  double relax = 1;     
	size_t iter = 0;
	Eigen::VectorXd Xhprev = Eigen::VectorXd::Zero(n_cell_dofs + n_edges_dofs);
	Eigen::VectorXd Xh = Eigen::VectorXd::Zero(n_cell_dofs + n_edges_dofs);

	// Vector of Dirichlet BC, either on u (if weight>0) or zeta(u) (if weight=0)
	Eigen::VectorXd DirBC = Eigen::VectorXd::Zero(n_cell_dofs + n_edges_dofs);
	for (size_t ibF = 0; ibF < mesh->n_b_edges(); ibF++){
    Edge* edge = mesh->b_edge(ibF);
    if (m_BC.type(*edge)=="dir"){
		  size_t iF = edge->global_index();
		  QuadratureRule quadF = generate_quadrature_rule(*edge, 5);
		  for (QuadratureNode quadrule : quadF){
			  if (_weight == 0){
				  DirBC(n_cell_dofs + iF) += quadrule.w * zeta(exact_solution(quadrule.x, quadrule.y), "fct");
			  }else{
				  DirBC(n_cell_dofs + iF) += quadrule.w * exact_solution(quadrule.x, quadrule.y);
			  }
		  }
		  // Take average
		  DirBC(n_cell_dofs + iF) /= mesh->b_edge(ibF)->measure();
    }
	}
	Xhprev = DirBC;


	// ---- NEWTON ITERATONS ------ //
	Eigen::VectorXd RES = residual_scheme(Xhprev); 		// Scheme is F(X)=C, then RES = C - F(u_prev)
  std::vector<double> residual(maxiter, 1.0);
	residual[iter] = RES.norm();

	while ( (iter < maxiter) && (residual[iter] > tol) ){
		iter++;

		// System matrix and initialisation of RHS (only on the edge dofs)
		Eigen::SparseMatrix<double> GlobMat(n_edges_dofs, n_edges_dofs);
		std::vector<Eigen::Triplet<double>> triplets_GlobMat;
		Eigen::VectorXd GlobRHS = RES.tail(n_edges_dofs);

		// Static condensation: matrix and source to recover cell unknowns
		Eigen::SparseMatrix<double> ScMat(n_cell_dofs, n_edges_dofs);
		std::vector<Eigen::Triplet<double>> triplets_ScMat;
		Eigen::VectorXd ScRHS = Eigen::VectorXd::Zero(n_cell_dofs);

		// Assemble local contributions
		for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
			Cell* cell = mesh->cell(iT);
			size_t n_edgesT = cell->n_edges();
		
			// Local static condensation of element unknowns
			Eigen::VectorXd XTprev = nc.nc_restr(Xhprev, iT);
			Eigen::VectorXd ZetaprimeXTprev = apply_nonlinearity(XTprev, "der");
			Eigen::MatrixXd MatZetaprimeXTprev = Eigen::MatrixXd::Zero(DimPoly<Cell>(1)+n_edgesT, DimPoly<Cell>(1)+n_edgesT);
			for (size_t i = 0; i < DimPoly<Cell>(1) + n_edgesT; i++){
				MatZetaprimeXTprev(i, i) = ZetaprimeXTprev(i);
			}
			Eigen::MatrixXd MatT = DiffT[iT] * MatZetaprimeXTprev + MassT[iT];
			Eigen::MatrixXd ATT = MatT.topLeftCorner(DimPoly<Cell>(1), DimPoly<Cell>(1));
			Eigen::MatrixXd ATF = MatT.topRightCorner(DimPoly<Cell>(1), n_edgesT);
			Eigen::MatrixXd AFT = MatT.bottomLeftCorner(n_edgesT, DimPoly<Cell>(1));
			Eigen::MatrixXd AFF = MatT.bottomRightCorner(n_edgesT, n_edgesT);

			Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
			invATT.compute(ATT);
				
			Eigen::MatrixXd invATT_ATF = invATT.solve(ATF);
			Eigen::VectorXd RES_cell = RES.segment(iT * DimPoly<Cell>(1), DimPoly<Cell>(1));
			Eigen::VectorXd invATT_RES_cell = invATT.solve(RES_cell);

			// Local matrix and right-hand side on the face unknowns
			Eigen::MatrixXd MatF = Eigen::MatrixXd::Zero(n_edgesT, n_edgesT);
			Eigen::VectorXd bF = Eigen::VectorXd::Zero(n_edgesT);
			MatF = AFF - AFT * invATT_ATF;
			bF = - AFT * invATT_RES_cell;		// only additional RHS coming from static condensation, beyond RES
							
			// Assemble static condensation operator
			ScRHS.segment(iT * DimPoly<Cell>(1), DimPoly<Cell>(1)) = invATT_RES_cell;
			for (size_t i = 0; i < DimPoly<Cell>(1); i++){
				size_t iGlobal = iT * DimPoly<Cell>(1) + i;
				for (size_t j = 0; j < n_edgesT; j++){
					size_t jGlobal = cell->edge(j)->global_index();
		      triplets_ScMat.emplace_back(iGlobal, jGlobal, invATT_ATF(i, j));
				}
			}

		  // GLOBAL MATRIX and RHS on the edge unknowns, using the matrices MatF obtained after static condensation
			for (size_t i = 0; i < n_edgesT; i++){
				size_t iGlobal = cell->edge(i)->global_index();
				for (size_t j = 0; j < n_edgesT; j++){
					size_t jGlobal = cell->edge(j)->global_index();
					triplets_GlobMat.emplace_back(iGlobal, jGlobal, MatF(i, j));
				}
		  	GlobRHS(iGlobal) += bF(i);
			}
		}

		// Assemble the global linear system (without BC), and matrix to recover statically-condensed cell dofs
		GlobMat.setFromTriplets(std::begin(triplets_GlobMat), std::end(triplets_GlobMat));
		ScMat.setFromTriplets(std::begin(triplets_ScMat), std::end(triplets_ScMat));

		//	Dirichlet boundary conditions are trivial. Newton requires to solve
		//			F'(X^k) (X^{k+1} - X^k) = B - F(X^k)
		//	Both X^k and X^{k+1} have fixed Dirichlet values on boundary edges, so the system we solve
    //  is Dirichlet homogeneous, and we just select the non-dirichlet edges (first ones in the list)
    size_t n_edge_unknowns = mesh->n_edges() - m_BC.n_dir_edges();
		Eigen::SparseMatrix<double> A = GlobMat.topLeftCorner(n_edge_unknowns, n_edge_unknowns);
//std::cout << "Number non-zero terms in matrix: " << A.nonZeros() << "\n";
		Eigen::VectorXd B = GlobRHS.head(n_edge_unknowns);

		// Solve condensed system and recover cell unknowns. dX = X^{k+1}-X^k//
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
		solver.compute(A);
		Eigen::VectorXd dX_edge_unknowns = solver.solve(B);

		Eigen::VectorXd dX = Eigen::VectorXd::Zero(n_cell_dofs + n_edges_dofs);
		dX.segment(n_cell_dofs, n_edge_unknowns) = dX_edge_unknowns;
		dX.head(n_cell_dofs) = ScRHS - ScMat * dX.tail(n_edges_dofs);

		// Recover the fixed boundary values  and cell unknowns (from static condensation and Newton).
		// We start by putting the Dirichlet BC at the end.
		Xh = DirBC;
    Xh.head(n_cell_dofs + n_edge_unknowns) = relax*dX.head(n_cell_dofs + n_edge_unknowns) + Xhprev.head(n_cell_dofs + n_edge_unknowns);

    // Compute new residual. We might not accept new Xh but start over reducing relax. Otherwise, accept Xh and increase relax
    // to a max of 1
		Eigen::VectorXd RES_temp = residual_scheme(Xh);
    double residual_temp = RES_temp.norm();
    if (iter > 1 && residual_temp > 10*residual[iter-1]){
      // We don't update the residual and Xh, but we will do another iteration with relax reduced
      relax = relax/2.0;
      iter--;
    } else {
      // Increase relax
      relax = std::min(1.0, 1.2*relax);
  		// Store Xh in Xhprev, compute new residual
  		Xhprev = Xh;
  		RES = RES_temp;
  		residual[iter] = residual_temp;
    }

		m_output << "     ...Iteration: " << iter << "; residual = " << residual[iter] << "; relax = "<< relax << "\n";

	} //-------- END NEWTON ITERATIONS -----------//


  return Xh;
}

//****************************************
//		apply nonlinearity zeta to a vector
//****************************************

Eigen::VectorXd LEPNC_StefanPME::apply_nonlinearity(const Eigen::VectorXd& Y, const std::string type) const {
	// Type="fct", "nlin" or "der" as per the nonlinear function zeta
	// The function applied to each coefficient of Y depends on its nature. If it's a coefficient appearing in
	// the mass-lumping (which depends on _weight), the unknown represents u and we apply the full nonlinear function zeta.
	// Otherwise, the unknown is zeta(u) itself and the function we apply is the identity s->s (with nlin/der: s-> 1)
	// 

	Eigen::VectorXd ZetaY;
	// Initialisation depends on "type". We initialise as if the nonlinearity was the identity, we'll change the
	// appropriate coefficients afterwards
	if (type == "fct"){
		ZetaY = Y;
	}else if (type == "der" || type == "nlin"){
		ZetaY = Eigen::VectorXd::Ones(Y.size());
	}else{
		m_output << "type unknown in apply_nonlinearity: " << type << "\n";
		exit(EXIT_FAILURE);
	}

	const Mesh* mesh = nc.get_mesh();
	// We check if Y is a local vector or a global one, this determines the number of cell and edge unknowns
	size_t n_cell_dofs = 0;
	if (size_t(Y.size()) == DimPoly<Cell>(1)*mesh->n_cells() + mesh->n_edges()){
		n_cell_dofs = DimPoly<Cell>(1)*mesh->n_cells();
	}else{
		n_cell_dofs = DimPoly<Cell>(1);
	}
	size_t n_edges_dofs = Y.size() - n_cell_dofs;

	// The type of nonlinearity we apply on each coefficient depends on _weight
	if (_weight == 0){
		for (size_t i = 0; i < n_cell_dofs; i++){
			ZetaY(i) = zeta(Y(i), type);
		}
	}else if (_weight == 1){
		for (size_t i = n_cell_dofs; i < n_cell_dofs + n_edges_dofs; i++){
			ZetaY(i) = zeta(Y(i), type);
		}
	}else{
		for (size_t i = 0; i < n_cell_dofs + n_edges_dofs; i++){
			ZetaY(i) = zeta(Y(i), type);
		}
	}


	return ZetaY;
}


//******************************** 
//		local diffusion matrix 
//********************************

Eigen::MatrixXd LEPNC_StefanPME::diffusion_operator(const size_t iT) const {

  const auto mesh = nc.get_mesh();
	Cell* cell = mesh->cell(iT);
	const size_t n_edgesT = cell->n_edges();

  // Total number of degrees of freedom local to this cell
  size_t local_dofs = n_edgesT + DimPoly<Cell>(1);
	
	// Diffusion in the cell.
	std::function<Eigen::Matrix2d(double,double)> kappaT = [&](double x, double y){
			return kappa(x, y, cell);
	};

	// Compute cell quadrature nodes and values of cell basis functions and gradients at these nodes
	size_t doe = 2 + _deg_kappa;
	QuadratureRule quadT = generate_quadrature_rule(*cell, doe, true);
	std::vector<Eigen::ArrayXXd> Dnc_phiT_quadT = nc.grad_nc_basis_quad(iT, quadT);
	// Diffusion tensor at the quadrature nodes
	std::vector<Eigen::Matrix2d> kappaT_quadT(quadT.size());
	std::transform(quadT.begin(), quadT.end(), kappaT_quadT.begin(),
			[&kappaT,&cell](QuadratureNode qr) -> Eigen::MatrixXd { return kappaT(qr.x, qr.y); });

	return nc.gram_matrix(Dnc_phiT_quadT, Dnc_phiT_quadT, local_dofs, local_dofs, quadT, true, kappaT_quadT);  
}


//******************************** 
//		local load term 
//********************************

Eigen::VectorXd LEPNC_StefanPME::load_operator(const size_t iT) const {
	// Load for the cell DOFs (first indices) and face DOFs (last indices)
  const auto mesh = nc.get_mesh();
	Cell* cell = mesh->cell(iT);
	size_t n_edgesT = cell->n_edges();
	size_t local_dofs = DimPoly<Cell>(1) + n_edgesT;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(local_dofs);

	// Quadrature nodes and values of cell basis functions at these nodes
	size_t doe = 5;
	QuadratureRule quadT = generate_quadrature_rule(*cell, doe, true);
	size_t nbq = quadT.size();
	std::vector<Eigen::ArrayXd> nc_phiT_quadT = nc.nc_basis_quad(iT, quadT);

	// Value of source times quadrature weights at the quadrature nodes
	Eigen::ArrayXd weight_source_quad = Eigen::ArrayXd::Zero(nbq);
	for (size_t iqn = 0; iqn < nbq; iqn++){
		weight_source_quad(iqn) = quadT[iqn].w * source(quadT[iqn].x, quadT[iqn].y, cell);
	}

	for (size_t i=0; i < local_dofs; i++){
		b(i) = (weight_source_quad * nc_phiT_quadT[i]).sum();
	}

	// Boundary values, if we have a boundary cell and Neumann edges
	if (cell->is_boundary()){
    // Boundary values only on boundary Neumann edges
    for (size_t ilF = 0; ilF < cell->n_edges(); ilF++) {
      Edge* edge = cell->edge(ilF);
      if (m_BC.type(*edge)=="neu"){
				// BC on boundary faces
				if (edge->is_boundary()){
				  const auto& nTF = cell->edge_normal(ilF);
			    const auto& phi_i = nc.nc_basis(iT, DimPoly<Cell>(1) + ilF);
					QuadratureRule quadF = generate_quadrature_rule(*edge, 5);
					std::function<double(VectorRd)> Kgrad_n = [&](VectorRd p){
              return zeta(exact_solution(p.x(), p.y()), "der") * nTF.dot(kappa(p.x(),p.y(),cell) * grad_exact_solution(p.x(),p.y(),cell)) * phi_i(p.x(),p.y());
						};
					for (QuadratureNode qr : quadF){
						b(DimPoly<Cell>(1) + ilF) += qr.w * Kgrad_n(qr.vector());
					}
				}
			}
		}
	}
	
  return b;
}

Eigen::VectorXd LEPNC_StefanPME::load_operator_ml(const size_t iT) const {
	// Load for the cell DOFs (first indices) and face DOFs (last indices, will be 0)
  const auto mesh = nc.get_mesh();
	Cell* cell = mesh->cell(iT);
	size_t n_edgesT = cell->n_edges();
	size_t local_dofs = DimPoly<Cell>(1) + n_edgesT;

	// Beware, MassT[iT] must have been created
	if (MassT.size() < iT || size_t(MassT[iT].rows()) != local_dofs){
		m_output << "Called load_operator_ml without creating MassT[iT]\n";
		exit(EXIT_FAILURE);
	}
	Eigen::VectorXd fvec = Eigen::VectorXd::Zero(DimPoly<Cell>(1) + n_edgesT);
	for (size_t i = 0; i < local_dofs; i++){
		VectorRd node;
		if (i < DimPoly<Cell>(1)){
			node = nc.ml_node(iT, i);
		}else{
			node = cell->edge(i-DimPoly<Cell>(1))->center_mass();
		}
		fvec(i) = source(node.x(), node.y(), cell);
	}
  // Contribution of source to loading term
  Eigen::VectorXd load = MassT[iT] * fvec;

  // Adding Neumann boundary conditions
  // The boundary conditions coming from Neumann edges are computed based on the piecewise constant
  // discrete trace: if g is the value of the outer flux,
  //      \int g v -> \int g Tv
  //  where (Tv)_e = average of v on e, for e Neumann edge.
  // Only edge-based basis functions therefore have a contribution, which will simply be \int_e g.
  if (cell->is_boundary()){
    for (size_t ilF = 0; ilF < cell->n_edges(); ilF++) {
      Edge* edge = cell->edge(ilF);
      if (m_BC.type(*edge)=="neu"){
				if (edge->is_boundary()){
				  const auto& nTF = cell->edge_normal(ilF);
					QuadratureRule quadF = generate_quadrature_rule(*edge, 5);
					std::function<double(VectorRd)> Kgrad_n = [&](VectorRd p){
								return zeta(exact_solution(p.x(), p.y()), "der") * nTF.dot(kappa(p.x(),p.y(),cell) * grad_exact_solution(p.x(),p.y(),cell));
						};
					for (QuadratureNode qr : quadF){
						load(DimPoly<Cell>(1) + ilF) += qr.w * Kgrad_n(qr.vector());
					}
				}
			}
		}
	}

  return load;
}

//****************************
// Residual non-linear scheme
//****************************

Eigen::VectorXd LEPNC_StefanPME::residual_scheme(const Eigen::VectorXd& Xh) const{

	const Mesh* mesh = nc.get_mesh();	
	Eigen::VectorXd RES = Eigen::VectorXd::Zero(mesh->n_cells() * DimPoly<Cell>(1) + mesh->n_edges());

	// Compute the residual: C - (DIFF * ZetaXh + MASS * Xh), where DIFF is the diffusion matrix, MASS is the mass matrix and C is the source term
	size_t n_cell_dofs = mesh->n_cells() * DimPoly<Cell>(1);
	for (size_t iT = 0; iT < mesh->n_cells(); iT++){
		Cell* cell =  mesh->cell(iT);
		size_t n_edgesT = cell->n_edges();
		// Local unknown in the cell, and zeta of these unknowns
		Eigen::VectorXd XT = nc.nc_restr(Xh, iT);
		Eigen::VectorXd ZetaXT = apply_nonlinearity(XT, "fct");

		Eigen::VectorXd REST = SourceT[iT] - (DiffT[iT] * ZetaXT + MassT[iT] * XT);
		RES.segment(iT*DimPoly<Cell>(1), DimPoly<Cell>(1)) += REST.head(DimPoly<Cell>(1));
		for (size_t ilE = 0; ilE < n_edgesT; ilE++){
			size_t iE = cell->edge(ilE)->global_index();
			RES(n_cell_dofs + iE) += REST(DimPoly<Cell>(1)+ilE);
		}
	}

	// No residual on Dirichlet edges (do not correspond to equations of the system)
	RES.tail(m_BC.n_dir_edges()) = Eigen::VectorXd::Zero(m_BC.n_dir_edges());

	return RES;
}

//*******************************************************
// Norms: L2 mass lumped, Energy (pure diffusion)
//*******************************************************

double LEPNC_StefanPME::L2_MassLumped(const Eigen::VectorXd Xh) const {
  const Mesh* mesh = nc.get_mesh();
	double value = 0.0;

  for (size_t iT = 0; iT < mesh-> n_cells(); iT++) {
		Eigen::VectorXd XTF = nc.nc_restr(Xh, iT);
		value += XTF.transpose() * MassT[iT] * XTF;
  }

  return sqrt(value);
}

double LEPNC_StefanPME::EnergyNorm(const Eigen::VectorXd Xh) const {
  const Mesh* mesh = nc.get_mesh();
	double value = 0.0;

  for (size_t iT = 0; iT < mesh-> n_cells(); iT++) {
		Eigen::VectorXd XTF = nc.nc_restr(Xh, iT);
		value += XTF.transpose() * DiffT[iT] * XTF;
  }

  return sqrt(value);
}


double LEPNC_StefanPME::get_assembly_time() const {
  return double(_assembly_time) * pow(10, -9);
}

double LEPNC_StefanPME::get_solving_time() const {
  return double(_solving_time) * pow(10, -9);
}

double LEPNC_StefanPME::get_itime(size_t idx) const {
  return double(_itime[idx]) * pow(10, -9);
}

double LEPNC_StefanPME::get_solving_error() const {
  return _solving_error;
}

//@}
} // end of namespace HArDCore2D

#endif //_LEPNC_STEFANPME_HPP
