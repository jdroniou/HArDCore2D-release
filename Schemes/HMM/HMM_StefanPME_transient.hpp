// Implementation of the HMM in 2D for the transient Stefan and PME models
//
//   { d_t u - div(K \grad(zeta(u))) = f,       inside Omega
//   { K \grad(zeta(u)) . nTF = g,       on GammaN
//   { 								u = g,			 on GammaD
//   { u(0) = u_ini
//
//	At the moment, only pure Neumann or pure Dirichlet
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef _HMM_STEFANPME_TRANSIENT_HPP
#define _HMM_STEFANPME_TRANSIENT_HPP

#include <functional>
#include <utility>
#include <iostream>

#include <boost/timer/timer.hpp>

// Matrices and linear solvers
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include "Eigen/MA41.cpp"

#include "mesh.hpp"
#include "hybridcore.hpp"
#include "quad2d.hpp"
#include "TestCase/TestCaseNonLinearity.hpp"
#include "TestCase/TestCaseStefanPME.hpp"
#include "BoundaryConditions/BoundaryConditions.hpp"

/*!
* @defgroup HMM 
* @brief Hybrid Mimetic Mixed method
*/

namespace HArDCore2D {

/*!
*	@addtogroup HMM
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------
/* The HMM_StefanPME_Transient class provides an implementation of the HMM method for the transient Stefan and PME problems.
*
*/
/// The vector Xh manipulated in the resolution has mixed components, corresponding either to the unknown u or to \f$\zeta(u)\f$, depending on the choice of weight of mass-lumping for the cell/edge unknowns. If no weight is put on the edges (resp. the cells), then the edge (resp. cell) unknowns represent \f$\zeta(u)\f$. Otherwise, they represent u.

class HMM_StefanPME_Transient { 

// Types
public:
  using solution_function_type = std::function<double(const double&, const VectorRd&)>;		///< type for solution
  using source_function_type = std::function<double(const double&, const VectorRd&, const Cell*)>;		///< type for source
  using grad_function_type = std::function<VectorRd(const double&, const VectorRd&, const Cell*)>;		///< type for gradient
  using tensor_function_type = std::function<Eigen::Matrix2d(const double, const double, const Cell*)>;		///< type for diffusion tensor (does not depend on time)

	///@brief Constructor of the class
  HMM_StefanPME_Transient(
		HybridCore &hmm,										///< instance of the hybridcore class (with basis functions, etc.)
		tensor_function_type kappa, 	///< diffusion tensor
		source_function_type source,  ///< source term
		BoundaryConditions BC,					///< type of boundary conditions
		solution_function_type exact_solution, 	///< exact solution
		grad_function_type grad_exact_solution, 	///< gradient of the exact solution
		TestCaseNonLinearity::nonlinearity_function_type zeta,				///< function describing the nonlinearity
		double weight,				///< proportion of weight (between 0 and 1) of mass-lumping on the edges
		std::string solver_type,		///< type of solver to use for the global system (bicgstab at the moment)
    std::ostream & output = std::cout    ///< optional argument for output of messages
		);

  /// Execute one time iteration
  UVector iterate(
      const double tps,     ///< current time (for computing source)
      const double dt,    ///< time step
      const UVector& Xn    ///< Unkown at previous time step
      );

  /// Compute non-linearity on vector (depends if weight=0, weight=1 or weight\in (0,1) )
  Eigen::VectorXd apply_nonlinearity(const Eigen::VectorXd& Y, const std::string type) const;
  UVector apply_nonlinearity(const UVector& Y, const std::string type) const;

	/// Mass-lumped L2 norm of a function given by a vector
  double L2_MassLumped(const UVector& Xh) const; 

	/// Mass-lumped Lp norm of a function given by a vector
  double Lp_MassLumped(const UVector& Xh, double p) const; 

	/// Discrete energy norm (associated to the diffusion operator)
  double EnergyNorm(const UVector& Xh) const; 

  /// cpu time to assemble the scheme
  inline double get_assembly_time() const {
    return double(_assembly_time) * pow(10, -9);
  }
  /// cpu time to solve the scheme
  inline double get_solving_time() const {
    return double(_solving_time) * pow(10, -9);
  }
  /// various intermediate assembly times
  inline double get_itime(size_t idx) const {
    return double(_itime[idx]) * pow(10, -9);
  }
  /// residual after solving the scheme
  inline double get_solving_error() const {
    return _solving_error; 
  }
  /// number of Newton iterations
  inline size_t get_nb_newton() const {
    return m_nb_newton;
  }
  /// Mass matrix in cell iT
  inline Eigen::MatrixXd get_MassT(size_t iT) const {
    return MassT[iT];
  }


private:
  /// Compute the local diffusion operator in the cell iT
  Eigen::MatrixXd diffusion_operator(const size_t iT) const;

  /// Compute the local load operator (only for the source term) at times tps in the cell iT, using mass-lumping
  Eigen::VectorXd load_operator(const double tps, const size_t iT) const;

  /// Compute residual of non-linear scheme on vector Xh: source - scheme(Xh)
	/// The vectors DiffT, MassT and RightHandSideT must have been calculated and stored before calling this function
	Eigen::VectorXd residual_scheme(const double dt, const UVector& Xh) const;

	const HybridCore& hmm;
  const tensor_function_type kappa;
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
	std::vector<Eigen::VectorXd> RightHandSideT;

  // Computation statistics
  size_t _assembly_time;
  size_t _solving_time;
  double _solving_error;
	mutable std::vector<size_t> _itime = std::vector<size_t>(10, 0);
  size_t m_nb_newton;

};

HMM_StefanPME_Transient::HMM_StefanPME_Transient(HybridCore &hmm, tensor_function_type kappa, source_function_type source, BoundaryConditions BC, solution_function_type exact_solution, grad_function_type grad_exact_solution, TestCaseNonLinearity::nonlinearity_function_type zeta, double weight, std::string solver_type, std::ostream & output)
  : hmm(hmm),
		kappa(kappa),
    source(source),
		m_BC(BC),
		exact_solution(exact_solution),
		grad_exact_solution(grad_exact_solution),
		zeta(zeta),
		_weight(weight),
		solver_type(solver_type),
    m_output(output) {

		//---- Compute diffusion and mass matrix (do not change during time stepping or Newton) ------//
		const Mesh* mesh = hmm.get_mesh();
		DiffT.resize(mesh->n_cells());
		MassT.resize(mesh->n_cells());
    // Source is pre-allocated but not computed as it changes with each time iteration
		RightHandSideT.resize(mesh->n_cells());
		for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
			size_t n_edgesT = mesh->cell(iT)->n_edges();
			size_t n_local_dofs = 1 + n_edgesT;
			double measT = mesh->cell(iT)->measure();
			DiffT[iT] = diffusion_operator(iT);
			MassT[iT] = Eigen::MatrixXd::Zero(n_local_dofs, n_local_dofs);
			MassT[iT](0,0) = (1-_weight) * measT;
			MassT[iT].bottomRightCorner(n_edgesT, n_edgesT) = _weight * (measT / n_edgesT) * Eigen::MatrixXd::Identity(n_edgesT, n_edgesT);
		}
}

// PERFORM ONE TIME ITERATION
UVector HMM_StefanPME_Transient::iterate(const double tps, const double dt, const UVector& Xn) {

	boost::timer::cpu_timer timer; 
	boost::timer::cpu_timer timerint; 
  const auto mesh = hmm.get_mesh();
	size_t n_edges_dofs = mesh->n_edges();
	size_t n_cell_dofs = mesh->n_cells();

  //----- COMPUTE SOURCE: load + previous time -----//
  for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
  	RightHandSideT[iT] = dt * load_operator(tps, iT) + MassT[iT]*Xn.restr(iT);
  }

	//----- PARAMETERS FOR NEWTON ------------//
	constexpr double tol = 1e-8;
	constexpr size_t maxiter = 400;
  // relaxation parameter
  double relax = 1;     
	size_t iter = 0;
	UVector Xhprev = UVector(Eigen::VectorXd::Zero(n_cell_dofs + n_edges_dofs), *hmm.get_mesh(), 0, 0);
	UVector Xh = UVector(Eigen::VectorXd::Zero(n_cell_dofs + n_edges_dofs), *hmm.get_mesh(), 0, 0);

	// Vector of Dirichlet BC, either on u (if weight>0) or zeta(u) (if weight=0)
	Eigen::VectorXd DirBC = Eigen::VectorXd::Zero(n_cell_dofs + n_edges_dofs);
	for (size_t ibF = 0; ibF < mesh->n_b_edges(); ibF++){
    Edge* edge = mesh->b_edge(ibF);
    if (m_BC.type(*edge)=="dir"){
		  size_t iF = edge->global_index();
		  QuadratureRule quadF = generate_quadrature_rule(*edge, 5);
		  for (QuadratureNode quadrule : quadF){
			  if (_weight == 0){
				  DirBC(n_cell_dofs + iF) += quadrule.w * zeta(exact_solution(tps, quadrule.vector()), "fct");
			  }else{
				  DirBC(n_cell_dofs + iF) += quadrule.w * exact_solution(tps, quadrule.vector());
			  }
		  }
		  // Take average
		  DirBC(n_cell_dofs + iF) /= mesh->b_edge(ibF)->measure();
    }
	}
	Xhprev.asVectorXd() = DirBC;

	// ---- NEWTON ITERATONS ------ //
	Eigen::VectorXd RES = residual_scheme(dt, Xhprev); 		// Scheme is F(X)=C, then RES = C - F(u_prev)
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
			Eigen::VectorXd XTprev = Xhprev.restr(iT);
			Eigen::VectorXd ZetaprimeXTprev = apply_nonlinearity(XTprev, "der");
			Eigen::MatrixXd MatZetaprimeXTprev = Eigen::MatrixXd::Zero(1+n_edgesT, 1+n_edgesT);

		  for (size_t i = 0; i < 1 + n_edgesT; i++){
				MatZetaprimeXTprev(i, i) = ZetaprimeXTprev(i);
			}

////// Can be simplified afterwards (using 1x1 matrices here)
			Eigen::MatrixXd MatT = dt * DiffT[iT] * MatZetaprimeXTprev + MassT[iT];
			Eigen::MatrixXd ATT = MatT.topLeftCorner(1, 1);
			Eigen::MatrixXd ATF = MatT.topRightCorner(1, n_edgesT);
			Eigen::MatrixXd AFT = MatT.bottomLeftCorner(n_edgesT, 1);
			Eigen::MatrixXd AFF = MatT.bottomRightCorner(n_edgesT, n_edgesT);

			Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
			invATT.compute(ATT);
				
			Eigen::MatrixXd invATT_ATF = invATT.solve(ATF);
			Eigen::VectorXd RES_cell = RES.segment(iT, 1);
			Eigen::VectorXd invATT_RES_cell = invATT.solve(RES_cell);

			// Local matrix and right-hand side on the face unknowns
			Eigen::MatrixXd MatF = Eigen::MatrixXd::Zero(n_edgesT, n_edgesT);
			Eigen::VectorXd bF = Eigen::VectorXd::Zero(n_edgesT);
			MatF = AFF - AFT * invATT_ATF;
			bF = - AFT * invATT_RES_cell;		// only additional RHS coming from static condensation, beyond RES
							
			// Assemble static condensation operator
			ScRHS.segment(iT, 1) = invATT_RES_cell;

///// TO SIMPLIFY
			for (size_t i = 0; i < 1; i++){
				size_t iGlobal = iT + i;
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
		Xh.asVectorXd() = DirBC;
    Xh.asVectorXd().head(n_cell_dofs + n_edge_unknowns) = relax*dX.head(n_cell_dofs + n_edge_unknowns) + Xhprev.asVectorXd().head(n_cell_dofs + n_edge_unknowns);

    // Compute new residual. We might not accept new Xh but start over reducing relax. Otherwise, accept Xh and increase relax
    // to a max of 1
		Eigen::VectorXd RES_temp = residual_scheme(dt, Xh);
    double residual_temp = RES_temp.norm();

    if (iter > 1 && residual_temp > 10*residual[iter-1]){
      // We don't update the residual and Xh, but we will do another iteration with relax reduced
      relax = relax/2.0;
      iter--;
    } else {
      // Increase relax
      relax = std::min(1.0, 1.2*relax);
  		// Store Xh in Xhprev, compute new residual
  		Xhprev.asVectorXd() = Xh.asVectorXd();
  		RES = RES_temp;
  		residual[iter] = residual_temp;
    }

		m_output << "     ...Iteration: " << iter << "; residual = " << residual[iter] << "; relax = "<< relax << "\n";

	} //-------- END NEWTON ITERATIONS -----------//
  m_nb_newton = iter;

  return Xh;
}

//****************************************
//		apply nonlinearity zeta to a vector
//****************************************

Eigen::VectorXd HMM_StefanPME_Transient::apply_nonlinearity(const Eigen::VectorXd& Y, const std::string type) const {
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

	const Mesh* mesh = hmm.get_mesh();
	// We check if Y is a local vector or a global one, this determines the number of cell and edge unknowns
	size_t n_cell_dofs = 0;
	if (size_t(Y.size()) == mesh->n_cells() + mesh->n_edges()){
		n_cell_dofs = mesh->n_cells();
	}else{
		n_cell_dofs = 1;
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

// Overloaded version: apply to the values of UVector, don't change the other parameters
UVector HMM_StefanPME_Transient::apply_nonlinearity(const UVector& Y, const std::string type) const {

  UVector ZetaY = Y;
  ZetaY.asVectorXd() = apply_nonlinearity(Y.asVectorXd(), type);
  return ZetaY;
}


//******************************** 
//		local diffusion matrix 
//********************************

Eigen::MatrixXd HMM_StefanPME_Transient::diffusion_operator(const size_t iT) const {

  const auto mesh = hmm.get_mesh();
  size_t dim = mesh->dim();
	Cell* T = mesh->cell(iT);
  double mT = T->measure();
	const size_t n_edgesT = T->n_edges();

  size_t local_dofs = 1 + n_edgesT;
  Eigen::MatrixXd StiffT = Eigen::MatrixXd::Zero(local_dofs, local_dofs);
	
	// Diffusion in the cell.
  Eigen::Vector2d xT = T->center_mass();
  Eigen::Matrix2d kappaT = kappa(xT.x(), xT.y(), T);

  // Matrix of \nabla_K
  Eigen::MatrixXd nablaK = Eigen::MatrixXd::Zero(dim, local_dofs);
  for (size_t iE=0; iE < n_edgesT; iE++){
    Edge* E = T->edge(iE);
    double mE = E->measure();
    nablaK.block(0, iE+1, dim, 1) = mE * T->edge_normal(iE);
  }
  nablaK /= mT;

  // Contribution per diamond
  for (size_t iE=0; iE < n_edgesT; iE++){
    Edge* E = T->edge(iE);
    Eigen::Vector2d xE = E->center_mass();
    Eigen::Vector2d nTE = T->edge_normal(iE);
    double mE = E->measure();
    // Distance center to edge and measure of diamond
    double dTE =  (xE-xT).dot(nTE);
    double mTE = mE * dTE / dim;

    // stabilisation (without scaling or normal: u_E-u_K-nablaK u . (xE-xK))
    Eigen::MatrixXd stabE = Eigen::MatrixXd::Zero(1, local_dofs);
    stabE(0) = -1.0;
    stabE(1+iE) = 1.0;
    stabE -= (xE-xT).transpose() * nablaK;

    // Complete gradient with stabilisation
    Eigen::MatrixXd GRAD = nablaK + (pow(dim, 0.5)/dTE) * nTE * stabE;

    // Contribution to stiffness
    StiffT += mTE * GRAD.transpose() * kappaT * GRAD;
  }

	return StiffT; 
}


//******************************** 
//		local load term 
//********************************

Eigen::VectorXd HMM_StefanPME_Transient::load_operator(const double tps, const size_t iT) const {

  const auto mesh = hmm.get_mesh();
	Cell* T = mesh->cell(iT);
	size_t n_edgesT = T->n_edges();
	size_t local_dofs = 1 + n_edgesT;

	// Beware, MassT[iT] must have been created
	if (MassT.size() < iT || size_t(MassT[iT].rows()) != local_dofs){
		m_output << "Called load_operator without creating MassT[iT]\n";
		exit(EXIT_FAILURE);
	}
	Eigen::VectorXd fvec = Eigen::VectorXd::Zero(1 + n_edgesT);
	for (size_t i = 0; i < local_dofs; i++){
		VectorRd node;
		if (i == 0){
			node = T->center_mass();
		}else{
			node = T->edge(i-1)->center_mass();
		}
		fvec(i) = source(tps, node, T);
	}
  // Contribution of source to loading term
  Eigen::VectorXd load = MassT[iT] * fvec;

  // Adding Neumann boundary conditions
  // The boundary conditions coming from Neumann edges are computed based on the piecewise constant
  // discrete trace: if g is the value of the outer flux,
  //      \int g v -> \int g Tv
  //  where (Tv)_e = average of v on e, for e Neumann edge.
  // Only edge-based basis functions therefore have a contribution, which will simply be \int_e g.
  if (T->is_boundary()){
    for (size_t ilF = 0; ilF < T->n_edges(); ilF++) {
      Edge* edge = T->edge(ilF);
      if (m_BC.type(*edge)=="neu"){
				if (edge->is_boundary()){
				  const auto& nTF = T->edge_normal(ilF);
					QuadratureRule quadF = generate_quadrature_rule(*edge, 5);
					std::function<double(VectorRd)> Kgrad_n = [&](VectorRd p){
								return zeta(exact_solution(tps, p), "der") * nTF.dot(kappa(p.x(),p.y(),T) * grad_exact_solution(tps, p, T));
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

Eigen::VectorXd HMM_StefanPME_Transient::residual_scheme(const double dt, const UVector& Xh) const{

	const Mesh* mesh = hmm.get_mesh();	
	Eigen::VectorXd RES = Eigen::VectorXd::Zero(mesh->n_cells() + mesh->n_edges());

	// Compute the residual: C - (dt * DIFF * ZetaXh + MASS * Xh), where DIFF is the diffusion matrix, MASS is the mass matrix and C is the source term
	size_t n_cell_dofs = mesh->n_cells();
	for (size_t iT = 0; iT < mesh->n_cells(); iT++){
		Cell* T =  mesh->cell(iT);
		size_t n_edgesT = T->n_edges();
		// Local unknown in the cell, and zeta of these unknowns
		Eigen::VectorXd XT = Xh.restr(iT);
		Eigen::VectorXd ZetaXT = apply_nonlinearity(XT, "fct");

		Eigen::VectorXd REST = RightHandSideT[iT] - (dt * DiffT[iT] * ZetaXT + MassT[iT] * XT);
		RES(iT) += REST(0);
		for (size_t ilE = 0; ilE < n_edgesT; ilE++){
			size_t iE = T->edge(ilE)->global_index();
			RES(n_cell_dofs + iE) += REST(1+ilE);
		}
	}

	// No residual on Dirichlet edges (do not correspond to equations of the system)
	RES.tail(m_BC.n_dir_edges()) = Eigen::VectorXd::Zero(m_BC.n_dir_edges());

	return RES;
}

//*******************************************************
// Norms: L2 mass lumped, Energy (pure diffusion)
//*******************************************************

double HMM_StefanPME_Transient::L2_MassLumped(const UVector& Xh) const {
  const Mesh* mesh = hmm.get_mesh();
	double value = 0.0;

  for (size_t iT = 0; iT < mesh-> n_cells(); iT++) {
		Eigen::VectorXd XTF = Xh.restr(iT);
		value += XTF.transpose() * MassT[iT] * XTF;
  }

  return sqrt(value);
}


double HMM_StefanPME_Transient::Lp_MassLumped(const UVector& Xh, const double p) const {
  const Mesh* mesh = hmm.get_mesh();
	double value = 0.0;

  for (size_t iT = 0; iT < mesh-> n_cells(); iT++) {
		Eigen::VectorXd XTF = Xh.restr(iT);
    // Coefficient-wise power p
    Eigen::VectorXd XTF_powerp = (XTF.array().abs()).pow(p);
    // Sum local contributions; this assumes that MassT is diagonal
		value += (MassT[iT] * XTF_powerp.matrix() ).sum();
  }

  return std::pow(value, 1.0/p);
}

double HMM_StefanPME_Transient::EnergyNorm(const UVector& Xh) const {
  const Mesh* mesh = hmm.get_mesh();
	double value = 0.0;

  for (size_t iT = 0; iT < mesh-> n_cells(); iT++) {
		Eigen::VectorXd XTF = Xh.restr(iT);
		value += XTF.transpose() * DiffT[iT] * XTF;
  }

  return sqrt(value);
}


//@}
} // end of namespace HArDCore2D

#endif //_HMM_STEFANPME_TRANSIENT_HPP
