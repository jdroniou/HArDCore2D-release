// Derived class from HybridCore, provides non-conforming basis functions on generic polygonal cells
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef LEPNCCORE_HPP
#define LEPNCCORE_HPP

#include <hybridcore.hpp>
#include <quadraturerule.hpp>

/*!	
* @defgroup LEPNC 
* @brief Locally Enriched Polytopal Non-Conforming method
*/

namespace HArDCore2D {


/*!
*	\addtogroup LEPNC
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/** The LEPNCCore class provides basis functions for non-conforming schemes on generic polygonal meshes
*
*
**/

class LEPNCCore : public HybridCore {

public:
	///@brief Class constructor: initialise the LEPNCCore class, and creates non-conforming basis functions (and gradients)
  LEPNCCore(
						const Mesh* mesh_ptr, ///< A pointer to the loaded mesh 
            const size_t K, ///< The degree of the edge polynomials 
						const size_t L, ///< The degree of the cell polynomials 
            std::ostream & output = std::cout ///< Optional argument for specifying outputs of messages.
             ); 

  /// Types for basis functions
  using cell_basis_type = std::function<double(double, double)>;    ///< type for cell basis
  using cell_gradient_type = std::function<VectorRd(double, double)>;   ///< type for gradients of cell basis


  /// Return a reference to the i'th non-conforming basis function of the cell iT
	/** For i=DimPoly<Cell>(1) to nedge-1, nc_basis(iT,i) is the basis function associated with edge i (that is,
		product of all distances to all other edges, scaled to have an integral over edge i equal to 1).
		For i=0, ..., DimPoly<Cell>(1)-1 the basis function corresponds to 1, x, y, to which we subtract
		a linear combination of the basis functions corresponding to the edge in order to obtain basis
		functions with a zero integral on each edge **/
  inline const cell_basis_type& nc_basis(
      size_t iT, /**< The global region number of the cell */
      size_t i   /**< The index of the desired basis function */
      ) const;
	
  /// Return a reference to the gradient of the i'th non-conforming basis function of the cell iT
  /** The gradient functions are indexed the same as the basis functions. */
  inline const cell_gradient_type& nc_gradient(
      size_t iT, /**< The global region number of the cell */
      size_t i   /**< The index of the desired basis function */
      ) const;

  /// Return the i'th nodal point in cell iT (for mass lumping)
  inline const VectorRd ml_node(
      size_t iT, /**< The global region number of the cell */
      size_t i   /**< The index of the desired nodal point */
	     ) const;

	/// Computes non-conforming basis functions at the given quadrature nodes
	const std::vector<Eigen::ArrayXd> nc_basis_quad(
			const size_t iT, 		///<	global index of the cell/edge
			const QuadratureRule quad 	///< quadrature nodes and weights on the cell/edge
	) const;	///< @returns nc_phi_quad[i] = array listing the nbq (=nb of quadrature nodes) values of the nonconforming \f$\phi^{nc}_i\f$ basis function at the quadrature nodes

	/// Compute \f$(\nabla \phi_i^{nc})_{i\in I}\f$ at the quadrature nodes, where \f$(\phi_i^{nc})_{i\in I}\f$ are the cell basis functions
	const std::vector<Eigen::ArrayXXd> grad_nc_basis_quad(
			const size_t iT, 											///< global index of the cell
			const QuadratureRule quad 				///< quadrature rules in the cell
	) const; ///< @returns Dnc_phi_quad[i]: array of size 2*nbq (where nbq=nb of quadrature nodes), with each column being \f$\nabla \phi_i^{nc}\f$ at the corresponding quadrature node 

  /// Create the matrix of L2 products of two families (f_i) and (g_j) of functions
  /// (this is not really a Gram matrix, unless the two families are the same)
  Eigen::MatrixXd gram_matrix(
      const std::vector<Eigen::ArrayXd>& f_quad, ///< Values of functions (f1,f2,...) at the quadrature nodes 
      const std::vector<Eigen::ArrayXd>& g_quad, ///< Values of functions (g1,g2,...) at the quadrature nodes 
      const size_t& nrows, ///< Number of rows of the matrix - typically number of functions f_i (but could be less) 
      const size_t& ncols, ///< Number of columns of the matrix - typically number of functions g_j (but could be less)
      const QuadratureRule& quad,   ///< Quadrature nodes for integration 
      const bool& sym,    ///< True if the matrix is pseudo-symmetric (that is, #f<=#g and f_i=g_i if i<=#f) 
      std::vector<double> L2weight = {}   ///< Optional weight for the L2 product. If provided, should be a std::vector<double> of the weight at the quadrature nodes
  ) const;  /**< @returns The matrix \f$(\int f_i g_j)_{i=1\ldots nrows; j=1\ldots ncols}\f$ */

  /// Overloaded version of the previous one for vector-valued functions: the functions (F_i) and (G_j) are vector-valued functions
  Eigen::MatrixXd gram_matrix(
      const std::vector<Eigen::ArrayXXd>& F_quad,    ///< Values of functions (F1,F2,...) at the quadrature nodes 
      const std::vector<Eigen::ArrayXXd>& G_quad,    ///< Values of functions (G1,G2,...) at the quadrature nodes 
      const size_t& nrows,    ///< Number of rows of the matrix - typically number of functions F_i (but could be less) 
      const size_t& ncols,    ///< Number of rows of the matrix - typically number of functions G_j (but could be less) 
      const QuadratureRule& quad,    ///< Quadrature nodes for integration 
      const bool& sym,    ///< True if the matrix is pseudo-symmetric (that is, #F<=#G and F_i=G_i if i<=#F) 
      std::vector<Eigen::Matrix2d> L2Weight = {}  ///< Optional weight for the L2 product. If provided, should be a std::vector<Eigen::Matrix2d> of the weight at the quadrature nodes
    ) const;  /**< @returns The matrix \f$(\int F_i \cdot G_j)_{i=1\ldots nrows; j=1\ldots ncols}\f$ */

	/// Interpolates a continuous function on the degrees of freedom, using the moments on the basis functions. The first ones are the cell DOFs (DimPoly<Cell>(1) for each cell), the last ones are the edge DOFs (one for each edge)
	template<typename Function>
	Eigen::VectorXd nc_interpolate_moments(const Function& f, size_t doe) const; /// @returns XTF = vector of coefficients on the basis functions; the first "DimPoly<Cell>(1)*nb cells" correspond to the first three basis functions in each cell (1, x, y with adjustments for the averages on the edges and mass-lumping points), and the last "nb edges" to the edge basis functions.

	/// Interpolates a continuous function on the degrees of freedom, using the moments on the basis functions associated to the edges and the nodal values (corresponding to mass-lumping) on the cell basis functions. The first ones are the cell DOFs (DimPoly<Cell>(1) for each cell), the last ones are the edge DOFs (one for each edge)
	template<typename Function>
	Eigen::VectorXd nc_interpolate_ml(const Function& f, size_t doe) const; /// @returns XTF = vector of coefficients on the basis functions; the first "DimPoly<Cell>(1)*nb cells" correspond to the first three basis functions in each cell (1, x, y with adjustments for the averages on the edges and mass-lumping points), and the last "nb edges" to the edge basis functions.

	Eigen::VectorXd nc_restr(const Eigen::VectorXd &Xh, size_t iT) const; ///< Extract from a global vector Xh of unknowns the non-conforming unknowns corresponding to cell iT

	double nc_L2norm(const Eigen::VectorXd &Xh) const; ///< Compute L2 norm of a discrete function (given by coefficients on the basis functions)

	double nc_L2norm_ml(const Eigen::VectorXd &Xh) const; ///< Compute L2 norm of the mass-lumped discrete function (given by coefficients on the basis functions)

	double nc_H1norm(const Eigen::VectorXd &Xh) const; ///< Compute broken H1 norm of a discrete function (given by coefficients on the basis functions)

	/// Evaluates a non-conforming discrete function in the cell iT at point (x,y)
	double nc_evaluate_in_cell(const Eigen::VectorXd XTF, size_t iT, double x, double y) const; 

	/// From a non-conforming discrete function, computes a vector of values at the vertices of the mesh
	Eigen::VectorXd nc_VertexValues(
		const Eigen::VectorXd Xh, 		///< non-conforming function (coefficients on basis)
		const double weight = 0 		///< weight put on the edge values when evaluating the functions at the vertices
		);


private:
	std::tuple<std::vector<cell_basis_type>, 
				std::vector<cell_gradient_type>, 
				Eigen::MatrixXd>
  create_nc_basis(size_t iT) const;			///< creates the non-conforming basis for cell iT.

	std::vector<VectorRd> create_ml_nodes(size_t iT) const;		///< creates the nodes for mass-lumping in each cell

	// Collections of non-conforming basis functions for each cell
	std::vector<std::vector<cell_basis_type> > _nc_bases;
	std::vector<std::vector<cell_gradient_type> > _nc_gradients;
	// Basis functions to represent bubble basis functions based on monomials/edge-based functions
  std::vector<Eigen::MatrixXd> _M_basis;

	// Collection of mass-lumping nodes
	std::vector<std::vector<VectorRd>> _ml_nodes;

	// Sign function
	double sign(const double s) const;

  // Output stream
  std::ostream & m_output;

};

// ************************************************************************
//                            Implementation
// ************************************************************************

//-----------------------------------------------------
//------ Create class
//-----------------------------------------------------

LEPNCCore::LEPNCCore(const Mesh* mesh_ptr,
					const size_t K,
					const size_t L,
          std::ostream & output):
						HybridCore (mesh_ptr, std::min(size_t(1),K), L, false, output),
						_nc_bases(0),
						_nc_gradients(0),
						_M_basis(0),
						_ml_nodes(0),
            m_output(output) {
			size_t ncells = get_mesh()->n_cells();
			// Resize vectors
			_nc_bases.reserve(ncells);
			_nc_gradients.reserve(ncells);
			_ml_nodes.reserve(ncells);
			_M_basis.reserve(ncells);
			// Create nodes for mass-lumping and initialise non-conforming basis functions on cells
			for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++) {
				std::vector<VectorRd> ml_nodes = create_ml_nodes(iT);
				_ml_nodes.push_back(ml_nodes);
				std::tuple<std::vector<cell_basis_type>, std::vector<cell_gradient_type>, Eigen::MatrixXd> nc_basis = create_nc_basis(iT);
				_nc_bases.push_back(std::get<0>(nc_basis));
				_nc_gradients.push_back(std::get<1>(nc_basis));
				_M_basis.push_back(std::get<2>(nc_basis));
			}

}

//-----------------------------------------------------
//------ Create mass-lumping nodes
//-----------------------------------------------------
std::vector<VectorRd> LEPNCCore::create_ml_nodes(size_t iT) const {
	const Mesh* mesh = get_mesh();
	Cell* cell = mesh->cell(iT);
	size_t nedges = cell->n_edges();
	std::vector<VectorRd> ml_cell(mesh->dim()+1,VectorRd::Zero());

	// We look for the triangle of vertices of the cell that maximises the area
	double area_max = 0;
	std::vector<VectorRd> v_cur(mesh->dim()+1,VectorRd::Zero());

	for (size_t i=0; i < nedges; i++){
		v_cur[0] = cell->vertex(i)->coords();
		for (size_t j=i+1; j < nedges; j++){
			v_cur[1] = cell->vertex(j)->coords();
			for (size_t k=j+1; k < nedges; k++){
				v_cur[2] = cell->vertex(k)->coords();
				Eigen::MatrixXd matv = Eigen::MatrixXd::Zero(mesh->dim(), mesh->dim());
				for (size_t ll=0; ll < mesh->dim(); ll++){
					matv.col(ll) = v_cur[ll+1]-v_cur[0];
				}
				double area_cur = std::abs(matv.determinant());
				if (area_cur > area_max){
					for (size_t z=0; z < mesh->dim()+1; z++){
						ml_cell[z] = v_cur[z];
					}
					area_max = area_cur;
				}
			}
		}
	}

	return ml_cell;
}

//-----------------------------------------------------
//------ Create and get basis functions and gradients
//-----------------------------------------------------

std::tuple<std::vector<LEPNCCore::cell_basis_type>,
			std::vector<LEPNCCore::cell_gradient_type>,
			Eigen::MatrixXd>
	LEPNCCore::create_nc_basis(size_t iT) const {

		Cell* cell = get_mesh()->cell(iT);
		size_t nedgesT = cell->n_edges();
		std::vector<cell_basis_type> nc_basis(DimPoly<Cell>(1)+nedgesT,[](double x, double y)->double {return 0;});
		std::vector<cell_gradient_type> nc_gradient(nedgesT+DimPoly<Cell>(1),[](double x, double y)->VectorRd {return VectorRd::Zero();});

		// Computation of basis functions: we start by the functions associated with the edges,
		// that we put at the end of the lists
		// For each edge, the function is constructed the following way: the triangle with the center of the cell
		// is created, and the function is the product of the positive part of the signed distance to the
		// edges of this triangle internal to the cell. It's therefore positive in the triangle and zero outside
		// It is then scaled to have an average of one on the edge
		VectorRd xT = cell->center_mass();
		for (size_t ilE = 0; ilE < nedgesT; ilE++){
			Edge* edge = cell->edge(ilE);
			// Compute normals to the internal edges of the triangle 
			Eigen::Matrix2d Rot;
			Rot << 0, 1, -1, 0;
			std::vector<VectorRd> Nint_edge(2, VectorRd::Zero());
			for (size_t ilV = 0; ilV < 2; ilV++){
				// Non-normalised, non-oriented vector orthogonal to edge
				VectorRd temp_norm = Rot * (edge->vertex(ilV)->coords() - xT);
				// Orient and normalise to be external to the triangle
				double orientation = sign( (xT - edge->center_mass()).dot(temp_norm) );
				Nint_edge[ilV] = orientation * temp_norm.normalized();
			}

			// Unscaled basis function
			double edge_length = edge->measure();
			cell_basis_type unscaled_phi = [Nint_edge, xT, edge_length](double x, double y)->double {
				VectorRd posX = VectorRd(x,y);
				double val = 1;

				for (size_t ilV = 0; ilV < 2; ilV++){
					double signed_dist = (xT - posX).dot(Nint_edge[ilV]);
					val *= std::max(double(0), signed_dist) / edge_length;
				}
				return val;
			};

			// Compute factor to scale basis function and get an average of one on its associated edge.
			double scaling_factor_phi = 0;
			QuadratureRule quadE = generate_quadrature_rule(*edge, 2);
			for (QuadratureNode quadrule : quadE){
				scaling_factor_phi += quadrule.w * unscaled_phi(quadrule.x, quadrule.y);
			}
			if (std::abs(scaling_factor_phi) < 1e-10){
				m_output << "Scaling factor of non-conforming basis function too small: " << std::abs(scaling_factor_phi) << "\n";
				exit(EXIT_FAILURE);
			}

			// Scale basis function
			scaling_factor_phi /= edge->measure();
			cell_basis_type phi = [unscaled_phi, scaling_factor_phi](double x, double y)->double {
				return unscaled_phi(x, y) / scaling_factor_phi;
			};

			// Gradient of basis function
			cell_gradient_type dphi = [Nint_edge, xT, scaling_factor_phi, edge_length](double x, double y)->VectorRd {
				VectorRd val = VectorRd::Zero();

				VectorRd posX = VectorRd(x,y);
				for (size_t ilV = 0; ilV < 2; ilV++){
					double coef_ilV = 1;
					if ( (xT - posX).dot(Nint_edge[ilV]) < 0){
						coef_ilV = 0;
					}
					for (size_t ilVp = 0; ilVp < 2; ilVp++){
						if (ilVp != ilV){
							double signed_dist_ilVp = (xT - posX).dot(Nint_edge[ilVp]);
							coef_ilV *= std::max(double(0), signed_dist_ilVp) / edge_length;
						}
					}
					val -= coef_ilV * Nint_edge[ilV] / edge_length;
				}
				return val / scaling_factor_phi;
			};

			// Store basis function and gradient
			nc_basis[ilE+DimPoly<Cell>(1)] = std::move(phi);
			nc_gradient[ilE+DimPoly<Cell>(1)] = std::move(dphi);
		}

		// We then put the basis functions corresponding to the first three HybridCore cell basis functions (basis of P^1)
		//	to which we subtract a combination of the previous basis functions in order to obtain functions 
		// 	with zero averages on all edges. In a second step, we also modify this basis function to make sure
		//	it's nodal with respect to the selected points in _ml_nodes[iT] 
		//  These transformations of basis are done by computing the matrix to pass from the 
		//  basis functions (basis P^1, edges basis functions) to the new basis functions.
		Eigen::MatrixXd M_basis = Eigen::MatrixXd::Zero(DimPoly<Cell>(1), DimPoly<Cell>(1) + nedgesT);
		M_basis.topLeftCorner(DimPoly<Cell>(1), DimPoly<Cell>(1)) = Eigen::MatrixXd::Identity(DimPoly<Cell>(1), DimPoly<Cell>(1));
		for (size_t i = 0; i < DimPoly<Cell>(1); i++){
			// Basis function 1, x or y, and gradient
//////			cell_basis_type phi0 = cell_basis(iT, i);
			cell_basis_type phi0 = [&](double x, double y)->double {
                return CellBasis(iT).function(i,VectorRd(x,y));
          };
			// Averages of basis function over edges, that we store in matrix M_basis
			// At this point, M_basis transforms the basis (1,x,y) into a basis of functions
			// that have zero averages on the edges
			for (size_t ilE = 0; ilE < nedgesT; ilE++){
        Edge* edge = cell->edge(ilE);
				QuadratureRule quadE = generate_quadrature_rule(*edge, 1);
				for (QuadratureNode quadrule : quadE){
					M_basis(i, DimPoly<Cell>(1)+ilE) -= quadrule.w * phi0(quadrule.x, quadrule.y);
				}
				M_basis(i, DimPoly<Cell>(1)+ilE) /= edge->measure();
			}
		}

		// We then modify M_basis to ensure that the transformed basis function of P^1 are nodal with
		// respect to _ml_nodes[iT]. 
		// M_nodal_values is (phi_i(x_k))_ik where phi_i are the original basis functions and x_k 
		// are the nodal points
		// Note that the edge basis functions all have zero value at the nodal points
		const Mesh* mesh = get_mesh();
		std::vector<VectorRd> nodal_points = _ml_nodes[iT];
		Eigen::MatrixXd M_nodal_values = Eigen::MatrixXd::Zero(DimPoly<Cell>(1), mesh->dim()+1);
		for (size_t i=0; i < DimPoly<Cell>(1); i++){
			for (size_t k=0; k < mesh->dim() + 1; k++){
//////				cell_basis_type phi = cell_basis(iT, i);
//////				M_nodal_values(i, k) = phi(nodal_points[k].x(), nodal_points[k].y() );
				M_nodal_values(i, k) = CellBasis(iT).function(i, nodal_points[k]);
			}
		}
		// Modification M_basis to get nodal basis functions
		M_basis = (M_nodal_values.inverse() ) * M_basis;
		
		// Store new basis functions replacing (1,x,y)
		for (size_t i=0; i < DimPoly<Cell>(1); i++){
			Eigen::VectorXd M_basis_rowi = M_basis.row(i);
			cell_basis_type phi = [iT, nedgesT, nc_basis, M_basis_rowi, this](double x, double y)->double {
				double val = 0;
				for (size_t j = 0; j < DimPoly<Cell>(1)+nedgesT; j++){
					if (j < DimPoly<Cell>(1)){
//////						cell_basis_type phi0 = cell_basis(iT, j);
						val += M_basis_rowi(j) * CellBasis(iT).function(j, VectorRd(x, y));
					}else{
						cell_basis_type phi0 = nc_basis[j];
						val += M_basis_rowi(j) * phi0(x, y);
					}
				}
				return val;			
			};
			cell_gradient_type dphi = [iT, nedgesT, nc_gradient, M_basis_rowi, this](double x, double y)->VectorRd {
				VectorRd val = VectorRd::Zero();
				for (size_t j = 0; j < DimPoly<Cell>(1)+nedgesT; j++){
					if (j < DimPoly<Cell>(1)){
//////						cell_gradient_type dphi0 = cell_gradient(iT, j);
						val += M_basis_rowi(j) * CellBasis(iT).gradient(j, VectorRd(x, y));
					}else{
						cell_gradient_type dphi0 = nc_gradient[j];
						val += M_basis_rowi(j) * dphi0(x, y);
					}
				}
				return val;
			};


			// Store basis function and gradient
			nc_basis[i] = std::move(phi);
			nc_gradient[i] = std::move(dphi);
		}

		return std::make_tuple(std::move(nc_basis), std::move(nc_gradient), std::move(M_basis));
}


// Get basis functions, nodal points

inline const LEPNCCore::cell_basis_type& LEPNCCore::nc_basis(size_t iT, size_t i) const {
  assert(iT < this->get_mesh()->n_cells());
  assert(i < _nc_bases[iT].size());
  return _nc_bases[iT][i];
}

inline const LEPNCCore::cell_gradient_type& LEPNCCore::nc_gradient(size_t iT, size_t i) const {
  assert(iT < this->get_mesh()->n_cells());
  assert(i < _nc_gradients[iT].size());
  return _nc_gradients[iT][i];
}

inline const VectorRd LEPNCCore::ml_node(size_t iT, size_t i) const {
  assert(iT < this->get_mesh()->n_cells());
  assert(i < get_mesh()->dim() + 1);
  return _ml_nodes[iT][i];
}

// ----------------------------------------------------------------
// ------- Non-conforming basis functions and gradients at quadrature nodes
// ----------------------------------------------------------------

const std::vector<Eigen::ArrayXd> LEPNCCore::nc_basis_quad(const size_t iT, const QuadratureRule quad) const {

	size_t nbq = quad.size();
	const Mesh* mesh = get_mesh();
	size_t nbasis = DimPoly<Cell>(1) + mesh->cell(iT)->n_edges();
	std::vector<Eigen::ArrayXd> nc_phi_quad(nbasis, Eigen::ArrayXd::Zero(nbq));	

	// Compute first on monomial basis functions and edge-based basis functions
	std::vector<Eigen::ArrayXd> nc_phi_quad_tmp(nbasis, Eigen::ArrayXd::Zero(nbq));	
	for (size_t i = 0; i < nbasis; i++){
		if (i < DimPoly<Cell>(1)){
//////			const auto &monomial = cell_basis(iT, i);

			for (size_t iqn = 0; iqn < nbq; iqn++){
//////				nc_phi_quad_tmp[i](iqn) = monomial(quad[iqn].x, quad[iqn].y);
				nc_phi_quad_tmp[i](iqn) = CellBasis(iT).function(i, quad[iqn].vector());
			}
		}else{
			const auto &nc_phi_i = nc_basis(iT, i);

			for (size_t iqn = 0; iqn < nbq; iqn++){
				nc_phi_quad_tmp[i](iqn) = nc_phi_i(quad[iqn].x, quad[iqn].y);
			}
		}
	}

	// Adjust for bubble cell basis functions
	for (size_t i = DimPoly<Cell>(1); i < nbasis; i++){
		nc_phi_quad[i] = nc_phi_quad_tmp[i];
	}
	for (size_t i = 0; i < DimPoly<Cell>(1); i++){
		for (size_t j = 0; j < nbasis; j++){
			nc_phi_quad[i] += _M_basis[iT](i,j) * nc_phi_quad_tmp[j];
		}
	}

	return nc_phi_quad;
}

const std::vector<Eigen::ArrayXXd> LEPNCCore::grad_nc_basis_quad(const size_t iT, const QuadratureRule quad) const {

	size_t nbq = quad.size();
	const Mesh* mesh = get_mesh();
	size_t nbasis = DimPoly<Cell>(1) + mesh->cell(iT)->n_edges();
	std::vector<Eigen::ArrayXXd> Dnc_phi_quad(nbasis, Eigen::ArrayXXd::Zero(get_mesh()->dim(), nbq));	

	// Compute first on monomial basis functions and edge-based basis functions
	std::vector<Eigen::ArrayXXd> Dnc_phi_quad_tmp(nbasis, Eigen::ArrayXXd::Zero(get_mesh()->dim(), nbq));	
	for (size_t i = 0; i < nbasis; i++){
		if (i < DimPoly<Cell>(1)){
//////			const auto &grad_mono = cell_gradient(iT, i);

			for (size_t iqn = 0; iqn < nbq; iqn++){
//////				Dnc_phi_quad_tmp[i].col(iqn) = grad_mono(quad[iqn].x, quad[iqn].y);
				Dnc_phi_quad_tmp[i].col(iqn) = CellBasis(iT).gradient(i, quad[iqn].vector());
			}
		}else{
			const auto &Dnc_phi_i = nc_gradient(iT, i);

			for (size_t iqn = 0; iqn < nbq; iqn++){
				Dnc_phi_quad_tmp[i].col(iqn) = Dnc_phi_i(quad[iqn].x, quad[iqn].y);
			}
		}
	}

	// Adjust for bubble cell basis functions
	for (size_t i = DimPoly<Cell>(1); i < nbasis; i++){
		Dnc_phi_quad[i] = Dnc_phi_quad_tmp[i];
	}
	for (size_t i = 0; i < DimPoly<Cell>(1); i++){
		for (size_t j = 0; j < nbasis; j++){
			Dnc_phi_quad[i] += _M_basis[iT](i,j) * Dnc_phi_quad_tmp[j];
		}
	}

	return Dnc_phi_quad;
}

//-----------------------------------------------------
//------ Old versions of gram matrices, to replace later perhaps
//-----------------------------------------------------


  /// Create the matrix of L2 products of two families (f_i) and (g_j) of functions
  /// (this is not really a Gram matrix, unless the two families are the same)
  Eigen::MatrixXd LEPNCCore::gram_matrix(const std::vector<Eigen::ArrayXd>& f_quad, const std::vector<Eigen::ArrayXd>& g_quad, const size_t& nrows, const size_t& ncols, const QuadratureRule& quad, const bool& sym, std::vector<double> L2weight) const {
    Eigen::MatrixXd GGM = Eigen::MatrixXd::Zero(nrows, ncols);

    size_t nbq = quad.size();

    // Recast product of quadrature and L2weight into an Eigen::ArrayXd
    Eigen::ArrayXd quad_L2_weights = Eigen::ArrayXd::Zero(nbq);
    if (L2weight.size() == 0){
      for (size_t iqn = 0; iqn < nbq; iqn++){
        quad_L2_weights(iqn) = quad[iqn].w;
      }
    }else{
      for (size_t iqn = 0; iqn < nbq; iqn++){
        quad_L2_weights(iqn) = quad[iqn].w * L2weight[iqn];
      }
    }

    for (size_t i = 0; i < nrows; i++){
      size_t jcut = 0;
      if (sym) jcut = i;
      for (size_t j = 0; j < jcut; j++){
          GGM(i, j) = GGM(j, i);
      }
      for (size_t j = jcut; j < ncols; j++){
        // Integrate f_i * g_j
        // The products here are component-wise since the terms are Eigen::ArrayXd
        GGM(i, j) = (quad_L2_weights * f_quad[i] * g_quad[j]).sum();
      }
    }

    return GGM;
  }

  /// Overloaded version of the previous one for vector-valued functions: the functions (F_i) and (G_j) are vector-valued functions
  Eigen::MatrixXd LEPNCCore::gram_matrix(const std::vector<Eigen::ArrayXXd>& F_quad, const std::vector<Eigen::ArrayXXd>& G_quad, const size_t& nrows, const size_t& ncols, const QuadratureRule& quad, const bool& sym, std::vector<Eigen::Matrix2d> L2Weight) const {
    Eigen::MatrixXd GSM = Eigen::MatrixXd::Zero(nrows, ncols);

    size_t nbq = quad.size();
    for (size_t i = 0; i < nrows; i++){
      size_t jcut = 0;
      if (sym) jcut = i;
      for (size_t j = 0; j < jcut; j++){
          GSM(i, j) = GSM(j, i);
      }
      // Multiply F_i by quadrature weights and matrix L2Weight, if required
      Eigen::ArrayXXd WeightsTimesF_quad = Eigen::ArrayXXd::Zero(get_mesh()->dim(), nbq);
      if (L2Weight.size() == 0){
        for (size_t iqn = 0; iqn < nbq; iqn++){
          WeightsTimesF_quad.col(iqn) = quad[iqn].w * F_quad[i].col(iqn);
        }
      }else{
        for (size_t iqn = 0; iqn < nbq; iqn++){
          WeightsTimesF_quad.col(iqn) = quad[iqn].w * (L2Weight[iqn] * F_quad[i].col(iqn).matrix());
        }
      }
      for (size_t j = jcut; j < ncols; j++){
        // Integrate F_i * G_j
        GSM(i, j) = (WeightsTimesF_quad * G_quad[j]).sum();
      }
    }


    return GSM;
  }

//-----------------------------------------------------
//------ Interpolate a continuous function
//-----------------------------------------------------

template<typename Function>
Eigen::VectorXd LEPNCCore::nc_interpolate_moments(const Function& f, size_t doe) const {
	const Mesh* mesh_ptr = get_mesh();
	size_t ncells = mesh_ptr->n_cells();
	size_t nedges = mesh_ptr->n_edges();
  Eigen::VectorXd Xh = Eigen::VectorXd::Zero(DimPoly<Cell>(1)*ncells + nedges);

	// Each edge components of the vector of DOFs is the integral of f over the edge
	size_t offset_edges = ncells * DimPoly<Cell>(1);
  for (size_t iF = 0; iF < nedges; iF++) {
		QuadratureRule quadF = generate_quadrature_rule(*mesh_ptr->edge(iF), doe);
		for (QuadratureNode quadrule : quadF){
			Xh(offset_edges + iF) += quadrule.w * f(quadrule.x, quadrule.y);
		}
		Xh(offset_edges + iF) /= mesh_ptr->edge(iF)->measure();
	}

	// For components corresponding to the first DimPoly<Cell>(1) basis functions in each cell,
	// we put the L2 projection of the function minus its interpolate on the other basis functions
	// of the cell (computed using the integrals above)

  for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++) {

    // Mass matrix of first basis functions in cell (not the ones associated with edges)
		QuadratureRule quadT = generate_quadrature_rule(*mesh_ptr->cell(iT), doe, true);
		std::vector<Eigen::ArrayXd> nc_phi_quadT = nc_basis_quad(iT, quadT);
		Eigen::MatrixXd MT = gram_matrix(nc_phi_quadT, nc_phi_quadT, DimPoly<Cell>(1), DimPoly<Cell>(1), quadT, true);

		// Vector of integrals of f minus combination of functions corresponding to the edges, against basis functions
		std::vector<Edge*> edgesT = mesh_ptr->cell(iT)->get_edges();
    Eigen::VectorXd bT = Eigen::VectorXd::Zero(DimPoly<Cell>(1));
    for (size_t i = 0; i < DimPoly<Cell>(1); i++) {
      const auto& phi_i = nc_basis(iT, i);
			std::function<double(double,double)> adjusted_f = [&f, &edgesT, &Xh, &offset_edges, &iT, this](double x, double y) {
				double val = f(x,y);
				for (size_t ilF = 0; ilF < edgesT.size(); ilF++){
					size_t iF = edgesT[ilF]->global_index();
					val -= Xh(offset_edges + iF) * nc_basis(iT, DimPoly<Cell>(1)+ilF)(x, y);
				}
        return val;
      };
			for (QuadratureNode quadrule : quadT){
				bT(i) += quadrule.w * phi_i(quadrule.x, quadrule.y) * adjusted_f(quadrule.x, quadrule.y);
			}

    }

		// Vector of coefficients (on non-conforming cell basis functions) of the L2(T) projection of f
    Eigen::VectorXd UT = MT.ldlt().solve(bT);

		// Fill in the complete vector of DOFs
    Xh.segment(iT*DimPoly<Cell>(1), DimPoly<Cell>(1)) = UT;
	}

  return Xh;
}


template<typename Function>
Eigen::VectorXd LEPNCCore::nc_interpolate_ml(const Function& f, size_t doe) const {
	const Mesh* mesh_ptr = get_mesh();
	size_t ncells = mesh_ptr->n_cells();
	size_t nedges = mesh_ptr->n_edges();
  Eigen::VectorXd Xh = Eigen::VectorXd::Zero(DimPoly<Cell>(1)*ncells + nedges);

	// Each edge components of the vector of DOFs is the integral of f over the edge
	size_t offset_edges = ncells * DimPoly<Cell>(1);
  for (size_t iF = 0; iF < nedges; iF++) {
		QuadratureRule quadF = generate_quadrature_rule(*mesh_ptr->edge(iF), doe);
		for (QuadratureNode quadrule : quadF){
			Xh(offset_edges + iF) += quadrule.w * f(quadrule.x, quadrule.y);
		}
		Xh(offset_edges + iF) /= mesh_ptr->edge(iF)->measure();
	}

	// For components corresponding to the first DimPoly<Cell>(1) basis functions in each cell,
	// it's just the values at the mass-lumping nodes

  for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++) {
		for (size_t i = 0; i < DimPoly<Cell>(1); i++){
			VectorRd node = ml_node(iT, i);
			Xh(iT*DimPoly<Cell>(1)+i) = f(node.x(), node.y());
		}
	}

  return Xh;
}
//-------------------------------------------
//------------ Norms of discrete functions
//-------------------------------------------


Eigen::VectorXd LEPNCCore::nc_restr(const Eigen::VectorXd &Xh, size_t iT) const {
	
	size_t nedgesT = get_mesh()->cell(iT)->n_edges();
	Eigen::VectorXd XTF = Eigen::VectorXd::Zero(DimPoly<Cell>(1)+nedgesT);

	XTF.head(DimPoly<Cell>(1)) = Xh.segment(iT*DimPoly<Cell>(1), DimPoly<Cell>(1));

	size_t offset_edges = get_mesh()->n_cells() * DimPoly<Cell>(1);
	for (size_t ilE = 0; ilE < nedgesT; ilE++){
		size_t iE = get_mesh()->cell(iT)->edge(ilE)->global_index();
		XTF(DimPoly<Cell>(1)+ilE) = Xh(offset_edges + iE);
	}

	return XTF;
}


double LEPNCCore::nc_L2norm(const Eigen::VectorXd &Xh) const {
  double value = 0.0;
	
	size_t ncells = get_mesh()->n_cells();

  for (size_t iT = 0; iT < ncells; iT++) {
		// L2 norm computed using the mass matrix
		// Compute cell quadrature nodes and values of cell basis functions at these nodes
    Cell* cell = get_mesh()->cell(iT);
		size_t nedgesT = cell->n_edges();
		size_t nlocal_dofs = DimPoly<Cell>(1)+nedgesT;
		size_t doe = 4;
		QuadratureRule quadT = generate_quadrature_rule(*cell, doe, true);
		std::vector<Eigen::ArrayXd> nc_phiT_quadT = nc_basis_quad(iT, quadT);

		Eigen::MatrixXd MTT = gram_matrix(nc_phiT_quadT, nc_phiT_quadT, nlocal_dofs, nlocal_dofs, quadT, true);
		Eigen::VectorXd XTF = nc_restr(Xh, iT);

		value += XTF.dot(MTT*XTF);

  }
  return sqrt(value);
}


double LEPNCCore::nc_L2norm_ml(const Eigen::VectorXd &Xh) const {
  double value = 0.0;
	
	size_t ncells = get_mesh()->n_cells();

  for (size_t iT = 0; iT < ncells; iT++) {
		Cell* cell = get_mesh()->cell(iT);

		// L2 norm computed using mass-lumped version, so only based on first DimPoly<Cell>(1) basis
		// functions, with very simple coefficients
		for (size_t i = 0; i < DimPoly<Cell>(1); i++){
			value += (cell->measure()/DimPoly<Cell>(1)) * std::pow(Xh(iT*DimPoly<Cell>(1)+i), 2);
		}
	}

  return sqrt(value);
}


double LEPNCCore::nc_H1norm(const Eigen::VectorXd &Xh) const {
  double value = 0.0;
	
	size_t ncells = get_mesh()->n_cells();

  for (size_t iT = 0; iT < ncells; iT++) {
    Cell* cell = get_mesh()->cell(iT);
		size_t nlocal_dofs = DimPoly<Cell>(1) + cell->n_edges();
		size_t doe = 2;
		QuadratureRule quadT = generate_quadrature_rule(*cell, doe, true);
		std::vector<Eigen::ArrayXXd> Dnc_phiT_quad = grad_nc_basis_quad(iT, quadT);
		Eigen::MatrixXd StiffT = gram_matrix(Dnc_phiT_quad, Dnc_phiT_quad, nlocal_dofs, nlocal_dofs, quadT, true);
		Eigen::VectorXd XTF = nc_restr(Xh, iT);
		value += XTF.dot(StiffT * XTF);

  }
  return sqrt(value);
}


// ----------------------------------------------------------------
// ------- Compute vertex values from a non-conforming function
// ----------------------------------------------------------------

double LEPNCCore::nc_evaluate_in_cell(const Eigen::VectorXd Xh, size_t iT, double x, double y) const {
  double value = 0.0;
	const Mesh* mesh_ptr = get_mesh();
  for (size_t i = 0; i < DimPoly<Cell>(1); i++) {
    const auto &phi_i = nc_basis(iT, i);
    const size_t index = iT * DimPoly<Cell>(1) + i;
    value += Xh(index) * phi_i(x,y);
  }

	size_t nedgesT = mesh_ptr->cell(iT)->n_edges();
	size_t offset_edges = mesh_ptr->n_cells() * DimPoly<Cell>(1);
  for (size_t ilE = 0; ilE < nedgesT; ilE++) {
    const auto &phi_i = nc_basis(iT, DimPoly<Cell>(1) + ilE);
    const size_t index = offset_edges + mesh_ptr->cell(iT)->edge(ilE)->global_index();
    value += Xh(index) * phi_i(x,y);
  }

  return value;
}

Eigen::VectorXd LEPNCCore::nc_VertexValues(const Eigen::VectorXd Xh, const double weight) {
	const Mesh* mesh_ptr = get_mesh();

	// Average of values, at the vertices, given by basis functions associated with the edges around each vertex
	Eigen::VectorXd function_vertex_cells = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
	for (size_t iV = 0; iV < mesh_ptr->n_vertices(); iV++){
		auto xV = mesh_ptr->vertex(iV)->coords();
		auto cList = mesh_ptr->vertex(iV)->get_cells();
		auto nb_cells = cList.size();
		for (size_t ilT = 0; ilT < nb_cells; ilT++){
			auto temp = nc_evaluate_in_cell(Xh, cList[ilT]->global_index(), xV.x(), xV.y());
			function_vertex_cells(iV) += temp;
		}
		function_vertex_cells(iV) /= nb_cells;
	}

	// Average of values, at the vertices, given by basis functions associated with the edges, considered at the edge
	//	midpoints (they vanish on the vertex)
	Eigen::VectorXd function_vertex_edges = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
	for (size_t iV = 0; iV < mesh_ptr->n_vertices(); iV++){
		auto eList = mesh_ptr->vertex(iV)->get_edges();
		auto nb_edges = eList.size();
		for (Edge* edge : eList){
			auto cList_edge = edge->get_cells();
			VectorRd xE = edge->center_mass();
			// We average the values considered from both cells around the edge
			auto nb_cells_edges = cList_edge.size();
			for (Cell* cell : cList_edge){
				function_vertex_edges(iV) += nc_evaluate_in_cell(Xh, cell->global_index(), xE.x(), xE.y()) / nb_cells_edges;
			}
		}
		function_vertex_edges(iV) /= nb_edges;
	}

	return (1-weight)*function_vertex_cells + weight*function_vertex_edges;
}


// Sign function
double LEPNCCore::sign(const double s) const {
	double val = 0;
  if (s>0){
		val = 1;
	}else if (s<0){
		val = -1;
	}

	return val;	
}

//@}

}  // end of namespace HArDCore2D

#endif /* LEPNCCORE_HPP */
