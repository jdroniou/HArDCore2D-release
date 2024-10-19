#include <iostream>

#include <basis.hpp>
#include <mesh.hpp>

namespace HArDCore2D
{

  //------------------------------------------------------------------------------
  // Scalar monomial basis on a cell
  //------------------------------------------------------------------------------

  MonomialScalarBasisCell::MonomialScalarBasisCell(const Cell & T, size_t degree)
    : m_degree(degree),
      m_xT(T.center_mass()),
      m_hT(T.diam()),
      m_powers(MonomialPowers<Cell>::complete(degree))
  {
    // create rotation pi/2 for curl
    m_rot.row(0) << 0., 1.;
    m_rot.row(1) << -1., 0.;
  }

  MonomialScalarBasisCell::FunctionValue MonomialScalarBasisCell::function(size_t i, const VectorRd & x) const
  {
    VectorRd y = _coordinate_transform(x);
    const VectorZd & powers = m_powers[i];
    return std::pow(y(0), powers(0)) * std::pow(y(1), powers(1));
  }

  MonomialScalarBasisCell::GradientValue MonomialScalarBasisCell::gradient(size_t i, const VectorRd & x) const
  {
    VectorRd y = _coordinate_transform(x);
    const VectorZd & powers = m_powers[i];
    VectorRd grad = VectorRd::Zero();
    grad(0) = (powers(0) == 0 ? 0. : powers(0) * std::pow(y(0), powers(0)-1) * std::pow(y(1), powers(1)) );
    grad(1) = (powers(1) == 0 ? 0. : std::pow(y(0), powers(0)) * powers(1) * std::pow(y(1), powers(1)-1) );
    return grad / m_hT;
  }  

  MonomialScalarBasisCell::CurlValue MonomialScalarBasisCell::curl(size_t i, const VectorRd & x) const
  {
    return m_rot * gradient(i, x);
  }

  MonomialScalarBasisCell::HessianValue MonomialScalarBasisCell::hessian(size_t i, const VectorRd & x) const
  {
    VectorRd y = _coordinate_transform(x);
    const VectorZd & powers = m_powers[i];
    Eigen::Matrix2d hess = Eigen::Matrix2d::Zero();
    hess(0, 0) = ( powers(0) < 2 ? 0. : powers(0) * (powers(0) - 1) * std::pow(y(0), powers(0)-2) * std::pow(y(1), powers(1)) );
    hess(1, 1) = ( powers(1) < 2 ? 0. : powers(1) * (powers(1) - 1) * std::pow(y(0), powers(0)) * std::pow(y(1), powers(1)-2) );
    hess(1, 0) = hess(0, 1) = ( powers(0) == 0 ? 0. : powers(0) * std::pow(y(0), powers(0) - 1) )
      * ( powers(1) == 0 ? 0. : powers(1) * std::pow(y(1), powers(1) - 1) );
    return hess / std::pow(m_hT, 2);
  }

  //------------------------------------------------------------------------------
  // Scalar monomial basis on an edge
  //------------------------------------------------------------------------------

  MonomialScalarBasisEdge::MonomialScalarBasisEdge(const Edge & E, size_t degree)
    : m_degree(degree),
      m_xE(E.center_mass()),
      m_hE(E.diam()),
      m_tE(E.tangent())
  {
    // Do nothing    
  }

  MonomialScalarBasisEdge::FunctionValue MonomialScalarBasisEdge::function(size_t i, const VectorRd & x) const
  {
    return std::pow(_coordinate_transform(x), i);
  }

  MonomialScalarBasisEdge::GradientValue MonomialScalarBasisEdge::gradient(size_t i, const VectorRd & x) const
  {
    return (i == 0 ? 0. : i * std::pow(_coordinate_transform(x), i - 1) / m_hE) * m_tE;
  }

  //------------------------------------------------------------------------------
  // Basis for R^{c,k}(T)
  //------------------------------------------------------------------------------

  RolyComplBasisCell::RolyComplBasisCell(const Cell &T, size_t degree)
    : m_degree(degree),
      m_xT(T.center_mass()),
      m_hT(T.diam())
  {
    // Monomial powers for P^{k-1}(T)
    if (m_degree >= 1){
      m_powers = MonomialPowers<Cell>::complete(m_degree-1);
    }else{
      std::cout << "Attempting to construct RckT with degree 0, stopping" << std::endl;
      exit(1);
    }
  }

  RolyComplBasisCell::FunctionValue RolyComplBasisCell::function(size_t i, const VectorRd &x) const
  {
    VectorRd y = _coordinate_transform(x);
    const VectorZd &powers = m_powers[i];
    return std::pow(y(0), powers(0)) * std::pow(y(1), powers(1)) * y;
  }
  
  RolyComplBasisCell::DivergenceValue RolyComplBasisCell::divergence(size_t i, const VectorRd &x) const
  {
    VectorRd y = _coordinate_transform(x);
    const VectorZd &powers = m_powers[i];
    return (powers(0)+powers(1)+2) * std::pow(y(0), powers(0)) * std::pow(y(1), powers(1)) / m_hT;
  }

  //------------------------------------------------------------------------------
  // Basis for G^{c,k}(T)
  //------------------------------------------------------------------------------
  
  GolyComplBasisCell::GolyComplBasisCell(const Cell &T, size_t degree)
    : m_degree(degree)
  {
    m_rot.row(0) << 0., 1.;
    m_rot.row(1) << -1., 0.;
    m_Rck_basis.reset(new RolyComplBasisCell(T, degree));
  }

  GolyComplBasisCell::FunctionValue GolyComplBasisCell::function(size_t i, const VectorRd &x) const
  {
    // The basis of Gck(T) is a simple rotation of the basis of Rck
    return m_rot * m_Rck_basis->function(i, x);
  }

  GolyComplBasisCell::RotorValue GolyComplBasisCell::rotor(size_t i, const VectorRd &x) const
  {
    // The basis of Gck(T) is a rotation of the basis of Rck of \f$-\frac{\pi}{2}\f$, so the rotor of the \f$i\f$-th basis function of Gck is opposite of the divergence of the corresponding function in Rck
    return -m_Rck_basis->divergence(i, x);
  }
  
  //------------------------------------------------------------------------------
  // Basis for H^{c,k}(T)
  //------------------------------------------------------------------------------

  HolyComplBasisCell::HolyComplBasisCell(const Cell &T, size_t degree)
    : m_degree(degree),
      m_xT(T.center_mass()),
      m_hT(T.diam()),
      m_Pkmo(MonomialScalarBasisCell(T, degree - 1))
  {   
    // create rotation pi/2
    m_rot.row(0) << 0., 1.;
    m_rot.row(1) << -1., 0.;
  }

  HolyComplBasisCell::FunctionValue HolyComplBasisCell::function(size_t i, const VectorRd &x) const
  {
    VectorRd y = m_rot * _coordinate_transform(x);
    Eigen::Matrix2d M = y * m_Pkmo.function(i, x).transpose();
    return 0.5 * (M + M.transpose());
  }
  
  //------------------------------------------------------------------------------
  // A common notion of scalar product for scalars and vectors
  //------------------------------------------------------------------------------
  
  double scalar_product(const double & x, const double & y) {
    return x * y;
  }

  double scalar_product(const VectorRd & x, const VectorRd & y) {
    return x.dot(y);
  }

  //------------------------------------------------------------------------------
  // Matrix-vector product
  //------------------------------------------------------------------------------

  boost::multi_array<VectorRd, 2> matrix_vector_product(
                                                        const boost::multi_array<MatrixRd, 2> & basis_quad, ///< The basis evaluation
                                                        const std::vector<VectorRd> & v_quad ///< The vector to take the product with
                                                        )
  {
    assert(basis_quad.shape()[1] == v_quad.size());
    boost::multi_array<VectorRd, 2> basis_v_quad(boost::extents[basis_quad.shape()[0]][basis_quad.shape()[1]]);
    for (std::size_t i = 0; i < basis_quad.shape()[0]; i++) {
      for (std::size_t iqn = 0; iqn < basis_quad.shape()[1]; iqn++) {
        basis_v_quad[i][iqn] = basis_quad[i][iqn] * v_quad[iqn];
      } // for iqn
    } // for i
    return basis_v_quad;
  }

  boost::multi_array<VectorRd, 2> matrix_vector_product(
                                                        const std::vector<MatrixRd> & m_quad, ///< The matrix to take the product with
                                                        const boost::multi_array<VectorRd, 2> & basis_quad ///< The basis evaluation
                                                        )
  {
    assert(basis_quad.shape()[1] == m_quad.size());
    boost::multi_array<VectorRd, 2> m_basis_quad(boost::extents[basis_quad.shape()[0]][basis_quad.shape()[1]]);
    for (std::size_t i = 0; i < basis_quad.shape()[0]; i++) {
      for (std::size_t iqn = 0; iqn < basis_quad.shape()[1]; iqn++) {
        m_basis_quad[i][iqn] = m_quad[iqn] * basis_quad[i][iqn];
      } // for iqn
    } // for i
    return m_basis_quad;
  }

  boost::multi_array<VectorRd, 2> vector_matrix_product(
                                                        const std::vector<VectorRd> & v_quad, ///< The vector to take the product with
                                                        const boost::multi_array<MatrixRd, 2> & basis_quad ///< The basis evaluation                                                
                                                        )
  {
    assert(basis_quad.shape()[1] == v_quad.size());
    boost::multi_array<VectorRd, 2> basis_v_quad(boost::extents[basis_quad.shape()[0]][basis_quad.shape()[1]]);
    for (std::size_t i = 0; i < basis_quad.shape()[0]; i++) {
      for (std::size_t iqn = 0; iqn < basis_quad.shape()[1]; iqn++) {
        basis_v_quad[i][iqn] = v_quad[iqn].transpose() * basis_quad[i][iqn];
      } // for iqn
    } // for i
    return basis_v_quad;
  }
  
  //------------------------------------------------------------------------------
  //      Gram matrices
  //------------------------------------------------------------------------------

  // Vector2d for B1, tensorised double for B2
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<VectorRd, 2> & B1,
                                      const boost::multi_array<double, 2> & B2,
                                      const QuadratureRule & qr)
  {
    // Check that the basis evaluation and quadrature rule are coherent
    assert ( qr.size() == B1.shape()[1] && qr.size() == B2.shape()[1] );

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(B1.shape()[0], dimspace * B2.shape()[0]);
    for (size_t i = 0; i < B1.shape()[0]; i++) {
      for (size_t k = 0; k < dimspace; k++) {
        VectorRd ek = VectorRd::Zero();
        ek(k) = 1.;
        for(size_t j = 0; j < B2.shape()[0]; j++) {
          for (size_t iqn = 0; iqn < qr.size(); iqn++) {
            M(i,k * B2.shape()[0] + j) += qr[iqn].w * B1[i][iqn].dot(ek) * B2[j][iqn];
          } // for iqn
        } // for j
      } // for k
    } // for i
    return M;
  }

  // Gramm matrix for double-valued B1, B2
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<double, 2> & B1,
                                      const boost::multi_array<double, 2> & B2,
                                      const QuadratureRule & qr,
                                      const size_t nrows,
                                      const size_t ncols,
                                      const std::string sym
                                      )
  {
    // Check that the basis evaluation and quadrature rule are coherent
    assert ( qr.size() == B1.shape()[1] && qr.size() == B2.shape()[1] );
    // Check that we don't ask for more members of family than available
    assert ( nrows <= B1.shape()[0] && ncols <= B2.shape()[0] );

    // Re-cast quadrature weights into ArrayXd to facilitate computations
    Eigen::ArrayXd qr_weights = Eigen::ArrayXd::Zero(qr.size());
    for (size_t iqn = 0; iqn < qr.size(); iqn++){
      qr_weights(iqn) = qr[iqn].w;
    }

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nrows, ncols);
    for (size_t i = 0; i < nrows; i++) {
      size_t jcut = 0;
      if ( sym=="sym" ) jcut = i;
      for (size_t j = 0; j < jcut; j++){
        M(i, j) = M(j, i);
      }
      for(size_t j = jcut; j < ncols; j++) {
        std::vector<double> tmp(B1.shape()[1]);
        // Extract values at quadrature nodes for elements i of B1 and j of B2
        boost::multi_array<double, 1> B1i = B1[ boost::indices[i][boost::multi_array_types::index_range(0, B1.shape()[1])] ];
        boost::multi_array<double, 1> B2j = B2[ boost::indices[j][boost::multi_array_types::index_range(0, B1.shape()[1])] ];
        double *p1 = &B1i[0];
        double *p2 = &B2j[0];
        Eigen::ArrayXd B1i_as_array = Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(p1, B1i.shape()[0]);
        Eigen::ArrayXd B2j_as_array = Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(p2, B2j.shape()[0]);

        // Multiply by quadrature weights and sum (using .sum() of ArrayXd makes this step faster than a loop
        M(i,j) = (qr_weights * B1i_as_array * B2j_as_array).sum();
      } // for j
    } // for i
    return M;
  }

  // Gram matrix for double-valued complete family. Do not make this inline, this slows down calculations.
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<double, 2> & B1,
                                      const boost::multi_array<double, 2> & B2,
                                      const QuadratureRule & qr,
                                      const std::string sym
                                      )
  {
    return compute_gram_matrix(B1, B2, qr, B1.shape()[0], B2.shape()[0], sym);
  }


  // Gram matrix for Vector2d-valued B1 and B2
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<VectorRd, 2> & B1,
                                      const boost::multi_array<VectorRd, 2> & B2,
                                      const QuadratureRule & qr,
                                      const size_t nrows,
                                      const size_t ncols,
                                      const std::string sym
                                      )
  {
    // Check that the basis evaluation and quadrature rule are coherent
    assert ( qr.size() == B1.shape()[1] && qr.size() == B2.shape()[1] );
    // Check that we don't ask for more members of family than available
    assert ( nrows <= B1.shape()[0] && ncols <= B2.shape()[0] );

    // Re-cast quadrature weights into ArrayXd to make computations faster
    Eigen::ArrayXd qr_weights = Eigen::ArrayXd::Zero(qr.size());
    for (size_t iqn = 0; iqn < qr.size(); iqn++){
      qr_weights(iqn) = qr[iqn].w;
    }

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nrows, ncols);
    for (size_t i = 0; i < nrows; i++) {
      size_t jcut = 0;
      if ( sym=="sym" ) jcut = i;
      for (size_t j = 0; j < jcut; j++){
        M(i, j) = M(j, i);
      }
      for(size_t j = jcut; j < ncols; j++) {
        // Array of scalar products
        Eigen::ArrayXd B1i_dot_B2j = Eigen::ArrayXd::Zero(qr.size());
        for (size_t iqn = 0; iqn < qr.size(); iqn++){
          B1i_dot_B2j(iqn) = B1[i][iqn].dot(B2[j][iqn]);
        }

        // Multiply component-wise by weights and sum
        M(i,j) = (qr_weights * B1i_dot_B2j).sum();
      } // for j
    } // for i
    return M;
  }

  // Gram matrix for Vector2d-valued complete family. Do not make this inline, this slows down actual calculations.
  Eigen::MatrixXd compute_gram_matrix(const boost::multi_array<VectorRd, 2> & B1, 
                                      const boost::multi_array<VectorRd, 2> & B2, 
                                      const QuadratureRule & qr,       
                                      const std::string sym
                                      )
  {
    return compute_gram_matrix(B1, B2, qr, B1.shape()[0], B2.shape()[0], sym);
  }


  Eigen::MatrixXd compute_weighted_gram_matrix(const FType<VectorRd> &f, const BasisQuad<VectorRd> &B1, const BasisQuad<double> &B2, const QuadratureRule &qr, size_t n_rows, size_t n_cols)
  {
    // If default, set n_rows and n_cols to size of families
    if (n_rows == 0 && n_cols == 0)
      {
        n_rows = B1.shape()[0];
        n_cols = B2.shape()[0];
      }

    // Number of quadrature nodes
    const size_t num_quads = qr.size();
    // Check number of quadrature nodes is compatible with B1 and B2
    assert(num_quads == B1.shape()[1] && num_quads == B2.shape()[1]);
    // Check that we don't ask for more members of family than available
    assert(n_rows <= B1.shape()[0] && n_cols <= B2.shape()[0]);

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n_rows, n_cols);
    for (size_t iqn = 0; iqn < num_quads; iqn++)
      {
        double qr_weight = qr[iqn].w;
        VectorRd f_on_qr = f(qr[iqn].vector());
        for (size_t i = 0; i < n_rows; i++)
          {
            double f_dot_B1 = f_on_qr.dot(B1[i][iqn]);
            for (size_t j = 0; j < n_cols; j++)
              {
                M(i, j) += qr_weight * f_dot_B1 * B2[j][iqn];
              }
          }
      }
    return M;
  }

  Eigen::MatrixXd compute_weighted_gram_matrix(const FType<VectorRd> &f, const BasisQuad<double> &B1, const BasisQuad<VectorRd> &B2, const QuadratureRule &qr, size_t n_rows, size_t n_cols)
  {
    return compute_weighted_gram_matrix(f, B2, B1, qr, n_cols, n_rows).transpose();
  }

} // end of namespace HArDCore2D
