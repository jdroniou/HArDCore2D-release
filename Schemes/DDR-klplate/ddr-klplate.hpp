// Implementation of the discrete de Rham sequence for the Kirchhoff-Love plate problem.
//
// Authors: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr) and Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef DDR_KLPLATE_HPP
#define DDR_KLPLATE_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <xdivdiv.hpp>

#include <local_static_condensation.hpp>

/*!
 * @defgroup DDR_klplate
 * @brief Implementation of the DDR scheme for the Kirchoff-Love plate problem
 */

namespace HArDCore2D
{

  /*!
   * @addtogroup DDR_klplate
   * @{
   */
      
  /// Assemble a Kirchhoff-Love plate problem 
  struct KirchhoffLove
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<double(const Eigen::Vector2d &)> ForcingTermType;
    typedef std::function<double(const Eigen::Vector2d &)> DeflectionType;
    typedef std::function<Eigen::Vector2d(const Eigen::Vector2d &)> GradientDeflectionType;
    typedef std::function<Eigen::Matrix2d(const Eigen::Vector2d &)> MomentTensorType;
    typedef std::function<double(const Eigen::Vector2d &, const Edge &)> MomentTensorEdgeDerivativeType;

    typedef typename XDivDiv::ConstitutiveLawType ConstitutiveLawType;

    /// Constructor
    KirchhoffLove(
                   const PlatesCore & platescore,       ///< Core for the DDR space sequence
                   const ConstitutiveLawType & law,     ///< Constitutive law: law(sigma)+hess u = 0
                   bool use_threads,                    ///< True for parallel execution, false for sequential execution
                   std::ostream & output = std::cout    ///< Output stream to print status messages
                   );

    /// Returns the dimension of the moment + deflection space
    inline size_t dimensionSpace() const
    {
      return m_xdivdiv.dimension() + m_Pkm2_Th.dimension();
    }

    /// Returns the number of statically condensed DOFs (here, the cell moments DOFs)
    inline size_t nbSCDOFs() const
    {
      return m_platescore.mesh().n_cells() * 
          (PolynomialSpaceDimension<Cell>::HolyCompl(m_platescore.degree() - 1) 
           + PolynomialSpaceDimension<Cell>::Holy(m_platescore.degree() - 4));
    }

    /// Returns the size of the statically condensed system
    inline size_t sizeSystem() const
    {
      return dimensionSpace() - nbSCDOFs();
    }

    /// Returns the space XDivDiv
    inline const XDivDiv & xDivDiv() const {
      return m_xdivdiv;
    }

    /// Returns the space \f$\mathbb{P}^{k-2}(\mathcal{T}_h)\f$
    inline const GlobalDOFSpace polykm2Th() const {
      return m_Pkm2_Th;
    }
    
    /// Returns the linear system matrix
    inline const SystemMatrixType & systemMatrix() const {
      return m_A;
    }

    /// Returns the linear system matrix
    inline SystemMatrixType & systemMatrix() {
      return m_A;
    }

    /// Returns the linear system right-hand side vector
    inline const Eigen::VectorXd & systemVector() const {
      return m_b;
    }

    /// Returns the linear system right-hand side vector
    inline Eigen::VectorXd & systemVector() {
      return m_b;
    }  
    
    /// Returns the static condensation recovery operator
    inline const SystemMatrixType & scMatrix() const {
      return m_sc_A;
    }

    /// Returns the static condensation rhs
    inline Eigen::VectorXd & scVector() {
      return m_sc_b;
    }

    /// Returns the stabilization parameter
    inline const double & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline double & stabilizationParameter() {
      return m_stab_par;
    }

    /// Assemble the global system
    void assembleLinearSystem(
                              const ForcingTermType & f, ///< Forcing term
                              const DeflectionType & u,   ///< Deflection (for BC)
                              const GradientDeflectionType & grad_u   ///< Gradient of deflection (for BC)
                              );

    /// Interpolate deflection
    Eigen::VectorXd interpolateDeflection(
                                          const DeflectionType & u, ///< The function to interpolate
                                          int deg_quad = -1 ///< The degree of the quadrature rules (equal to \f$2(k+1)\f$ with \f$k\f$ equal to the degree of the complex by default)
                                          ) const;
    
    /// Compute the discrete norm
    double computeNorm(
                       const Eigen::VectorXd & v ///< The vector
                      ) const;

    
  private:
    const PlatesCore & m_platescore;
    ConstitutiveLawType m_law;
    bool m_use_threads;
    std::ostream & m_output;

    XDivDiv m_xdivdiv;
    GlobalDOFSpace m_Pkm2_Th;
    
    SystemMatrixType m_A;   // Matrix and RHS of statically condensed system
    Eigen::VectorXd m_b;
    SystemMatrixType m_sc_A; // Static condensation operator and RHS (to recover statically condensed DOFs)
    Eigen::VectorXd m_sc_b;

    double m_stab_par;

    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_local_contribution(
                                size_t iT, ///< Element index
                                const ForcingTermType & f, ///< Orthogonal load
                                const DeflectionType & u, ///< Deflection (for BC)
                                const GradientDeflectionType & grad_u ///< Gradient of deflection (for BC)
                                );

    /// Creates the permutation matrix and the global DOFs for the local static condensation
    LocalStaticCondensation _compute_static_condensaton(const size_t & iT) const;

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
                                      size_t iT,                                               ///< Element index
                                      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
                                      std::list<Eigen::Triplet<double> > & Ah_sys,             ///< List of triplets for system
                                      Eigen::VectorXd & bh_sys,                                ///< Vector for RHS for sysem
                                      std::list<Eigen::Triplet<double> > & Ah_sc,              ///< List of triplets for SC recovery
                                      Eigen::VectorXd & bh_sc                                  ///< Vector for RHS for SC recovery
                                      );

  };

  //------------------------------------------------------------------------------
  // Exact solutions
  //------------------------------------------------------------------------------

  static const double PI = boost::math::constants::pi<double>();
  using std::sin;
  using std::cos;
  using std::exp;
  using std::pow;

  //------------------------------------------------------------------------------
  // Trigonometric
  
  static KirchhoffLove::DeflectionType
  trigonometric_u = [](const VectorRd & x) -> double {
    return sin(PI*x(0))*sin(PI*x(1));
  };

  static KirchhoffLove::GradientDeflectionType
  trigonometric_grad_u = [](const VectorRd & x) -> VectorRd {
    return PI * VectorRd(sin(PI*x(1))*cos(PI*x(0)), sin(PI*x(0))*cos(PI*x(1)));
  };

  static KirchhoffLove::MomentTensorType
  trigonometric_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
    M.row(0) << -sin(PI*x(0))*sin(PI*x(1)), cos(PI*x(0))*cos(PI*x(1));
    M.row(1) << cos(PI*x(0))*cos(PI*x(1)), -sin(PI*x(0))*sin(PI*x(1));

    return pow(PI, 2) * M;
  };

  static KirchhoffLove::MomentTensorType
  trigonometric_dx_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
    M.row(0) << -cos(PI*x(0))*sin(PI*x(1)), -sin(PI*x(0))*cos(PI*x(1));
    M.row(1) << -sin(PI*x(0))*cos(PI*x(1)), -cos(PI*x(0))*sin(PI*x(1));

    return pow(PI, 3) * M;
  };


  static KirchhoffLove::MomentTensorType
  trigonometric_dy_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
    M.row(0) << -sin(PI*x(0))*cos(PI*x(1)), -cos(PI*x(0))*sin(PI*x(1));
    M.row(1) << -cos(PI*x(0))*sin(PI*x(1)), -sin(PI*x(0))*cos(PI*x(1));

    return pow(PI, 3) * M;
  };

  static KirchhoffLove::GradientDeflectionType
  trigonometric_grad_Delta_u = [](const VectorRd & x) -> VectorRd {
  
    Eigen::Vector2d G {
      trigonometric_dx_hess_u(x)(0,0)+trigonometric_dx_hess_u(x)(1,1), 
      trigonometric_dy_hess_u(x)(0,0)+trigonometric_dy_hess_u(x)(1,1)
      };
    
    return G;
  };

  static KirchhoffLove::MomentTensorEdgeDerivativeType
  trigonometric_hess_u_DE = [](const Vector2d & x, const Edge & E) -> double {

    Eigen::Vector2d tE = E.tangent();
    Eigen::Vector2d nE = E.normal();

    Eigen::Matrix2d dtE_hess = trigonometric_dx_hess_u(x) * tE.x() + trigonometric_dy_hess_u(x) * tE.y();

    Eigen::Vector2d div_hess {
      trigonometric_dx_hess_u(x)(0,0) + trigonometric_dy_hess_u(x)(0,1),
      trigonometric_dx_hess_u(x)(1,0) + trigonometric_dy_hess_u(x)(1,1),
    };

    return (dtE_hess*nE).dot(tE) + div_hess.dot(nE);
  };

  static KirchhoffLove::ForcingTermType
  trigonometric_divdiv_hess_u = [](const VectorRd & x) -> double {
    return 4.*pow(PI, 4)*sin(PI*x(0))*sin(PI*x(1));
  };

  //------------------------------------------------------------------------------
  // Quartic
  
  static KirchhoffLove::DeflectionType
  quartic_u = [](const VectorRd & x) -> double {
    return pow(x(0),2)*pow(1.-x(0),2) + pow(x(1),2)*pow(1.-x(1),2);
  };
  
  static KirchhoffLove::GradientDeflectionType
  quartic_grad_u = [](const VectorRd & x) -> VectorRd {
    return VectorRd(pow(x(0), 2)*(2*x(0) - 2) + 2*x(0)*pow(1 - x(0), 2), pow(x(1), 2)*(2*x(1) - 2) + 2*x(1)*pow(1 - x(1), 2));
  };

  static KirchhoffLove::MomentTensorType
  quartic_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {

    Eigen::Matrix2d M;
    M <<
      2.*(pow(x(0),2)+4.*x(0)*(x(0)-1.)+pow(1.-x(0),2)), 0.,
      0., 2.*(pow(x(1),2)+4.*x(1)*(x(1)-1.)+pow(1.-x(1),2));

    return M;
  };

  static KirchhoffLove::MomentTensorType
  quartic_dx_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
    M.row(0) << 24*x(0) - 12, 0;
    M.row(1) << 0., 0.;

    return M;
  };


  static KirchhoffLove::MomentTensorType
  quartic_dy_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
    M.row(0) << 0., 0.;
    M.row(1) << 0., 24*x(1) - 12;

    return M;
  };

  static KirchhoffLove::GradientDeflectionType
  quartic_grad_Delta_u = [](const VectorRd & x) -> VectorRd {
  
    Eigen::Vector2d G {
      quartic_dx_hess_u(x)(0,0)+quartic_dx_hess_u(x)(1,1), 
      quartic_dy_hess_u(x)(0,0)+quartic_dy_hess_u(x)(1,1)
      };
    
    return G;
  };

  static KirchhoffLove::MomentTensorEdgeDerivativeType
  quartic_hess_u_DE = [](const Vector2d & x, const Edge & E) -> double {

    Eigen::Vector2d tE = E.tangent();
    Eigen::Vector2d nE = E.normal();

    Eigen::Matrix2d dtE_hess = quartic_dx_hess_u(x) * tE.x() + quartic_dy_hess_u(x) * tE.y();

    Eigen::Vector2d div_hess {
      quartic_dx_hess_u(x)(0,0) + quartic_dy_hess_u(x)(0,1),
      quartic_dx_hess_u(x)(1,0) + quartic_dy_hess_u(x)(1,1),
    };

    return (dtE_hess*nE).dot(tE) + div_hess.dot(nE);
  };


  static KirchhoffLove::ForcingTermType
  quartic_divdiv_hess_u = [](const VectorRd & x) -> double {
    return 48.;
  };

  //------------------------------------------------------------------------------
  // Biquartic
  
  static KirchhoffLove::DeflectionType
  biquartic_u = [](const VectorRd & x) -> double {
    return pow(x(0),2)*pow(1.-x(0),2)*pow(x(1),2)*pow(1.-x(1),2);
  };
  
  static KirchhoffLove::GradientDeflectionType
  biquartic_grad_u = [](const VectorRd & x) -> VectorRd {
    return VectorRd(pow(x(0), 2)*pow(x(1), 2)*pow(1 - x(1), 2)*(2*x(0) - 2) + 2*x(0)*pow(x(1), 2)*pow(1 - x(0), 2)*pow(1 - x(1), 2), pow(x(0), 2)*pow(x(1), 2)*pow(1 - x(0), 2)*(2*x(1) - 2) + 2*pow(x(0), 2)*x(1)*pow(1 - x(0), 2)*pow(1 - x(1), 2));
  };

  static KirchhoffLove::MomentTensorType
  biquartic_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
      M.row(0) << 2*pow(x(1), 2)*pow(x(1) - 1, 2)*(pow(x(0), 2) + 4*x(0)*(x(0) - 1) + pow(x(0) - 1, 2)), 4*x(0)*x(1)*(x(0) - 1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1));      
      M.row(1) << 4*x(0)*x(1)*(x(0) - 1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)), 2*pow(x(0), 2)*pow(x(0) - 1, 2)*(pow(x(1), 2) + 4*x(1)*(x(1) - 1) + pow(x(1) - 1, 2));
    return M;
  };

  static KirchhoffLove::MomentTensorType
  biquartic_dx_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
    M.row(0) << 2*pow(x(1), 2)*(12*x(0) - 6)*pow(x(1) - 1, 2), 4*x(0)*x(1)*(x(0) - 1)*(x(1) - 1)*(4*x(1) - 2) + 4*x(0)*x(1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)) + 4*x(1)*(x(0) - 1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1));
    M.row(1) << 4*x(0)*x(1)*(x(0) - 1)*(x(1) - 1)*(4*x(1) - 2) + 4*x(0)*x(1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)) + 4*x(1)*(x(0) - 1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)), 2*pow(x(0), 2)*(2*x(0) - 2)*(pow(x(1), 2) + 4*x(1)*(x(1) - 1) + pow(x(1) - 1, 2)) + 4*x(0)*pow(x(0) - 1, 2)*(pow(x(1), 2) + 4*x(1)*(x(1) - 1) + pow(x(1) - 1, 2));

    return M;
  };


  static KirchhoffLove::MomentTensorType
  biquartic_dy_hess_u = [](const VectorRd & x) -> Eigen::Matrix2d {
    Eigen::Matrix2d M;
    M.row(0) << 2*pow(x(1), 2)*(2*x(1) - 2)*(pow(x(0), 2) + 4*x(0)*(x(0) - 1) + pow(x(0) - 1, 2)) + 4*x(1)*pow(x(1) - 1, 2)*(pow(x(0), 2) + 4*x(0)*(x(0) - 1) + pow(x(0) - 1, 2)), 4*x(0)*x(1)*(x(0) - 1)*(4*x(0) - 2)*(x(1) - 1) + 4*x(0)*x(1)*(x(0) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)) + 4*x(0)*(x(0) - 1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1));
    M.row(1) << 4*x(0)*x(1)*(x(0) - 1)*(4*x(0) - 2)*(x(1) - 1) + 4*x(0)*x(1)*(x(0) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)) + 4*x(0)*(x(0) - 1)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)), 2*pow(x(0), 2)*pow(x(0) - 1, 2)*(12*x(1) - 6);

    return M;
  };

  static KirchhoffLove::GradientDeflectionType
  biquartic_grad_Delta_u = [](const VectorRd & x) -> VectorRd {
  
    Eigen::Vector2d G {
      biquartic_dx_hess_u(x)(0,0)+biquartic_dx_hess_u(x)(1,1), 
      biquartic_dy_hess_u(x)(0,0)+biquartic_dy_hess_u(x)(1,1)
      };
    
    return G;
  };

  static KirchhoffLove::MomentTensorEdgeDerivativeType
  biquartic_hess_u_DE = [](const Vector2d & x, const Edge & E) -> double {

    Eigen::Vector2d tE = E.tangent();
    Eigen::Vector2d nE = E.normal();

    Eigen::Matrix2d dtE_hess = biquartic_dx_hess_u(x) * tE.x() + biquartic_dy_hess_u(x) * tE.y();

    Eigen::Vector2d div_hess {
      biquartic_dx_hess_u(x)(0,0) + biquartic_dy_hess_u(x)(0,1),
      biquartic_dx_hess_u(x)(1,0) + biquartic_dy_hess_u(x)(1,1),
    };

    return (dtE_hess*nE).dot(tE) + div_hess.dot(nE);
  };
  
  static KirchhoffLove::ForcingTermType
  biquartic_divdiv_hess_u = [](const VectorRd & x) -> double {
    return 24*pow(x(0), 2)*pow(x(0) - 1, 2) + 32*x(0)*x(1)*(x(0) - 1)*(x(1) - 1) + 8*x(0)*x(1)*(x(0) - 1)*(4*x(1) - 2) + 8*x(0)*x(1)*(4*x(0) - 2)*(x(1) - 1) + 8*x(0)*x(1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)) + 8*x(0)*(x(0) - 1)*(x(1) - 1)*(4*x(1) - 2) + 8*x(0)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)) + 24*pow(x(1), 2)*pow(x(1) - 1, 2) + 8*x(1)*(x(0) - 1)*(4*x(0) - 2)*(x(1) - 1) + 8*x(1)*(x(0) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1)) + 2*(4*x(0) - 4)*(x(1) - 1)*(x(0)*x(1) + x(0)*(x(1) - 1) + x(1)*(x(0) - 1) + (x(0) - 1)*(x(1) - 1));
  };

  
} // end of namespace HArDCore2D

#endif
