// Implementation of the discrete de Rham sequence for the Reissner-Mindlin plate problem.
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
 * The DDR method is described in
 *
 **  An arbitrary-order method for magnetostatics on polyhedral meshes based on a discrete de Rham
 sequence. 
 **   D. A. Di Pietro and J. Droniou, 31p, 2020. url: https://arxiv.org/abs/2005.06890.
 *
 *
 * The design and analysis of the DDR scheme for the Reissner-Mindlin plate implemented here is described in
 *  
 ** A DDR method for the Reissner–Mindlin plate bending problem on polygonal meshes
 ** D. A. Di Pietro and J. Droniou, 23p, 2021. url: https://arxiv.org/abs/2105.11773
 *
 * If you use this code in a scientific publication, please mention the above articles.
 *
 */

#ifndef MAGNETOSTATICS_HPP
#define MAGNETOSTATICS_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <BoundaryConditions/BoundaryConditions.hpp>   // To re-order the boundary edges and vertices
#include <BoundaryConditions/BChandlers.hpp>   

#include <xgrad.hpp>
#include <excurl.hpp>

/*!
 * @defgroup DDR_rmplate
 * @brief Implementation of the DDR scheme for the Reissner-Mindlin plate problem
 */

namespace HArDCore2D
{

  /*!
   * @addtogroup DDR_rmplate
   * @{
   */
   
  /// Structure to store model data
  struct RMParameters
  {
    /// Constructor
    RMParameters( const double thickness, 
                  const double young_modulus,
                  const double poisson_ratio
                ):
             t(thickness),
             E(young_modulus),
             nu(poisson_ratio)
             {
                // Compute derived quantities
                beta0 = E/(12*(1.+nu));
                beta1 = E*nu/(12*(1.-std::pow(nu,2)));
                kappa = 5*E/(12*(1.+nu));
             }
       
    double t; // plate thickness
    double E; // Young modulus
    double nu; // Poisson ratio
    double beta0;
    double beta1;
    double kappa;
  };

  /// Assemble a RM problem 
  struct ReissnerMindlin
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<double(const RMParameters &, const Eigen::Vector2d &)> ForcingTermType;
    typedef std::function<Eigen::Vector2d(const RMParameters &, const Eigen::Vector2d &)> SolutionRotationType;
    typedef std::function<Eigen::Matrix2d(const RMParameters &, const Eigen::Vector2d &)> GradientRotationType;
    typedef std::function<double(const RMParameters &, const Eigen::Vector2d &)> SolutionDisplacementType;

    /// Constructor
    ReissnerMindlin(
                   const DDRCore & ddrcore,               ///< Core for the DDR space sequence
                   const RMParameters & para,             ///< Physical parameters
                   const BoundaryConditions & BC_theta,   ///< Boundary conditions for rotation
                   const BoundaryConditions & BC_u,       ///< Boundary conditions for displacement
                   bool use_threads,                      ///< True for parallel execution, false for sequential execution
                   std::ostream & output = std::cout      ///< Output stream to print status messages
                   );

    /// Assemble the global system    
    void assembleLinearSystem(
                              const ForcingTermType & f,      ///< Forcing term
                              const SolutionRotationType & theta,  ///< Boundary value
                              const GradientRotationType & grad_theta,  ///< Boundary value
                              const SolutionDisplacementType & u  ///< Boundary value
                              );


    /// Returns the dimension of the rotation + displacement space (with BC)
    inline size_t dimensionSpace() const
    {
      return m_excurl.dimension() + m_xgrad.dimension();
    }

    /// Returns the nb of DOFs for BC
    inline size_t nb_bdryDOFs() const
    {
      return m_BC_theta.n_dir_edges() * m_excurl.numLocalDofsEdge() 
             + m_BC_u.n_dir_edges() * m_xgrad.numLocalDofsEdge()
             + m_BC_u.n_dir_vertices() * m_xgrad.numLocalDofsVertex();
    }

    /// Returns the size of the system without BC
    inline size_t sizeSystem() const
    {
      return dimensionSpace() - nb_bdryDOFs();
    }

    /// Returns the location of the unknowns among the DOFs
    inline const std::vector<std::pair<size_t,size_t>> & locUKN() const {
      return m_locUKN;
    }

    /// Create the vector of DOF indices for cell T, which combines the DOFs for the spaces EXcurl and Xgrad
    std::vector<size_t> globalDOFIndices(const Cell &T) const
    {
      std::vector<size_t> I_excurl_T = m_excurl.globalDOFIndices(T);
      std::vector<size_t> I_xgrad_T = m_xgrad.globalDOFIndices(T);
      size_t dim_T = m_excurl.dimensionCell(T.global_index()) + m_xgrad.dimensionCell(T.global_index());
      std::vector<size_t> I_T(dim_T);
      size_t dim_excurl = m_excurl.dimension();
      auto it_I_T = std::copy(I_excurl_T.begin(), I_excurl_T.end(), I_T.begin());
      std::transform(I_xgrad_T.begin(), I_xgrad_T.end(), it_I_T, [&dim_excurl](const size_t & index) { return index + dim_excurl; });
      return I_T;
    }

    /// Returns the parameters
    inline const RMParameters & para() const
    {
      return m_para;
    }
    
    /// Returns the space EXCurl
    inline const EXCurl & exCurl() const
    {
      return m_excurl;
    }

    /// Returns the space XGrad
    inline const XGrad & xGrad() const
    {
      return m_xgrad;
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

    /// Returns the Matrix for BC
    inline const SystemMatrixType & bdryMatrix() const {
      return m_bdryMatrix;
    }

    /// Returns the boundary values
    inline const Eigen::VectorXd & bdryValues() const {
      return m_bdryValues;
    }

    /// Returns the stabilization parameter
    inline const double & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline double & stabilizationParameter() {
      return m_stab_par;
    }

    /// Compute the discrete norm
    double computeNorm(
                       const Eigen::VectorXd & v ///< The vector
                      ) const;

    /// Takes a function dependent on RMParameter and a position x, and returns a function depending only on x (using the parameters of this class)
    template<typename outValue, typename Fct>
    std::function<outValue(const Eigen::Vector2d &)> contractPara(const Fct &F) const 
    {
      std::function<outValue(const Eigen::Vector2d &)> f = [this, &F](const Eigen::Vector2d &x)->outValue { return F(m_para,x);};
      return f;
    }
    

  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_local_contribution(
                                size_t iT,                      ///< Element index
                                const GradientRotationType & grad_theta,      ///< Gradient of rotation (for simply supported BC)
                                const ForcingTermType & f       ///< Forcing term
                                );

    /// Compute the local contribution for low-order stabilisation, jump around edge iE
    Eigen::MatrixXd
    _compute_local_jump_stab(size_t iE);

    /// Assemble the local contribution for the element of index iT into the global system
    /** The strategy to assemble the matrix is as follows. The global matrix has DOFs corresponding to Dirichlet BC, and other DOFs (for short "internal" and "boundary"). The rows and columns corresponding to internal DOFs form the system matrix m_A. The columns, in the global matrix, corresponding to boundary DOFs are used to adjust the RHS of the system to account for BC. During the assembly, we create two sets of triplets: one for the system matrix (of size nb internal DOFs x nb internal DOFs), one for the boundary matrix (of size nb internal DOFs x total number of DOFs); the vector m_bdryValues will be of size "total number of DOFs" so that the BC are already well positioned in the whole set of DOFs. 
    */
    void _assemble_local_contribution(
                                      size_t iT,                                               ///< Element index
                                      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
                                      std::list<Eigen::Triplet<double> > & A1,                 ///< List of triplets for system
                                      Eigen::VectorXd & b1,                                    ///< Vector for RHS for system
                                      std::list<Eigen::Triplet<double> > & A2                  ///< List of triplets for BC
                                      );

    /// Assemble the low-order jump stabilisation (for degree=0)
    void _assemble_local_jump_stab(
                          size_t iE,                                  ///< Edge index
                          const Eigen::MatrixXd & locJ,               ///< Local jumps bilinear form
                          std::list<Eigen::Triplet<double> > & A1,    ///< List of triplets for lo stab
                          std::list<Eigen::Triplet<double> > & A2     ///< List of triplets for BC
                          );

    const DDRCore & m_ddrcore;
    const RMParameters m_para;
    bool m_use_threads;
    std::ostream & m_output;
    EXCurl m_excurl;
    XGrad m_xgrad;
    double m_t;
    double m_E;
    double m_nu;
    BoundaryConditions m_BC_theta;
    BoundaryConditions m_BC_u;
    SystemMatrixType m_bdryMatrix;  // Matrix with columns corresponding to Dirichlet DOFs; multiplied by the boundary values, yields the modification to the system RHS for non-homogeneous BC
    Eigen::VectorXd m_bdryValues;   // Boundary values
    SystemMatrixType m_A;           // Matrix and RHS for system (after removing Dirichlet DOFs)
    Eigen::VectorXd m_b;
    double m_stab_par;
    std::vector<std::pair<size_t,size_t>> m_locUKN;    // Indicates the location, among the DOFs, of the unknowns after the Dirichlet DOFs are removed (the unknowns are between m_locUCK[i].first and m_locUCK[i].first+m_locUKN[i].second for all i)
    Eigen::VectorXi m_DOFtoUKN;       // Maps a dof i to its unknown in the system without BC, or to -1 if the DOF is a Dirichlet DOF
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
  //    Constant

  static ReissnerMindlin::SolutionRotationType
  constant_theta = [](const RMParameters & para, const VectorRd & x) -> VectorRd {
                 return VectorRd(1., 2.);
               };

  static ReissnerMindlin::GradientRotationType
  constant_grad_theta = [](const RMParameters & para, const VectorRd & x) -> Eigen::Matrix2d {
                 return Eigen::Matrix2d::Zero();
               };

  static ReissnerMindlin::SolutionDisplacementType  
  constant_u = [](const RMParameters & para, const VectorRd & x) -> double {
                     return 1.;
                   };

  static ReissnerMindlin::ForcingTermType
  constant_f = [](const RMParameters & para, const VectorRd & x) -> double {
                 return 0.;
               };
 
  //------------------------------------------------------------------------------
  // Polynomial
  
  static ReissnerMindlin::SolutionRotationType
  polynomial_theta = [](const RMParameters & para, const VectorRd & x) -> VectorRd {
                 return VectorRd(
                          pow(x(1),3)*pow(x(1)-1,3)*pow(x(0),2)*pow(x(0)-1,2)*(2*x(0)-1),
                          pow(x(0),3)*pow(x(0)-1,3)*pow(x(1),2)*pow(x(1)-1,2)*(2*x(1)-1)
                          );
               };
               
  static ReissnerMindlin::GradientRotationType
  polynomial_grad_theta = [](const RMParameters & para, const VectorRd & x) -> Eigen::Matrix2d {
                 Eigen::Matrix2d G = Eigen::Matrix2d::Zero();
                 G(0,0) = 2.*pow(x(0),2)*pow(x(1),3)*pow(x(0)-1.0,2)*pow(x(1)-1.0,3)+2.*x(0)*pow(x(1),3)*(2.*x(0)-1.)*pow(x(0)-1.,2)*pow(x(1)-1.0,3)+pow(x(0),2)*pow(x(1),3)*(2.*x(0)-1.)*(2.*x(0)-2.)*pow(x(1)-1.,3);
                 G(0,1) = 3.*pow(x(0)*x(1),2)*(2.*x(0)-1.)*pow(x(0)-1.,2)*pow(x(1)-1.,3)+3.*pow(x(0),2)*pow(x(1),3)*(2.*x(0)-1.)*pow(x(0)-1.,2)*pow(x(1)-1.,2);
                 G(1,0) = 2.*pow(x(1),2)*pow(x(0),3)*pow(x(1)-1.0,2)*pow(x(0)-1.0,3)+2.*x(1)*pow(x(0),3)*(2.*x(1)-1.)*pow(x(1)-1.,2)*pow(x(0)-1.0,3)+pow(x(1),2)*pow(x(0),3)*(2.*x(1)-1.)*(2.*x(1)-2.)*pow(x(0)-1.,3);
                 G(1,1) = 3.*pow(x(1)*x(0),2)*(2.*x(1)-1.)*pow(x(1)-1.,2)*pow(x(0)-1.,3)+3.*pow(x(1),2)*pow(x(0),3)*(2.*x(1)-1.)*pow(x(1)-1.,2)*pow(x(0)-1.,2);

                 return G;
               };               

  static ReissnerMindlin::SolutionDisplacementType  
  polynomial_u = [](const RMParameters & para, const VectorRd & x) -> double {
                 double val = 0;
                 
                 val += (1./3.)*pow(x(0),3)*pow(x(0)-1,3)*pow(x(1),3)*pow(x(1)-1,3);
                 
                 val -= ( 2*pow(para.t,2) / (5*(1.-para.nu)) ) *
                          ( pow(x(1),3)*pow(x(1)-1,3)*x(0)*(x(0)-1)*(5*x(0)*x(0)-5*x(0)+1)
                            + pow(x(0),3)*pow(x(0)-1,3)*x(1)*(x(1)-1)*(5*x(1)*x(1)-5*x(1)+1) );
                 
                 return val;
                 };

  static ReissnerMindlin::ForcingTermType
  polynomial_f = [](const RMParameters & para, const VectorRd & x) -> double {
                double val = 0;
                
                val += 12*x(1)*(x(1)-1)*(5*x(0)*x(0)-5*x(0)+1) * 
                        ( 2*x(1)*x(1)*(x(1)-1)*(x(1)-1) + x(0)*(x(0)-1)*(5*x(1)*x(1)-5*x(1)+1) );
                val += 12*x(0)*(x(0)-1)*(5*x(1)*x(1)-5*x(1)+1) * 
                        ( 2*x(0)*x(0)*(x(0)-1)*(x(0)-1) + x(1)*(x(1)-1)*(5*x(0)*x(0)-5*x(0)+1) );
                
                return val * para.E / (12* (1.-para.nu*para.nu) );
               };

  //------------------------------------------------------------------------------
  // Analytical solution (see paper)
  //  (BC are not zero here)  

  static std::function<double(const VectorRd &)> 
          an_g = [](const VectorRd &x)->double { return sin(PI*x(0))*sin(PI*x(1)); };
          
  static std::function<VectorRd(const VectorRd &)>
          an_GRAD_g = [](const VectorRd & x) -> VectorRd {
                 return VectorRd( PI*cos(PI*x(0))*sin(PI*x(1)), PI*sin(PI*x(0))*cos(PI*x(1)) ); 
               };
               
  static std::function<MatrixRd(const VectorRd &)>
          an_HESS_g = [](const VectorRd & x) -> MatrixRd {
                MatrixRd H = MatrixRd::Zero();
                H.row(0) << -PI*PI*sin(PI*x(0))*sin(PI*x(1)), PI*PI*cos(PI*x(0))*cos(PI*x(1));
                H.row(1) << PI*PI*cos(PI*x(0))*cos(PI*x(1)), -PI*PI*sin(PI*x(0))*sin(PI*x(1));
                return H;
              };              

  static std::function<double(const VectorRd &)> 
          an_LAPL_g = [](const VectorRd &x)->double { return -2*pow(PI,2)*sin(PI*x(0))*sin(PI*x(1)); };

  static std::function<double(const VectorRd &)> 
          an_V = [](const VectorRd &x)->double { return x(0) * exp(-x(0))*cos(x(1)); };

  static std::function<VectorRd(const VectorRd &)>
          an_GRAD_V = [](const VectorRd & x) -> VectorRd {
                 return VectorRd((1.-x(0))*exp(-x(0))*cos(x(1)), -x(0)*exp(-x(0))*sin(x(1))); 
               };
            
  static std::function<MatrixRd(const VectorRd &)>
          an_HESS_V = [](const VectorRd & x) -> MatrixRd {
                MatrixRd H = MatrixRd::Zero();
                H.row(0) << (x(0)-2.)*exp(-x(0))*cos(x(1)), -(1.-x(0))*exp(-x(0))*sin(x(1));
                H.row(1) << -(1.-x(0))*exp(-x(0))*sin(x(1)), -x(0)*exp(-x(0))*cos(x(1));
                return H;
              };              

  static std::function<double(const VectorRd &)> 
          an_LAPL_V = [](const VectorRd &x)->double { return -2.* exp(-x(0))*cos(x(1)); };


  static ReissnerMindlin::SolutionDisplacementType  
          an_v = [](const RMParameters & para, const VectorRd & x) -> double {
                     return pow(para.t,3) * an_V(x/para.t) + an_g(x);
                   };

  static ReissnerMindlin::SolutionRotationType
          an_GRAD_v = [](const RMParameters & para, const VectorRd & x) -> VectorRd {
                 return pow(para.t, 2) * an_GRAD_V(x/para.t) + an_GRAD_g(x);
               };
  
  static ReissnerMindlin::GradientRotationType
          an_HESS_v = [](const RMParameters & para, const VectorRd & x) -> MatrixRd {
                  return para.t * an_HESS_V(x/para.t) + an_HESS_g(x);
              };
          
  static ReissnerMindlin::SolutionDisplacementType  
          an_LAPL_v = [](const RMParameters & para, const VectorRd & x) -> double {
                     return para.t * (an_LAPL_V)(x/para.t) + an_LAPL_g(x);
                   };

  static ReissnerMindlin::SolutionRotationType analytical_theta = an_GRAD_v;
  
  static ReissnerMindlin::GradientRotationType analytical_grad_theta = an_HESS_v;

  static ReissnerMindlin::SolutionDisplacementType  
  analytical_u = [](const RMParameters & para, const VectorRd & x) -> double {
                     return an_v(para, x) - pow(para.t,2) * ( (para.beta0+para.beta1)/para.kappa ) * an_LAPL_v(para, x);
                   };

  static ReissnerMindlin::ForcingTermType
  analytical_f = [](const RMParameters & para, const VectorRd & x) -> double {
                 return (para.beta0+para.beta1) * 4 * pow(PI,4)*sin(PI*x(0))*sin(PI*x(1));
               };

  //------------------------------------------------------------------------------
  // Unkown exact solution, just useful to fix BC and load and plot the approximation solution
  
  static ReissnerMindlin::SolutionRotationType
  ukn_theta = [](const RMParameters & para, const VectorRd & x) -> VectorRd {
                 return VectorRd(0., 0.);
               };

  static ReissnerMindlin::GradientRotationType
  ukn_grad_theta = [](const RMParameters & para, const VectorRd & x) -> Eigen::Matrix2d {
                 return Eigen::Matrix2d::Zero();
               };               

  static ReissnerMindlin::SolutionDisplacementType  
  ukn_u = [](const RMParameters & para, const VectorRd & x) -> double {
                     return x(0);
                   };

  static ReissnerMindlin::ForcingTermType
  ukn_f = [](const RMParameters & para, const VectorRd & x) -> double {
                double val = 1.;
                if (x(0)<.5){
                  val = -1.;
                }
                return val;
               };

  //------------------------------------------------------------------------------
  // Kirchoff limit: solution of E/(12(1-nu^2)) Delta^2 u = g
  //  Useful to evaluate convergence as t->0
  
  static ReissnerMindlin::SolutionRotationType
  kir_theta = [](const RMParameters & para, const VectorRd & x) -> VectorRd {
                double coef = 12.*(1-pow(para.nu,2))/ (para.E * 4.*pow(PI, 4));
                return coef * PI * VectorRd(cos(PI*x(0))*sin(PI*x(1)), sin(PI*x(0))*cos(PI*x(1)));
               };

  static ReissnerMindlin::GradientRotationType
  kir_grad_theta = [](const RMParameters & para, const VectorRd & x) -> Eigen::Matrix2d {
                 double coef = 12.*(1-pow(para.nu,2))/ (para.E * 4.*pow(PI, 4));
                 Eigen::Matrix2d H = Eigen::Matrix2d::Zero();
                 H.row(0) << -PI*PI*sin(PI*x(0))*sin(PI*x(1)), PI*PI*cos(PI*x(0))*cos(PI*x(1));
                 H.row(1) <<  PI*PI*cos(PI*x(0))*cos(PI*x(1)), -PI*PI*sin(PI*x(0))*sin(PI*x(1));
                 return coef * H;
               };               

  static ReissnerMindlin::SolutionDisplacementType  
  kir_u = [](const RMParameters & para, const VectorRd & x) -> double {
                  double coef = 12.*(1-pow(para.nu,2))/ (para.E * 4.*pow(PI, 4));
                  return coef  * sin(PI*x(0))*sin(PI*x(1));
               };

  static ReissnerMindlin::ForcingTermType
  kir_f = [](const RMParameters & para, const VectorRd & x) -> double {
                return sin(PI*x(0))*sin(PI*x(1));
               };


} // end of namespace HArDCore2D

#endif
