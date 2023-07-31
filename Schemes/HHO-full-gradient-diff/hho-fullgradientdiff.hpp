// Implementation of the HHO scheme, with full gradient, for -div(K nabla u)=f with K scalar piecewise function.
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


#ifndef HHOFULLGRADIENTDIFF_HPP
#define HHOFULLGRADIENTDIFF_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <local_static_condensation.hpp>

#include <hhospace.hpp>
#include <BoundaryConditions/BoundaryConditions.hpp>
#include <integralweight.hpp>

/*!
 * @defgroup HHO_fullgradientdiff
 * @brief Implementation of the HHO scheme with full gradient for variable diffusion equations
 */

namespace HArDCore2D
{

  /*!
   * @addtogroup HHO_fullgradientdiff
   * @{
   */

  /// Assemble a diffusion problem 
  struct FullGradientDiffusion
  {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<double(const VectorRd &)> ForcingTermType;
    typedef std::function<double(const VectorRd &)> SolutionType;
    typedef std::function<VectorRd(const VectorRd &)> SolutionGradientType;
    typedef IntegralWeight PermeabilityType;

    /// Constructor
    FullGradientDiffusion(
                   const HHOSpace & hho_space,           ///< HHO space and operators
                   const BoundaryConditions & BC,     ///< Boundary conditions
                   bool use_threads,                  ///< True for parallel execution, false for sequential execution
                   std::ostream & output = std::cout  ///< Output stream to print status messages
                   );

    /// Assemble the global system    
    void assembleLinearSystem(
                              const ForcingTermType & f,      ///< Forcing term
                              const PermeabilityType & kappa,    ///< Permeability
                              const SolutionType & u,          ///< Exact solution for boundary condition
                              Eigen::VectorXd & UDir          ///< Vector filled in by Dirichlets BCs
                              );

    /// Returns the number of statically condensed DOFs (the cell DOFs)
    inline size_t numSCDOFs() const
    {
      return m_hhospace.mesh().n_cells() * m_nloc_sc;
    }

    /// Returns the number of Dirichlet DOFs
    inline size_t numDirDOFs() const
    {
      return m_BC.n_dir_edges()*m_hhospace.numLocalDofsEdge();
    }

    /// Returns the number of DOFs after SC but before eliminating Dirichlet DOFs
    inline size_t numSkeletalDOFs() const
    {
      return m_hhospace.dimension() - numSCDOFs();
    }

    /// Returns the size of the final system, after application of SC and removal of Dirichlet BCs
    inline size_t sizeSystem() const
    {
      return m_hhospace.dimension() - numSCDOFs() - numDirDOFs();
    }

    /// Returns the HHO space
    inline const HHOSpace & hhospace() const
    {
      return m_hhospace;
    }
    
    /// Returns the mesh
    inline const Mesh & mesh() const
    {
      return m_hhospace.mesh();
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

    /// Returns the stabilization parameter (scaling)
    inline const double & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline double & stabilizationParameter() {
      return m_stab_par;
    }

    /// Returns the static condensation recovery operator
    inline const SystemMatrixType & scMatrix() const {
      return m_sc_A;
    }

    /// Returns the static condensation rhs
    inline Eigen::VectorXd & scVector() {
      return m_sc_b;
    }

    /// Compute the discrete energy norm of a family of vectors representing the dofs
    std::vector<double> computeEnergyNorms(
                       const std::vector<Eigen::VectorXd> & list_dofs   ///< The list of vectors representing the dofs
                      ) const;

  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
          _compute_local_contribution(
                                size_t iT,                      ///< Element index
                                const ForcingTermType & f,      ///< Forcing term
                                const PermeabilityType & kappa     ///< Permeability
                                );

    /// Creates the permutation matrix and the global DOFs for the local static condensation
    LocalStaticCondensation _compute_static_condensation(const size_t & iT) const;

    /// Assemble the local contribution for the element of index iT into the global system
    void _assemble_local_contribution(
                                      size_t iT,                                               ///< Element index
                                      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
                                      std::list<Eigen::Triplet<double> > & A1,                 ///< List of triplets for system
                                      Eigen::VectorXd & b1,                                    ///< Vector for RHS for sysem
                                      std::list<Eigen::Triplet<double> > & A2,                 ///< List of triplets for sc
                                      Eigen::VectorXd & b2                                     ///< Vector for RHS for sc
                                      );
    
    const HHOSpace & m_hhospace;
    const BoundaryConditions & m_BC;
    bool m_use_threads;
    std::ostream & m_output;
    const size_t m_nloc_sc;   // Number of statically condensed DOFs in each cell
    SystemMatrixType m_A;   // Matrix and RHS of statically condensed system
    Eigen::VectorXd m_b;
    SystemMatrixType m_sc_A; // Static condensation operator and RHS (to recover statically condensed DOFs)
    Eigen::VectorXd m_sc_b;
    double m_stab_par;
  };

  //------------------------------------------------------------------------------
  // Exact solutions
  //------------------------------------------------------------------------------

  static const double PI = boost::math::constants::pi<double>();
  using std::sin;
  using std::cos;

  //------------------------------------------------------------------------------
  
  static const VectorRd vec_a = VectorRd(1.,2.);
  
  static FullGradientDiffusion::SolutionType
  linear_u = [](const VectorRd & x) -> double {
                 return vec_a.dot(x) + 4.;
               };

  static FullGradientDiffusion::SolutionGradientType  
  linear_gradu = [](const VectorRd & x) -> VectorRd {
                     return vec_a;
                   };

  static FullGradientDiffusion::ForcingTermType
  linear_f = [](const VectorRd & x) -> double {
                 return 0.;
               };

  static FullGradientDiffusion::PermeabilityType
  linear_kappa = FullGradientDiffusion::PermeabilityType(1.);
  

  //------------------------------------------------------------------------------

  static FullGradientDiffusion::SolutionType
  trigonometric_u = [](const VectorRd & x) -> double {
                      return sin(PI*x(0)) * sin(PI*x(1));
                    };

  static FullGradientDiffusion::SolutionGradientType  
  trigonometric_gradu = [](const VectorRd & x) -> VectorRd {
                          return PI * VectorRd(
                                             cos(PI*x(0)) * sin(PI*x(1)) * sin(PI*x(2)),
                                             sin(PI*x(0)) * cos(PI*x(1)) * sin(PI*x(2))
                                             );
                        };

  static FullGradientDiffusion::ForcingTermType
  trigonometric_f = [](const VectorRd & x) -> double {
                      return 2. * std::pow(PI, 2) * trigonometric_u(x);
                    };

  static FullGradientDiffusion::PermeabilityType
  trigonometric_kappa = FullGradientDiffusion::PermeabilityType(1.);


} // end of namespace HArDCore2D

#endif
