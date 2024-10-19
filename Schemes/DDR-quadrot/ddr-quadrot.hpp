// 
// Solver for a quad-rot problem
// Author: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
//


#ifndef QUADROT_HPP
#define QUADROT_HPP

#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <mesh.hpp>
#include <BoundaryConditions/BoundaryConditions.hpp>   // To re-order the boundary edges and vertices
#include <BoundaryConditions/BChandlers.hpp>   

#include <xgrad.hpp>
#include <xrotrot.hpp>
#include <xrot.hpp>

/*!
 * @defgroup DDR_quadrot
 * @brief Implementation of the DDR scheme for the quad-rot problem
 */

namespace HArDCore2D
{

  /*!
   * @addtogroup DDR_quadrot
   * @{
   */

  struct QuadRot {
    typedef Eigen::SparseMatrix<double> SystemMatrixType;
    
    typedef std::function<VectorRd(const VectorRd &)> PotentialType;
    typedef std::function<double(const VectorRd &)> RotorType;
    typedef std::function<VectorRd(const VectorRd &)> RotRotType;
    typedef std::function<double(const VectorRd &)> LagrangeMultiplierType;
    typedef std::function<VectorRd(const VectorRd &)> ForcingTermType;
    
    /// Constructor
    QuadRot(
	    const DDRCore & ddrcore,          ///< Core for the DDR spaces sequence
            const BoundaryConditions & bc_u,  ///< Boundary condition for the potential
            const BoundaryConditions & bc_p,  ///< Boundary condition for the Lagrange multiplier
	    bool use_threads = true,          ///< True for parallel execution, false for sequential execution
	    std::ostream & output = std::cout ///< Output stream to print status messages
	    );

    /// Returns the space XGrad
    inline const XGrad & xGrad() const
    {
      return m_xgrad;
    }

    /// Returns the space XRotRot
    inline const XRotRot & xRotRot() const
    {
      return m_xrotrot;
    }

    /// Returns the space XGrad
    inline const XRot & xRot() const
    {
      return m_xrot;
    }
    
    /// Returns the stabilization parameter
    inline const double & stabilizationParameter() const {
      return m_stab_par;
    }

    /// Returns the stabilization parameter
    inline double & stabilizationParameter() {
      return m_stab_par;
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

    /// Create the vector of DOF indices for the cell T, which combines the DOFs for the spaces XRotRot and XGrad
    std::vector<size_t> globalDOFIndices(const Cell & T) const;

    /// Returns the dimension of the local potential + Lagrange multiplier space on the element of index iT
    inline size_t dimensionLocalSpace(size_t iT) const
    {
      return m_xrotrot.dimensionCell(iT) + m_xgrad.dimensionCell(iT);
    }

    /// Returns the dimension of the potential + Lagrange multiplier space
    inline size_t dimensionSpace() const
    {
      return m_xrotrot.dimension() + m_xgrad.dimension();
    }

    /// Returns the number of boundary DOFs
    inline size_t numBoundaryDofs() const
    {
      return m_bc_u.n_dir_vertices() * m_xrotrot.numLocalDofsVertex()
        + m_bc_u.n_dir_edges() * m_xrotrot.numLocalDofsEdge()
        + m_bc_p.n_dir_vertices() * m_xgrad.numLocalDofsVertex()
        + m_bc_p.n_dir_edges() * m_xgrad.numLocalDofsEdge();
    }
    
    /// Returns the dimension of the linear system
    inline size_t dimensionLinearSystem() const
    {
      return dimensionSpace() - numBoundaryDofs();
    }

    /// DOF-to-unknown mapping
    inline int dofToUnknown(size_t i) const
    {
      return m_dof_to_unknown(i);
    }

    /// Assemble linear system and returns a lifting of the boundary condition
    /// (which coincides with the interpolate of the exact solution when the
    /// latter is passed as an argument)
    Eigen::VectorXd assembleLinearSystem(
                                         const ForcingTermType & f,       ///< Forcing term
                                         const PotentialType & u,         ///< Solution
                                         const RotorType & rot_u,         ///< Rotor of the exact solution
                                         const LagrangeMultiplierType & p ///< Lagrange multiplier
                                         );
  private:
    /// Compute the local contribution for the element of index iT
    std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    _compute_local_contribution(
				size_t iT,                ///< Element index
				const ForcingTermType & f ///< Forcing term
				);

    /// Assemble the local contribution for the element of index iT
    void _assemble_local_contribution(
                                      size_t iT,                                               ///< Element index
                                      const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT, ///< Local contribution
                                      std::list<Eigen::Triplet<double> > & triplets_A,         ///< List of triplets for the system matrix
                                      Eigen::VectorXd & rsh_b,                                 ///< System vector
                                      std::list<Eigen::Triplet<double> > & triplets_boundary_A ///< List of triplets for the portion of the matrix associated with boundary values
                                      );
    
    const DDRCore & m_ddrcore;
    BoundaryConditions m_bc_u;
    BoundaryConditions m_bc_p;
    bool m_use_threads;
    std::ostream & m_output;
    double m_stab_par;
    XGrad m_xgrad;
    XRotRot m_xrotrot;
    XRot m_xrot;
    SystemMatrixType m_A;
    SystemMatrixType m_B;
    Eigen::VectorXd m_b;
    // Indicates the location, among the DOFs, of the unknowns after the Dirichlet DOFs are removed
    // (m_unknowns_stride[i].first = starting index, m_unknowns_stride[i].second = number of unknowns)
    std::vector<std::pair<size_t, size_t>> m_unknowns_stride;
    Eigen::VectorXi m_dof_to_unknown;
  };

  /*!
   * @}
   */

  //------------------------------------------------------------------------------
  // Exact solutions
  //------------------------------------------------------------------------------

  static const double PI = boost::math::constants::pi<double>();
  using std::sin;
  using std::cos;

  //------------------------------------------------------------------------------
  // Linear

  static QuadRot::PotentialType
  linear_u = [](const VectorRd & x) -> VectorRd {
    VectorRd u_x;
    u_x << -(x(1) - 0.5), (x(0) - 0.5);
    return u_x;
  };

  static QuadRot::RotorType
  linear_rot_u = [](const VectorRd & x) -> double {
    return 2.;
  };

  static QuadRot::RotRotType
  linear_rotrot_u = [](const VectorRd & x)->VectorRd {
    return VectorRd::Zero();
  };
  
  static QuadRot::LagrangeMultiplierType
  linear_p = [](const VectorRd & x) -> double {
    return 0.;
  };

  static QuadRot::ForcingTermType
  linear_f = [](const VectorRd & x) -> VectorRd {
    VectorRd f_x;
    f_x << 0, 0;
    return f_x;
  };

  //------------------------------------------------------------------------------
  // Quadratic

  static QuadRot::PotentialType
  quadratic_u = [](const VectorRd & x) -> VectorRd {
    VectorRd u_x;
    u_x << -pow(x(1) - 0.5, 2), pow(x(0) - 0.5, 2);
    return u_x;
  };

  static QuadRot::RotorType
  quadratic_rot_u = [](const VectorRd & x) -> double {
    return 2. * ( x(0) + x(1) - 1 );
  };

  static QuadRot::RotRotType
  quadratic_rotrot_u = [](const VectorRd & x)->VectorRd {
    VectorRd rotrot_u_x;
    rotrot_u_x << 2., -2.;
    return rotrot_u_x;
  };
  
  static QuadRot::LagrangeMultiplierType
  quadratic_p = [](const VectorRd & x) -> double {
    return 0.;
  };

  static QuadRot::ForcingTermType
  quadratic_f = [](const VectorRd & x) -> VectorRd {
    VectorRd f_x;
    f_x << 0., 0.;
    return f_x;
  };

  //------------------------------------------------------------------------------
  // Cubic
  
  static QuadRot::PotentialType
  cubic_u = [](const VectorRd & x) -> VectorRd {
    VectorRd u_x;
    u_x << -pow(x(1) - 0.5, 3), pow(x(0) - 0.5, 3);
    return u_x;
  };

  static QuadRot::RotorType
  cubic_rot_u = [](const VectorRd & x) -> double {
    return 3. * ( pow(x(0) - 0.5, 2) + pow(x(1) - 0.5, 2) );
  };

  static QuadRot::RotRotType
  cubic_rotrot_u = [](const VectorRd & x)->VectorRd {
    VectorRd rotrot_u_x;
    rotrot_u_x << 6. * (x(1) - 0.5), -6. * (x(0) - 0.5);
    return rotrot_u_x;
  };
  
  static QuadRot::LagrangeMultiplierType
  cubic_p = [](const VectorRd & x) -> double {
    return 0.;
  };

  static QuadRot::ForcingTermType
  cubic_f = [](const VectorRd & x) -> VectorRd {
    VectorRd f_x;
    f_x << 0., 0.;
    return f_x;
  };
  
  //------------------------------------------------------------------------------
  // Quartic

  static QuadRot::PotentialType
  quartic_u = [](const VectorRd & x) -> VectorRd {
    VectorRd u_x;
    u_x << -pow(x(1) - 0.5, 4), pow(x(0) - 0.5, 4);
    return u_x;
  };

  static QuadRot::RotorType
  quartic_rot_u = [](const VectorRd & x) -> double {
    return 4. * ( pow(x(0)- 0.5, 3) + pow(x(1) - 0.5, 3) );
  };

  static QuadRot::RotRotType
  quartic_rotrot_u = [](const VectorRd & x)->VectorRd {
    VectorRd rotrot_u_x;
    rotrot_u_x << 12. * pow(x(1) - 0.5, 2),
      -12. * pow(x(0) - 0.5 , 2);
    return rotrot_u_x;
  };

  // rot rot rot (u) = -24. * ( x(0) + x(1) - 1 )
  
  static QuadRot::LagrangeMultiplierType
  quartic_p = [](const VectorRd & x) -> double {
    return x(0) * (1. - x(0)) * x(1) * (1. - x(1));
  };

  static QuadRot::ForcingTermType
  quartic_f = [](const VectorRd & x) -> VectorRd {
    VectorRd f_x;
    f_x << -24. + (1. - 2. * x(0)) * x(1) * (1. - x(1)),
      24. + x(0) * (1. - x(0)) * (1. - 2. * x(1));
    return f_x;
  }; 
  
  //------------------------------------------------------------------------------
  // Trigonometric

  static QuadRot::PotentialType
  trigonometric_u = [](const VectorRd & x) -> VectorRd {
    VectorRd u_x;
    u_x << -sin(PI * (x(0) + x(1))), sin(PI * (x(0) + x(1)));
    return u_x;
  };

  static QuadRot::RotorType
  trigonometric_rot_u = [](const VectorRd & x) -> double {
    return 2. * PI * cos(PI * (x(0) + x(1)));
  };

  static QuadRot::RotRotType
  trigonometric_rotrot_u = [](const VectorRd & x) -> VectorRd {
    VectorRd rotrot_u_x;
    rotrot_u_x << -sin(PI * (x(0) + x(1))), sin(PI * (x(0) + x(1)));
    return pow(PI, 2) * rotrot_u_x;
  };

  static QuadRot::RotorType
  trigonometric_rotrotrot_u = [](const VectorRd & x) -> double {
    return 4. * pow(PI, 3) * cos(PI * (x(0) + x(1)));
  };
  
  static QuadRot::LagrangeMultiplierType
  trigonometric_p = [](const VectorRd & x) -> double {
    return sin(PI * x(0)) * sin(PI * x(1));
  };

  static QuadRot::ForcingTermType
  trigonometric_f = [](const VectorRd & x) -> VectorRd {
    VectorRd f_x;
    f_x << -4. * pow(PI, 4) * sin(PI * (x(0) + x(1))) + PI * cos(PI * x(0)) * sin(PI * x(1)),
      4. * pow(PI, 4) * sin(PI * (x(0) + x(1))) + PI * sin(PI * x(0)) * cos(PI * x(1));
    return f_x;
  };  

} // end of namespace HArDCore2D

#endif
