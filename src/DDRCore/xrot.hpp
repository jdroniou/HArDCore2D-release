#ifndef XROT_HPP
#define XROT_HPP

#include <globaldofspace.hpp>
#include <ddrcore.hpp>
#include <integralweight.hpp>
#include <xrotrot.hpp>

namespace HArDCore2D
{
  /*!
   *  \addtogroup DDRCore
   * @{
   */

  /// Discrete H1 space: local operators, L2 product and global interpolator
  class XRot : public GlobalDOFSpace
  {
  public:
    typedef std::function<double(const Eigen::Vector2d &)> FunctionType;

    /// A structure to store local operators (vector rotor and potential)
    struct LocalOperators
    {
      LocalOperators(
                     const Eigen::MatrixXd & _rotor,     ///< Vector rotor operator
                     const Eigen::MatrixXd & _rotor_rhs, ///< Vector rotor right-hand side
                     const Eigen::MatrixXd & _potential  ///< Potential operator
                     )
        : rotor(_rotor),
	  rotor_rhs(_rotor_rhs),
          potential(_potential)
      {
        // Do nothing
      }
      
      Eigen::MatrixXd rotor;
      Eigen::MatrixXd rotor_rhs;
      Eigen::MatrixXd potential;
    };
    
    /// Constructor
    XRot(const DDRCore & ddr_core, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_ddr_core.mesh();
    }

    /// Return true if we use thread-based parallelism
    const bool useThreads() const
    {
      return m_use_threads;
    }
    
    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_ddr_core.degree();
    }
    
    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
                                const FunctionType & q, ///< The function to interpolate
                                const int deg_quad = -1 ///< The optional degre of quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
                                ) const;

    /// Return edge potential for the edge of index iE
    inline const Eigen::MatrixXd & edgePotential(size_t iE) const
    {
      return *m_edge_potentials[iE];
    }

    /// Return edge potential for the edge E
    inline const Eigen::MatrixXd & edgePotential(const Edge& E) const
    {
      return *m_edge_potentials[E.global_index()];
    }

    /// Return cell operators for the cell of index iT
    inline const LocalOperators & cellOperators(size_t iT) const
    {
      return *m_cell_operators[iT];
    }

    /// Return cell operators for cell T
    inline const LocalOperators & cellOperators(const Cell & T) const
    {
      return *m_cell_operators[T.global_index()];
    }

    /// Return cell bases for the cell of index iT
    inline const DDRCore::CellBases & cellBases(size_t iT) const
    {
      return m_ddr_core.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const DDRCore::CellBases & cellBases(const Cell & T) const
    {
      return m_ddr_core.cellBases(T.global_index());
    }
    
    /// Return edge bases for the edge of index iE
    inline const DDRCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_ddr_core.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const DDRCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_ddr_core.edgeBases(E.global_index());
    }

    /// Compute rotor L2-product
    Eigen::MatrixXd computeRotorL2Product(
					  const size_t iT,                         ///< Index of the cell
					  const double & penalty_factor_cell = 1., ///< Coefficient for the cell-based term
                                          const double & penalty_factor_edge = 1.  ///< Coefficient for the edge-based term
                                          
					  ) const;

    /// Compute rotor L2-norm
    double computeRotorL2Norm(const Eigen::VectorXd & v) const;
    
  private:    
    Eigen::MatrixXd _compute_edge_potential(size_t iE);
    LocalOperators  _compute_cell_operators(size_t iT);

    const DDRCore & m_ddr_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Containers for local operators
    std::vector<std::unique_ptr<Eigen::MatrixXd> > m_edge_potentials;
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
  };
   
} // end of namespace HArDCore2D

#endif
