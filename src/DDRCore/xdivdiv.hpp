#ifndef XDIVDIV_HPP
#define XDIVDIV_HPP

#include <globaldofspace.hpp>
#include <platescore.hpp>
#include <integralweight.hpp>

namespace HArDCore2D
{
  /*!
   *	\addtogroup PlatesCore
   * @{
   */

  /// Discrete Hdivdiv space
    /** At each vertex, the DOFs are the coefficients on the basis "m_basis" (accessible via m_symm_vertex) of Symm.
  On each edge \f$E\f$, the DOFs are \f$\tau_E\f$ on \f$P^{k-3}(E)\f$, then \f$D_{E,\tau}\f$ on \f$P^{k-2}(E)\f$.
  In each element \f$T\f$, the DOFs are \f$\tau_{\mathcal H,T}^c\f$ in \f$\mathcal H^{c,k-1}(T)\f$, and then (if it exists) \f$\tau_{\mathcal H,T}\f$ in \f$\mathcal H^{k-4}(T)\f$. */
  class XDivDiv : public GlobalDOFSpace
  {
  public:
    typedef std::function<Eigen::Matrix2d(const Eigen::Vector2d &)> FunctionType;
    typedef std::function<double(const Eigen::Vector2d &, const Edge &)> EdgeFunctionType;
    typedef std::function<Eigen::Matrix2d(const Eigen::Matrix2d &)> ConstitutiveLawType;

    /// A structure to store the local operators (divdiv and potential)
    struct LocalOperators
    {
      LocalOperators(
		     const Eigen::MatrixXd & _divdiv,     ///< Div-div operator
                     const Eigen::MatrixXd & _divdiv_rhs, ///< Linear form corresponding to the div-div operator
		     const Eigen::MatrixXd & _potential   ///< Potential operator
		     )
	    : divdiv(_divdiv),
              divdiv_rhs(_divdiv_rhs),
	      potential(_potential)
      {
      	// Do nothing
      }

      Eigen::MatrixXd divdiv;
      Eigen::MatrixXd divdiv_rhs;
      Eigen::MatrixXd potential;
    };
    
    /// Basis for the space of symmetric matrices
    struct SymmetricMatrixBasisVertex {
      typedef MatrixRd FunctionValue;

      SymmetricMatrixBasisVertex();      

      inline FunctionValue function(size_t i) const
      {
        assert(i >= 0 && i <= 2);
        return m_basis[i];
      }

      inline size_t size() const {
        return 3;
      }
      
    private:
      std::array<MatrixRd, 3> m_basis;
    };

    /// Constructor
    XDivDiv(const PlatesCore & plates_core, bool use_threads = true, std::ostream & output = std::cout);

    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_plates_core.mesh();
    }
    
    /// Return the polynomial degree. Notice that this is the degree of the complex,
    /// not of the space (the latter being \f$k-1\f$)
    const size_t & degree() const
    {
      return m_plates_core.degree();
    }

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
                                const FunctionType & tau,      ///< The function to interpolate \f$\boldsymbol{\tau}\f$
                                const EdgeFunctionType & Dtau, ///< Given a mesh edge \f$E\in\mathcal{E}_h\f$, this function should return $\partial_{\boldsymbol{t}_E}(\boldsymbol{\tau}_{|E}\boldsymbol{n}_E\cdot\boldsymbol{t}_E) + ({\bf div}\boldsymbol{\tau})_{|E}\cdot\normal_E\f$
                                const int deg_quad = -1        ///< The optional degree of quadrature rules to compute the interpolate. If negative a default value will be used
                                ) const;

    /// Matrix the L2-product for the cell of index iT, weighted by a tensor A.
    Eigen::MatrixXd computeL2Product(
                                     const size_t iT, ///< index of the cell
                                     const ConstitutiveLawType & A = [](const Eigen::Matrix2d tau)->Eigen::Matrix2d { return tau; }, ///< Tensor to transform the potentials, defaults to identity
                                     const double & penalty_factor = 1. ///< pre-factor for stabilisation term
                                     ) const;


    /// Return cell operators for the cell of index iT
    inline const LocalOperators & cellOperators(size_t iT) const
    {
      return *m_cell_operators[iT];
    }

    /// Return edge potential for the edge of index iE
    inline const Eigen::MatrixXd & edgePotential(size_t iE) const
    {
      return m_edge_potentials[iE];
    }

    /// Return edge potential for the edge E
    inline const Eigen::MatrixXd & edgePotential(const Edge & E) const
    {
      return m_edge_potentials[E.global_index()];
    }

    /// Return cell bases for the cell of index iT
    inline const PlatesCore::CellBases & cellBases(size_t iT) const
    {
      return m_plates_core.cellBases(iT);
    }

    /// Return cell bases for cell T
    inline const PlatesCore::CellBases & cellBases(const Cell & T) const
    {
      return m_plates_core.cellBases(T.global_index());
    }
        
    /// Return edge bases for the edge of index iE
    inline const PlatesCore::EdgeBases & edgeBases(size_t iE) const
    {
      return m_plates_core.edgeBases(iE);
    }

    /// Return edge bases for edge E
    inline const PlatesCore::EdgeBases & edgeBases(const Edge & E) const
    {
      return m_plates_core.edgeBases(E.global_index());
    }
    
    /// Return the structure of basis for the symmetric matrices at the vertices
    inline const SymmetricMatrixBasisVertex & SymBasisVertex() const
    {
      return m_symm_vertex;
    }
    
  private:
    LocalOperators _compute_cell_divdiv_potential(size_t iT);
    Eigen::MatrixXd _compute_edge_potential(size_t iE);

    const PlatesCore & m_plates_core;
    bool m_use_threads;
    std::ostream & m_output;

    // Basis for the space of symmetric tensor-valued fields on a vertex
    SymmetricMatrixBasisVertex m_symm_vertex;

    // Containers for local operators
    std::vector<std::unique_ptr<LocalOperators> > m_cell_operators;
    std::vector<Eigen::MatrixXd> m_edge_potentials;
  };
  
} // end of namespace HArDCore2D

#endif
