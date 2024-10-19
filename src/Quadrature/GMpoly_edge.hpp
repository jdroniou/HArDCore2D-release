# ifndef _GMPOLY_EDGE_HPP
# define _GMPOLY_EDGE_HPP

/*

  Classes to implement the quadrature-free polynomial integration rules
  
*/

#include <cassert>
#include <cmath>

#include <mesh.hpp>
#include <quadraturerule.hpp>
#include <basis.hpp>
#include <GMpoly_cell.hpp>

namespace HArDCore2D {

/*!
*	@addtogroup Quadratures
* @{
*/


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// INTEGRALS OF MONOMIALS
//----------------------------------------------------------------------
//----------------------------------------------------------------------


/// Type for list of edge integrals of monomials
typedef std::vector<double> MonomialEdgeIntegralsType;

/// Compute all integrals of edge monomials up to a total degree
MonomialEdgeIntegralsType IntegrateEdgeMonomials
          (const Edge & E,      ///< Edge
          const size_t maxdeg   ///< Maximal total degree
          );


/// Checks if the degree of an existing list of monomial integrals is sufficient, other re-compute and return a proper list
MonomialEdgeIntegralsType CheckIntegralsDegree
          (const Edge & E,              ///< Edge
           const size_t degree,         ///< Expected degree
           const MonomialEdgeIntegralsType & mono_int_map = {}    ///< Existing list, optional
          );

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  FUNCTIONS TO COMPUTE GRAM MATRICES
//----------------------------------------------------------------------
//----------------------------------------------------------------------

//
/* SCALAR bases */

/// Computes the Gram Matrix of a pair of local scalar monomial bases
Eigen::MatrixXd GramMatrix(
                    const Edge & E,                         ///< Edge to which the basis corresponds
                    const MonomialScalarBasisEdge & basis1, ///< First basis
                    const MonomialScalarBasisEdge & basis2, ///< Second basis
                    MonomialEdgeIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Generic template to compute the Gram Matrix of any pair of bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(const Edge& E,         ///< Edge to which the basis corresponds
                     const BasisType1 & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,   ///< Second basis (columns of the Gram matrix)
                     MonomialEdgeIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrix(E, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrix(E, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrix(E, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }

  };

/// GramMatrix of two tensorized families on the edge
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
            const Edge& E, ///< Edge to which the basis corresponds
            const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix)
            const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix)
            MonomialEdgeIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
    )
    {
        Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(basis1.dimension(), basis2.dimension());

        Eigen::MatrixXd anc_gm = GramMatrix(E, basis1.ancestor(), basis2.ancestor(), mono_int_map);
        size_t dim1 = anc_gm.rows();
        size_t dim2 = anc_gm.cols();

        for (size_t i=0; i<N; i++){
            gm.block(i*dim1, i*dim2, dim1, dim2) = anc_gm;
        }
        return gm;
    }

/// This overload to simplify the call to GramMatrix in case the two bases are the same
template<typename BasisType>
Eigen::MatrixXd GramMatrix(const Edge& E, const BasisType & basis, MonomialEdgeIntegralsType mono_int_map = {})
  {
    return GramMatrix(E, basis, basis, mono_int_map);
  };
 
/* GRADIENT-SCALAR bases */ 
 
/// Computes the Gram Matrix of the derivative of a monomial basis with another monomial basis
Eigen::MatrixXd GMDer(
                  const Edge & E,  ///< Edge to which the basis corresponds
                  const MonomialScalarBasisEdge & basis1, ///< First basis, to be differentiated
                  const MonomialScalarBasisEdge & basis2, ///< Second basis
                  MonomialEdgeIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                  );

/// Generic template for GMDer with derived bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GMDer(const Edge& E,         ///< Edge to which the basis corresponds
                     const BasisType1 & basis1,   ///< First basis (rows of the Gram matrix), to be differentiated
                     const BasisType2 & basis2,   ///< Second basis (columns of the Gram matrix)
                     MonomialEdgeIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GMDer(E, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GMDer(E, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GMDer(E, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }

  };


/// Computes the Gram Matrix of a gradient basis (considering the tangential gradient as a scalar) and a scalar basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                    const Edge & E,                         ///< Edge to which the basis corresponds
                    const GradientBasis<BasisType1> & basis1, ///< Gradient basis
                    const BasisType2 & basis2, ///< Scalar basis
                    MonomialEdgeIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GMDer(E, basis1.ancestor(), basis2, mono_int_map)/E.diam();
  };
  
/// Computes the Gram Matrix of a scalar basis and a gradient basis (considering the tangential gradient as a scalar)
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                    const Edge & E,                         ///< Edge to which the basis corresponds
                    const BasisType1 & basis1, ///< Scalar basis
                    const GradientBasis<BasisType2> & basis2, ///< Gradient basis
                    MonomialEdgeIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return (GMDer(E, basis2.ancestor(), basis1, mono_int_map).transpose())/E.diam();
  };




/*@}*/
} // end of namespace HArDCore2D

#endif // end of _GMPOLY_EDGE_HPP
