# ifndef _GMPOLY_CELL_HPP
# define _GMPOLY_CELL_HPP

/*

  Classes to implement the quadrature-free polynomial integration rules
  
*/

#include <cassert>
#include <cmath>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
#include <unordered_map>

#include <mesh.hpp>
#include <cell.hpp>
#include <edge.hpp>
#include <vertex.hpp>
#include <quadraturerule.hpp>
#include <basis.hpp>
#include <typeinfo>

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


/// Hash function for VectorZd type
struct VecHash
{
  size_t operator()(const VectorZd & p) const
  {
    constexpr size_t b = 100; // max degree must be less than the basis b
    return p(0) + p(1) * b;
  }
};

/// Type for list of integrals of monomials
typedef std::unordered_map<VectorZd, double, VecHash> MonomialCellIntegralsType;
    
 
/// Compute all integrals on a cell of monomials up to a total degree, using vertex values
MonomialCellIntegralsType IntegrateCellMonomials
          (const Cell & T,      ///< Cell
          const size_t maxdeg           ///< Maximal total degree
          );

/// Checks if the degree of an existing list of monomial integrals is sufficient, other re-compute and return a proper list
MonomialCellIntegralsType CheckIntegralsDegree
         (const Cell & T,              ///< Cell
           const size_t degree,         ///< Expected degree
           const MonomialCellIntegralsType & mono_int_map = {}    ///< Existing list, optional
          );

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  TRANSFORMATIONS OF GRAM MATRICES FOR DERIVED BASES
//----------------------------------------------------------------------
//----------------------------------------------------------------------

/// Transforms a Gram Matrix from an ancestor to a family basis
template<typename BasisType>
inline Eigen::MatrixXd transformGM(const Family<BasisType> & family_basis,     ///< Family 
                               const char RC,                ///< R if transformation applied on rows (left), C if applied on columns (right)
                               const Eigen::MatrixXd & anc_GM          ///< Gram matrix of the ancestor basis
                               )
  {
    if (RC=='R'){
      return family_basis.matrix() * anc_GM;
    }else{
      return anc_GM * (family_basis.matrix()).transpose();
    }
  }

/// Transforms a Gram Matrix from an ancestor to a restricted basis
template<typename BasisType>
inline Eigen::MatrixXd transformGM(const RestrictedBasis<BasisType> & restr_basis,     ///< Restricted basis 
                               const char RC,                ///< R if transformation applied on rows (left), C if applied on columns (right)
                               const Eigen::MatrixXd & anc_GM              ///< Gram matrix of the ancestor basis
                               )
  {
    if (RC=='R'){
      return anc_GM.topRows(restr_basis.dimension());
    }else{
      return anc_GM.leftCols(restr_basis.dimension());
    }
  }

/// Transforms a Gram Matrix from an ancestor to a shifted basis
template<typename BasisType>
inline Eigen::MatrixXd transformGM(const ShiftedBasis<BasisType> & shifted_basis,     ///< Shifted basis 
                               const char RC,                ///< R if transformation applied on rows (left), C if applied on columns (right)
                               const Eigen::MatrixXd & anc_GM              ///< Gram matrix of the ancestor basis
                               )
  {
    if (RC=='R'){
      return anc_GM.bottomRows(shifted_basis.dimension());
    }else{
      return anc_GM.rightCols(shifted_basis.dimension());
    }
  }


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  FUNCTIONS TO COMPUTE GRAM MATRICES
//----------------------------------------------------------------------
//----------------------------------------------------------------------


/* BASIC ONES: Monomial, tensorized, and generic template */

/// Computes the Gram Matrix of a pair of local scalar monomial bases
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const MonomialScalarBasisCell & basis1, ///< First basis
                    const MonomialScalarBasisCell & basis2, ///< Second basis
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
/// This overload to simplify the call to GramMatrix in case the two bases are the same
template<typename BasisType>
Eigen::MatrixXd GramMatrix(const Cell& T, const BasisType & basis, MonomialCellIntegralsType mono_int_map = {})
  {
    return GramMatrix(T, basis, basis, mono_int_map);
  };

/// Template to compute the Gram Matrix of any pair of tensorized scalar bases
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const TensorizedVectorFamily<BasisType1, N> & basis1, ///< First basis (rows of the Gram matrix)
                     const TensorizedVectorFamily<BasisType2, N> & basis2,  ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(basis1.dimension(), basis2.dimension());

    Eigen::MatrixXd anc_gm = GramMatrix(T, basis1.ancestor(), basis2.ancestor(), mono_int_map);
    size_t dim1 = anc_gm.rows();
    size_t dim2 = anc_gm.cols();

    for (size_t i=0; i<N; i++){
      gm.block(i*dim1, i*dim2, dim1, dim2) = anc_gm;    
    }
    return gm;
  };

/// Computes the Gram Matrix of a pair of RolyCompl bases
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & basis1,      ///< First basis
                    const RolyComplBasisCell & basis2,      ///< Second basis
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Template to compute the Gram Matrix of a RolyCompl basis and a tensorized scalar basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & rolycompl_basis,             ///< First basis (RolyCompl basis)
                    const TensorizedVectorFamily<BasisType1, N> & tens_family,  ///< Second basis (tensorized basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = rolycompl_basis.dimension();
    size_t dim2 = tens_family.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = rolycompl_basis.max_degree()+tens_family.max_degree();
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t m=0; m<N; m++){
      gm.block(0, m*dim2/N, dim1, dim2/N) = GMRolyComplScalar(T, rolycompl_basis, tens_family.ancestor(), m, intmap);
    }
    return gm;
  };
  
/// Template to compute the Gram Matrix of a tensorized scalar basis and a RolyCompl basis
template<typename BasisType1, size_t N>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const TensorizedVectorFamily<BasisType1, N> & tens_family,  ///< First basis (tensorized basis)
                    const RolyComplBasisCell & rolycompl_basis,             ///< Second basis (RolyCompl basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrix(T, rolycompl_basis, tens_family, mono_int_map).transpose();
  };
  
/// Gram Matrix of a pair of GolyCompl bases
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const GolyComplBasisCell & basis1,      ///< First basis
                    const GolyComplBasisCell & basis2,      ///< Second basis
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
  
/// Computes the Gram Matrix of the mth component of a RolyCompl Basis and a monomial basis
Eigen::MatrixXd GMRolyComplScalar(
                    const Cell & T, ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & rolycompl_basis, ///< First basis
                    const MonomialScalarBasisCell & mono_basis, ///< Second basis
                    const size_t m, ///< Add one to the power of the mth variable
                    MonomialCellIntegralsType mono_int_map ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );

/// Generic template to compute the Gram Matrix of the mth component of a RolyCompl Basis and any basis
template<typename BasisType>
Eigen::MatrixXd GMRolyComplScalar(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const RolyComplBasisCell & basis1, ///< First basis
                     const BasisType & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     MonomialCellIntegralsType mono_int_map ///< list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType::hasAncestor, "No method to compute this Gram matrix of derivatives");
    return transformGM(basis2, 'C', GMRolyComplScalar(T, basis1, basis2.ancestor(), m, mono_int_map) );
  };  
  
/// Determines if the ancestor of a basis will be used to compute a Gram matrix for this basis
template<typename BasisType>
constexpr bool useAncestor()
  {
    // For a given basis, we only use the ancestor if it has an ancestor and has the same rank as its
    //  ancestor. If it has a different rank than its ancestor, it must be dealt with a specific overload

    if constexpr(BasisType::hasAncestor){
      return (BasisType::tensorRank == BasisType::AncestorType::tensorRank);
    } else {
      return false;
    }
  };

/// Generic template to compute the Gram Matrix of any pair of bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(const Cell& T,         ///< Cell to which the basis corresponds
                     const BasisType1 & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,   ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrix(T, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrix(T, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrix(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }

  };


/* Gram matrix of scalar basis with one or two derivatives */

/// Computes the Gram Matrix of a pair of local scalar monomial bases, taking a partial derivative of the first (w.r.t. the homogeneous coordinates, without scaling)
Eigen::MatrixXd GMScalarDerivative(
                  const Cell & T,  ///< Cell to which the basis corresponds
                  const MonomialScalarBasisCell & basis1, ///< First basis
                  const MonomialScalarBasisCell & basis2, ///< Second basis
                  const size_t m, ///< Differentiate basis1 with respect to the mth variable
                  MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                  );
  
/// Computes the Gram Matrix of a pair of local scalar monomial bases, taking partial derivatives of each of them  (w.r.t. the homogeneous coordinates, without scaling)
Eigen::MatrixXd GMScalarDerivative(
                  const Cell & T,  ///< Cell to which the basis corresponds
                  const MonomialScalarBasisCell & basis1, ///< First basis
                  const MonomialScalarBasisCell & basis2, ///< Second basis
                  const size_t m, ///< Differentiate basis1 with respect to the mth variable
                  const size_t l, ///< Differentiate basis2 with respect to the lth variable
                  MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                  );

/// Generic template to compute the Gram Matrix of any pair of scalar bases, taking a partial derivative of the first  (w.r.t. the homogeneous coordinates, without scaling)
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GMScalarDerivative(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                             )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType1::hasAncestor || BasisType2::hasAncestor, "No method to compute this Gram matrix of derivatives");
      
    if constexpr (!BasisType1::hasAncestor && BasisType2::hasAncestor) {
      return transformGM(basis2, 'C', GMScalarDerivative(T, basis1, basis2.ancestor(), m, mono_int_map) );
    } else if constexpr (BasisType1::hasAncestor && !BasisType2::hasAncestor) {
      return transformGM(basis1, 'R', GMScalarDerivative(T, basis1.ancestor(), basis2, m, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), m, mono_int_map) ) );
    }
  };
  
/// Generic template to compute the Gram Matrix of any pair of scalar bases, taking partial derivatives of each of them  (w.r.t. the homogeneous coordinates, without scaling)
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GMScalarDerivative(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,  ///< Second basis (columns of the Gram matrix)
                     const size_t m, ///< Differentiate basis1 with respect to the mth variable
                     const size_t l, ///< Differentiate basis2 with respect to the lth variable
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                             )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(BasisType1::hasAncestor || BasisType2::hasAncestor, "No method to compute this Gram matrix of derivatives");
      
    if constexpr (!BasisType1::hasAncestor && BasisType2::hasAncestor) {
      return transformGM(basis2, 'C', GMScalarDerivative(T, basis1, basis2.ancestor(), m, l, mono_int_map) );
    } else if constexpr (BasisType1::hasAncestor && !BasisType2::hasAncestor) {
      return transformGM(basis1, 'R', GMScalarDerivative(T, basis1.ancestor(), basis2, m, l, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), m, l, mono_int_map) ) );
    }
    
  };
  

/* Gram matrices of vector-valued basis with gradient, curl... */

/// Template to compute the Gram Matrix of a gradient basis and a tensorized scalar basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GradientBasis<BasisType1> & grad_basis, ///< First basis (rows of the Gram matrix) - gradient basis
                     const TensorizedVectorFamily<BasisType2, N> & tens_family,  ///< Second basis (columns of the Gram matrix) - tensorized basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = grad_basis.dimension();
    size_t dim2 = tens_family.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = grad_basis.ancestor().max_degree()+tens_family.max_degree()-1;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t m=0; m<N; m++){
      gm.block(0, m*dim2/N, dim1, dim2/N) = GMScalarDerivative(T, grad_basis.ancestor(), tens_family.ancestor(), m, intmap);
    }
    return gm/T.diam();
  };


/// Template to compute the Gram Matrix of a tensorized scalar basis and a gradient basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
             const Cell& T, ///< Cell to which the basis corresponds
             const TensorizedVectorFamily<BasisType1, N> & tens_family, ///< First basis (rows of the Gram matrix) - gradient basis
             const GradientBasis<BasisType2> & grad_basis,  ///< Second basis (columns of the Gram matrix) - tensorized basis
             MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
             )
  {
    return GramMatrix(T, grad_basis, tens_family, mono_int_map).transpose();
  };


/// Template to compute the Gram Matrix of a gradient basis and another gradient basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const GradientBasis<BasisType1> & grad_basis1, ///< First basis (rows of the Gram matrix) - gradient basis
                     const GradientBasis<BasisType2> & grad_basis2,  ///< Second basis (columns of the Gram matrix) - gradient basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = grad_basis1.dimension();
    size_t dim2 = grad_basis2.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = grad_basis1.ancestor().max_degree()+grad_basis2.ancestor().max_degree()-2;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t m=0; m<2; m++){
      gm += GMScalarDerivative(T, grad_basis1.ancestor(), grad_basis2.ancestor(), m, m, intmap);
    }
    return gm/std::pow(T.diam(), 2);
  };
  
  

/// Generic template to compute the Gram Matrix of a pair of Curl bases
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                             ///< Cell to which the basis corresponds
                    const CurlBasis<BasisType1> & basis1,       ///< First basis
                    const CurlBasis<BasisType2> & basis2,       ///< Second basis
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    // Dimension of the gram matrix
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    size_t totaldegree = basis1.ancestor().max_degree()+basis2.ancestor().max_degree()-2;
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
    
    // Obtain integration data from IntegrateCellMonomials
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
    
    for (size_t m=0; m<2; m++){
      gm += GMScalarDerivative(T, basis1.ancestor(), basis2.ancestor(), m, m, intmap);
    }
    
    return gm/std::pow(T.diam(), 2);
  };
  

/// Template to compute the Gram Matrix of a curl basis and a tensorized scalar basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrix(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const CurlBasis<BasisType1> & curl_basis, ///< First basis (rows of the Gram matrix) - curl basis
                     const TensorizedVectorFamily<BasisType2, 2> & tens_family,  ///< Second basis (columns of the Gram matrix) - tensorized basis
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    size_t dim1 = curl_basis.dimension();
    size_t dim2 = tens_family.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = curl_basis.ancestor().max_degree()+tens_family.max_degree()-1;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    gm.block(0, 0, dim1, dim2/2) = GMScalarDerivative(T, curl_basis.ancestor(), tens_family.ancestor(), 1, intmap);
    gm.block(0, dim2/2, dim1, dim2/2) = -GMScalarDerivative(T, curl_basis.ancestor(), tens_family.ancestor(), 0, intmap);

    return gm/T.diam();
  };


/// Template to compute the Gram Matrix of a tensorized scalar basis and a gradient basis
template<typename BasisType1, typename BasisType2, size_t N>
Eigen::MatrixXd GramMatrix(
             const Cell& T, ///< Cell to which the basis corresponds
             const TensorizedVectorFamily<BasisType1, N> & tens_family, ///< First basis (rows of the Gram matrix) - gradient basis
             const CurlBasis<BasisType2> & curl_basis,  ///< Second basis (columns of the Gram matrix) - tensorized basis
             MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
             )
  {
    return GramMatrix(T, curl_basis, tens_family, mono_int_map).transpose();
  };


/// Generic template to compute the Gram Matrix of a Divergence basis and any other basis
template<typename BasisType1, typename BasisType2>
typename boost::disable_if<boost::is_same<BasisType2, MonomialCellIntegralsType>, Eigen::MatrixXd>::type GramMatrix(
                     const Cell& T,                          ///< Cell to which the basis corresponds
                     const DivergenceBasis<BasisType1> & basis1,   ///< First basis (rows of the Gram matrix)
                     const BasisType2 & basis2,              ///< Second basis (columns of the Gram matrix)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                     )
  {
    return GramMatrixDiv(T, basis1.ancestor(), basis2, mono_int_map);
  };
  
/// Template to compute the Gram Matrix of any basis and a Divergence basis
template<typename BasisType1, typename Basis2>
Eigen::MatrixXd GramMatrix(
                    const Cell & T,                         ///< Cell to which the basis corresponds
                    const BasisType1 & basis1,                  ///< First basis (tensorized basis)
                    const DivergenceBasis<Basis2> & basis2,       ///< Second basis (divergence basis)
                    MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    return GramMatrixDiv(T, basis2.ancestor(), basis1, mono_int_map).transpose();
  };
  
/// Template to compute the Gram Matrix of a Divergence<Tensorized> basis and a monomial scalar basis
template<typename BasisType1>
Eigen::MatrixXd GramMatrixDiv(
                    const Cell & T,                                         ///< Cell to which the basis corresponds
                    const TensorizedVectorFamily<BasisType1, 2> & basis1,   ///< First basis (divergence basis)
                    const MonomialScalarBasisCell & basis2,                 ///< Second basis (monomial scalar basis)
                    MonomialCellIntegralsType mono_int_map = {}                 ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    )
  {
    size_t dim1 = basis1.dimension();
    size_t dim2 = basis2.dimension();
    Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1, dim2);

    // Integrals of monomials
    size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
    MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

    for (size_t i = 0; i < 2; i++) {
      gm.block(i*dim1/2, 0, dim1/2, dim2) = GMScalarDerivative(T, basis1.ancestor(), basis2, i, intmap);
    }
    
    return gm/T.diam();
  };
  
/// Computes the Gram Matrix of a Divergence<RolyCompl> basis and a monomial scalar basis
Eigen::MatrixXd GramMatrixDiv(
                    const Cell & T,                           ///< Cell to which the basis corresponds
                    const RolyComplBasisCell & basis1,        ///< First basis (RolyCompl basis)
                    const MonomialScalarBasisCell & basis2,   ///< Second basis (tensorized basis)
                    MonomialCellIntegralsType mono_int_map = {}   ///< Optional list of integrals of monomials up to the sum of max degree of basis1 and basis2
                    );
  
/// Template to compute the Gram Matrix of the divergence of any basis and any other basis
template<typename BasisType1, typename BasisType2>
Eigen::MatrixXd GramMatrixDiv(
                     const Cell& T, ///< Cell to which the basis corresponds
                     const BasisType1 & basis1, ///< First basis (vector basis)
                     const BasisType2 & basis2,  ///< Second basis (scalar basis)
                     MonomialCellIntegralsType mono_int_map = {} ///< Optional list of integrals of monomials up to the sum of max degree of the bases
                     )
  {
    // If no ancestor is to be used, we shouldn't be in this overload
    static_assert(useAncestor<BasisType1>() || useAncestor<BasisType2>(), "No method to compute this Gram matrix");
      
    if constexpr (!useAncestor<BasisType1>() && useAncestor<BasisType2>()) {
      return transformGM(basis2, 'C', GramMatrixDiv(T, basis1, basis2.ancestor(), mono_int_map) );
    } else if constexpr (useAncestor<BasisType1>() && !useAncestor<BasisType2>()) {
      return transformGM(basis1, 'R', GramMatrixDiv(T, basis1.ancestor(), basis2, mono_int_map) );
    } else {
      return transformGM(basis1, 'R', transformGM(basis2, 'C', GramMatrixDiv(T, basis1.ancestor(), basis2.ancestor(), mono_int_map) ) );
    }
  };


/*@}*/
} // end of namespace HArDCore2D

#endif // end of _GMPOLY_CELL_HPP
