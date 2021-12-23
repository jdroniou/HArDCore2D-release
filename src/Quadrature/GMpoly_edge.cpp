#include <GMpoly_edge.hpp>

using namespace HArDCore2D;

// IntegrateEdgeMonomials
MonomialEdgeIntegralsType HArDCore2D::IntegrateEdgeMonomials(const Edge & E, const size_t maxdeg) 
  {
  MonomialEdgeIntegralsType integrals(maxdeg+1, 0.);
  
  // Loop over the monomials
  for (size_t m = 0; m<maxdeg+1; m+=2){
    integrals[m] = E.diam()/((m+1)*std::pow(2,m));     
  }
    
  return integrals;
}

// Check integrals list
MonomialEdgeIntegralsType HArDCore2D::CheckIntegralsDegree(const Edge & E, const size_t degree, const MonomialEdgeIntegralsType & mono_int_map){
  if (mono_int_map.size() >= PolynomialSpaceDimension<Edge>::Poly(degree)){
    return mono_int_map;
  }else{
    return IntegrateEdgeMonomials(E, degree);
  }
}

// GramMatrix
Eigen::MatrixXd HArDCore2D::GramMatrix(const Edge & E, const MonomialScalarBasisEdge & basis1, const MonomialScalarBasisEdge & basis2, MonomialEdgeIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateEdgeMonomials
  MonomialEdgeIntegralsType intmap = CheckIntegralsDegree(E, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = intmap[i+j];
    } // for j
  }   // for i
  
  return gm;
}

// GMDer
Eigen::MatrixXd HArDCore2D::GMDer(const Edge & E, const MonomialScalarBasisEdge & basis1, const MonomialScalarBasisEdge & basis2, MonomialEdgeIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateEdgeMonomials
  MonomialEdgeIntegralsType intmap = CheckIntegralsDegree(E, totaldegree, mono_int_map);
  
  for (size_t i = 1; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = i*intmap[i+j-1];
    } // for j
  }   // for i
  
  return gm;
}


