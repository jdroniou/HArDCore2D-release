
#include <GMpoly_cell.hpp>

using namespace HArDCore2D;

// IntegrateCellMonomials
MonomialCellIntegralsType HArDCore2D::IntegrateCellMonomials(const Cell & T, const size_t maxdeg) {
  // Number of monomials
  const size_t n = PolynomialSpaceDimension<Cell>::Poly(maxdeg);
  
  std::vector<double> ints (n,0.);
  MonomialCellIntegralsType powersmap;
  
  // Dimension assumed
  const size_t d = 2;
  
  // Coordinate transformation data
  VectorRd xT = T.center_mass();
  double hT = T.diam();
  
  // Create powers and degrees of all monomials
  std::vector<VectorZd> powers = MonomialPowers<Cell>::complete(maxdeg);
  std::vector<size_t> degrees; degrees.reserve(n);
  for (size_t m = 0; m < n; m++) {
    degrees.push_back(powers[m].sum());
  } // m (monomial)
  
  // Loop over each edge
  for (size_t iE = 0; iE < T.n_edges(); iE++)
  {
    Edge *E = T.edge(iE);
    // Outer unit normal and arbitrary point on face
    VectorRd nTE = T.edge_normal(iE);
    VectorRd xE = E->center_mass();
        
    // The integrals of each monomial over just this edge
    std::vector<double> ints_e_i (n,0.);
    
    // Loop over all the monomials
    for (size_t m = 0; m < n; m++)
    {
      // Loop over this edge's vertices
      for (size_t iV = 0; iV < 2; iV++)
      {
        VectorRd xV = E->vertex(iV)->coords();
        VectorRd var = (xV - xT)/hT;
        double val_V = std::pow(var.x(),powers[m](0)) * std::pow(var.y(),powers[m](1));
        ints_e_i[m] += abs( (E->tangent()).dot(xV-xE) ) * val_V; // abs(...) computes the norm of xV-xE
      } // for iV
      // Gradient correction: linear combination of up to three face integrals of one lower degree
      for (size_t ip = 0; ip < d; ip++)
      {
        if (powers[m](ip) > 0){
          VectorZd powers_diff = VectorZd::Zero(d);
          powers_diff(ip) = -1;
          int im = std::find(powers.begin(),powers.end(), powers[m]+powers_diff) - powers.begin();
          ints_e_i[m] += (xE(ip)-xT(ip))/hT * powers[m](ip) * ints_e_i[im];
        }
      }
      ints_e_i[m] /= d-1+degrees[m];
      ints[m] += nTE.dot(xE-xT)*ints_e_i[m];
    } // for m
    
  } // for iE
  
  for (size_t m = 0; m < n; m++) {
    ints[m] /= d+degrees[m];
    powersmap[powers[m]] = ints[m];
  }
  
  return powersmap;
}


// Check integrals list
MonomialCellIntegralsType HArDCore2D::CheckIntegralsDegree(const Cell & T, const size_t degree, const MonomialCellIntegralsType & mono_int_map){
  if (mono_int_map.size() >= PolynomialSpaceDimension<Cell>::Poly(degree)){
    return mono_int_map;
  }else{
    return IntegrateCellMonomials(T, degree);
  }
}


// GramMatrix
Eigen::MatrixXd HArDCore2D::GramMatrix(const Cell & T, const MonomialScalarBasisCell & basis1, const MonomialScalarBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = intmap.at(basis1.powers(i) + basis2.powers(j));
    } // for j
  }   // for i
  
  return gm;
}

// GramMatrix for RolyCompl
Eigen::MatrixXd HArDCore2D::GramMatrix(const Cell & T, const RolyComplBasisCell & basis1, const RolyComplBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t k = 0; k < 2; k++) {
    VectorZd powers = VectorZd::Zero(2);
    powers(k) = 2;
    for (size_t i = 0; i < dim1; i++) {
      for (size_t j = 0; j < dim2; j++) {
          gm(i,j) += intmap.at(powers + basis1.powers(i) + basis2.powers(j));
      } // for j
    }   // for i
  }     // for k
  
  return gm;
}

// Computes the Gram Matrix of a pair of GolyCompl bases
Eigen::MatrixXd HArDCore2D::GramMatrix(const Cell & T, const GolyComplBasisCell & basis1, const GolyComplBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  return GramMatrix(T, *basis1.rck(), *basis2.rck(), mono_int_map);
}


// GMScalarDerivative, one derivative
Eigen::MatrixXd HArDCore2D::GMScalarDerivative(const Cell & T, const MonomialScalarBasisCell & basis1, const MonomialScalarBasisCell & basis2, const size_t m, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    if(basis1.powers(i)(m) > 0) {
      VectorZd powers1 = basis1.powers(i);
      powers1(m) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        gm(i,j) = basis1.powers(i)(m) * intmap.at(powers1 + basis2.powers(j));
      } // for j
    }   // if
  }     // for i
  
  return gm;
}

// GMScalarDerivative, two derivatives
Eigen::MatrixXd HArDCore2D::GMScalarDerivative(const Cell & T, const MonomialScalarBasisCell & basis1, const MonomialScalarBasisCell & basis2, const size_t m, const size_t l, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-2;
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);
  
  // Obtain integration data from IntegrateCellMonomials
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    if(basis1.powers(i)(m) > 0) {
      VectorZd powers1 = basis1.powers(i);
      powers1(m) -= 1;
      for (size_t j = 0; j < dim2; j++) {
        if (basis2.powers(j)(l) > 0){
          VectorZd powers2 = basis2.powers(j);
          powers2(l) -= 1;
          gm(i,j) = basis1.powers(i)(m) * basis2.powers(j)(l) * intmap.at(powers1 + powers2);
        }
      } // for j
    }   // if
  }     // for i
  
  return gm;
}
 
// GMRolyScalar, the mth component of RolyCompl and the scalar ancestor of a tensorized basis
Eigen::MatrixXd HArDCore2D::GMRolyComplScalar(const Cell & T, const RolyComplBasisCell & rolycompl_basis, const MonomialScalarBasisCell & mono_basis, const size_t m, MonomialCellIntegralsType mono_int_map)
{
  // Dimension of the gram matrix
  size_t dim1 = rolycompl_basis.dimension();
  size_t dim2 = mono_basis.dimension();
  Eigen::MatrixXd gm = Eigen::MatrixXd::Zero(dim1,dim2);

  // Integrals of monomials
  size_t totaldegree = rolycompl_basis.max_degree()+mono_basis.max_degree();
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);

  VectorZd powers = VectorZd::Zero(2);
  powers(m) = 1;
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = intmap.at(powers + rolycompl_basis.powers(i) + mono_basis.powers(j));
    }   // for j
  }     // for i
  
  return gm;
}


/// Computes the Gram Matrix of a Divergence<RolyCompl> basis and a monomial scalar basis
Eigen::MatrixXd HArDCore2D::GramMatrixDiv(const Cell & T, const RolyComplBasisCell & basis1, const MonomialScalarBasisCell & basis2, MonomialCellIntegralsType mono_int_map)
{
  size_t dim1 = basis1.dimension();
  size_t dim2 = basis2.dimension();
  Eigen::MatrixXd gm(dim1, dim2);
  
  // Integrals of monomials
  size_t totaldegree = basis1.max_degree()+basis2.max_degree()-1;
  MonomialCellIntegralsType intmap = CheckIntegralsDegree(T, totaldegree, mono_int_map);
  
  for (size_t i = 0; i < dim1; i++) {
    for (size_t j = 0; j < dim2; j++) {
      gm(i,j) = (2 + basis1.powers(i).sum()) * intmap.at(basis1.powers(i) + basis2.powers(j));
    } // for j
  }   // for i
  
  return gm / T.diam();
};
   

