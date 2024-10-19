#include <operatorfillers.hpp>

using namespace HArDCore2D;

std::vector<Eigen::MatrixXd>
XRotRotDetail::_fill_potential_operators(size_t iT, const XRotRot * x_rotrot)
{
  const Cell & T = *x_rotrot->mesh().cell(iT);
  
  std::vector<Eigen::MatrixXd> potentialOp(T.n_vertices() + 2 * T.n_edges() + 2);

  // Vertex operators
  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    const Vertex & V = *T.vertex(iV);
    
    Eigen::MatrixXd CvV = Eigen::MatrixXd::Zero(1, x_rotrot->dimensionCell(iT));
    CvV(0, x_rotrot->localOffset(T, V)) = 1.;
    
    potentialOp[iV] = CvV;
  } // for iV
  
  // Edge operators
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    size_t offset_PE = T.n_vertices() + 2 * iE;
    size_t offset_RE = offset_PE + 1;

    size_t dim_PkE = PolynomialSpaceDimension<Edge>::Poly(x_rotrot->degree());
    size_t dim_PkmoE = PolynomialSpaceDimension<Edge>::Poly(x_rotrot->degree() - 1);
    
    Eigen::MatrixXd vE = Eigen::MatrixXd::Zero(dim_PkE, x_rotrot->dimensionCell(iT));
    vE.block(0, x_rotrot->localOffset(T, E), dim_PkE, dim_PkE)
      = Eigen::MatrixXd::Identity(dim_PkE, dim_PkE); 
    potentialOp[offset_PE] = vE;

    Eigen::MatrixXd CvE = Eigen::MatrixXd::Zero(dim_PkmoE, x_rotrot->dimensionCell(iT));
    CvE.block(0, x_rotrot->localOffset(T, E) + dim_PkE, dim_PkmoE, dim_PkmoE)
      = Eigen::MatrixXd::Identity(dim_PkmoE, dim_PkmoE);
    potentialOp[offset_RE] = CvE;
  } // for iE

  size_t offset_PT = T.n_vertices() + 2 * T.n_edges();
  size_t offset_RT = offset_PT + 1;
  
  potentialOp[offset_PT] = x_rotrot->cellOperators(iT).potential;
  potentialOp[offset_RT] = x_rotrot->cellRotor(iT);

  return potentialOp;
}

//------------------------------------------------------------------------------

std::vector<Eigen::MatrixXd>
XRotRotDetail::_fill_gradient_operators(size_t iT, const XGrad * x_grad)
{
  const Cell & T = *x_grad->mesh().cell(iT);
  
  std::vector<Eigen::MatrixXd> gradientOp(T.n_vertices() + 2 * T.n_edges() + 2);

  // Vertex components
  for (size_t iV = 0; iV < T.n_vertices(); iV++) {
    gradientOp[iV] = Eigen::MatrixXd::Zero(1, x_grad->dimensionCell(iT));
  } // for iV
  
  for (size_t iE = 0; iE < T.n_edges(); iE++) {
    const Edge & E = *T.edge(iE);
    
    size_t offset_PE = T.n_vertices() + 2 * iE;
    size_t offset_RE = offset_PE + 1;
    
    gradientOp[offset_PE] = x_grad->extendOperator(T, E, x_grad->edgeOperators(E).gradient);
    gradientOp[offset_RE] = Eigen::MatrixXd::Zero(
                                                  PolynomialSpaceDimension<Edge>::Poly(x_grad->degree() - 1),
                                                  x_grad->dimensionCell(iT)
                                                  );
  } // for iE

  size_t offset_PT = T.n_vertices() + 2 * T.n_edges();
  size_t offset_RT = offset_PT + 1;
  
  gradientOp[offset_PT] = x_grad->cellOperators(iT).gradient;
  gradientOp[offset_RT]
    = Eigen::MatrixXd::Zero(
                            PolynomialSpaceDimension<Cell>::Poly(x_grad->degree()),
                            x_grad->dimensionCell(iT)
                            );

  return gradientOp;
}
