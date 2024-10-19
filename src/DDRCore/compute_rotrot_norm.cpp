#include <compute_rotrot_norm.hpp>

//------------------------------------------------------------------------------

double HArDCore2D::compute_squared_rotrot_norm(
                                               size_t iT,
                                               const Eigen::VectorXd & vT,
                                               const XRot & x_rot,
                                               const XRotRot & x_rotrot
                                               )
{
  auto RvT = x_rotrot.cellOperators(iT).rotor * vT;  
  return RvT.transpose() * x_rot.computeRotorL2Product(iT, 1., 1.) * RvT;
}

//------------------------------------------------------------------------------

double HArDCore2D::compute_rotrot_norm(
                                       const Eigen::VectorXd & v,
                                       const XRot & x_rot,
                                       const XRotRot & x_rotrot
                                       )
{
  Eigen::VectorXd local_rotrot_norms = Eigen::VectorXd::Zero(x_rot.mesh().n_cells());
  
  std::function<void(size_t, size_t)> compute_squared_rotrot_norm_all
    = [&v, &local_rotrot_norms, &x_rot, &x_rotrot](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++) {
        local_rotrot_norms(iT) = compute_squared_rotrot_norm(
                                                             iT,
                                                             x_rotrot.restrictCell(iT, v),
                                                             x_rot,
                                                             x_rotrot
                                                             );
      } // for iT
    };

  parallel_for(x_rot.mesh().n_cells(), compute_squared_rotrot_norm_all, x_rot.useThreads());
  
  return std::sqrt(local_rotrot_norms.sum());
}

