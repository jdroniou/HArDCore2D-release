#ifndef COMPUTE_ROTROT_NORM
#define COMPUTE_ROTROT_NORM

#include <xrot.hpp>
#include <xrotrot.hpp>
#include <parallel_for.hpp>

namespace HArDCore2D {
  /// Compute the squared rotrot norm on the element of index iT
  double compute_squared_rotrot_norm(
                                     size_t iT,
                                     const Eigen::VectorXd & vT,
                                     const XRot & x_rot,
                                     const XRotRot & x_rotrot
                                     );
  
  /// Compute rotrot norm
  double compute_rotrot_norm(
                             const Eigen::VectorXd & v,
                             const XRot & x_rot,
                             const XRotRot & x_rotrot
                             );

} // end of namespace HArDCore2D

#endif
