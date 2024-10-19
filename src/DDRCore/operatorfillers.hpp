#ifndef OPERATORFILLERS_HPP
#define OPERATORFILLERS_HPP

#include <xgrad.hpp>
#include <xrotrot.hpp>

namespace HArDCore2D {
  namespace XRotRotDetail {
  std::vector<Eigen::MatrixXd> _fill_potential_operators(size_t iT, const XRotRot * x_rotrot);
  std::vector<Eigen::MatrixXd> _fill_gradient_operators(size_t iT, const XGrad * x_grad);
  } // namespace XRotRotDetail
} // namespace HArDCore2D
#endif
