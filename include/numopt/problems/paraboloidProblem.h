#pragma once

#include "problem.h"

// paraboloid problem with analytic df

namespace numopt {
namespace problems {

class ParaboloidProblem : public ProblemBase {
public:
  ParaboloidProblem() {}
  double operator()(const VectorX &in) const override {
    return (in.squaredNorm() + (5));
  }
  
  std::pair<VectorX, double> gradient(const VectorX &in) const override {
    const VectorX grad = 2 * in.transpose();
    return std::make_pair(grad, this->operator()(in));
  }

  std::pair<SparseMatrixX, VectorX> hessian(const VectorX &in) const {
    SparseMatrixX H(in.rows(), in.rows());
    H.setIdentity();
    H = 2 * H;
    return std::make_pair(H, gradient(in).first);
  }
};
} // namespace problems
} // namespace numopt