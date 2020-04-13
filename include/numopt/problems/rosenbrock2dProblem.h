#pragma once

#include "problem.h"

// Rosenbrock problem of 2 vars with analytic df

namespace numopt {
namespace problems {

class Rosenbrock2dProblem : public ProblemBase {
public:
  Rosenbrock2dProblem() {}
  double operator()(const VectorX &in) const override {
    const auto p1 = (in.x() * in.x() - in.y());
    const auto p2 = (in.x() - 1);
    return 100 * p1 * p1 + p2 * p2;
  }

  std::pair<VectorX, double> gradient(const VectorX &in) const override {
    const double dx = 2 * (in.x() - 1) + 400 * in.x() * (in.x() * in.x() - in.y());
    const double dy = -200 * (in.x() * in.x() - in.y());
    const VectorX grad = (Eigen::Matrix<double, 2, 1>() << dx, dy).finished();
    return std::make_pair(grad, this->operator()(in));
  }

  std::pair<SparseMatrixX, VectorX> hessian(const VectorX &in) const {
    MatrixX H(in.rows(), in.rows());
    H(0, 0) = 1200 * in.x() * in.x() - 400 * in.y() + 2;
    H(0, 1) = -400 * in.x();
    H(1, 0) = -400 * in.x();
    H(1, 1) = 200;
    return std::make_pair(H.sparseView(), gradient(in).first);
  }
};
} // namespace problems
} // namespace numopt