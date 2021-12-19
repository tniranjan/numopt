#pragma once

#include "problem.h"


namespace numopt {
namespace problems {
/**
 * Describes a simple paraboloid problem with analytic gradients and hessian.
 * The minima is at T(5)
 */ 
class ParaboloidProblem : public ProblemBase {

public:
  ParaboloidProblem() {}
  /** \fn operator()(const VectorX &in) const
   * \brief returns function value at \a in
   * \param in The vector to evaluate
   */
  double operator()(const VectorX &in) const override {
    return (in.squaredNorm() + (5));
  }
  /** 
   * \brief returns the gradient vector at \a in
   * \param in The vector at which gradient is computed
   */ 
  std::pair<VectorX, double> gradient(const VectorX &in) const override {
    const VectorX grad = 2 * in.transpose();
    return std::make_pair(grad, this->operator()(in));
  }
  /**
   * \brief returns hessian matrix at \a in
   * \param in The vector at which the hessian is computed
   */ 
  std::pair<SparseMatrixX, VectorX> hessian(const VectorX &in) const {
    SparseMatrixX H(in.rows(), in.rows());
    H.setIdentity();
    H = 2 * H;
    return std::make_pair(H, gradient(in).first);
  }
};
} // namespace problems
} // namespace numopt