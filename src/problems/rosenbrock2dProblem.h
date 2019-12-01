#pragma once

#include "problem.h"

// Rosenbrock problem of 2 vars with analytic df

namespace numopt {
namespace problems {
template <typename _Functor>
class Rosenbrock2dProblem : public Problem<_Functor> {

public:
  typedef _Functor TFunctor;
  typedef typename TFunctor::InputType InputType;
  typedef typename TFunctor::ValueType ValueType;
  typedef typename TFunctor::JacobianType JacobianType;
  typedef typename TFunctor::HessianType HessianType;
  Rosenbrock2dProblem() : Problem<TFunctor>::Problem() {}
  Rosenbrock2dProblem(TFunctor &functor)
      : Problem<TFunctor>::Problem(functor) {}

  std::pair<JacobianType, ValueType> jacobian(const InputType &in) const {
    const ValueType out = this->operator()(in);
    const double dx = 2 * (in.x() - 1) + 400 * in.x() * (in.x() * in.x() - in.y());
    const double dy = -200 * (in.x() * in.x() - in.y());
    const JacobianType jac = JacobianType(dx, dy);
    return std::make_pair(jac, out);
  }

  HessianType hessian(const InputType &in) const {
    HessianType H;
    H(0, 0) = 1200 * in.x() * in.x() - 400 * in.y() + 2;
    H(0, 1) = -400 * in.x();
    H(1, 0) = -400 * in.x();
    H(1, 1) = 200;
    return H;
  }
};
} // namespace problems
} // namespace numopt