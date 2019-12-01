#pragma once

#include "problem.h"

// paraboloid problem with analytic df

namespace numopt {
namespace problems {
template <typename _Functor>
class ParaboloidProblem : public Problem<_Functor> {
public:
  typedef _Functor TFunctor;
  typedef typename TFunctor::InputType InputType;
  typedef typename TFunctor::ValueType ValueType;
  typedef typename TFunctor::JacobianType JacobianType;
  typedef typename TFunctor::HessianType HessianType;

  ParaboloidProblem() : Problem<TFunctor>::Problem() {}
  ParaboloidProblem(TFunctor &functor) : Problem<TFunctor>::Problem(functor) {}

  std::pair<JacobianType, ValueType>  jacobian(const InputType &in) const {
    ValueType out;
    out = this->operator()(in);
    JacobianType jac = (2 * in.transpose());
    return std::make_pair(jac,out);
  }

  HessianType hessian(const InputType &in) const{
    HessianType H;
    H.setIdentity();
    H = 2*H;
    return H;
  }
};
} // namespace problems
} // namespace numopt