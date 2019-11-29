#pragma once

#include <base/types.h>
#include <unsupported/Eigen/AutoDiff>

namespace numopt {
template <typename _Functor> class Problem {
public:
  typedef _Functor TFunctor;
  typedef Eigen::AutoDiffJacobian<TFunctor> TJacobianAD;
  typedef typename TFunctor::InputType InputType;
  typedef typename TFunctor::ValueType ValueType;
  typedef typename TFunctor::JacobianType JacobianType;

  Problem() : functor_(TFunctor()) {}
  Problem(TFunctor &functor) : functor_(functor) {}

  ValueType operator()(const InputType &in) const {
    ValueType out;
    functor_(in, &out);
    return out;
  }

  JacobianType jacobian(const InputType &in, ValueType &out) const {
    JacobianType jac;
    jacFunctor_(in, &out, &jac);
    return jac;
  }

private:
  TFunctor functor_;
  TJacobianAD jacFunctor_;
};
} // namespace numopt