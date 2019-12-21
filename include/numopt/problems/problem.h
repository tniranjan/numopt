#pragma once

#include <base/types.h>
#include <unsupported/Eigen/AutoDiff>

namespace numopt {
template <typename _Functor> class Problem {
public:
  typedef _Functor TFunctor;
  typedef typename TFunctor::InputType InputType;
  typedef typename TFunctor::ValueType ValueType;
  typedef typename TFunctor::JacobianType JacobianType;
  typedef typename TFunctor::HessianType HessianType;

  // use with AutoDiffScalar
  typedef InputType DerivativeType;
  typedef Eigen::AutoDiffScalar<DerivativeType> ActiveScalar;
  typedef Eigen::Matrix<ActiveScalar, TFunctor::InputsAtCompileTime, 1>
      ActiveInput;
  typedef Eigen::Matrix<ActiveScalar, TFunctor::InputsAtCompileTime, 1>
      OuterDerivativeType;
  typedef Eigen::AutoDiffScalar<OuterDerivativeType> OuterActiveScalar;
  typedef Eigen::Matrix<OuterActiveScalar, TFunctor::InputsAtCompileTime, 1>
      OuterActiveInput;
  typedef Eigen::Matrix<ActiveScalar, 1, 1> ActiveValue;
  typedef Eigen::Matrix<OuterActiveScalar, 1, 1> OuterActiveValue;

  Problem() : functor_(TFunctor()) {}
  Problem(TFunctor &functor) : functor_(functor) {}

  ValueType operator()(const InputType &in) const {
    ValueType out;
    functor_(in, &out);
    return out;
  }

  std::pair<JacobianType, ValueType> jacobian(const InputType &in) const {
    ActiveInput ax = in.template cast<ActiveScalar>();
    ActiveValue av(1);
    ValueType out;
    av[0].derivatives().resize(in.rows());

    for (Index i = 0; i < in.rows(); i++)
      ax[i].derivatives() = DerivativeType::Unit(in.rows(), i);

    functor_(ax, &av);

    out.x() = av[0].value();
    const JacobianType jac = av[0].derivatives();
    return std::make_pair(jac, out);
  }

  HessianType hessian(const InputType &in) const {
    OuterActiveInput axx = in.template cast<OuterActiveScalar>();
    for (int i = 0; i < axx.rows(); i++) {
      // initialize derivative direction in value field of outer active variable
      axx(i).value().derivatives() =
          OuterActiveScalar::DerType::Scalar::DerType::Unit(in.rows(), i);
      // initialize derivatives direction of the variable
      axx(i).derivatives() = DerivativeType::Unit(in.rows(), i);
      // initialize Hessian matrix of variable to zero
      for (int idx = 0; idx < in.rows(); idx++) {
        axx(i).derivatives()(idx).derivatives() =
            OuterActiveScalar::DerType::Scalar::DerType::Zero(in.rows());
      }
    }
    OuterActiveValue ayy(1);
    functor_(axx, &ayy);
    HessianType H;
    for (Index i = 0; i < in.rows(); i++)
      H.middleRows(i, 1) =
          ayy.value().derivatives()(i).derivatives().transpose();
    return H;
  }

private:
  TFunctor functor_;
};
} // namespace numopt