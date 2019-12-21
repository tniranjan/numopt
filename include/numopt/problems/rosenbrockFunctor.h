#pragma once

#include "functor.h"

namespace numopt {
namespace functors {
template <Index nIN> struct RosenbrockFunctor : Functor<nIN> {
  typedef Functor<nIN> Base;
  using typename Base::InputType;
  using typename Base::JacobianType;
  using typename Base::ValueType;

  template <typename T>
  void operator()(const Eigen::Matrix<T, nIN, 1> &in,
                  Eigen::Matrix<T, 1, 1> *out) const {
    auto inEven = in(Eigen::seq(0, Eigen::last, 2), Eigen::all);
    auto inOdd = in(Eigen::seq(1, Eigen::last, 2), Eigen::all);
    auto onesT = inEven.eval();
    onesT.setConstant(T(1.0));
    auto p1 = (inEven.cwiseProduct(inEven) - inOdd).eval();
    auto p2 = (inEven - onesT).eval();
    T val = (T(100) * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
    *out << val;
  }
};
} // namespace functors
} // namespace numopt