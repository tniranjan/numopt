#pragma once

#include "functor.h"

namespace numopt {
namespace functors {
template <Index nIN> struct ParaboloidFunctor : Functor<nIN> {
  typedef Functor<nIN> Base;
  using typename Base::InputType;
  using typename Base::JacobianType;
  using typename Base::ValueType;

  template <typename T>
  void operator()(const Eigen::Matrix<T, nIN, 1> &in,
                  Eigen::Matrix<T, 1, 1> *out) const {
    *out << (in.squaredNorm() + T(5));
  }
};
} // namespace functors
} // namespace numopt
