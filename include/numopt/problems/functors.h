#pragma once

#include <base/types.h>

namespace numopt {
/**
 * Define test functors here
 */

namespace functors {
/**
 * returns the paraboloid with an offset of 5.
 */
template <typename T> T paraboloid(const VectorS<T> &in) {
  return (in.squaredNorm() + T(5));
}

/**
 * returns the rosenbrock function
 */
template <typename T> T rosenbrock(const VectorS<T> &in) {
  const VectorS<T> inEven =
      in(Eigen::seq(0, Eigen::last, 2), Eigen::all).eval();
  const VectorS<T> inOdd = in(Eigen::seq(1, Eigen::last, 2), Eigen::all).eval();
  const VectorS<T> onesT = VectorS<T>::Ones(inEven.rows());
  const VectorS<T> p1 = (inEven.cwiseProduct(inEven) - inOdd).eval();
  const VectorS<T> p2 = (inEven - onesT).eval();
  return (T(100) * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}
} // namespace functors
} // namespace numopt