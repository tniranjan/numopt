#pragma once

#include <base/types.h>

namespace numopt {
// A Generic Functor type, that can be used in combination with Eigen AD.
// From GGael's example https://forum.kde.org/viewtopic.php?f=74&t=111409
template <Index nIN> struct Functor {
  typedef double Scalar;
  enum { InputsAtCompileTime = nIN};
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, 1, 1> ValueType;
  typedef Eigen::Matrix<Scalar, 1, InputsAtCompileTime> JacobianType;
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, InputsAtCompileTime> HessianType;

  const int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(1) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  // you should define that in the subclass :

  //   virtual void operator() (const InputType& x, ValueType* v) const = 0;
};
} // namespace numopt
