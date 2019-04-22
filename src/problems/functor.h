#pragma once

#include <base/types.h>

namespace numopt {
// A Generic Functor type, that can be used in combination with Eigen AD.
// From GGael's example https://forum.kde.org/viewtopic.php?f=74&t=111409
template <Index nIN, Index nOUT> 
struct Functor {
  typedef double Scalar;
  enum {
    InputsAtCompileTime = nIN,
    ValuesAtCompileTime = nOUT
  };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  const int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  // you should define that in the subclass :

  //   virtual void operator() (const InputType& x, ValueType* v) const = 0;
};
} // namespace NumOpt
