#pragma once

namespace numopt {
namespace solver {
template <typename InputType, typename ValueType> struct SolverData {
  InputType argmin; // solution, x
  double min;    // function value, f(x)
  double paramNorm; // alpha * param.norm()
  int nIter;
};
} // namespace solver
} // namespace numopt