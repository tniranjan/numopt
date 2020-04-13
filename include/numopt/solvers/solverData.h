#pragma once

namespace numopt {
namespace solver {
struct SolverData {
  VectorX argmin;   // solution, x
  double min;       // function value, f(x)
  double paramNorm; // alpha * param.norm()
  int nIter;        // #iterations
};
} // namespace solver
} // namespace numopt