#pragma once

namespace numopt {
namespace solver {
struct SolverData {
  VectorX argmin;   // solution, x
  double min;       // function value, f(x)
  double paramNorm; // alpha * param.norm()
  int nIter;        // #iterations
  void printSummary() const {
    std::cout << "Achieved a function minimum of : " << min
              << " after : " << nIter
              << " iterations. Final parameter change : " << paramNorm
              << " Argmin was : " << argmin.transpose() << std::endl;
              }
  };
} // namespace solver
} // namespace numopt