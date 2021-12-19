#pragma once

namespace numopt {
namespace solver {
struct SolverSettings {
  enum class LineSearchType { BackTracking, Interpolation };
  double functionTolerance;
  double parameterTolerance;
  unsigned maxSolverIterations;
  unsigned maxLineSearchIterations;
  int verbosity;
  LineSearchType linesearchtype;
  SolverSettings()
      : functionTolerance(1e-6), parameterTolerance(1e-6),
        maxSolverIterations(100), maxLineSearchIterations(50), verbosity(0),
        linesearchtype(LineSearchType::BackTracking) {}
};


} // namespace solver
} // namespace numopt
