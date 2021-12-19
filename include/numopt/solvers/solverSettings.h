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

constexpr static double LS_InitalAlpha = 0.8;
constexpr static double LS_Rho = 0.5;
constexpr static double LS_C = 0.0001;
} // namespace solver
} // namespace numopt
