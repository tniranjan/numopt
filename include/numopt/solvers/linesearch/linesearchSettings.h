#pragma once

namespace numopt::solver::linesearch {
constexpr static double LS_MaxAlpha = 1e5;
constexpr static double LS_InitAlpha = 0.8;
constexpr static double LS_Rho = 0.5;
constexpr static double LS_C1 = 1e-4;
constexpr static double LS_C2 = 0.9;
/**
 * @brief  Structure that holds tolerances & max # iterations.
 */
struct LinesearchSettings {
  /** Tolerance for change in f(x) */
  double functionTolerance;
  /** Tolerance for change in x*/
  double parameterTolerance;
  /**  max no. of linesearch iterations. */
  unsigned maxIterations;
  LinesearchSettings(double funcTol, double paramTol, unsigned maxIter)
      : functionTolerance(funcTol), parameterTolerance(paramTol),
        maxIterations(maxIter) {}
};
} // namespace numopt::solver::linesearch