#pragma once

#include "solverData.h"
#include "solverSettings.h"
#include <base/types.h>
#include <problems/problem.h>

#include <string>

namespace numopt {
template <typename _Problem> class Solver {
  typedef _Problem TProblem;

public:
  Solver(const TProblem &problem, const solver::SolverSettings &settings)
      : problem_(problem),
        settings_(settings) {};
  virtual solver::SolverData minimize(const VectorX &initialValue) = 0;
  void printSummary(const solver::SolverData &solverData) const {
    std::cout << "Achieved a function minimum of : " << solverData.min
              << " after : " << solverData.nIter
              << " . Final parameter change : " << solverData.paramNorm
              << std::endl;
  }

protected:
  virtual VectorX direction(const VectorX &in) const = 0;
  int verbosity() const { return settings_.verbosity; }
  bool hasConverged(const double delFunc, const double delParam,
                    const unsigned curIter) const {
    return (
        checkCondition(delFunc < settings_.functionTolerance,
                       "Achieved function tolerance at iteration " +
                           std::to_string(curIter)) ||
        checkCondition(delParam < settings_.parameterTolerance,
                       "Achieved parameter tolerance at iteration " +
                           std::to_string(curIter)) ||
        checkCondition(curIter == (settings_.maxSolverIterations - 1),
                       "Reached maximum no.of iterations without convergence"));
  }
  const TProblem &problem() const { return problem_; }
  const solver::SolverSettings &settings() const { return settings_; }

private:
  bool checkCondition(const bool condition, const std::string message) const {
    if (settings_.verbosity && condition)
      std::cout << message << std::endl;

    return condition;
  }

  const TProblem &problem_;
  const solver::SolverSettings &settings_;
};
} // namespace numopt