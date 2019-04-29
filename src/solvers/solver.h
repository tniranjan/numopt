#pragma once

#include "solverData.h"
#include "solverSettings.h"
#include <base/types.h>
#include <problems/problem.h>

#include <string>

namespace numopt {
template <typename _Problem> class Solver {
  typedef _Problem TProblem;
  typedef typename TProblem::InputType InputType;
  typedef typename TProblem::ValueType ValueType;
  typedef typename TProblem::JacobianType JacobianType;
  typedef solver::SolverData<InputType, ValueType> TSolverData;

public:
  Solver(TProblem &problem, solver::SolverSettings &settings)
      : problem_(problem), settings_(settings){};
  virtual TSolverData minimize(const InputType &initialValue) = 0;
  TProblem problem_;
  solver::SolverSettings settings_;

protected:
  virtual InputType direction(const InputType &in) const = 0;
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

private:
  bool checkCondition(const bool condition, const std::string message) const {
    if (settings_.verbosity)
      std::cout << message << std::endl;

    return condition;
  }

};
} // namespace numopt