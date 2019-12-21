#pragma once
#include <solvers/armijoLineSearch.h>
#include <solvers/solver.h>
#include <solvers/solverSettings.h>

namespace numopt {
template <typename _Problem>
class GradientDescentSolver : public Solver<_Problem> {
  typedef _Problem TProblem;
  typedef typename TProblem::InputType InputType;
  typedef typename TProblem::ValueType ValueType;
  typedef typename TProblem::JacobianType JacobianType;
  typedef solver::SolverData<InputType, ValueType> TSolverData;

public:
  GradientDescentSolver(const TProblem &problem,
                        const solver::SolverSettings &settings)
      : Solver<TProblem>::Solver(problem, settings) {}
  using Solver<_Problem>::problem;
  using Solver<_Problem>::settings;
  using Solver<_Problem>::hasConverged;
  TSolverData minimize(const InputType &initValue) {
    TSolverData solverData;
    solverData.argmin = initValue;
    double prevNorm;
    solverData.min = problem()(initValue).x();
    double alpha(1);
    for (unsigned iter = 0; iter < settings().maxSolverIterations; iter++) {
      prevNorm = solverData.min;
      const InputType dir = direction(solverData.argmin);
      const InputType xCur = solverData.argmin;
      solverData.min = //(settings().linesearchtype == LineSearchType::Armijo) ?
          solver::ArmijoLineSearch(xCur, dir, problem(), solverData.argmin,
                                   &alpha, settings().verbosity,
                                   settings().functionTolerance,
                                   settings().parameterTolerance,
                                   settings().maxLineSearchIterations);      
      solverData.nIter = iter + 1;
      solverData.paramNorm = alpha * dir.norm();
      if (hasConverged(std::abs(prevNorm - solverData.min), alpha * dir.norm(), iter)) {
        break;
      }
    }
    return solverData;
  }

protected:
  InputType direction(const InputType &in) const override {
    // gradient = jacobianT here
    ValueType out;
    const JacobianType grad = -problem().jacobian(in).first;
    return grad.transpose();
  }

private:
};
} // namespace numopt