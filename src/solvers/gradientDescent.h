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
  GradientDescentSolver(TProblem &problem, solver::SolverSettings &settings)
      : Solver<TProblem>::Solver(problem, settings) {}
  using Solver<_Problem>::problem_;
  using Solver<_Problem>::settings_;
  using Solver<_Problem>::hasConverged;
  TSolverData minimize(const InputType &initValue) {
    TSolverData solverData;
    solverData.argmin = initValue;
    const double prevNorm = problem_(initValue).x();
    double alpha(1);
    for (unsigned iter = 0; iter < settings_.maxSolverIterations; iter++) {
      const InputType dir = - direction(solverData.argmin);
      const auto
          curNorm = //(settings_.linesearchtype == LineSearchType::Armijo) ?
          solver::ArmijoLineSearch(
              solverData.argmin, dir, problem_, solverData.argmin, &alpha,
              settings_.verbosity, settings_.functionTolerance,
              settings_.parameterTolerance, settings_.maxLineSearchIterations);
      if (hasConverged(prevNorm - curNorm, alpha * dir.norm(), iter)) {
        solverData.nIter = iter;
        solverData.paramNorm = alpha * dir.norm();
        solverData.min = curNorm;
        /*if (settings_.verbosity) */{
          std::cout << "Gradient Descent minimum of : " << solverData.min
                    << " after : " << solverData.nIter
                    << " . Final parameter change : " << solverData.paramNorm;
        }
        break;
      }
    }
    return solverData;
  }

protected:
  InputType direction(const InputType &in) const override {
    // gradient = jacobianT here
    ValueType out;
    const JacobianType grad = problem_.jacobian(in, out);
    return grad.transpose();
  }

private:
};
} // namespace numopt