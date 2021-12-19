#pragma once
#include <solvers/linesearch/backtrackingLineSearch.h>
#include <solvers/solver.h>
#include <solvers/solverSettings.h>

namespace numopt {
template <typename _Problem>
class GradientDescentSolver : public Solver<_Problem> {
  typedef _Problem TProblem;

public:
  GradientDescentSolver(const TProblem &problem,
                        const solver::SolverSettings &settings)
      : Solver<TProblem>::Solver(problem, settings) {}
  using Solver<_Problem>::problem;
  using Solver<_Problem>::settings;
  using Solver<_Problem>::hasConverged;
  solver::SolverData minimize(const VectorX &initValue) {
    solver::SolverData solverData;
    solverData.argmin = initValue;
    double prevNorm;
    solverData.min = problem()(initValue);
    double alpha(1);
    for (unsigned iter = 0; iter < settings().maxSolverIterations; iter++) {
      prevNorm = solverData.min;
      const VectorX dir = direction(solverData.argmin);
      const VectorX xCur = solverData.argmin;
      solverData.min = //(settings().linesearchtype == solver::SolverSettings::LineSearchType::BackTracking) ?
          solver::BackTrackingLineSearch(xCur, dir, problem(), solverData.argmin,
                                   &alpha, settings().verbosity,
                                   settings().functionTolerance,
                                   settings().parameterTolerance,
                                   settings().maxLineSearchIterations);
      solverData.nIter = iter + 1;
      solverData.paramNorm = alpha * dir.norm();
      if (hasConverged(std::abs(prevNorm - solverData.min), alpha * dir.norm(),
                       iter)) {
        break;
      }
    }
    return solverData;
  }

protected:
  VectorX direction(const VectorX &in) const override {
    const VectorX grad = -problem().gradient(in).first;
    return grad;
  }

private:
};
} // namespace numopt