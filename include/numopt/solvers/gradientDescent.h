#pragma once
#include <solvers/linesearch/backtrackingLineSearch.h>
#include <solvers/linesearch/strongWolfeLineSearch.h>
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
    double alpha(1), dphi(0), dphi_im1(0);
    for (unsigned iter = 0; iter < settings().maxSolverIterations; iter++) {
      prevNorm = solverData.min;
      const VectorX dir = direction(solverData.argmin);
      const VectorX xCur = solverData.argmin;
      dphi = dir.dot(dir);
      const double alphaCur = iter ? alpha * dphi_im1 / dphi : 1e-5;
      dphi_im1 = dphi;
      solverData
          .min = (settings().linesearchtype ==
                  solver::SolverSettings::LineSearchType::BackTracking) ?
          solver::linesearch::BackTrackingLinesearch<_Problem>().run(
              xCur, dir, problem(), solverData.argmin, &alpha, settings(),
              settings().verbosity):
          solver::linesearch::StrongWolfeLinesearch<_Problem>(problem()).run(
              xCur, dir, solverData.argmin, alphaCur, &alpha, settings(),
              settings().verbosity);
      if (settings().verbosity) {
        std::cout << "Gradient Descent iteration : " << iter
                  << " f(x):" << solverData.min
                  << " at x_o:" << solverData.argmin << std::endl;
      }
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