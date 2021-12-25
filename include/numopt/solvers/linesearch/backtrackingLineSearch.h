#pragma once

#include "linesearchSettings.h"
#include "solvers/solverSettings.h"
#include <base/types.h>

namespace numopt::solver::linesearch {
/**
 * Implements the Backtracking algorithm. (Algorithm 3.1)
 * Only ensures sufficient decrease (weak Wolfe Conditions)
 */

template <typename Evaluator> class BackTrackingLinesearch {
public:
  double run(const VectorX &xcur, const VectorX &dir, const Evaluator func,
             VectorX &xnext, double *palpha, const SolverSettings &settings,
             const unsigned verbose) {
    double alphak = LS_InitAlpha;
    const double initNorm = func(xcur);
    double prevNorm = initNorm;

    if (verbose > 1)
      std::cout << "BackTracking - Init Norm : " << initNorm << std::endl;

    for (unsigned iter = 0; iter < settings.maxLineSearchIterations; iter++) {
      xnext = (xcur + alphak * dir).eval();
      const double curNorm = func(xnext);

      if (verbose > 1)
        std::cout << "BackTracking - Cur Norm : " << curNorm
                  << " Step Size : " << alphak << std::endl;

      if (curNorm <= initNorm - LS_C1 * alphak * dir.norm()) {
        if (palpha)
          *palpha = alphak;

        return curNorm;
      }

      if ((iter == (settings.maxLineSearchIterations - 1)) ||
          (std::fabs(curNorm - prevNorm) < settings.functionTolerance) ||
          (alphak * dir.norm() < settings.parameterTolerance))
        break;

      alphak = alphak * LS_Rho;
      prevNorm = curNorm;
    }

    if (palpha)
      *palpha = 0.0;

    xnext = xcur;

    return initNorm;
  }
};

} // namespace numopt::solver::linesearch
