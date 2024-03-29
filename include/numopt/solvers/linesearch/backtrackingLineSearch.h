#pragma once

#include <base/types.h>

namespace numopt {
namespace solver {
constexpr static double LS_InitalAlpha = 0.8;
constexpr static double LS_Rho = 0.5;
constexpr static double LS_C = 0.0001;
/**
 * Implements the Backtracking algorithm. (3.1)
 * Only ensures sufficient decrease (weak Wolfe Conditions)
 */

template <typename DerivedX, typename Evaluator>
double BackTrackingLineSearch(const Eigen::MatrixBase<DerivedX> &xcur,
                              const Eigen::MatrixBase<DerivedX> &dir,
                              const Evaluator func,
                              Eigen::MatrixBase<DerivedX> &xnext,
                              double *palpha, const unsigned verbose,
                              const double funcTol, const double paramTol,
                              const unsigned maxIterations) {
  double alphak = LS_InitalAlpha;
  const double initNorm = func(xcur);
  double prevNorm = initNorm;

  if (verbose > 1)
    std::cout << "BackTracking - Init Norm : " << initNorm << std::endl;

  for (unsigned iter = 0; iter < maxIterations; iter++) {
    xnext = (xcur + alphak * dir).eval();
    const double curNorm = func(xnext);

    if (verbose > 1)
      std::cout << "BackTracking - Cur Norm : " << curNorm
                << " Step Size : " << alphak << std::endl;

    if (curNorm <= initNorm - LS_C * alphak * dir.norm()) {
      if (palpha)
        *palpha = alphak;

      return curNorm;
    }

    if ((iter == (maxIterations - 1)) ||
        (std::fabs(curNorm - prevNorm) < funcTol) ||
        (alphak * dir.norm() < paramTol))
      break;

    alphak = alphak * LS_Rho;
    prevNorm = curNorm;
  }

  if (palpha)
    *palpha = 0.0;

  xnext = xcur;

  return initNorm;
}
} // namespace solver
} // namespace numopt