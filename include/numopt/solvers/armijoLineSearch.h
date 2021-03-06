#pragma once

#include <base/types.h>
#include <solvers/solverSettings.h>

namespace numopt {
namespace solver {
template <typename DerivedX, typename Evaluator>
double ArmijoLineSearch(const Eigen::MatrixBase<DerivedX> &xcur,
                        const Eigen::MatrixBase<DerivedX> &dir,
                        const Evaluator func,
                        Eigen::MatrixBase<DerivedX> &xnext, double *palpha,
                        const unsigned verbose, const double funcTol,
                        const double paramTol, const unsigned maxIterations) {
  double alphak = LS_InitalAlpha;
  const double initNorm = func(xcur);
  double prevNorm = initNorm;

  if (verbose > 1)
    std::cout << "Armijo - Init Norm : " << initNorm << std::endl;

  for (unsigned iter = 0; iter < maxIterations; iter++) {
    xnext = (xcur + alphak * dir).eval();
    const double curNorm = func(xnext);

    if (verbose > 1)
      std::cout << "Armijo - Cur Norm : " << curNorm
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