#pragma once

#include "linesearchSettings.h"
#include <base/types.h>

namespace numopt::solver::linesearch {
/**
 * Implements interpolation based linesearch algorithm (Algorithm 3.5 & 3.6)
 * Ensures Strong Wolfe conditions are satisfied
 */
template <typename DerivedX, typename Evaluator>
double StrongWolfeLineSearch(const Eigen::MatrixBase<DerivedX> &xcur,
                             const Eigen::MatrixBase<DerivedX> &dir,
                             const Evaluator func,
                             Eigen::MatrixBase<DerivedX> &xnext,
                             const double alpha_0, double *palpha,
                             const LinesearchSettings &settings,
                             const unsigned verbose) {
  double alpha_im1 = min(1.0, alpha_0 * 1.01);
  double alpha_i = std::min(alpha_im1 / LS_Rho, LS_MaxAlpha);
  const auto [grad_0, phi_0] = func.gradient(xcur);
  const auto dphi_0 = grad_0.dot(dir);
  double phi_im1 = phi_0;
  for (unsigned iter = 0; iter < settings.maxIterations; iter++) {
    const auto [grad_i, phi_i] = func.gradient(xcur + alpha_i * dir);
    const auto dphi_i = grad_i.dot(dir);
    if ((phi_i > (phi_0 + LS_C1 * alpha_i * dphi_0)) ||
        (iter && (phi_i >= phi_im1))) {
      *palpha = zoom(alpha_im1, alpha_i, xcur, dir, func);
      xnext = xcur + *palpha + dir;
      return func(xnext);
    }
    if (abs(dphi_i) <= -LS_C2 * dphi_0) {
      *palpha = alpha_i;
      xnext = xcur + *palpha * dir;
      return phi_i;
    }
    if (dphi_i >= 0) {
      *palpha = zoom(alpha_i, alpha_im1, xcur, dir, func);
      xnext = xcur + *palpha + dir;
      return func(xnext);
    }
    alpha_i = std::min(alpha_i / LS_Rho, LS_MaxAlpha);
  }
  return phi_0;
}
double interpolate_cubic(const double alpha_im1, const double phi_im1,
                         const double dphi_im1, const double alpha_i,
                         const double phi_i, const double dphi_i) {
  const auto d1 =
      dphi_im1 + dphi_i - 3 * (phi_im1 - phi_i) / (alpha_im1 - alpha_i);
  const auto d2 = std::copysign(1.0, alpha_i - alpha_im1) *
                  std::sqrt(d1 * d1 - dphi_i * dphi_im1);
  return alpha_i - (alpha_i - alpha_im1) * (dphi_i + d2 - d1) /
                       (dphi_i - dphi_im1 + 2 * d2);
}
template <typename DerivedX, typename Evaluator>
double zoom(double alpha_l, double alpha_h,
            const Eigen::MatrixBase<DerivedX> &xcur,
            const Eigen::MatrixBase<DerivedX> &dir, Evaluator func) {
  for (int j = 0; j < 10; j++) {      
      double alpha_j = interpolate_cubic()
  
  }
}
} // namespace numopt::solver::linesearch