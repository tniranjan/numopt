#pragma once

#include "linesearchSettings.h"
#include "solvers/solverSettings.h"
#include <base/types.h>

namespace numopt::solver::linesearch {
/**
 * Implements interpolation based linesearch algorithm (Algorithm 3.5 & 3.6)
 * Ensures Strong Wolfe conditions are satisfied
 */
template <typename Evaluator> class StrongWolfeLinesearch {
public:
  StrongWolfeLinesearch(const Evaluator &functor) : functor_(functor) {}
  double run(const VectorX &xcur, const VectorX &dir, VectorX &xnext,
             const double alpha_0, double *palpha,
             const SolverSettings &settings, const unsigned verbose) {
    xcur_ = xcur;
    dir_ = dir;
    auto phiInfo_i =getPhiInfo(std::min(1.0, alpha_0 * 1.01));
    phiInfo_0_ = getPhiInfo(0);
    PhiInfo phiInfo_im1 = phiInfo_0_;
    PhiInfo phiInfo_amax = getPhiInfo(LS_MaxAlpha);
    for (unsigned iter = 0; iter < settings.maxLineSearchIterations; iter++) {
      if (check_armijo(phiInfo_i) || (iter && (phiInfo_i.phi >= phiInfo_im1.phi))) {
        auto phiOptimal = zoom(phiInfo_im1, phiInfo_i);
        *palpha = phiOptimal.alpha;
        xnext = xcur_ + *palpha * dir_;
        return functor_(xnext);
      }
      if (abs(phiInfo_i.dphi) <= -LS_C2 * phiInfo_0_.dphi) {
        *palpha = phiInfo_i.alpha;
        xnext = xcur_ + *palpha * dir_;
        return phiInfo_i.phi;
      }
      if (phiInfo_i.dphi >= 0) {
        auto phiOptimal = zoom(phiInfo_i, phiInfo_im1);
        *palpha = phiOptimal.alpha;
        xnext = xcur_ + *palpha * dir_;
        return functor_(xnext);
      }
      phiInfo_im1 = phiInfo_i;
      phiInfo_i = //getphiinfo(std::min(phiinfo_im1.alpha / ls_rho, ls_maxalpha));
      getPhiInfo(interpolate_cubic(phiInfo_i,phiInfo_amax ));
      if (verbose > 1) {
        std::cout << "StrongWolfe iteration " << iter
                  << " f(x) : " << phiInfo_i.phi << " alpha " << phiInfo_i.alpha
                  << std::endl;
      }
    }
    return phiInfo_0_.phi;
  }

private:
  struct PhiInfo {
    double alpha;
    double phi;
    double dphi;
  };

  PhiInfo getPhiInfo(const double alpha) {
    const auto [grad, phi] = functor_.gradient(xcur_ + alpha * dir_);
    const auto dphi = grad.dot(dir_);
    return {alpha, phi, dphi};
  }

  double interpolate_cubic(const PhiInfo low, const PhiInfo high) {
    const auto d1 = low.dphi + high.dphi -
                    3 * (low.phi - high.phi) / (low.alpha - high.alpha);
    const auto d2 = std::copysign(1.0, high.alpha - low.alpha) *
                    std::sqrt(d1 * d1 - high.dphi * low.dphi);
    return high.alpha - (high.alpha - low.alpha) * (high.dphi + d2 - d1) /
                            (high.dphi - low.dphi + 2 * d2);
  }

  PhiInfo zoom(PhiInfo phiInfo_l, PhiInfo phiInfo_h) {
    for (int j = 0; j < 10; j++) {
      double alpha_j = interpolate_cubic(phiInfo_l, phiInfo_h);
      PhiInfo phiInfo_j = getPhiInfo(alpha_j);
    //   std::cout << "Z " << phiInfo_j.phi << " " << phiInfo_j.alpha << std::endl;
      if (check_armijo(phiInfo_j) || (phiInfo_j.phi >= phiInfo_l.phi)) {
        phiInfo_h = phiInfo_j;
      } else {
        if (std::abs(phiInfo_j.dphi) <= -LS_C2 * phiInfo_0_.dphi) {
          return phiInfo_j;
        }
        if (phiInfo_j.dphi * (phiInfo_h.alpha - phiInfo_l.alpha) >= 0) {
          phiInfo_h = phiInfo_l;
        }
        phiInfo_l = phiInfo_j;
      }
    }
    return phiInfo_l;
  }

  bool check_armijo(const PhiInfo phiInfo) const {
    return (phiInfo.phi >
            (phiInfo_0_.phi + LS_C1 * phiInfo.alpha * phiInfo_0_.dphi));
  }

  PhiInfo phiInfo_0_;
  const Evaluator &functor_;
  VectorX xcur_;
  VectorX dir_;
};
} // namespace numopt::solver::linesearch