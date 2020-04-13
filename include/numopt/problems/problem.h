#pragma once

#include <autodiff/reverse.hpp>
#include <autodiff/reverse/eigen.hpp>
#include <base/types.h>
#include <functional>

namespace numopt {
typedef autodiff::var ScalarAD;
typedef VectorS<ScalarAD> VectorAD;

class ProblemBase {
public:
  virtual double operator()(const VectorX &in) const = 0;
  virtual std::pair<VectorX, double> gradient(const VectorX &in) const = 0;
};

class ProblemAD : public ProblemBase {
public:
  ProblemAD(const FunctorAD<double> functorX,
            const FunctorAD<ScalarAD> functorAD)
      : functorX_(functorX), functorAD_(functorAD) {}

  double operator()(const VectorX &in) const override { return functorX_(in); }

  std::pair<VectorX, double> gradient(const VectorX &in) const override {
    VectorAD inAD = in.cast<ScalarAD>();
    ScalarAD fValue = functorAD_(inAD);
    const VectorX grad = autodiff::gradient(fValue, inAD);
    return std::make_pair(grad, double(fValue));
  }

  std::pair<SparseMatrixX, VectorX> hessian(const VectorX &in) const {
    VectorAD inAD = in.cast<ScalarAD>();
    ScalarAD fValue = functorAD_(inAD);
    VectorX gradient;
    const SparseMatrixX H =
        autodiff::hessian(fValue, inAD, gradient).sparseView();
    return std::make_pair(H, gradient);
  }

private:
  FunctorAD<double> functorX_;
  FunctorAD<ScalarAD> functorAD_;
};

} // namespace numopt