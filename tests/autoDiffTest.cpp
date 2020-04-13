#include <base/constants.h>
#include <gtest/gtest.h>
#include <problems/functors.h>
#include <problems/problem.h>
#include <problems/paraboloidProblem.h>
#include <problems/rosenbrock2dProblem.h>

TEST(autodiff_parab_test_case, autodiff_test) {
  const VectorX x = VectorX::Random(256, 1);
  numopt::ProblemAD parabProbAD(numopt::functors::paraboloid<double>,
                                numopt::functors::paraboloid<numopt::ScalarAD>);
  numopt::problems::ParaboloidProblem parabProb;
  EXPECT_TRUE(parabProbAD.gradient(x).first.isApprox(
      parabProb.gradient(x).first, numopt::constants::s_eps10));
  EXPECT_TRUE(
      std::abs(parabProbAD.gradient(x).second - parabProb.gradient(x).second) <
      numopt::constants::s_eps10);
  EXPECT_TRUE(parabProbAD.hessian(x).first.isApprox(
      parabProb.hessian(x).first, numopt::constants::s_eps10));
}

TEST(autodiff_rosenbrock_test_case, autodiff_test) {
  const VectorX x = VectorX::Random(2, 1);
  numopt::ProblemAD rosenbrockProbAD(
      numopt::functors::rosenbrock<double>,
      numopt::functors::rosenbrock<numopt::ScalarAD>);
  numopt::problems::Rosenbrock2dProblem rosenbrockProb;
  EXPECT_TRUE(rosenbrockProbAD.gradient(x).first.isApprox(
      rosenbrockProb.gradient(x).first, numopt::constants::s_eps10));
  EXPECT_TRUE(std::abs(rosenbrockProbAD.gradient(x).second -
                       rosenbrockProb.gradient(x).second) <
              numopt::constants::s_eps10);
  EXPECT_TRUE(rosenbrockProbAD.hessian(x).first.isApprox(
      rosenbrockProb.hessian(x).first, numopt::constants::s_eps10));
}