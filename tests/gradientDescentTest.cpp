#include <base/constants.h>
#include <gtest/gtest.h>
#include <problems/functors.h>
#include <problems/paraboloidProblem.h>
#include <problems/problem.h>
#include <problems/rosenbrock2dProblem.h>
#include <solvers/gradientDescent.h>

TEST(gradientdescent_rosenrock_sameAD_test_case, gradientdescent_test) {
  const VectorX x = VectorX::Random(2, 1);
  numopt::ProblemAD rosenbrockProbAD(
      numopt::functors::rosenbrock<double>,
      numopt::functors::rosenbrock<numopt::ScalarAD>);
  numopt::problems::Rosenbrock2dProblem rosenbrockProb;

  numopt::solver::SolverSettings settings;

  numopt::GradientDescentSolver<numopt::problems::Rosenbrock2dProblem> gdSolver(
      rosenbrockProb, settings);
  const auto data = gdSolver.minimize(x);

  numopt::GradientDescentSolver<numopt::ProblemAD> gdSolverAD(rosenbrockProbAD,
                                                              settings);
  const auto dataAD = gdSolverAD.minimize(x);

  EXPECT_TRUE(data.argmin.isApprox(dataAD.argmin, numopt::constants::s_eps10));
  EXPECT_TRUE(std::abs(data.min - dataAD.min) < numopt::constants::s_eps10);
}

TEST(gradientdescent_parab_sameAD_test_case, gradientdescent_test) {
  const VectorX x = VectorX::Random(2048, 1);
  numopt::ProblemAD parabProbAD(numopt::functors::paraboloid<double>,
                                numopt::functors::paraboloid<numopt::ScalarAD>);
  numopt::problems::ParaboloidProblem parabProb;

  numopt::solver::SolverSettings settings;

  numopt::GradientDescentSolver<numopt::problems::ParaboloidProblem> gdSolver(
      parabProb, settings);
  const auto data = gdSolver.minimize(x);

  numopt::GradientDescentSolver<numopt::ProblemAD> gdSolverAD(parabProbAD,
                                                              settings);
  const auto dataAD = gdSolverAD.minimize(x);

  EXPECT_TRUE(data.argmin.isApprox(dataAD.argmin, numopt::constants::s_eps10));
  EXPECT_TRUE(std::abs(data.min - dataAD.min) < numopt::constants::s_eps10);
}

TEST(gradientdescent_parab_test_case, gradientdescent_test) {
  const VectorX x = VectorX::Random(2048, 1);
  numopt::problems::ParaboloidProblem parabProb;
  numopt::solver::SolverSettings settings;
  settings.functionTolerance = numopt::constants::s_eps10;
  settings.parameterTolerance = numopt::constants::s_eps10;

  numopt::GradientDescentSolver<numopt::problems::ParaboloidProblem> gdSolver(
      parabProb, settings);
  const auto data = gdSolver.minimize(x);

  EXPECT_TRUE(data.argmin.squaredNorm() <
              numopt::constants::s_eps6); // can be better if properly scaled!
  EXPECT_TRUE(std::abs(data.min - 5) < numopt::constants::s_eps6);
}