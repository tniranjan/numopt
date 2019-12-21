#include <gtest/gtest.h>
#include <base/constants.h>
#include <problems/problem.h>
#include <problems/paraboloidFunctor.h>
#include <problems/paraboloidProblem.h>
#include <problems/rosenbrock2dProblem.h>
#include <problems/rosenbrockFunctor.h>
#include <solvers/gradientDescent.h>

TEST(gradientdescent_rosenrock_sameAD_test_case, gradientdescent_test) {
  const Index nin = 2;
  const Index nout = 1;
  typedef Eigen::Matrix<double, nin, 1> IType;

  IType x;
  x.setRandom();

  numopt::functors::RosenbrockFunctor<nin> rosenFunctor;
  const numopt::problems::Rosenbrock2dProblem<
      numopt::functors::RosenbrockFunctor<nin>>
      rosenbrock = numopt::problems::Rosenbrock2dProblem<
          numopt::functors::RosenbrockFunctor<nin>>();
  
  numopt::solver::SolverSettings settings;

  numopt::GradientDescentSolver<numopt::problems::Rosenbrock2dProblem<
      numopt::functors::RosenbrockFunctor<nin>>>
      gdSolver(rosenbrock, settings);
  const auto data = gdSolver.minimize(x);

  numopt::Problem<numopt::functors::RosenbrockFunctor<nin>> rosenbrockAD(
      rosenFunctor);
  numopt::GradientDescentSolver<
      numopt::Problem<numopt::functors::RosenbrockFunctor<nin>>>
      gdSolverAD(rosenbrockAD, settings);
  const auto dataAD = gdSolverAD.minimize(x);

  EXPECT_TRUE(data.argmin.isApprox(dataAD.argmin, numopt::constants::s_eps10));
  EXPECT_TRUE(std::abs(data.min-dataAD.min) < numopt::constants::s_eps10);
}

TEST(gradientdescent_parab_sameAD_test_case, gradientdescent_test) {
  const Index nin = 32;
  const Index nout = 1;
  typedef Eigen::Matrix<double, nin, 1> IType;

  IType x;
  x.setRandom();

  numopt::functors::ParaboloidFunctor<nin> parabFunctor;
  const numopt::problems::ParaboloidProblem<
      numopt::functors::ParaboloidFunctor<nin>>
      paraboloid = numopt::problems::ParaboloidProblem<
          numopt::functors::ParaboloidFunctor<nin>>();
  
  numopt::solver::SolverSettings settings;

  numopt::GradientDescentSolver<numopt::problems::ParaboloidProblem<
      numopt::functors::ParaboloidFunctor<nin>>>
      gdSolver(paraboloid, settings);
  const auto data = gdSolver.minimize(x);

  numopt::Problem<numopt::functors::ParaboloidFunctor<nin>> paraboloidAD(
      parabFunctor);
  numopt::GradientDescentSolver<
      numopt::Problem<numopt::functors::ParaboloidFunctor<nin>>>
      gdSolverAD(paraboloidAD, settings);
  const auto dataAD = gdSolverAD.minimize(x);

  EXPECT_TRUE(data.argmin.isApprox(dataAD.argmin, numopt::constants::s_eps10));
  EXPECT_TRUE(std::abs(data.min-dataAD.min) < numopt::constants::s_eps10);
  }

  TEST(gradientdescent_parab_test_case, gradientdescent_test) {
  const Index nin = 32;
  const Index nout = 1;
  typedef Eigen::Matrix<double, nin, 1> IType;

  IType x;
  x.setRandom();

  const numopt::problems::ParaboloidProblem<
      numopt::functors::ParaboloidFunctor<nin>>
      paraboloid = numopt::problems::ParaboloidProblem<
          numopt::functors::ParaboloidFunctor<nin>>();
  
  numopt::solver::SolverSettings settings;
  settings.functionTolerance = numopt::constants::s_eps10;
  settings.parameterTolerance = numopt::constants::s_eps10;

  numopt::GradientDescentSolver<numopt::problems::ParaboloidProblem<
      numopt::functors::ParaboloidFunctor<nin>>>
      gdSolver(paraboloid, settings);
  const auto data = gdSolver.minimize(x);

  EXPECT_TRUE(data.argmin.squaredNorm() < numopt::constants::s_eps6); //can be better if properly scaled!
  EXPECT_TRUE(std::abs(data.min-5) < numopt::constants::s_eps6);
  }