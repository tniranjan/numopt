#include <gtest/gtest.h>
#include <base/constants.h>
#include <problems/paraboloidFunctor.h>
#include <problems/paraboloidProblem.h>
#include <problems/problem.h>
#include <problems/rosenbrock2dProblem.h>
#include <problems/rosenbrockFunctor.h>
#include <solvers/gradientDescent.h>

TEST(autodiff_parab_test_case, autodiff_test) {
  const Index nin = 2;
  typedef Eigen::Matrix<double, nin, 1> IType;

  numopt::functors::ParaboloidFunctor<nin> parabFunctor;
  numopt::Problem<numopt::functors::ParaboloidFunctor<nin>> paraboloidAD(
      parabFunctor);
  const auto paraboloid = numopt::problems::ParaboloidProblem<
      numopt::functors::ParaboloidFunctor<nin>>();
  IType x;
  x.setRandom();
  EXPECT_TRUE(paraboloidAD.jacobian(x).first.isApprox(paraboloid.jacobian(x).first, numopt::constants::s_eps10));
  EXPECT_TRUE(paraboloidAD.jacobian(x).second.isApprox(paraboloid.jacobian(x).second, numopt::constants::s_eps10));
  EXPECT_TRUE(paraboloidAD.hessian(x).isApprox(paraboloid.hessian(x), numopt::constants::s_eps10));
}

TEST(autodiff_rosenbrock_test_case, autodiff_test) {
  const Index nin = 2;
  typedef Eigen::Matrix<double, nin, 1> IType;

  numopt::functors::RosenbrockFunctor<nin> rosenFunctor;
  numopt::Problem<numopt::functors::RosenbrockFunctor<nin>> rosenbrockAD(
      rosenFunctor);
  const numopt::problems::Rosenbrock2dProblem<
      numopt::functors::RosenbrockFunctor<nin>>
      rosenbrock = numopt::problems::Rosenbrock2dProblem<
          numopt::functors::RosenbrockFunctor<nin>>();
  IType x;
  x.setRandom();
  EXPECT_TRUE(rosenbrockAD.jacobian(x).first.isApprox(rosenbrock.jacobian(x).first, numopt::constants::s_eps10));
  EXPECT_TRUE(rosenbrockAD.jacobian(x).second.isApprox(rosenbrock.jacobian(x).second, numopt::constants::s_eps10));
  EXPECT_TRUE(rosenbrockAD.hessian(x).isApprox(rosenbrock.hessian(x), numopt::constants::s_eps10));
}
