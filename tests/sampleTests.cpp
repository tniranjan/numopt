#include <iostream>

#include "gtest/gtest.h"
#include <problems/paraboloidFunctor.h>
#include <problems/paraboloidProblem.h>
#include <problems/problem.h>
#include <problems/rosenbrock2dProblem.h>
#include <problems/rosenbrockFunctor.h>
#include <solvers/gradientDescent.h>

TEST(sample_test_case, sample_test)
{
    EXPECT_EQ(1, 1);
}
TEST(sample_numopt_test_case, sample_numopt_test)
{
const Index nin = 2;
  const Index nout = 1;
  typedef Eigen::Matrix<double, nin, 1> IType;
  typedef Eigen::Matrix<double, nout, 1> OType;
  typedef Eigen::Matrix<double, nout, nin> JType;

  numopt::functors::ParaboloidFunctor<nin> parabFunctor;
  numopt::Problem<numopt::functors::ParaboloidFunctor<nin>> paraboloidAD(
      parabFunctor);
  const auto paraboloid = numopt::problems::ParaboloidProblem<
      numopt::functors::ParaboloidFunctor<nin>>();
  IType x;
  x.setRandom();
  OType out;
  
    EXPECT_EQ(paraboloidAD.jacobian(x).first, paraboloid.jacobian(x).first);
}

// #include <iostream>
// #include <problems/paraboloidFunctor.h>
// #include <problems/paraboloidProblem.h>
// #include <problems/problem.h>
// #include <problems/rosenbrock2dProblem.h>
// #include <problems/rosenbrockFunctor.h>
// #include <solvers/gradientDescent.h>

// int main() {
//   const Index nin = 2;
//   const Index nout = 1;
//   typedef Eigen::Matrix<double, nin, 1> IType;
//   typedef Eigen::Matrix<double, nout, 1> OType;
//   typedef Eigen::Matrix<double, nout, nin> JType;

//   numopt::functors::ParaboloidFunctor<nin> parabFunctor;
//   numopt::Problem<numopt::functors::ParaboloidFunctor<nin>> paraboloidAD(
//       parabFunctor);
//   const auto paraboloid = numopt::problems::ParaboloidProblem<
//       numopt::functors::ParaboloidFunctor<nin>>();
//   IType x;
//   x.setRandom();
//   OType out;
//   std::cout << "\n Paraboloid: \n" << std::endl;
//   std::cout << "AD: \n"
//             << "In : " << x.transpose() << " Out : " << paraboloidAD(x)
//             << " J : " << paraboloidAD.jacobian(x).first << "\n H: \n "
//             << paraboloidAD.hessian(x) << std::endl;
//   std::cout << "AnD: \n"
//             << "In : " << x.transpose() << " Out : " << paraboloid(x)
//             << " J : " << paraboloid.jacobian(x).first << "\n H: \n "
//             << paraboloid.hessian(x) << std::endl;
//   numopt::functors::RosenbrockFunctor<nin> rosenFunctor;
//   numopt::Problem<numopt::functors::RosenbrockFunctor<nin>> rosenbrockAD(
//       rosenFunctor);
//   const numopt::problems::Rosenbrock2dProblem<
//       numopt::functors::RosenbrockFunctor<nin>>
//       rosenbrock = numopt::problems::Rosenbrock2dProblem<
//           numopt::functors::RosenbrockFunctor<nin>>();
//   std::cout << "\n Rosenbrock: \n" << std::endl;
//   std::cout << "AD: \n"
//             << "In : " << x.transpose() << " Out : " << rosenbrockAD(x)
//             << " J : " << rosenbrockAD.jacobian(x).first << "\n H: \n "
//             << rosenbrockAD.hessian(x) << std::endl;
//   std::cout << "AnD: \n"
//             << "In : " << x.transpose() << " Out : " << rosenbrock(x)
//             << " J : " << rosenbrock.jacobian(x).first<< "\n H: \n "
//             << rosenbrock.hessian(x) << std::endl;

//   numopt::solver::SolverSettings settings;
//   numopt::GradientDescentSolver<numopt::problems::Rosenbrock2dProblem<
//       numopt::functors::RosenbrockFunctor<nin>>>
//       gdSolver(rosenbrock, settings);
//   const auto data = gdSolver.minimize(x);
//   gdSolver.printSummary(data);

//   return 0;
// }
