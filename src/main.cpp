#include <iostream>
#include <problems/paraboloidFunctor.h>
#include <problems/paraboloidProblem.h>
#include <problems/problem.h>
#include <problems/rosenbrock2dProblem.h>
#include <problems/rosenbrockFunctor.h>
#include <solvers/gradientDescent.h>

int main() {
  const Index nin = 2;
  const Index nout = 1;
  typedef Eigen::Matrix<double, nin, 1> IType;
  typedef Eigen::Matrix<double, nout, 1> OType;
  typedef Eigen::Matrix<double, nout, nin> JType;

  numopt::functors::ParaboloidFunctor<nin, nout> parabFunctor;
  numopt::Problem<numopt::functors::ParaboloidFunctor<nin, nout>> paraboloidAD(
      parabFunctor);
  const auto paraboloid = numopt::problems::ParaboloidProblem<
      numopt::functors::ParaboloidFunctor<nin, nout>>();
  IType x;
  x.setRandom();
  OType out;
  std::cout << "\n Paraboloid: \n" << std::endl;
  std::cout << "AD: \n"
            << "In : " << x.transpose() << " Out : " << paraboloidAD(x)
            << " J : " << paraboloidAD.jacobian(x, out) << std::endl;
  std::cout << "AnD: \n"
            << "In : " << x.transpose() << " Out : " << paraboloid(x)
            << " J : " << paraboloid.jacobian(x, out) << std::endl;
  numopt::functors::RosenbrockFunctor<nin, nout> rosenFunctor;
  numopt::Problem<numopt::functors::RosenbrockFunctor<nin, nout>> rosenbrockAD(
      rosenFunctor);
  const numopt::problems::Rosenbrock2dProblem<
      numopt::functors::RosenbrockFunctor<nin, nout>>
      rosenbrock = numopt::problems::Rosenbrock2dProblem<
          numopt::functors::RosenbrockFunctor<nin, nout>>();
  std::cout << "\n Rosenbrock: \n" << std::endl;
  std::cout << "AD: \n"
            << "In : " << x.transpose() << " Out : " << rosenbrockAD(x)
            << " J : " << rosenbrockAD.jacobian(x, out) << std::endl;
  std::cout << "AnD: \n"
            << "In : " << x.transpose() << " Out : " << rosenbrock(x)
            << " J : " << rosenbrock.jacobian(x, out) << std::endl;

  numopt::solver::SolverSettings settings;
  numopt::GradientDescentSolver<numopt::problems::Rosenbrock2dProblem<
      numopt::functors::RosenbrockFunctor<nin, nout>>>
      gdSolver(rosenbrock, settings);
  const auto data = gdSolver.minimize(x);
  gdSolver.printSummary(data);

  return 0;
}
