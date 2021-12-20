numopt::
========

The goal is to implement all algorithms from **Numerical Optimization, Jorge Nocedal & Stephen J. Wright**, and in the process better my software engineering skills.
A homepage will eventually reside here_.

.. _here: https://tniranjan.github.io/numopt/

*Here be dragons* 🐉🐉🐉
---------------------
Or a yali_. Whatever sinks your boat.

.. _yali: https://4.bp.blogspot.com/-rCFOkVcUoO0/WjqOLLDjY-I/AAAAAAAAAeM/Jjvq5OBQ37UuCddKNhAbR1loJTEYm2VegCLcBGAs/s1600/yazhi%2BWIP.JPG

Docs
====

Problems
--------
A problem defines the objective funtion that is to be minimized (unconstrained as of now). It must atleast define the evaluation function and a gradient function. 

.. doxygenclass:: numopt::ProblemBase
   :members:
.. doxygenclass:: numopt::problems::Rosenbrock2dProblem
   :members:
.. doxygenclass:: numopt::problems::ParaboloidProblem
   :members:



.. doxygenclass:: numopt::ProblemAD
   :members:

Functors
--------
It might not always be easy to derive your gradient vector, and even more hard to get your hessian matrix. In such situations, we can use autodiff to obtain the gradient and Hessian. The class :class:`numopt::ProblemAD` implements this functionality, as long as the desired function is defined as a template function.
The following functors serve as examples.

.. doxygenfunction:: numopt::functors::paraboloid
.. doxygenfunction:: numopt::functors::rosenbrock

Solvers
-------
As of now this is an unconstrained optimization solver, that takes a problem, and finds its minima.

.. doxygenstruct:: numopt::solver::SolverData
   :members:
.. doxygenstruct:: numopt::solver::SolverSettings
   :members:
.. doxygenclass:: numopt::Solver
   :members:
.. doxygenfunction:: numopt::solver::BackTrackingLineSearch
.. doxygenclass:: numopt::GradientDescentSolver
   :members:

