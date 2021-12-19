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

.. doxygenfunction:: numopt::functors::paraboloid
.. doxygenfunction:: numopt::functors::rosenbrock
.. doxygenclass:: numopt::ProblemBase
   :members:
.. doxygenclass:: numopt::ProblemAD
   :members:

.. doxygenclass:: numopt::problems::Rosenbrock2dProblem
   :members:
.. doxygenclass:: numopt::problems::ParaboloidProblem
   :members:

.. doxygenclass:: numopt::Solver
   :members:
.. doxygenstruct:: numopt::solver::SolverData
   :members:
.. doxygenstruct:: numopt::solver::SolverSettings
   :members:
.. doxygenfunction:: numopt::solver::BackTrackingLineSearch
.. doxygenclass:: numopt::GradientDescentSolver
   :members:

