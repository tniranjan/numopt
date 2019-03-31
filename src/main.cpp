#include <iostream>
#include <solvers/solver.h>
#include <problems/problem.h>

int main()
{
    std::cout<<"Hello World"<<std::endl;
    VectorX x(10); x.setRandom();
    std::cout<<x<<std::endl;
    return 0;
}
