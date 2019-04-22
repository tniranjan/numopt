#pragma once

#include "problem.h"

// Rosenbrock problem of 2 vars with analytic df

namespace numopt
{
    namespace problems
    {
        template <typename _Functor>
        class Rosenbrock2dProblem : public Problem<_Functor>
        {
            typedef _Functor TFunctor;
            typedef typename TFunctor::InputType InputType;
            typedef typename TFunctor::ValueType ValueType;
            typedef typename TFunctor::JacobianType JacobianType;

        public:
            Rosenbrock2dProblem(TFunctor &functor) : Problem<TFunctor>::Problem(functor) {}
    
            JacobianType jacobian(const InputType &in, ValueType &out) const
            {
                out = this->operator()(in);
                double dx = 2 * (in.x() - 1) + 400 * in.x() * (in.x() * in.x() -in.y());
                double dy = - 200 * (in.x() * in.x() - in.y());
                JacobianType jac = JacobianType(dx, dy);
                return jac;
            }
        };
    }
}