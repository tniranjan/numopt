#pragma once

#include "problem.h"

// paraboloid problem with analytic df

namespace numopt
{
    namespace problems
    {
        template <typename _Functor>
        class ParaboloidProblem : public Problem<_Functor>
        {
            typedef _Functor TFunctor;
            typedef typename TFunctor::InputType InputType;
            typedef typename TFunctor::ValueType ValueType;
            typedef typename TFunctor::JacobianType JacobianType;

        public:
            ParaboloidProblem(TFunctor &functor) : Problem<TFunctor>::Problem(functor) {}
    
            JacobianType jacobian(const InputType &in, ValueType &out) const
            {
                out = this->operator()(in);
                JacobianType jac = (2*in.transpose());
                return jac;
            }
        };
    }
}