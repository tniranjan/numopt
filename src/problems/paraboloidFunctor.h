#pragma once

#include "functor.h"

namespace numopt
{
namespace functors
{
    template <Index nIN, Index nOUT> 
    struct ParaboloidFunctor : Functor<nIN, nOUT>
    {
        typedef Functor<nIN, nOUT> Base;
        using typename Base::InputType;
        using typename Base::ValueType;
        using typename Base::JacobianType;
        
        template <typename T>
        void operator()(const Eigen::Matrix<T, nIN, 1>& in, Eigen::Matrix<T, nOUT, 1>* out) const
        {
              *out<<(in.cwiseProduct(in).sum() + 5);
        }
    };
}
}
