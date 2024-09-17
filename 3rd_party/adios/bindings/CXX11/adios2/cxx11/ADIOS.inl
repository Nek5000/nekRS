/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOS.inl :
 *
 *  Created on: Jun 7, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_CXX11_CXX11_ADIOS_INL_
#define ADIOS2_BINDINGS_CXX11_CXX11_ADIOS_INL_
#ifndef ADIOS2_BINDINGS_CXX11_CXX11_ADIOS_H_
#error "Inline file should only be included from it's header, never on it's own"
#endif

#include <stdexcept>

namespace adios2
{

template <class R, class... Args>
Operator ADIOS::DefineOperator(const std::string name,
                               const std::function<R(Args...)> &function,
                               const Params &parameters)
{
    CheckPointer("for operator name " + name +
                 ", in call to ADIOS::DefineOperator");
    return Operator(DefineCallBack(name, function, parameters));
}

} // end namespace adios2

#endif /* ADIOS2_BINDINGS_CXX11_CXX11_ADIOS_INL_ */
