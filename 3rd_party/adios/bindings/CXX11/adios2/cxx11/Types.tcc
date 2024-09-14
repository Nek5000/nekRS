/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Types.tcc : implement template functions
 *
 *  Created on: Feb 11, 2019
 *      Author: William F Godoy
 */

#ifndef ADIOS2_BINDINGS_CXX11_TYPES_TCC_
#define ADIOS2_BINDINGS_CXX11_TYPES_TCC_

#include "Types.h"

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{

template <class T>
std::string GetType() noexcept
{
    return ToString(helper::GetDataType<typename TypeInfo<T>::IOType>());
}

} // end namespace adios2

#endif /* ADIOS2_BINDINGS_CXX11_TYPES_TCC_ */
