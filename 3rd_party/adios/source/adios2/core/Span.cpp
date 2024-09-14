/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Span.cpp
 *
 *  Created on: Apr 17, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "Span.tcc"

namespace adios2
{
namespace core
{

#define declare_template_instantiation(T) template class Span<T>;
ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_type

} // end namespace core
} // end namespace adios2
