/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Attribute.cpp : needed for template separation using Attribute.tcc
 *
 *  Created on: Aug 3, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "Attribute.h"
#include "Attribute.tcc"

#include "adios2/common/ADIOSMacros.h"

namespace adios2
{
namespace core
{

#define declare_template_instantiation(T) template class Attribute<T>;
ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace core
} // end namespace adios2
