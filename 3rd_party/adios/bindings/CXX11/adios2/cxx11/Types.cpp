/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Types.cpp : implementation of Types.h
 *
 *  Created on: Feb 12, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include "Types.h"
#include "Types.tcc"

namespace adios2
{

#define declare_template_instantiation(T) template std::string GetType<T>() noexcept;

ADIOS2_FOREACH_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace adios2
