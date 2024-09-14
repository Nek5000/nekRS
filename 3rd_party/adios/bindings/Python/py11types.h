/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11types.h
 *
 *  Created on: Nov 7, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_PYTHON_PY11TYPES_H_
#define ADIOS2_BINDINGS_PYTHON_PY11TYPES_H_

#include <adios2.h>

namespace adios2
{
namespace py11
{

#if ADIOS2_USE_MPI

/**
 * MPI4PY_Comm provides automatic conversion of Python mpi4py communicators to
 * the C++ MPI4PY_Comm type, which in itself implicitly converts to a MPI_Comm.
 *
 * The actual work is done by the caster in py11glue.cpp
 */
struct MPI4PY_Comm
{
    MPI_Comm comm;

    // allow implicit conversion to MPI_Comm
    operator MPI_Comm() { return comm; }
};

#endif

#define ADIOS2_FOREACH_PYTHON_TYPE_1ARG(MACRO) ADIOS2_FOREACH_STDTYPE_1ARG(MACRO)

#define ADIOS2_FOREACH_NUMPY_TYPE_1ARG(MACRO) ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(MACRO)

#define ADIOS2_FOREACH_NUMPY_ATTRIBUTE_TYPE_1ARG(MACRO) ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(MACRO)

} // end namespace py11
} // end namespace adios2

#endif
