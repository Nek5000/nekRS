/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11IOMPI.cpp
 */

#include "py11IO.h"

#include "adios2/helper/adiosCommMPI.h"
#include <mpi4py/mpi4py.h>

#include "py11types.h"

namespace adios2
{
namespace py11
{

Engine IO::Open(const std::string &name, const int mode, MPI4PY_Comm comm)
{
    helper::CheckForNullptr(m_IO, "for engine " + name + ", in call to IO::Open");

    return Engine(&m_IO->Open(name, static_cast<adios2::Mode>(mode), helper::CommDupMPI(comm)));
}

} // end namespace py11
} // end namespace adios2
