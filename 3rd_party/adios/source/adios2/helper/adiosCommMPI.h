/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosCommMPI.h : MPI-specific Comm APIs
 */

#ifndef ADIOS2_HELPER_ADIOSCOMMMPI_H_
#define ADIOS2_HELPER_ADIOSCOMMMPI_H_

#include "adios2/common/ADIOSConfig.h"
#if !ADIOS2_USE_MPI
#error "Do not include adiosCommMPI.h without ADIOS2_USE_MPI."
#endif

#include "adiosComm.h"

#include <mpi.h>

namespace adios2
{
namespace helper
{

/**
 * @brief Create a communicator by taking ownership of a MPI communicator.
 */
Comm CommWithMPI(MPI_Comm mpiComm);

/**
 * @brief Create a communicator by duplicating a MPI communicator.
 */
Comm CommDupMPI(MPI_Comm mpiComm);

/**
 * @brief Get the underlying raw MPI communicator.
 *
 * Returns MPI_COMM_NULL if the communicator is empty or is not
 * backed by MPI.
 */
MPI_Comm CommAsMPI(Comm const &comm);

} // end namespace helper
} // end namespace adios2

#endif // ADIOS2_HELPER_ADIOSCOMMMPI_H_
