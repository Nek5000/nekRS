/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_f2c_adios_mpi.cpp
 */

#include "adios2_f2c_common.h"

#include <mpi.h>

extern "C" {

// this function is not exposed in the public APIs
extern adios2_adios *adios2_init_config_glue_mpi(const char *config_file, MPI_Comm comm,
                                                 const char *host_language);

void FC_GLOBAL(adios2_init_config_mpi_f2c, ADIOS2_INIT_CONFIG_MPI_F2C)(adios2_adios **adios,
                                                                       const char *config_file,
                                                                       MPI_Fint *comm, int *ierr)
{
    *adios = adios2_init_config_glue_mpi(config_file, MPI_Comm_f2c(*comm), "Fortran");
    *ierr = (*adios == NULL) ? static_cast<int>(adios2_error_exception)
                             : static_cast<int>(adios2_error_none);
}

void FC_GLOBAL(adios2_init_mpi_f2c, ADIOS2_INIT_MPI_F2C)(adios2_adios **adios, MPI_Fint *comm,
                                                         int *ierr)
{
    FC_GLOBAL(adios2_init_config_mpi_f2c, ADIOS2_INIT_CONFIG_MPI_F2C)
    (adios, "", comm, ierr);
}
}
