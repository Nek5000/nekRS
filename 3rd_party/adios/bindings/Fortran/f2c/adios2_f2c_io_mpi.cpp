/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_f2c_io_mpi.cpp
 */

#include "adios2_f2c_common.h"

#include "adios2/common/ADIOSTypes.h"
#include "adios2/helper/adiosFunctions.h"

extern "C" {

extern adios2_engine *adios2_open_new_comm(adios2_io *io, const char *name, const adios2_mode mode,
                                           MPI_Comm comm);

void FC_GLOBAL(adios2_open_new_comm_f2c,
               ADIOS2_OPEN_NEW_COMM_F2C)(adios2_engine **engine, adios2_io **io, const char *name,
                                         const int *open_mode, MPI_Fint *comm, int *ierr)
{
    *engine =
        adios2_open_new_comm(*io, name, static_cast<adios2_mode>(*open_mode), MPI_Comm_f2c(*comm));
    *ierr = (*engine == NULL) ? static_cast<int>(adios2_error_exception)
                              : static_cast<int>(adios2_error_none);
}
}
