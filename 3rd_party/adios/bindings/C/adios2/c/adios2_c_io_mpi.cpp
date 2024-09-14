/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_io_mpi.cpp : MPI-specific C bindings for IO
 */

#include "adios2_c_io.h"
#include "adios2_c_io.tcc"

#include <vector>

#include "adios2/core/IO.h"
#include "adios2/helper/adiosFunctions.h" //GetDataType<T>
#include "adios2_c_internal.h"

#include "adios2/helper/adiosCommMPI.h"

extern "C" {

adios2_engine *adios2_open_new_comm(adios2_io *io, const char *name, const adios2_mode mode,
                                    MPI_Comm comm)
{
    adios2_engine *engine = nullptr;
    try
    {
        adios2::helper::CheckForNullptr(io, "for adios2_io, in call to adios2_open");
        engine = reinterpret_cast<adios2_engine *>(&reinterpret_cast<adios2::core::IO *>(io)->Open(
            name, adios2_ToOpenMode(mode), adios2::helper::CommDupMPI(comm)));
    }
    catch (...)
    {
        adios2::helper::ExceptionToError("adios2_open_new_comm");
    }
    return engine;
}

} // end extern C
