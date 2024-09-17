/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_adios_mpi.cpp : MPI-specific C bindings
 */

#include "adios2_c_adios.h"

#include "adios2/core/ADIOS.h"
#include "adios2/helper/adiosFunctions.h"

#include "adios2/helper/adiosCommMPI.h"

extern "C" {

// to be called from other languages, hidden from the public apis
adios2_adios *adios2_init_config_glue_mpi(const char *config_file, MPI_Comm comm,
                                          const char *host_language)
{
    adios2_adios *adios = nullptr;

    try
    {
        adios2::helper::CheckForNullptr(
            config_file, "for config_file, in call to adios2_init or adios2_init_config");
        adios = reinterpret_cast<adios2_adios *>(
            new adios2::core::ADIOS(config_file, adios2::helper::CommDupMPI(comm), host_language));
    }
    catch (...)
    {
        adios2::helper::ExceptionToError("adios2_init or adios2_init_config");
    }
    return adios;
}

adios2_adios *adios2_init_mpi(MPI_Comm comm) { return adios2_init_config_glue_mpi("", comm, "C"); }

adios2_adios *adios2_init_config_mpi(const char *config_file, MPI_Comm comm)
{
    return adios2_init_config_glue_mpi(config_file, comm, "C");
}

} // end extern C
