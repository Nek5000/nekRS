/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <adios2.h>

#if !ADIOS2_USE_MPI
#error "ADIOS2_USE_MPI is not true for source using ADIOS2 MPI bindings"
#endif

#include <mpi.h>

int main(int argc, char **argv)
{
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    adios2::ADIOS adios(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
