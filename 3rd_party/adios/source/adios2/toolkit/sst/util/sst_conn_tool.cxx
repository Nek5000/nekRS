/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include "sst_comm_fwd.h"

#if ADIOS2_USE_MPI
#include "adios2/helper/adiosCommMPI.h"
static adios2::helper::Comm CommWorld =
    adios2::helper::CommWithMPI(MPI_COMM_WORLD);
#else
#include "adios2/helper/adiosCommDummy.h"
static adios2::helper::Comm CommWorld = adios2::helper::CommDummy();
#endif

SMPI_Comm SMPI_COMM_WORLD = &CommWorld;

extern "C" {
int SMPI_Init(int *argc, char ***argv)
{
#if ADIOS2_USE_MPI
    return MPI_Init(argc, argv);
#else
    static_cast<void>(argc);
    static_cast<void>(argv);
    return 0;
#endif
}
}
