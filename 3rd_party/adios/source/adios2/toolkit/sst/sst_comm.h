/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#ifndef ADIOS2_TOOLKIT_SST_SST_COMM_H_
#define ADIOS2_TOOLKIT_SST_SST_COMM_H_

// SMPI, short for "SST's MPI", provides a MPI-like C interface wrapping
// around the ADIOS multi-process communcation API in adios2::helper::Comm.
// The SST implementation uses the wrapper to access the C++ Comm API
// from C.

#include "sst_comm_fwd.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum
{
    SMPI_INT,
    SMPI_LONG,
    SMPI_SIZE_T,
    SMPI_CHAR,
    SMPI_BYTE
} SMPI_Datatype;

typedef enum
{
    SMPI_MAX,
    SMPI_LAND
} SMPI_Op;

int SMPI_Comm_rank(SMPI_Comm comm, int *rank);
int SMPI_Comm_size(SMPI_Comm comm, int *size);
int SMPI_Barrier(SMPI_Comm comm);
int SMPI_Bcast(void *buffer, int count, SMPI_Datatype datatype, int root, SMPI_Comm comm);
int SMPI_Gather(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                int recvcount, SMPI_Datatype recvtype, int root, SMPI_Comm comm);
int SMPI_Gatherv(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                 const size_t *recvcounts, const size_t *displs, SMPI_Datatype recvtype, int root,
                 SMPI_Comm comm);
int SMPI_Allgather(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                   int recvcount, SMPI_Datatype recvtype, SMPI_Comm comm);
int SMPI_Allgatherv(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                    const size_t *recvcounts, const size_t *displs, SMPI_Datatype recvtype,
                    SMPI_Comm comm);
int SMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, SMPI_Datatype datatype,
                   SMPI_Op op, SMPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif /* ADIOS2_TOOLKIT_SST_SST_COMM_H_ */
