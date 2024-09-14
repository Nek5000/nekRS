
/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 */

#include "sst_comm.h"

#include "adios2/helper/adiosComm.h"

#include <algorithm>
#include <vector>

#define CASE_FOR_EACH_TYPE(F)                                                                      \
    case SMPI_INT:                                                                                 \
        F(int);                                                                                    \
        break;                                                                                     \
    case SMPI_CHAR:                                                                                \
        F(char);                                                                                   \
        break;                                                                                     \
    case SMPI_LONG:                                                                                \
        F(long);                                                                                   \
        break;                                                                                     \
    case SMPI_SIZE_T:                                                                              \
        F(size_t);                                                                                 \
        break;                                                                                     \
    case SMPI_BYTE:                                                                                \
        F(unsigned char);                                                                          \
        break

#define CASE_FOR_EACH_OP(F)                                                                        \
    case SMPI_MAX:                                                                                 \
        F(adios2::helper::Comm::Op::Max);                                                          \
        break;                                                                                     \
    case SMPI_LAND:                                                                                \
        F(adios2::helper::Comm::Op::LogicalAnd);                                                   \
        break

int SMPI_Comm_rank(SMPI_Comm comm, int *rank)
{
    *rank = comm->Rank();
    return 0;
}

int SMPI_Comm_size(SMPI_Comm comm, int *size)
{
    *size = comm->Size();
    return 0;
}

int SMPI_Barrier(SMPI_Comm comm)
{
    comm->Barrier();
    return 0;
}

int SMPI_Bcast(void *buffer, int count, SMPI_Datatype datatype, int root, SMPI_Comm comm)
{
    switch (datatype)
    {
#define F(type) comm->Bcast(static_cast<type *>(buffer), count, root)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return 0;
}

namespace
{

template <typename TSend>
int SMPI_Gather_Impl(const TSend *sendbuf, int sendcount, void *recvbuf, int recvcount,
                     SMPI_Datatype recvtype, int root, SMPI_Comm comm)
{
    switch (recvtype)
    {
#define F(type) comm->Gather(sendbuf, sendcount, static_cast<type *>(recvbuf), recvcount, root)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return 0;
}

template <typename TSend>
int SMPI_Gatherv_Impl(const TSend *sendbuf, int sendcount, void *recvbuf, const size_t *recvcounts,
                      const size_t *displs, SMPI_Datatype recvtype, int root, SMPI_Comm comm)
{
    switch (recvtype)
    {
#define F(type)                                                                                    \
    comm->Gatherv(sendbuf, sendcount, static_cast<type *>(recvbuf), recvcounts, displs, root)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return 0;
}

template <typename TSend>
int SMPI_Allgather_Impl(const TSend *sendbuf, int sendcount, void *recvbuf, int recvcount,
                        SMPI_Datatype recvtype, SMPI_Comm comm)
{
    switch (recvtype)
    {
#define F(type) comm->Allgather(sendbuf, sendcount, static_cast<type *>(recvbuf), recvcount)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return 0;
}

template <typename TSend>
int SMPI_Allgatherv_Impl(const TSend *sendbuf, int sendcount, void *recvbuf,
                         const size_t *recvcounts, const size_t *displs, SMPI_Datatype recvtype,
                         SMPI_Comm comm)
{
    switch (recvtype)
    {
#define F(type)                                                                                    \
    comm->Allgatherv(sendbuf, sendcount, static_cast<type *>(recvbuf), recvcounts, displs)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return 0;
}

template <typename T>
int SMPI_Allreduce_Impl(const T *sendbuf, T *recvbuf, int count, SMPI_Op op, SMPI_Comm comm)
{
    switch (op)
    {
#define F(op) comm->Allreduce(sendbuf, recvbuf, count, op)
        CASE_FOR_EACH_OP(F);
#undef F
    }
    return 0;
}
}

int SMPI_Gather(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                int recvcount, SMPI_Datatype recvtype, int root, SMPI_Comm comm)
{
    int ret = 0;
    switch (sendtype)
    {
#define F(type)                                                                                    \
    ret = SMPI_Gather_Impl(static_cast<const type *>(sendbuf), sendcount, recvbuf, recvcount,      \
                           recvtype, root, comm)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return ret;
}

int SMPI_Gatherv(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                 const size_t *recvcounts, const size_t *displs, SMPI_Datatype recvtype, int root,
                 SMPI_Comm comm)
{
    int ret = 0;
    switch (sendtype)
    {
#define F(type)                                                                                    \
    ret = SMPI_Gatherv_Impl(static_cast<const type *>(sendbuf), sendcount, recvbuf, recvcounts,    \
                            displs, recvtype, root, comm)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return ret;
}

int SMPI_Allgather(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                   int recvcount, SMPI_Datatype recvtype, SMPI_Comm comm)
{
    int ret = 0;
    switch (sendtype)
    {
#define F(type)                                                                                    \
    ret = SMPI_Allgather_Impl(static_cast<const type *>(sendbuf), sendcount, recvbuf, recvcount,   \
                              recvtype, comm)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return ret;
}

int SMPI_Allgatherv(const void *sendbuf, int sendcount, SMPI_Datatype sendtype, void *recvbuf,
                    const size_t *recvcounts, const size_t *displs, SMPI_Datatype recvtype,
                    SMPI_Comm comm)
{
    int ret = 0;
    switch (sendtype)
    {
#define F(type)                                                                                    \
    ret = SMPI_Allgatherv_Impl(static_cast<const type *>(sendbuf), sendcount, recvbuf, recvcounts, \
                               displs, recvtype, comm)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return ret;
}

int SMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, SMPI_Datatype datatype,
                   SMPI_Op op, SMPI_Comm comm)
{
    switch (datatype)
    {
#define F(type)                                                                                    \
    SMPI_Allreduce_Impl(static_cast<const type *>(sendbuf), static_cast<type *>(recvbuf), count,   \
                        op, comm)
        CASE_FOR_EACH_TYPE(F);
#undef F
    }
    return 0;
}
