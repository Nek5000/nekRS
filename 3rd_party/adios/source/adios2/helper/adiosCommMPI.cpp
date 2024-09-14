/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosCommMPI.cpp
 */

#include "adiosCommMPI.h"
#include "adiosLog.h"

#include <algorithm>
#include <ios> //std::ios_base::failure
#include <iterator>
#include <utility>

#include "adiosComm.h"
#include "adiosCommDummy.h"

#include "adios2/common/ADIOSTypes.h"

#include <mpi.h>

namespace adios2
{
namespace core
{
void RegisterMPIEngines();
}
namespace helper
{

namespace
{
struct InitMPI
{
    InitMPI() { core::RegisterMPIEngines(); }
};

const MPI_Op OpToMPI[] = {
    MPI_OP_NULL, MPI_MAX,  MPI_MIN,  MPI_SUM,    MPI_PROD,   MPI_LAND,    MPI_BAND,  MPI_LOR,
    MPI_BOR,     MPI_LXOR, MPI_BXOR, MPI_MAXLOC, MPI_MINLOC, MPI_REPLACE, MPI_NO_OP,
};

MPI_Op ToMPI(Comm::Op op) { return OpToMPI[int(op)]; }

const int LockTypeToMPI[] = {MPI_LOCK_EXCLUSIVE, MPI_LOCK_SHARED};

int ToMPI(Comm::LockType lock_type) { return LockTypeToMPI[int(lock_type)]; }

const MPI_Datatype DatatypeToMPI[] = {
    MPI_SIGNED_CHAR,
    MPI_CHAR,
    MPI_SHORT,
    MPI_INT,
    MPI_LONG,
    MPI_UNSIGNED_CHAR,
    MPI_UNSIGNED_SHORT,
    MPI_UNSIGNED,
    MPI_UNSIGNED_LONG,
    MPI_UNSIGNED_LONG_LONG,
    MPI_LONG_LONG_INT,
    MPI_DOUBLE,
    MPI_LONG_DOUBLE,
    MPI_2INT,
    MPI_FLOAT_INT,
    MPI_DOUBLE_INT,
    MPI_LONG_DOUBLE_INT,
    MPI_SHORT_INT,
};

MPI_Datatype ToMPI(CommImpl::Datatype dt) { return DatatypeToMPI[int(dt)]; }

void CheckMPIReturn(const int value, const std::string &hint)
{
    if (value == MPI_SUCCESS)
    {
        return;
    }

    std::string error;
    if (value == MPI_ERR_COMM)
    {
        error = "MPI_ERR_COMM";
    }
    else if (value == MPI_ERR_INTERN)
    {
        error = "MPI_ERR_INTERN";
    }
    else
    {
        error = "MPI_ERR number: " + std::to_string(value);
    }

    helper::Throw<std::runtime_error>("Helper", "AdiosCommMPI", "CheckMPIReturn",
                                      "ADIOS2 detected " + error + ", " + hint);
}
}

class CommReqImplMPI : public CommReqImpl
{
public:
    CommReqImplMPI(MPI_Datatype datatype) : m_MPIDatatype(datatype) {}
    ~CommReqImplMPI() override;

    Comm::Status Wait(const std::string &hint) override;

    /** Encapsulated MPI datatype of the requested operation.  */
    MPI_Datatype m_MPIDatatype = MPI_DATATYPE_NULL;

    /** Encapsulated MPI request instances.  There may be more than
     *  one when we batch requests too large for MPI interfaces.  */
    std::vector<MPI_Request> m_MPIReqs;
};

CommReqImplMPI::~CommReqImplMPI() = default;

class CommWinImplMPI : public CommWinImpl
{
public:
    CommWinImplMPI() {}
    ~CommWinImplMPI() override;

    int Free(const std::string &hint) override;

    MPI_Win m_Win;
};

CommWinImplMPI::~CommWinImplMPI() = default;

class CommImplMPI : public CommImpl
{
public:
    CommImplMPI(MPI_Comm mpiComm) : m_MPIComm(mpiComm) {}

    MPI_Comm m_MPIComm;

    ~CommImplMPI() override;

    void Free(const std::string &hint) override;
    std::unique_ptr<CommImpl> Duplicate(const std::string &hint) const override;
    std::unique_ptr<CommImpl> Split(int color, int key, const std::string &hint) const override;
    std::unique_ptr<CommImpl> World(const std::string &hint) const override;
    virtual std::unique_ptr<CommImpl> GroupByShm(const std::string &hint) const override;

    int Rank() const override;
    int Size() const override;
    bool IsMPI() const override;
    void Barrier(const std::string &hint) const override;

    void Allgather(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                   size_t recvcount, Datatype recvtype, const std::string &hint) const override;

    void Allgatherv(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                    const size_t *recvcounts, const size_t *displs, Datatype recvtype,
                    const std::string &hint) const override;

    void Allreduce(const void *sendbuf, void *recvbuf, size_t count, Datatype datatype, Comm::Op op,
                   const std::string &hint) const override;

    void Bcast(void *buffer, size_t count, Datatype datatype, int root,
               const std::string &hint) const override;

    void Gather(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                size_t recvcount, Datatype recvtype, int root,
                const std::string &hint) const override;

    void Gatherv(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                 const size_t *recvcounts, const size_t *displs, Datatype recvtype, int root,
                 const std::string &hint) const override;

    void Reduce(const void *sendbuf, void *recvbuf, size_t count, Datatype datatype, Comm::Op op,
                int root, const std::string &hint) const override;

    void ReduceInPlace(void *buf, size_t count, Datatype datatype, Comm::Op op, int root,
                       const std::string &hint) const override;

    void Send(const void *buf, size_t count, Datatype datatype, int dest, int tag,
              const std::string &hint) const override;

    Comm::Status Recv(void *buf, size_t count, Datatype datatype, int source, int tag,
                      const std::string &hint) const override;

    void Scatter(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                 size_t recvcount, Datatype recvtype, int root,
                 const std::string &hint) const override;

    Comm::Req Isend(const void *buffer, size_t count, Datatype datatype, int dest, int tag,
                    const std::string &hint) const override;

    Comm::Req Irecv(void *buffer, size_t count, Datatype datatype, int source, int tag,
                    const std::string &hint) const override;

    Comm::Win Win_allocate_shared(size_t size, int disp_unit, void *baseptr,
                                  const std::string &hint) const override;
    int Win_shared_query(Comm::Win &win, int rank, size_t *size, int *disp_unit, void *baseptr,
                         const std::string &hint) const override;
    int Win_free(Comm::Win &win, const std::string &hint) const override;
    int Win_lock(Comm::LockType lock_type, int rank, int assert, Comm::Win &win,
                 const std::string &hint) const override;
    int Win_unlock(int rank, Comm::Win &win, const std::string &hint) const override;
    int Win_lock_all(int assert, Comm::Win &win, const std::string &hint) const override;
    int Win_unlock_all(Comm::Win &win, const std::string &hint) const override;
};

CommImplMPI::~CommImplMPI()
{
    // Handle the case where MPI is finalized before the ADIOS destructor is
    // called, which happens, e.g., with global / static ADIOS objects
    int flag;
    MPI_Finalized(&flag);
    if (!flag)
    {
        if (m_MPIComm != MPI_COMM_NULL && m_MPIComm != MPI_COMM_WORLD && m_MPIComm != MPI_COMM_SELF)
        {
            MPI_Comm_free(&m_MPIComm);
        }
    }
}

void CommImplMPI::Free(const std::string &hint)
{
    if (m_MPIComm != MPI_COMM_NULL && m_MPIComm != MPI_COMM_WORLD && m_MPIComm != MPI_COMM_SELF)
    {
        MPI_Comm mpiComm = m_MPIComm;
        m_MPIComm = MPI_COMM_NULL; // prevent freeing a second time
        CheckMPIReturn(MPI_Comm_free(&mpiComm), hint);
    }
}

std::unique_ptr<CommImpl> CommImplMPI::Duplicate(const std::string &hint) const
{
    MPI_Comm newComm;
    CheckMPIReturn(MPI_Comm_dup(m_MPIComm, &newComm), hint);
    return std::unique_ptr<CommImpl>(new CommImplMPI(newComm));
}

std::unique_ptr<CommImpl> CommImplMPI::Split(int color, int key, const std::string &hint) const
{
    MPI_Comm newComm;
    CheckMPIReturn(MPI_Comm_split(m_MPIComm, color, key, &newComm), hint);
    return std::unique_ptr<CommImpl>(new CommImplMPI(newComm));
}

std::unique_ptr<CommImpl> CommImplMPI::World(const std::string &) const
{
    return std::unique_ptr<CommImpl>(new CommImplMPI(MPI_COMM_WORLD));
}

std::unique_ptr<CommImpl> CommImplMPI::GroupByShm(const std::string &hint) const
{
    MPI_Comm nodeComm;
    MPI_Info info;
    MPI_Info_create(&info);
    CheckMPIReturn(MPI_Comm_split_type(m_MPIComm, MPI_COMM_TYPE_SHARED, 0, info, &nodeComm), hint);
    return std::unique_ptr<CommImpl>(new CommImplMPI(nodeComm));
}

int CommImplMPI::Rank() const
{
    int rank;
    CheckMPIReturn(MPI_Comm_rank(m_MPIComm, &rank), {});
    return rank;
}

int CommImplMPI::Size() const
{
    int size;
    CheckMPIReturn(MPI_Comm_size(m_MPIComm, &size), {});
    return size;
}

bool CommImplMPI::IsMPI() const { return true; }

void CommImplMPI::Barrier(const std::string &hint) const
{
    CheckMPIReturn(MPI_Barrier(m_MPIComm), hint);
}

void CommImplMPI::Allgather(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                            size_t recvcount, Datatype recvtype, const std::string &hint) const
{
    CheckMPIReturn(MPI_Allgather(sendbuf, static_cast<int>(sendcount), ToMPI(sendtype), recvbuf,
                                 static_cast<int>(recvcount), ToMPI(recvtype), m_MPIComm),
                   hint);
}

void CommImplMPI::Allgatherv(const void *sendbuf, size_t sendcount, Datatype sendtype,
                             void *recvbuf, const size_t *recvcounts, const size_t *displs,
                             Datatype recvtype, const std::string &hint) const
{
    std::vector<int> countsInt;
    std::vector<int> displsInt;
    {
        auto cast = [](size_t sz) -> int { return int(sz); };
        const int size = this->Size();
        countsInt.reserve(size);
        std::transform(recvcounts, recvcounts + size, std::back_inserter(countsInt), cast);
        displsInt.reserve(size);
        std::transform(displs, displs + size, std::back_inserter(displsInt), cast);
    }
    CheckMPIReturn(MPI_Allgatherv(sendbuf, static_cast<int>(sendcount), ToMPI(sendtype), recvbuf,
                                  countsInt.data(), displsInt.data(), ToMPI(recvtype), m_MPIComm),
                   hint);
}

void CommImplMPI::Allreduce(const void *sendbuf, void *recvbuf, size_t count, Datatype datatype,
                            Comm::Op op, const std::string &hint) const
{
    CheckMPIReturn(MPI_Allreduce(sendbuf, recvbuf, static_cast<int>(count), ToMPI(datatype),
                                 ToMPI(op), m_MPIComm),
                   hint);
}

void CommImplMPI::Bcast(void *buffer, size_t count, Datatype datatype, int root,
                        const std::string &hint) const
{
    size_t inputSize = count;
    const int MAXBCASTSIZE = 1073741824;
    size_t blockSize = (inputSize > MAXBCASTSIZE ? MAXBCASTSIZE : inputSize);
    unsigned char *blockBuf = static_cast<unsigned char *>(buffer);
    while (inputSize > 0)
    {
        CheckMPIReturn(
            MPI_Bcast(blockBuf, static_cast<int>(blockSize), ToMPI(datatype), root, m_MPIComm),
            hint);
        blockBuf += blockSize * CommImpl::SizeOf(datatype);
        inputSize -= blockSize;
        blockSize = (inputSize > MAXBCASTSIZE ? MAXBCASTSIZE : inputSize);
    }
}

void CommImplMPI::Gather(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                         size_t recvcount, Datatype recvtype, int root,
                         const std::string &hint) const
{
    CheckMPIReturn(MPI_Gather(sendbuf, static_cast<int>(sendcount), ToMPI(sendtype), recvbuf,
                              static_cast<int>(recvcount), ToMPI(recvtype), root, m_MPIComm),
                   hint);
}

void CommImplMPI::Gatherv(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                          const size_t *recvcounts, const size_t *displs, Datatype recvtype,
                          int root, const std::string &hint) const
{
    std::vector<int> countsInt;
    std::vector<int> displsInt;
    if (root == this->Rank())
    {
        auto cast = [](size_t sz) -> int { return int(sz); };
        const int size = this->Size();
        countsInt.reserve(size);
        std::transform(recvcounts, recvcounts + size, std::back_inserter(countsInt), cast);
        displsInt.reserve(size);
        std::transform(displs, displs + size, std::back_inserter(displsInt), cast);
    }
    CheckMPIReturn(MPI_Gatherv(sendbuf, static_cast<int>(sendcount), ToMPI(sendtype), recvbuf,
                               countsInt.data(), displsInt.data(), ToMPI(recvtype), root,
                               m_MPIComm),
                   hint);
}

void CommImplMPI::Reduce(const void *sendbuf, void *recvbuf, size_t count, Datatype datatype,
                         Comm::Op op, int root, const std::string &hint) const
{
    CheckMPIReturn(MPI_Reduce(sendbuf, recvbuf, static_cast<int>(count), ToMPI(datatype), ToMPI(op),
                              root, m_MPIComm),
                   hint);
}

void CommImplMPI::ReduceInPlace(void *buf, size_t count, Datatype datatype, Comm::Op op, int root,
                                const std::string &hint) const
{
    CheckMPIReturn(MPI_Reduce(MPI_IN_PLACE, buf, static_cast<int>(count), ToMPI(datatype),
                              ToMPI(op), root, m_MPIComm),
                   hint);
}

void CommImplMPI::Send(const void *buf, size_t count, Datatype datatype, int dest, int tag,
                       const std::string &hint) const
{
    CheckMPIReturn(MPI_Send(buf, static_cast<int>(count), ToMPI(datatype), dest, tag, m_MPIComm),
                   hint);
}

Comm::Status CommImplMPI::Recv(void *buf, size_t count, Datatype datatype, int source, int tag,
                               const std::string &hint) const
{
    MPI_Status mpiStatus;
    CheckMPIReturn(
        MPI_Recv(buf, static_cast<int>(count), ToMPI(datatype), source, tag, m_MPIComm, &mpiStatus),
        hint);

    Comm::Status status;
    status.Source = mpiStatus.MPI_SOURCE;
    status.Tag = mpiStatus.MPI_TAG;
    {
        int mpiCount = 0;
        CheckMPIReturn(MPI_Get_count(&mpiStatus, ToMPI(datatype), &mpiCount), hint);
        status.Count = mpiCount;
    }
    return status;
}

void CommImplMPI::Scatter(const void *sendbuf, size_t sendcount, Datatype sendtype, void *recvbuf,
                          size_t recvcount, Datatype recvtype, int root,
                          const std::string &hint) const
{
    CheckMPIReturn(MPI_Scatter(sendbuf, static_cast<int>(sendcount), ToMPI(sendtype), recvbuf,
                               static_cast<int>(recvcount), ToMPI(recvtype), root, m_MPIComm),
                   hint);
}

Comm::Req CommImplMPI::Isend(const void *buffer, size_t count, Datatype datatype, int dest, int tag,
                             const std::string &hint) const
{
    auto req = std::unique_ptr<CommReqImplMPI>(new CommReqImplMPI(ToMPI(datatype)));

    if (count > DefaultMaxFileBatchSize)
    {
        const size_t batches = count / DefaultMaxFileBatchSize;

        size_t position = 0;
        for (size_t b = 0; b < batches; ++b)
        {
            int batchSize = static_cast<int>(DefaultMaxFileBatchSize);
            MPI_Request mpiReq;
            CheckMPIReturn(MPI_Isend(static_cast<char *>(const_cast<void *>(buffer)) + position,
                                     batchSize, ToMPI(datatype), dest, tag, m_MPIComm, &mpiReq),
                           "in call to Isend batch " + std::to_string(b) + " " + hint + "\n");
            req->m_MPIReqs.emplace_back(mpiReq);

            position += DefaultMaxFileBatchSize;
        }
        const size_t remainder = count % DefaultMaxFileBatchSize;
        if (remainder > 0)
        {
            int batchSize = static_cast<int>(remainder);
            MPI_Request mpiReq;
            CheckMPIReturn(MPI_Isend(static_cast<char *>(const_cast<void *>(buffer)) + position,
                                     batchSize, ToMPI(datatype), dest, tag, m_MPIComm, &mpiReq),
                           "in call to Isend remainder batch " + hint + "\n");
            req->m_MPIReqs.emplace_back(mpiReq);
        }
    }
    else
    {
        int batchSize = static_cast<int>(count);
        MPI_Request mpiReq;
        CheckMPIReturn(MPI_Isend(static_cast<char *>(const_cast<void *>(buffer)), batchSize,
                                 ToMPI(datatype), dest, tag, m_MPIComm, &mpiReq),
                       " in call to Isend with single batch " + hint + "\n");
        req->m_MPIReqs.emplace_back(mpiReq);
    }

    return MakeReq(std::move(req));
}

Comm::Req CommImplMPI::Irecv(void *buffer, size_t count, Datatype datatype, int source, int tag,
                             const std::string &hint) const
{
    auto req = std::unique_ptr<CommReqImplMPI>(new CommReqImplMPI(ToMPI(datatype)));

    if (count > DefaultMaxFileBatchSize)
    {
        const size_t batches = count / DefaultMaxFileBatchSize;
        size_t position = 0;
        for (size_t b = 0; b < batches; ++b)
        {
            int batchSize = static_cast<int>(DefaultMaxFileBatchSize);
            MPI_Request mpiReq;
            CheckMPIReturn(MPI_Irecv(static_cast<char *>(buffer) + position, batchSize,
                                     ToMPI(datatype), source, tag, m_MPIComm, &mpiReq),
                           "in call to Irecv batch " + std::to_string(b) + " " + hint + "\n");
            req->m_MPIReqs.emplace_back(mpiReq);

            position += DefaultMaxFileBatchSize;
        }

        const size_t remainder = count % DefaultMaxFileBatchSize;
        if (remainder > 0)
        {
            int batchSize = static_cast<int>(remainder);
            MPI_Request mpiReq;
            CheckMPIReturn(MPI_Irecv(static_cast<char *>(buffer) + position, batchSize,
                                     ToMPI(datatype), source, tag, m_MPIComm, &mpiReq),
                           "in call to Irecv remainder batch " + hint + "\n");
            req->m_MPIReqs.emplace_back(mpiReq);
        }
    }
    else
    {
        int batchSize = static_cast<int>(count);
        MPI_Request mpiReq;
        CheckMPIReturn(
            MPI_Irecv(buffer, batchSize, ToMPI(datatype), source, tag, m_MPIComm, &mpiReq),
            " in call to Isend with single batch " + hint + "\n");
        req->m_MPIReqs.emplace_back(mpiReq);
    }

    return MakeReq(std::move(req));
}

Comm::Win CommImplMPI::Win_allocate_shared(size_t size, int disp_unit, void *baseptr,
                                           const std::string &hint) const
{
    auto w = std::unique_ptr<CommWinImplMPI>(new CommWinImplMPI());
    MPI_Aint asize = static_cast<MPI_Aint>(size);
    CheckMPIReturn(
        MPI_Win_allocate_shared(asize, disp_unit, MPI_INFO_NULL, m_MPIComm, baseptr, &w->m_Win),
        "in call to Win_allocate_shared " + hint + "\n");
    return MakeWin(std::move(w));
}

int CommImplMPI::Win_shared_query(Comm::Win &win, int rank, size_t *size, int *disp_unit,
                                  void *baseptr, const std::string &hint) const
{
    CommWinImplMPI *w = dynamic_cast<CommWinImplMPI *>(CommWinImpl::Get(win));
    MPI_Aint asize;
    int ret = MPI_Win_shared_query(w->m_Win, rank, &asize, disp_unit, baseptr);
    CheckMPIReturn(ret, "in call to Win_shared_query " + hint + "\n");
    *size = static_cast<size_t>(asize);
    return ret;
}

int CommImplMPI::Win_free(Comm::Win &win, const std::string &hint) const
{
    CommWinImplMPI *w = dynamic_cast<CommWinImplMPI *>(CommWinImpl::Get(win));
    int ret = MPI_Win_free(&w->m_Win);
    CheckMPIReturn(ret, "in call to Win_free " + hint + "\n");
    return ret;
}

int CommImplMPI::Win_lock(Comm::LockType lock_type, int rank, int assert, Comm::Win &win,
                          const std::string &hint) const
{
    CommWinImplMPI *w = dynamic_cast<CommWinImplMPI *>(CommWinImpl::Get(win));
    int mpi_lock_type = ToMPI(lock_type);
    int ret = MPI_Win_lock(mpi_lock_type, rank, assert, w->m_Win);
    CheckMPIReturn(ret, "in call to Win_Lock " + hint + "\n");
    return ret;
}
int CommImplMPI::Win_unlock(int rank, Comm::Win &win, const std::string &hint) const
{
    CommWinImplMPI *w = dynamic_cast<CommWinImplMPI *>(CommWinImpl::Get(win));
    int ret = MPI_Win_unlock(rank, w->m_Win);
    CheckMPIReturn(ret, "in call to Win_Lock " + hint + "\n");
    return ret;
}

int CommImplMPI::Win_lock_all(int assert, Comm::Win &win, const std::string &hint) const
{
    CommWinImplMPI *w = dynamic_cast<CommWinImplMPI *>(CommWinImpl::Get(win));
    int ret = MPI_Win_lock_all(assert, w->m_Win);
    CheckMPIReturn(ret, "in call to Win_Lock_all " + hint + "\n");
    return ret;
}
int CommImplMPI::Win_unlock_all(Comm::Win &win, const std::string &hint) const
{
    CommWinImplMPI *w = dynamic_cast<CommWinImplMPI *>(CommWinImpl::Get(win));
    int ret = MPI_Win_unlock_all(w->m_Win);
    CheckMPIReturn(ret, "in call to Win_Lock " + hint + "\n");
    return ret;
}

Comm::Status CommReqImplMPI::Wait(const std::string &hint)
{
    Comm::Status status;
    if (m_MPIReqs.empty())
    {
        return status;
    }

    std::vector<MPI_Request> mpiRequests = std::move(m_MPIReqs);
    std::vector<MPI_Status> mpiStatuses(mpiRequests.size());

    if (mpiRequests.size() > 1)
    {
        int mpiReturn = MPI_Waitall(static_cast<int>(mpiRequests.size()), mpiRequests.data(),
                                    mpiStatuses.data());
        if (mpiReturn == MPI_ERR_IN_STATUS)
        {
            for (auto &mpiStatus : mpiStatuses)
            {
                if (mpiStatus.MPI_ERROR != MPI_SUCCESS)
                {
                    mpiReturn = mpiStatus.MPI_ERROR;
                    break;
                }
            }
        }
        CheckMPIReturn(mpiReturn, hint);
    }
    else
    {
        CheckMPIReturn(MPI_Wait(mpiRequests.data(), mpiStatuses.data()), hint);
    }

    // Our batched operation should be from only one source and have one tag.
    status.Source = mpiStatuses.front().MPI_SOURCE;
    status.Tag = mpiStatuses.front().MPI_TAG;

    // Accumulate the total count of our batched operation.
    for (auto &mpiStatus : mpiStatuses)
    {
        int mpiCount = 0;
        CheckMPIReturn(MPI_Get_count(&mpiStatus, m_MPIDatatype, &mpiCount), hint);
        status.Count += mpiCount;
    }

    // Our batched operation was cancelled if any member was cancelled.
    for (auto &mpiStatus : mpiStatuses)
    {
        int mpiCancelled = 0;
        MPI_Test_cancelled(&mpiStatus, &mpiCancelled);
        if (mpiCancelled)
        {
            status.Cancelled = true;
            break;
        }
    }

    return status;
}

int CommWinImplMPI::Free(const std::string &hint) { return MPI_Win_free(&m_Win); }

Comm CommWithMPI(MPI_Comm mpiComm)
{
    static InitMPI const initMPI;
    if (mpiComm == MPI_COMM_NULL)
    {
        return CommDummy();
    }
    auto comm = std::unique_ptr<CommImpl>(new CommImplMPI(mpiComm));
    return CommImpl::MakeComm(std::move(comm));
}

Comm CommDupMPI(MPI_Comm mpiComm)
{
    MPI_Comm newComm;
    if (mpiComm != MPI_COMM_NULL)
    {
        MPI_Comm_dup(mpiComm, &newComm);
    }
    else
    {
        newComm = MPI_COMM_NULL;
    }
    return CommWithMPI(newComm);
}

MPI_Comm CommAsMPI(Comm const &comm)
{
    if (CommImplMPI *mpi = dynamic_cast<CommImplMPI *>(CommImpl::Get(comm)))
    {
        return mpi->m_MPIComm;
    }
    return MPI_COMM_NULL;
}

} // end namespace helper
} // end namespace adios2
