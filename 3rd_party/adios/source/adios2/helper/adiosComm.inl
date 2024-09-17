/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosComm.inl
 */

#ifndef ADIOS2_HELPER_ADIOSCOMM_INL_
#define ADIOS2_HELPER_ADIOSCOMM_INL_
#ifndef ADIOS2_HELPER_ADIOSCOMM_H_
#error "Inline file should only be included from it's header, never on it's own"
#endif

#include "adios2/helper/adiosLog.h"

#include <numeric>   //std::accumulate
#include <stdexcept> //std::runtime_error
#include <string>    //std::to_string
#include <utility>   //std::pair

namespace adios2
{
namespace helper
{

template <class T>
void Comm::GatherArrays(const T *source, size_t sourceCount, T *destination,
                        int rankDestination) const
{
    this->Gather(source, sourceCount, destination, sourceCount,
                 rankDestination);
}

template <class T>
std::vector<T> Comm::GatherValues(T source, int rankDestination) const
{
    int rank = this->Rank();
    int size = this->Size();

    std::vector<T> output;

    if (rank == rankDestination) // pre-allocate in destination rank
    {
        output.resize(size);
    }

    T sourceCopy = source; // so we can have an address for rvalues
    this->GatherArrays(&sourceCopy, 1, output.data(), rankDestination);

    return output;
}

template <class T>
void Comm::GathervArrays(const T *source, size_t sourceCount,
                         const size_t *counts, size_t countsSize,
                         T *destination, int rankDestination) const
{
    std::vector<size_t> displs;
    if (rankDestination == this->Rank())
    {
        displs = GetGathervDisplacements(counts, countsSize);
        const size_t totalElements =
            displs[countsSize - 1] + counts[countsSize - 1];
        if (totalElements > 2147483648)
        {
            helper::ThrowNested<std::runtime_error>(
                "Helper", "adiosComm", "GathervVectors",
                "ERROR: GathervArrays does not support gathering more than "
                "2^31 elements. Here it was tasked with " +
                    std::to_string(totalElements) + " elements\n");
        }
    }
    this->Gatherv(source, sourceCount, destination, counts, displs.data(),
                  rankDestination);
}

template <class T>
void Comm::GathervVectors(const std::vector<T> &in, std::vector<T> &out,
                          size_t &position, int rankDestination) const
{
    const size_t inSize = in.size();
    const std::vector<size_t> counts =
        this->GatherValues(inSize, rankDestination);

    size_t gatheredSize = 0;

    int rank = this->Rank();

    if (rank == rankDestination) // pre-allocate vector
    {
        gatheredSize = std::accumulate(counts.begin(), counts.end(), size_t(0));

        const size_t newSize = position + gatheredSize;
        try
        {
            out.reserve(newSize); // to avoid power of 2 growth
            out.resize(newSize);
        }
        catch (...)
        {
            helper::ThrowNested<std::runtime_error>(
                "Helper", "adiosComm", "GathervVectors",
                "buffer overflow when resizing to " + std::to_string(newSize));
        }
    }

    this->GathervArrays(in.data(), in.size(), counts.data(), counts.size(),
                        out.data() + position, rankDestination);
    position += gatheredSize;
}

template <class T>
std::vector<T> Comm::AllGatherValues(const T source) const
{
    int size = this->Size();
    std::vector<T> output(size);

    T sourceCopy = source; // so we can have an address for rvalues
    this->Allgather(&sourceCopy, 1, output.data(), 1);
    return output;
}

template <class T>
T Comm::ReduceValues(const T source, Op op, const int rankDestination) const
{
    T sourceLocal = source;
    T reduceValue = 0;
    this->Reduce(&sourceLocal, &reduceValue, 1, op, rankDestination);
    return reduceValue;
}

template <class T>
T Comm::BroadcastValue(const T &input, const int rankSource) const
{
    T output = 0;
    if (rankSource == this->Rank())
    {
        output = input;
    }

    this->Bcast(&output, 1, rankSource);

    return output;
}

// BroadcastValue full specializations implemented in 'adiosComm.tcc'.
template <>
std::string Comm::BroadcastValue(const std::string &input,
                                 const int rankSource) const;

template <class T>
void Comm::BroadcastVector(std::vector<T> &vector, const int rankSource) const
{
    if (this->Size() == 1)
    {
        return;
    }

    // First Broadcast the size, then the contents
    size_t inputSize = this->BroadcastValue(vector.size(), rankSource);

    if (rankSource != this->Rank())
    {
        vector.resize(inputSize);
    }

    if (inputSize > 0)
    {
        this->Bcast(vector.data(), inputSize, rankSource);
    }
}

template <typename TSend, typename TRecv>
void Comm::Allgather(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf,
                     size_t recvcount, const std::string &hint) const
{
    return m_Impl->Allgather(sendbuf, sendcount, CommImpl::GetDatatype<TSend>(),
                             recvbuf, recvcount, CommImpl::GetDatatype<TRecv>(),
                             hint);
}

template <typename TSend, typename TRecv>
void Comm::Allgatherv(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf,
                      const size_t *recvcounts, const size_t *displs,
                      const std::string &hint) const
{
    return m_Impl->Allgatherv(
        sendbuf, sendcount, CommImpl::GetDatatype<TSend>(), recvbuf, recvcounts,
        displs, CommImpl::GetDatatype<TRecv>(), hint);
}

template <typename T>
void Comm::Allreduce(const T *sendbuf, T *recvbuf, size_t count, Op op,
                     const std::string &hint) const
{
    return m_Impl->Allreduce(sendbuf, recvbuf, count,
                             CommImpl::GetDatatype<T>(), op, hint);
}

template <typename T>
void Comm::Bcast(T *buffer, const size_t count, int root,
                 const std::string &hint) const
{
    return m_Impl->Bcast(buffer, count, CommImpl::GetDatatype<T>(), root, hint);
}

template <typename TSend, typename TRecv>
void Comm::Gather(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf,
                  size_t recvcount, int root, const std::string &hint) const
{
    return m_Impl->Gather(sendbuf, sendcount, CommImpl::GetDatatype<TSend>(),
                          recvbuf, recvcount, CommImpl::GetDatatype<TRecv>(),
                          root, hint);
}

template <typename TSend, typename TRecv>
void Comm::Gatherv(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf,
                   const size_t *recvcounts, const size_t *displs, int root,
                   const std::string &hint) const
{
    return m_Impl->Gatherv(sendbuf, sendcount, CommImpl::GetDatatype<TSend>(),
                           recvbuf, recvcounts, displs,
                           CommImpl::GetDatatype<TRecv>(), root, hint);
}

template <typename T>
void Comm::Reduce(const T *sendbuf, T *recvbuf, size_t count, Op op, int root,
                  const std::string &hint) const
{
    return m_Impl->Reduce(sendbuf, recvbuf, count, CommImpl::GetDatatype<T>(),
                          op, root, hint);
}

template <typename T>
void Comm::ReduceInPlace(T *buf, size_t count, Op op, int root,
                         const std::string &hint) const
{
    return m_Impl->ReduceInPlace(buf, count, CommImpl::GetDatatype<T>(), op,
                                 root, hint);
}

template <typename T>
void Comm::Send(const T *buf, size_t count, int dest, int tag,
                const std::string &hint) const
{
    if (dest < 0 || dest > m_Impl->Size() - 1)
    {
        throw std::runtime_error(
            "Invalid MPI dest rank in Send: " + std::to_string(dest) +
            " for a communicator of size " + std::to_string(m_Impl->Size()));
    }
    return m_Impl->Send(buf, count, CommImpl::GetDatatype<T>(), dest, tag,
                        hint);
}

template <typename T>
Comm::Status Comm::Recv(T *buf, size_t count, int source, int tag,
                        const std::string &hint) const
{
    if (source < 0 || source > m_Impl->Size() - 1)
    {
        throw std::runtime_error(
            "Invalid MPI source rank in Recv: " + std::to_string(source) +
            " for a communicator of size " + std::to_string(m_Impl->Size()));
    }
    return m_Impl->Recv(buf, count, CommImpl::GetDatatype<T>(), source, tag,
                        hint);
}

template <typename TSend, typename TRecv>
void Comm::Scatter(const TSend *sendbuf, size_t sendcount, TRecv *recvbuf,
                   size_t recvcount, int root, const std::string &hint) const
{
    return m_Impl->Scatter(sendbuf, sendcount, CommImpl::GetDatatype<TSend>(),
                           recvbuf, recvcount, CommImpl::GetDatatype<TRecv>(),
                           root, hint);
}

template <typename T>
Comm::Req Comm::Isend(const T *buffer, const size_t count, int dest, int tag,
                      const std::string &hint) const
{
    if (dest < 0 || dest > m_Impl->Size() - 1)
    {
        throw std::runtime_error(
            "Invalid MPI dest rank in Isend: " + std::to_string(dest) +
            " for a communicator of size " + std::to_string(m_Impl->Size()));
    }
    return m_Impl->Isend(buffer, count, CommImpl::GetDatatype<T>(), dest, tag,
                         hint);
}

template <typename T>
Comm::Req Comm::Irecv(T *buffer, const size_t count, int source, int tag,
                      const std::string &hint) const
{
    if (source < 0 || source > m_Impl->Size() - 1)
    {
        throw std::runtime_error(
            "Invalid MPI source rank in Irecv: " + std::to_string(source) +
            " for a communicator of size " + std::to_string(m_Impl->Size()));
    }
    return m_Impl->Irecv(buffer, count, CommImpl::GetDatatype<T>(), source, tag,
                         hint);
}

// CommImpl::GetDatatype full specializations implemented in 'adiosComm.tcc'.
template <>
CommImpl::Datatype CommImpl::GetDatatype<signed char>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<char>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<short>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<int>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<long>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned char>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned short>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned int>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned long>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned long long>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<long long>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<double>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<long double>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<int, int>>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<float, int>>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<double, int>>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<long double, int>>();
template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<short, int>>();

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSCOMM_INL_ */
