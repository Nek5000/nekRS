/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscHelper.h
 *
 *  Created on: Sep 30, 2019
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCHELPER_H_
#define ADIOS2_ENGINE_SSCHELPER_H_

#include "adios2/core/IO.h"
#include <mpi.h>

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

class Buffer
{
public:
    Buffer(const size_t capacity = 1)
    {
        m_Buffer = reinterpret_cast<uint8_t *>(malloc(capacity));
        m_Capacity = capacity;
        m_Size = 0;
    }
    ~Buffer()
    {
        if (m_Buffer)
        {
            free(m_Buffer);
        }
    }
    void clear()
    {
        m_Buffer = reinterpret_cast<uint8_t *>(realloc(m_Buffer, 1));
        m_Capacity = 1;
        m_Size = 0;
    }
    void resize(const size_t size)
    {
        m_Size = size;
        if (size > m_Capacity)
        {
            m_Capacity = size * 2;
            m_Buffer = reinterpret_cast<uint8_t *>(realloc(m_Buffer, m_Capacity));
        }
        if (m_Buffer == nullptr)
        {
            helper::Throw<std::runtime_error>("Engine", "SscHelper", "resize",
                                              "ssc buffer realloc failed");
        }
    }
    template <typename T>
    T &value(const size_t pos = 0)
    {
        return *reinterpret_cast<T *>(m_Buffer + pos);
    }
    uint8_t &value(const size_t pos = 0) { return *reinterpret_cast<uint8_t *>(m_Buffer + pos); }
    template <typename T>
    T value(const size_t pos = 0) const
    {
        return *reinterpret_cast<T *>(m_Buffer + pos);
    }
    uint8_t value(const size_t pos = 0) const
    {
        return *reinterpret_cast<uint8_t *>(m_Buffer + pos);
    }
    template <typename T>
    T *data(const size_t pos = 0)
    {
        return reinterpret_cast<T *>(m_Buffer + pos);
    }
    template <typename T>
    const T *data(const size_t pos = 0) const
    {
        return reinterpret_cast<const T *>(m_Buffer + pos);
    }
    uint8_t *data(const size_t pos = 0) { return reinterpret_cast<uint8_t *>(m_Buffer + pos); }
    const uint8_t *data(const size_t pos = 0) const
    {
        return reinterpret_cast<const uint8_t *>(m_Buffer + pos);
    }
    size_t size() const { return m_Size; }
    uint8_t &operator[](const size_t pos) { return *(m_Buffer + pos); }
    const uint8_t &operator[](const size_t pos) const { return *(m_Buffer + pos); }

private:
    size_t m_Capacity = 0;
    size_t m_Size = 0;
    uint8_t *m_Buffer = nullptr;
};

struct BlockInfo
{
    std::string name;
    std::string structDef;
    DataType type;
    size_t elementSize;
    ShapeID shapeId;
    Dims shape;
    Dims start;
    Dims count;
    Dims memStart;
    Dims memCount;
    size_t bufferStart;
    size_t bufferCount;
    std::vector<char> value;
    void *data;
    bool performed;
};
using BlockVec = std::vector<BlockInfo>;
using BlockVecVec = std::vector<BlockVec>;
using RankPosMap = std::unordered_map<int, std::pair<size_t, size_t>>;
using MpiInfo = std::vector<std::vector<int>>;

void PrintDims(const Dims &dims, const std::string &label = std::string());
void PrintBlock(const BlockInfo &b, const std::string &label = std::string());
void PrintBlockVec(const BlockVec &bv, const std::string &label = std::string());
void PrintBlockVecVec(const BlockVecVec &bvv, const std::string &label = std::string());
void PrintRankPosMap(const RankPosMap &m, const std::string &label = std::string());
void PrintMpiInfo(const MpiInfo &writersInfo, const MpiInfo &readersInfo);

size_t TotalDataSize(const Dims &dims, const size_t elementSize, const ShapeID &shapeId);
size_t TotalDataSize(const BlockVec &bv);

RankPosMap CalculateOverlap(BlockVecVec &globalPattern, const BlockVec &localPattern);

void SerializeVariables(const BlockVec &input, Buffer &output, const int rank);
void SerializeAttributes(IO &input, Buffer &output);
void SerializeStructDefinitions(
    const std::unordered_multimap<std::string, StructDefinition> &definitions, Buffer &output);
void DeserializeVariable(const Buffer &input, const ShapeID shapeId, uint64_t &pos, BlockInfo &b,
                         IO &io, const bool regIO,
                         std::unordered_map<std::string, StructDefinition> &StructDefs);
void DeserializeAttribute(const Buffer &input, uint64_t &pos, IO &io, const bool regIO);
void DeserializeStructDefinitions(const Buffer &input, uint64_t &pos, IO &io, const bool regIO,
                                  std::unordered_map<std::string, StructDefinition> &StructDefs);
void Deserialize(const Buffer &input, BlockVecVec &output, IO &io, const bool regVars,
                 const bool regAttrs, const bool regDefs,
                 std::unordered_map<std::string, StructDefinition> &StructDefs);
void AggregateMetadata(const Buffer &localBuffer, Buffer &globalBuffer, MPI_Comm comm,
                       const bool finalStep, const bool locked);
void BroadcastMetadata(Buffer &globalBuffer, const int root, MPI_Comm comm);

void MPI_Gatherv64OneSidedPush(const void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const uint64_t *recvcounts, const uint64_t *displs,
                               MPI_Datatype recvtype, int root, MPI_Comm comm,
                               const int chunksize = std::numeric_limits<int>::max());
void MPI_Gatherv64OneSidedPull(const void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const uint64_t *recvcounts, const uint64_t *displs,
                               MPI_Datatype recvtype, int root, MPI_Comm comm,
                               const int chunksize = std::numeric_limits<int>::max());
void MPI_Gatherv64(const void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const uint64_t *recvcounts, const uint64_t *displs, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, const int chunksize = std::numeric_limits<int>::max());

bool AreSameDims(const Dims &a, const Dims &b);

} // end namespace ssc
} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif // ADIOS2_ENGINE_SSCHELPER_H_
