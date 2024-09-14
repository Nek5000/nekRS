/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscHelper.cpp
 *
 *  Created on: Sep 30, 2019
 *      Author: Jason Wang
 */

#include "SscHelper.h"

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

size_t TotalDataSize(const Dims &dims, const size_t elementSize, const ShapeID &shapeId)
{
    if (shapeId == ShapeID::GlobalArray || shapeId == ShapeID::LocalArray)
    {
        return std::accumulate(dims.begin(), dims.end(), elementSize, std::multiplies<size_t>());
    }
    else if (shapeId == ShapeID::GlobalValue || shapeId == ShapeID::LocalValue)
    {
        return elementSize;
    }
    helper::Throw<std::runtime_error>("Engine", "SscHelper", "TotalDataSize",
                                      "ShapeID not supported");
    return 0;
}

size_t TotalDataSize(const BlockVec &bv)
{
    size_t s = 0;
    for (const auto &b : bv)
    {
        if (b.type == DataType::String)
        {
            s += b.bufferCount;
        }
        else
        {
            s += TotalDataSize(b.count, b.elementSize, b.shapeId);
        }
    }
    return s;
}

RankPosMap CalculateOverlap(BlockVecVec &globalVecVec, const BlockVec &localVec)
{
    RankPosMap ret;
    int rank = 0;
    for (auto &rankBlockVec : globalVecVec)
    {
        for (auto &gBlock : rankBlockVec)
        {
            for (auto &lBlock : localVec)
            {
                if (lBlock.name == gBlock.name)
                {
                    if (gBlock.shapeId == ShapeID::GlobalValue)
                    {
                        ret[rank].first = 0;
                    }
                    else if (gBlock.shapeId == ShapeID::GlobalArray)
                    {
                        bool hasOverlap = true;
                        for (size_t i = 0; i < gBlock.start.size(); ++i)
                        {
                            if (gBlock.start[i] + gBlock.count[i] <= lBlock.start[i] ||
                                lBlock.start[i] + lBlock.count[i] <= gBlock.start[i])
                            {
                                hasOverlap = false;
                                break;
                            }
                        }
                        if (hasOverlap)
                        {
                            ret[rank].first = 0;
                        }
                    }
                    else if (gBlock.shapeId == ShapeID::LocalValue)
                    {
                    }
                    else if (gBlock.shapeId == ShapeID::LocalArray)
                    {
                    }
                }
            }
        }
        ++rank;
    }
    return ret;
}

void SerializeVariables(const BlockVec &input, Buffer &output, const int rank)
{
    for (const auto &b : input)
    {
        uint64_t pos = output.value<uint64_t>();
        output.resize(pos + 1024);

        output.value<uint8_t>(pos) = static_cast<uint8_t>(b.shapeId);
        ++pos;

        output.value<int>(pos) = rank;
        pos += 4;

        output.value(pos) = static_cast<uint8_t>(b.name.size());
        ++pos;

        std::memcpy(output.data(pos), b.name.data(), b.name.size());
        pos += b.name.size();

        output.value(pos) = static_cast<uint8_t>(b.type);
        ++pos;

        if (b.type == DataType::Struct)
        {
            output.value(pos) = static_cast<uint8_t>(b.structDef.size());
            ++pos;
            std::memcpy(output.data(pos), b.structDef.data(), b.structDef.size());
            pos += b.structDef.size();
        }

        output.value<uint64_t>(pos) = static_cast<uint64_t>(b.elementSize);
        pos += 8;

        output.value(pos) = static_cast<uint8_t>(b.shape.size());
        ++pos;

        for (const auto &s : b.shape)
        {
            output.value<uint64_t>(pos) = s;
            pos += 8;
        }

        for (const auto &s : b.start)
        {
            output.value<uint64_t>(pos) = s;
            pos += 8;
        }

        for (const auto &s : b.count)
        {
            output.value<uint64_t>(pos) = s;
            pos += 8;
        }

        output.value<uint64_t>(pos) = b.bufferStart;
        pos += 8;

        output.value<uint64_t>(pos) = b.bufferCount;
        pos += 8;

        output.value(pos) = static_cast<uint8_t>(b.value.size());
        ++pos;

        std::memcpy(output.data(pos), b.value.data(), b.value.size());
        pos += b.value.size();

        output.value<uint64_t>() = pos;
    }
}

void DeserializeVariable(const Buffer &input, const ShapeID shapeId, uint64_t &pos, BlockInfo &b,
                         IO &io, const bool regIO,
                         std::unordered_map<std::string, StructDefinition> &StructDefs)
{
    b.shapeId = static_cast<ShapeID>(shapeId);

    uint8_t nameSize = input[pos];
    ++pos;

    std::vector<char> name(nameSize);
    std::memcpy(name.data(), input.data(pos), nameSize);
    b.name = std::string(name.begin(), name.end());
    pos += nameSize;

    b.type = static_cast<DataType>(input[pos]);
    ++pos;

    if (b.type == DataType::Struct)
    {
        uint8_t structDefSize = input.value<uint8_t>(pos);
        ++pos;
        std::vector<char> structDefChar(structDefSize);
        std::memcpy(structDefChar.data(), input.data(pos), structDefSize);
        pos += structDefSize;
        b.structDef = std::string(structDefChar.begin(), structDefChar.end());
    }

    b.elementSize = input.value<uint64_t>(pos);
    pos += 8;

    uint8_t shapeSize = input[pos];
    ++pos;
    b.shape.resize(shapeSize);
    b.start.resize(shapeSize);
    b.count.resize(shapeSize);

    std::memcpy(b.shape.data(), input.data(pos), 8 * shapeSize);
    pos += 8 * shapeSize;

    std::memcpy(b.start.data(), input.data(pos), 8 * shapeSize);
    pos += 8 * shapeSize;

    std::memcpy(b.count.data(), input.data(pos), 8 * shapeSize);
    pos += 8 * shapeSize;

    b.bufferStart = input.value<uint64_t>(pos);
    pos += 8;

    b.bufferCount = input.value<uint64_t>(pos);
    pos += 8;

    uint8_t valueSize = input[pos];
    pos++;
    b.value.resize(valueSize);
    if (valueSize > 0)
    {
        std::memcpy(b.value.data(), input.data() + pos, valueSize);
        pos += valueSize;
    }

    if (regIO)
    {
        if (b.type == DataType::Struct)
        {
            auto v = io.InquireStructVariable(b.name);
            if (!v)
            {
                Dims vStart = b.start;
                Dims vShape = b.shape;
                if (io.m_ArrayOrder != ArrayOrdering::RowMajor)
                {
                    std::reverse(vStart.begin(), vStart.end());
                    std::reverse(vShape.begin(), vShape.end());
                }
                if (b.shapeId == ShapeID::GlobalArray)
                {
                    StructDefinition *def = nullptr;
                    auto it = StructDefs.find(b.structDef);
                    if (it != StructDefs.end())
                    {
                        def = &it->second;
                    }
                    if (def == nullptr)
                    {
                        helper::Throw<std::runtime_error>(
                            "Engine", "SscHelper", "DeserializeVariable",
                            "struct type " + b.structDef + " not defined");
                    }
                    else
                    {
                        io.DefineStructVariable(b.name, *def, vShape, vStart, vShape);
                    }
                }
            }
        }
#define declare_type(T)                                                                            \
    else if (b.type == helper::GetDataType<T>())                                                   \
    {                                                                                              \
        auto v = io.InquireVariable<T>(b.name);                                                    \
        if (!v)                                                                                    \
        {                                                                                          \
            Dims vStart = b.start;                                                                 \
            Dims vShape = b.shape;                                                                 \
            if (io.m_ArrayOrder != ArrayOrdering::RowMajor)                                        \
            {                                                                                      \
                std::reverse(vStart.begin(), vStart.end());                                        \
                std::reverse(vShape.begin(), vShape.end());                                        \
            }                                                                                      \
            if (b.shapeId == ShapeID::GlobalValue)                                                 \
            {                                                                                      \
                io.DefineVariable<T>(b.name);                                                      \
            }                                                                                      \
            else if (b.shapeId == ShapeID::GlobalArray)                                            \
            {                                                                                      \
                io.DefineVariable<T>(b.name, vShape, vStart, vShape);                              \
            }                                                                                      \
            else if (b.shapeId == ShapeID::LocalValue)                                             \
            {                                                                                      \
                io.DefineVariable<T>(b.name, {adios2::LocalValueDim});                             \
            }                                                                                      \
            else if (b.shapeId == ShapeID::LocalArray)                                             \
            {                                                                                      \
                io.DefineVariable<T>(b.name, {}, {}, vShape);                                      \
            }                                                                                      \
        }                                                                                          \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
        else
        {
            helper::Throw<std::runtime_error>("Engine", "SscHelper", "DeserializeVariable",
                                              "unknown variable data type");
        }
    }
}

void SerializeAttributes(IO &input, Buffer &output)
{
    const auto &attributeMap = input.GetAttributes();
    for (const auto &attributePair : attributeMap)
    {
        uint64_t pos = output.value<uint64_t>();
        output.resize(pos + 1024);

        if (attributePair.second->m_Type == DataType::String)
        {
            const auto &attribute = input.InquireAttribute<std::string>(attributePair.first);
            output[pos] = 66;
            ++pos;
            output[pos] = static_cast<uint8_t>(attribute->m_Type);
            ++pos;
            output[pos] = static_cast<uint8_t>(attribute->m_Name.size());
            ++pos;
            std::memcpy(output.data(pos), attribute->m_Name.data(), attribute->m_Name.size());
            pos += attribute->m_Name.size();
            output.value<uint64_t>(pos) = attribute->m_DataSingleValue.size();
            pos += 8;
            std::memcpy(output.data(pos), attribute->m_DataSingleValue.data(),
                        attribute->m_DataSingleValue.size());
            pos += attribute->m_DataSingleValue.size();
        }
#define declare_type(T)                                                                            \
    else if (attributePair.second->m_Type == helper::GetDataType<T>())                             \
    {                                                                                              \
        const auto &attribute = input.InquireAttribute<T>(attributePair.first);                    \
        output[pos] = 66;                                                                          \
        ++pos;                                                                                     \
        output[pos] = static_cast<uint8_t>(attribute->m_Type);                                     \
        ++pos;                                                                                     \
        output[pos] = static_cast<uint8_t>(attribute->m_Name.size());                              \
        ++pos;                                                                                     \
        std::memcpy(output.data(pos), attribute->m_Name.data(), attribute->m_Name.size());         \
        pos += attribute->m_Name.size();                                                           \
        if (attribute->m_IsSingleValue)                                                            \
        {                                                                                          \
            output.value<uint64_t>(pos) = sizeof(T);                                               \
            pos += 8;                                                                              \
            output.value<T>(pos) = attribute->m_DataSingleValue;                                   \
            pos += sizeof(T);                                                                      \
        }                                                                                          \
        else                                                                                       \
        {                                                                                          \
            uint64_t size = sizeof(T) * attribute->m_DataArray.size();                             \
            output.value<uint64_t>(pos) = size;                                                    \
            pos += 8;                                                                              \
            std::memcpy(output.data(pos), attribute->m_DataArray.data(), size);                    \
            pos += size;                                                                           \
        }                                                                                          \
    }
        ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_type)
#undef declare_type
        output.value<uint64_t>() = pos;
    }
}

void DeserializeAttribute(const Buffer &input, uint64_t &pos, IO &io, const bool regIO)
{
    const DataType type = static_cast<DataType>(input[pos]);
    ++pos;

    uint8_t nameSize = input[pos];
    ++pos;

    std::vector<char> namev(nameSize);
    std::memcpy(namev.data(), input.data(pos), nameSize);
    std::string name = std::string(namev.begin(), namev.end());
    pos += nameSize;

    uint64_t size = input.value<uint64_t>(pos);
    pos += 8;

    if (regIO)
    {
        const auto &attributes = io.GetAttributes();
        auto it = attributes.find(name);
        if (it == attributes.end())
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (type == DataType::String)
            {
                io.DefineAttribute<std::string>(name, std::string(input.data<char>(pos), size));
            }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        if (size == sizeof(T))                                                                     \
        {                                                                                          \
            io.DefineAttribute<T>(name, input.value<T>(pos));                                      \
        }                                                                                          \
        else                                                                                       \
        {                                                                                          \
            io.DefineAttribute<T>(name, input.data<T>(pos), size / sizeof(T));                     \
        }                                                                                          \
    }
            ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_type)
#undef declare_type
            else
            {
                helper::Throw<std::runtime_error>("Engine", "SscHelper", "DeserializeAttribute",
                                                  "unknown attribute data type");
            }
        }
    }
    pos += size;
}

void SerializeStructDefinitions(
    const std::unordered_multimap<std::string, StructDefinition> &definitions, Buffer &output)
{
    if (definitions.empty())
    {
        return;
    }
    uint64_t pos = output.value<uint64_t>();
    output[pos] = 65;
    ++pos;
    output[pos] = static_cast<uint8_t>(definitions.size());
    ++pos;
    for (const auto &p : definitions)
    {
        output.resize(pos + 1024);
        output[pos] = static_cast<uint8_t>(p.first.size());
        ++pos;
        std::memcpy(output.data(pos), p.first.data(), p.first.size());
        pos += p.first.size();
        output.value<uint64_t>(pos) = p.second.StructSize();
        pos += 8;
        output[pos] = static_cast<uint8_t>(p.second.Fields());
        ++pos;
        for (size_t i = 0; i < p.second.Fields(); ++i)
        {
            output[pos] = static_cast<uint8_t>(p.second.Name(i).size());
            ++pos;
            std::memcpy(output.data(pos), p.second.Name(i).data(), p.second.Name(i).size());
            pos += p.second.Name(i).size();
            output.value<uint64_t>(pos) = p.second.Offset(i);
            pos += 8;
            output.value<DataType>(pos) = p.second.Type(i);
            pos += sizeof(DataType);
            output.value<uint64_t>(pos) = p.second.ElementCount(i);
            pos += 8;
        }
    }
    output.value<uint64_t>() = pos;
}

void DeserializeStructDefinitions(const Buffer &input, uint64_t &pos, IO &io, const bool regIO,
                                  std::unordered_map<std::string, StructDefinition> &StructDefs)
{
    uint8_t defs = input.value<uint8_t>(pos);
    ++pos;
    for (uint8_t i = 0; i < defs; ++i)
    {
        uint8_t nameSize = input.value<uint8_t>(pos);
        ++pos;
        std::vector<char> defNameChar(nameSize);
        std::memcpy(defNameChar.data(), input.data(pos), nameSize);
        std::string defName = std::string(defNameChar.begin(), defNameChar.end());
        pos += nameSize;
        size_t structSize = input.value<uint64_t>(pos);
        pos += 8;
        uint8_t items = input.value<uint8_t>(pos);
        ++pos;
        StructDefinition *structDefinition = nullptr;
        if (StructDefs.find(defName) == StructDefs.end())
        {
            structDefinition =
                &StructDefs.emplace(defName, StructDefinition(defName, structSize)).first->second;
        }
        for (uint8_t j = 0; j < items; ++j)
        {
            nameSize = input.value<uint8_t>(pos);
            ++pos;
            defNameChar.resize(nameSize);
            std::memcpy(defNameChar.data(), input.data(pos), nameSize);
            std::string itemName = std::string(defNameChar.begin(), defNameChar.end());
            pos += nameSize;
            size_t itemOffset = input.value<uint64_t>(pos);
            pos += 8;
            DataType itemType = input.value<DataType>(pos);
            pos += sizeof(DataType);
            size_t itemSize = input.value<uint64_t>(pos);
            pos += 8;
            if (structDefinition != nullptr)
            {
                structDefinition->AddField(itemName, itemOffset, itemType, itemSize);
            }
        }
    }
}

void Deserialize(const Buffer &input, BlockVecVec &output, IO &io, const bool regVars,
                 const bool regAttrs, const bool regDefs,
                 std::unordered_map<std::string, StructDefinition> &StructDefs)
{
    for (auto &i : output)
    {
        i.clear();
    }

    uint64_t pos = 2;

    uint64_t blockSize = input.value<uint64_t>(pos);

    pos += 8;

    while (pos < blockSize)
    {

        uint8_t shapeId = input[pos];
        ++pos;

        if (shapeId == 65)
        {
            DeserializeStructDefinitions(input, pos, io, regDefs, StructDefs);
        }
        else if (shapeId == 66)
        {
            DeserializeAttribute(input, pos, io, regAttrs);
        }
        else
        {
            int rank = input.value<int>(pos);
            pos += 4;
            output[rank].emplace_back();
            auto &b = output[rank].back();

            DeserializeVariable(input, static_cast<ShapeID>(shapeId), pos, b, io, regVars,
                                StructDefs);
        }
    }
}

void AggregateMetadata(const Buffer &localBuffer, Buffer &globalBuffer, MPI_Comm comm,
                       const bool finalStep, const bool locked)
{
    int mpiSize;
    MPI_Comm_size(comm, &mpiSize);
    int localSize = static_cast<int>(localBuffer.value<uint64_t>()) - 8;
    std::vector<int> localSizes(mpiSize);
    MPI_Gather(&localSize, 1, MPI_INT, localSizes.data(), 1, MPI_INT, 0, comm);
    int globalSize = std::accumulate(localSizes.begin(), localSizes.end(), 0);
    globalBuffer.resize(globalSize + 10);

    std::vector<int> displs(mpiSize);
    for (size_t i = 1; i < static_cast<size_t>(mpiSize); ++i)
    {
        displs[i] = displs[i - 1] + localSizes[i - 1];
    }

    MPI_Gatherv(localBuffer.data() + 8, localSize, MPI_CHAR, globalBuffer.data() + 10,
                localSizes.data(), displs.data(), MPI_CHAR, 0, comm);
    globalBuffer[0] = finalStep;
    globalBuffer[1] = locked;
    globalBuffer.value<uint64_t>(2) = globalSize;
}

void BroadcastMetadata(Buffer &globalBuffer, const int root, MPI_Comm comm)
{
    int globalBufferSize = static_cast<int>(globalBuffer.size());
    MPI_Bcast(&globalBufferSize, 1, MPI_INT, root, comm);
    if (globalBuffer.size() < static_cast<size_t>(globalBufferSize))
    {
        globalBuffer.resize(globalBufferSize);
    }
    MPI_Bcast(globalBuffer.data(), globalBufferSize, MPI_CHAR, root, comm);
}

bool AreSameDims(const Dims &a, const Dims &b)
{
    if (a.size() != b.size())
    {
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i)
    {
        if (a[i] != b[i])
        {
            return false;
        }
    }
    return true;
}

void MPI_Gatherv64(const void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const uint64_t *recvcounts, const uint64_t *displs, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, const int chunksize)
{

    int mpiSize;
    int mpiRank;
    MPI_Comm_size(comm, &mpiSize);
    MPI_Comm_rank(comm, &mpiRank);

    int recvTypeSize;
    int sendTypeSize;

    MPI_Type_size(recvtype, &recvTypeSize);
    MPI_Type_size(sendtype, &sendTypeSize);

    std::vector<MPI_Request> requests;
    if (mpiRank == root)
    {
        for (int i = 0; i < mpiSize; ++i)
        {
            uint64_t recvcount = recvcounts[i];
            while (recvcount > 0)
            {
                requests.emplace_back();
                if (recvcount > static_cast<uint64_t>(chunksize))
                {
                    MPI_Irecv(reinterpret_cast<char *>(recvbuf) +
                                  (displs[i] + recvcounts[i] - recvcount) * recvTypeSize,
                              chunksize, recvtype, i, 0, comm, &requests.back());
                    recvcount -= static_cast<uint64_t>(chunksize);
                }
                else
                {
                    MPI_Irecv(reinterpret_cast<char *>(recvbuf) +
                                  (displs[i] + recvcounts[i] - recvcount) * recvTypeSize,
                              static_cast<int>(recvcount), recvtype, i, 0, comm, &requests.back());
                    recvcount = 0;
                }
            }
        }
    }

    uint64_t sendcountvar = sendcount;

    while (sendcountvar > 0)
    {
        requests.emplace_back();
        if (sendcountvar > static_cast<uint64_t>(chunksize))
        {
            MPI_Isend(reinterpret_cast<const char *>(sendbuf) +
                          (sendcount - sendcountvar) * sendTypeSize,
                      chunksize, sendtype, root, 0, comm, &requests.back());
            sendcountvar -= static_cast<uint64_t>(chunksize);
        }
        else
        {
            MPI_Isend(reinterpret_cast<const char *>(sendbuf) +
                          (sendcount - sendcountvar) * sendTypeSize,
                      static_cast<int>(sendcountvar), sendtype, root, 0, comm, &requests.back());
            sendcountvar = 0;
        }
    }

    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
}

void MPI_Gatherv64OneSidedPull(const void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const uint64_t *recvcounts, const uint64_t *displs,
                               MPI_Datatype recvtype, int root, MPI_Comm comm, const int chunksize)
{

    int mpiSize;
    int mpiRank;
    MPI_Comm_size(comm, &mpiSize);
    MPI_Comm_rank(comm, &mpiRank);

    int recvTypeSize;
    int sendTypeSize;

    MPI_Type_size(recvtype, &recvTypeSize);
    MPI_Type_size(sendtype, &sendTypeSize);

    MPI_Win win;
    MPI_Win_create(const_cast<void *>(sendbuf), sendcount * sendTypeSize, sendTypeSize,
                   MPI_INFO_NULL, comm, &win);

    if (mpiRank == root)
    {
        for (int i = 0; i < mpiSize; ++i)
        {
            uint64_t recvcount = recvcounts[i];
            while (recvcount > 0)
            {
                if (recvcount > static_cast<uint64_t>(chunksize))
                {
                    MPI_Get(reinterpret_cast<char *>(recvbuf) +
                                (displs[i] + recvcounts[i] - recvcount) * recvTypeSize,
                            chunksize, recvtype, i, recvcounts[i] - recvcount, chunksize, recvtype,
                            win);
                    recvcount -= static_cast<uint64_t>(chunksize);
                }
                else
                {
                    MPI_Get(reinterpret_cast<char *>(recvbuf) +
                                (displs[i] + recvcounts[i] - recvcount) * recvTypeSize,
                            static_cast<int>(recvcount), recvtype, i, recvcounts[i] - recvcount,
                            static_cast<int>(recvcount), recvtype, win);
                    recvcount = 0;
                }
            }
        }
    }

    MPI_Win_free(&win);
}

void MPI_Gatherv64OneSidedPush(const void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const uint64_t *recvcounts, const uint64_t *displs,
                               MPI_Datatype recvtype, int root, MPI_Comm comm, const int chunksize)
{

    int mpiSize;
    int mpiRank;
    MPI_Comm_size(comm, &mpiSize);
    MPI_Comm_rank(comm, &mpiRank);

    int recvTypeSize;
    int sendTypeSize;

    MPI_Type_size(recvtype, &recvTypeSize);
    MPI_Type_size(sendtype, &sendTypeSize);

    uint64_t recvsize = displs[mpiSize - 1] + recvcounts[mpiSize - 1];

    MPI_Win win;
    MPI_Win_create(recvbuf, recvsize * recvTypeSize, recvTypeSize, MPI_INFO_NULL, comm, &win);

    uint64_t sendcountvar = sendcount;

    while (sendcountvar > 0)
    {
        if (sendcountvar > static_cast<uint64_t>(chunksize))
        {
            MPI_Put(reinterpret_cast<const char *>(sendbuf) +
                        (sendcount - sendcountvar) * sendTypeSize,
                    chunksize, sendtype, root, displs[mpiRank] + sendcount - sendcountvar,
                    chunksize, sendtype, win);
            sendcountvar -= static_cast<uint64_t>(chunksize);
        }
        else
        {
            MPI_Put(reinterpret_cast<const char *>(sendbuf) +
                        (sendcount - sendcountvar) * sendTypeSize,
                    static_cast<int>(sendcountvar), sendtype, root,
                    displs[mpiRank] + sendcount - sendcountvar, static_cast<int>(sendcountvar),
                    sendtype, win);
            sendcountvar = 0;
        }
    }

    MPI_Win_free(&win);
}

void PrintDims(const Dims &dims, const std::string &label)
{
    std::cout << label;
    for (const auto &i : dims)
    {
        std::cout << i << ", ";
    }
    std::cout << std::endl;
}

void PrintBlock(const BlockInfo &b, const std::string &label)
{
    std::cout << label << std::endl;
    std::cout << b.name << std::endl;
    std::cout << "    DataType : " << b.type << std::endl;
    std::cout << "    ElementSize : " << b.elementSize << std::endl;
    PrintDims(b.shape, "    Shape : ");
    PrintDims(b.start, "    Start : ");
    PrintDims(b.count, "    Count : ");
    std::cout << "    Position Start : " << b.bufferStart << std::endl;
    std::cout << "    Position Count : " << b.bufferCount << std::endl;
}

void PrintBlockVec(const BlockVec &bv, const std::string &label)
{
    std::cout << label << std::endl;
    for (const auto &i : bv)
    {
        std::cout << i.name << std::endl;
        std::cout << "    DataType : " << i.type << std::endl;
        std::cout << "    ElementSize : " << i.elementSize << std::endl;
        PrintDims(i.shape, "    Shape : ");
        PrintDims(i.start, "    Start : ");
        PrintDims(i.count, "    Count : ");
        std::cout << "    Position Start : " << i.bufferStart << std::endl;
        std::cout << "    Position Count : " << i.bufferCount << std::endl;
    }
}

void PrintBlockVecVec(const BlockVecVec &bvv, const std::string &label)
{
    std::cout << label << std::endl;
    size_t rank = 0;
    for (const auto &bv : bvv)
    {
        std::cout << "Rank " << rank << std::endl;
        for (const auto &i : bv)
        {
            std::cout << "    " << i.name << std::endl;
            std::cout << "        DataType : " << i.type << std::endl;
            std::cout << "        ElementSize : " << i.elementSize << std::endl;
            PrintDims(i.shape, "        Shape : ");
            PrintDims(i.start, "        Start : ");
            PrintDims(i.count, "        Count : ");
            std::cout << "        Position Start : " << i.bufferStart << std::endl;
            std::cout << "        Position Count : " << i.bufferCount << std::endl;
        }
        ++rank;
    }
}

void PrintRankPosMap(const RankPosMap &m, const std::string &label)
{
    std::cout << label << std::endl;
    for (const auto &rank : m)
    {
        std::cout << "Rank = " << rank.first << ", bufferStart = " << rank.second.first
                  << ", bufferCount = " << rank.second.second << std::endl;
    }
}

void PrintMpiInfo(const MpiInfo &writersInfo, const MpiInfo &readersInfo)
{
    int s = 0;
    for (size_t i = 0; i < writersInfo.size(); ++i)
    {
        std::cout << "App " << s << " Writer App " << i << " Wrold Ranks : ";
        for (size_t j = 0; j < writersInfo[i].size(); ++j)
        {
            std::cout << writersInfo[i][j] << "  ";
        }
        std::cout << std::endl;
        ++s;
    }
    for (size_t i = 0; i < readersInfo.size(); ++i)
    {
        std::cout << "App " << s << " Reader App " << i << " Wrold Ranks : ";
        for (size_t j = 0; j < readersInfo[i].size(); ++j)
        {
            std::cout << readersInfo[i][j] << "  ";
        }
        std::cout << std::endl;
        ++s;
    }
    std::cout << std::endl;
}

} // end namespace ssc
} // end namespace engine
} // end namespace core
} // end namespace adios2
