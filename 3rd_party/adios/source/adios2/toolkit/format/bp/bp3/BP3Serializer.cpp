/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP3Serializer.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BP3Serializer.h"
#include "BP3Serializer.tcc"

#include <chrono>
#include <future>
#include <string>
#include <vector>

#include "adios2/helper/adiosFunctions.h" //helper::GetDataType<T>, helper::ReadValue<T>
#include <adios2-perfstubs-interface.h>

#ifdef _WIN32
#pragma warning(disable : 4503) // Windows complains about SubFileInfoMap levels
#endif

namespace adios2
{
namespace format
{

BP3Serializer::BP3Serializer(helper::Comm const &comm)
: BPBase(comm), BP3Base(comm), BPSerializer(comm, 3)
{
}

void BP3Serializer::PutProcessGroupIndex(const std::string &ioName, const std::string hostLanguage,
                                         const std::vector<std::string> &transportsTypes) noexcept
{
    m_Profiler.Start("buffering");
    std::vector<char> &metadataBuffer = m_MetadataSet.PGIndex.Buffer;

    std::vector<char> &dataBuffer = m_Data.m_Buffer;
    size_t &dataPosition = m_Data.m_Position;

    m_MetadataSet.DataPGLengthPosition = dataPosition;
    dataPosition += 8; // skip pg length (8)

    const std::size_t metadataPGLengthPosition = metadataBuffer.size();
    metadataBuffer.insert(metadataBuffer.end(), 2, '\0'); // skip pg length (2)

    // write name to metadata
    PutNameRecord(ioName, metadataBuffer);

    // write if data is column major in metadata and data
    const char columnMajor = (helper::IsRowMajor(hostLanguage) == false) ? 'y' : 'n';
    helper::InsertToBuffer(metadataBuffer, &columnMajor);
    helper::CopyToBuffer(dataBuffer, dataPosition, &columnMajor);

    // write name in data
    PutNameRecord(ioName, dataBuffer, dataPosition);

    // processID in metadata,
    const uint32_t processID = static_cast<uint32_t>(m_RankMPI);
    helper::InsertToBuffer(metadataBuffer, &processID);
    // skip coordination var in data ....what is coordination var?
    dataPosition += 4;

    // time step name to metadata and data
    const std::string timeStepName(std::to_string(m_MetadataSet.TimeStep));
    PutNameRecord(timeStepName, metadataBuffer);
    PutNameRecord(timeStepName, dataBuffer, dataPosition);

    // time step to metadata and data
    helper::InsertToBuffer(metadataBuffer, &m_MetadataSet.TimeStep);
    helper::CopyToBuffer(dataBuffer, dataPosition, &m_MetadataSet.TimeStep);

    // offset to pg in data in metadata which is the current absolute position
    helper::InsertU64(metadataBuffer, m_Data.m_AbsolutePosition);

    // Back to writing metadata pg index length (length of group)
    const uint16_t metadataPGIndexLength =
        static_cast<uint16_t>(metadataBuffer.size() - metadataPGLengthPosition - 2);

    size_t backPosition = metadataPGLengthPosition;
    helper::CopyToBuffer(metadataBuffer, backPosition, &metadataPGIndexLength);
    // DONE With metadataBuffer

    // here write method in data
    const std::vector<uint8_t> methodIDs = GetTransportIDs(transportsTypes);
    const uint8_t methodsCount = static_cast<uint8_t>(methodIDs.size());
    helper::CopyToBuffer(dataBuffer, dataPosition, &methodsCount); // count
    // methodID (1) + method params length(2), no parameters for now
    const uint16_t methodsLength = static_cast<uint16_t>(methodsCount * 3);

    helper::CopyToBuffer(dataBuffer, dataPosition, &methodsLength); // length

    for (const auto methodID : methodIDs)
    {
        helper::CopyToBuffer(dataBuffer, dataPosition, &methodID); // method ID,
        dataPosition += 2; // skip method params length = 0 (2 bytes) for now
    }

    // update absolute position
    m_Data.m_AbsolutePosition += dataPosition - m_MetadataSet.DataPGLengthPosition;
    // pg vars count and position
    m_MetadataSet.DataPGVarsCount = 0;
    m_MetadataSet.DataPGVarsCountPosition = dataPosition;
    // add vars count and length
    dataPosition += 12;
    m_Data.m_AbsolutePosition += 12; // add vars count and length

    ++m_MetadataSet.DataPGCount;
    m_MetadataSet.DataPGIsOpen = true;

    m_Profiler.Stop("buffering");
}

void BP3Serializer::CloseData(core::IO &io)
{
    m_Profiler.Start("buffering");

    if (!m_IsClosed)
    {
        if (m_MetadataSet.DataPGIsOpen)
        {
            SerializeDataBuffer(io);
        }

        SerializeMetadataInData();

        if (m_Profiler.m_IsActive)
        {
            m_Profiler.m_Bytes.at("buffering") = m_Data.m_AbsolutePosition;
        }

        m_Aggregator.Close();
        m_IsClosed = true;
    }

    m_Profiler.Stop("buffering");
}

void BP3Serializer::CloseStream(core::IO &io, const bool addMetadata)
{
    m_Profiler.Start("buffering");
    if (m_MetadataSet.DataPGIsOpen)
    {
        SerializeDataBuffer(io);
    }

    SerializeMetadataInData(false, addMetadata);

    if (m_Profiler.m_IsActive)
    {
        m_Profiler.m_Bytes.at("buffering") += m_Data.m_Position;
    }
    m_Profiler.Stop("buffering");
}

void BP3Serializer::CloseStream(core::IO &io, size_t &metadataStart, size_t &metadataCount,
                                const bool addMetadata)
{

    m_Profiler.Start("buffering");
    if (m_MetadataSet.DataPGIsOpen)
    {
        SerializeDataBuffer(io);
    }

    metadataStart = m_Data.m_Position;

    SerializeMetadataInData(false, addMetadata);

    metadataCount = m_Data.m_Position - metadataStart;

    if (m_Profiler.m_IsActive)
    {
        m_Profiler.m_Bytes.at("buffering") += m_Data.m_Position;
    }
    m_Profiler.Stop("buffering");
}

void BP3Serializer::ResetIndices()
{
    m_MetadataSet.PGIndex.Buffer.clear();
    m_MetadataSet.AttributesIndices.clear();
    m_MetadataSet.VarsIndices.clear();
}

void BP3Serializer::AggregateCollectiveMetadata(helper::Comm const &comm, BufferSTL &bufferSTL,
                                                const bool inMetadataBuffer)
{
    m_Profiler.Start("buffering");
    m_Profiler.Start("meta_sort_merge");

    const std::vector<size_t> indicesPosition = AggregateCollectiveMetadataIndices(comm, bufferSTL);

    int rank = comm.Rank();
    if (rank == 0)
    {
        PutMinifooter(static_cast<uint64_t>(indicesPosition[0]),
                      static_cast<uint64_t>(indicesPosition[1]),
                      static_cast<uint64_t>(indicesPosition[2]), bufferSTL.m_Buffer,
                      bufferSTL.m_Position, inMetadataBuffer);

        if (inMetadataBuffer)
        {
            bufferSTL.m_AbsolutePosition = bufferSTL.m_Position;
        }
        else
        {
            bufferSTL.m_AbsolutePosition += bufferSTL.m_Position;
        }
    }

    bufferSTL.Resize(bufferSTL.m_Position, "after collective metadata is done");

    m_Profiler.Stop("meta_sort_merge");
    m_Profiler.Stop("buffering");
}

// PRIVATE FUNCTIONS
void BP3Serializer::SerializeDataBuffer(core::IO &io) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    // vars count and Length (only for PG)
    helper::CopyToBuffer(buffer, m_MetadataSet.DataPGVarsCountPosition,
                         &m_MetadataSet.DataPGVarsCount);
    // without record itself and vars count
    const uint64_t varsLength = position - m_MetadataSet.DataPGVarsCountPosition - 8 - 4;
    helper::CopyToBuffer(buffer, m_MetadataSet.DataPGVarsCountPosition, &varsLength);

    // each attribute is only written to output once
    size_t attributesSizeInData = GetAttributesSizeInData(io);
    if (attributesSizeInData)
    {
        attributesSizeInData += 12; // count + length
        m_Data.Resize(position + attributesSizeInData, "when writing Attributes in rank=0\n");

        PutAttributes(io);
    }
    else
    {
        m_Data.Resize(position + 12, "for empty Attributes\n");
        // Attribute index header for zero attributes: 0, 0LL
        // Resize() already takes care of this
        position += 12;
        absolutePosition += 12;
    }

    // Finish writing pg group length without record itself
    const uint64_t dataPGLength = position - m_MetadataSet.DataPGLengthPosition - 8;
    helper::CopyToBuffer(buffer, m_MetadataSet.DataPGLengthPosition, &dataPGLength);

    m_MetadataSet.DataPGIsOpen = false;
}

void BP3Serializer::SerializeMetadataInData(const bool updateAbsolutePosition, const bool inData)
{
    auto lf_SetIndexCountLength = [](std::unordered_map<std::string, SerialElementIndex> &indices,
                                     uint32_t &count, uint64_t &length) {
        count = static_cast<uint32_t>(indices.size());
        length = 0;
        for (auto &indexPair : indices) // set each index length
        {
            auto &indexBuffer = indexPair.second.Buffer;
            const uint32_t indexLength = static_cast<uint32_t>(indexBuffer.size() - 4);
            size_t indexLengthPosition = 0;
            helper::CopyToBuffer(indexBuffer, indexLengthPosition, &indexLength);

            length += indexBuffer.size(); // overall length
        }
    };

    auto lf_FlattenIndices = [](const uint32_t count, const uint64_t length,
                                const std::unordered_map<std::string, SerialElementIndex> &indices,
                                std::vector<char> &buffer, size_t &position) {
        helper::CopyToBuffer(buffer, position, &count);
        helper::CopyToBuffer(buffer, position, &length);

        for (const auto &indexPair : indices) // set each index length
        {
            const auto &indexBuffer = indexPair.second.Buffer;
            helper::CopyToBuffer(buffer, position, indexBuffer.data(), indexBuffer.size());
        }
    };

    // Finish writing metadata counts and lengths
    // PG Index
    const uint64_t pgCount = m_MetadataSet.DataPGCount;
    const uint64_t pgLength = m_MetadataSet.PGIndex.Buffer.size();

    // var index count and length (total), and each index length
    uint32_t varsCount = 0;
    uint64_t varsLength = 0;
    lf_SetIndexCountLength(m_MetadataSet.VarsIndices, varsCount, varsLength);

    // attribute index count and length, and each index length
    uint32_t attributesCount = 0;
    uint64_t attributesLength = 0;
    lf_SetIndexCountLength(m_MetadataSet.AttributesIndices, attributesCount, attributesLength);

    if (!inData)
    {
        return;
    }

    const size_t footerSize =
        static_cast<size_t>((pgLength + 16) + (varsLength + 12) + (attributesLength + 12) +
                            m_MetadataSet.MiniFooterSize);

    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    // reserve data to fit metadata,
    // must replace with growth buffer strategy?
    m_Data.Resize(position + footerSize, " when writing metadata in bp data buffer");

    // write pg index
    helper::CopyToBuffer(buffer, position, &pgCount);
    helper::CopyToBuffer(buffer, position, &pgLength);
    helper::CopyToBuffer(buffer, position, m_MetadataSet.PGIndex.Buffer.data(),
                         static_cast<size_t>(pgLength));

    // Vars indices
    lf_FlattenIndices(varsCount, varsLength, m_MetadataSet.VarsIndices, buffer, position);
    // Attribute indices
    lf_FlattenIndices(attributesCount, attributesLength, m_MetadataSet.AttributesIndices, buffer,
                      position);

    // getting absolute offset start, minifooter is 28 bytes for now
    const uint64_t pgIndexStart = static_cast<uint64_t>(absolutePosition);
    const uint64_t variablesIndexStart = static_cast<uint64_t>(pgIndexStart + (pgLength + 16));
    const uint64_t attributesIndexStart =
        static_cast<uint64_t>(variablesIndexStart + (varsLength + 12));

    PutMinifooter(pgIndexStart, variablesIndexStart, attributesIndexStart, buffer, position);

    if (updateAbsolutePosition)
    {
        absolutePosition += footerSize;
    }

    if (m_Profiler.m_IsActive)
    {
        m_Profiler.m_Bytes.emplace("buffering", absolutePosition);
    }
}

std::vector<size_t> BP3Serializer::AggregateCollectiveMetadataIndices(helper::Comm const &comm,
                                                                      BufferSTL &bufferSTL)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    int rank = comm.Rank();
    int size = comm.Size();

    // pre-allocate with rank 0 data
    size_t pgCount = 0; //< tracks global PG count
    if (rank == 0)
    {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        // assumes that things are more or less balanced
        m_PGRankIndices.reserve(m_MetadataSet.PGIndex.Buffer.size() * static_cast<size_t>(size));
        m_PGRankIndices.resize(0);

        m_VariableRankIndices.clear();
        m_VariableRankIndices.reserve(m_MetadataSet.VarsIndices.size());

        m_AttributesRankIndices.clear();
        m_AttributesRankIndices.reserve(m_MetadataSet.AttributesIndices.size());
    }

    auto lf_IndicesSize =
        [&](const std::unordered_map<std::string, SerialElementIndex> &indices) -> size_t

    {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        size_t indicesSize = 0;
        for (const auto &indexPair : indices)
        {
            indicesSize += indexPair.second.Buffer.size();
        }
        return indicesSize;
    };

    auto lf_SerializeIndices =
        [&](const std::unordered_map<std::string, SerialElementIndex> &indices, size_t &position)

    {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        for (const auto &indexPair : indices)
        {
            const auto &buffer = indexPair.second.Buffer;
            helper::CopyToBuffer(m_SerializedIndices, position, buffer.data(), buffer.size());
        }
    };

    auto lf_SerializeAllIndices = [&](helper::Comm const &comm, const int rank) {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        const size_t pgIndicesSize = m_MetadataSet.PGIndex.Buffer.size();
        const size_t variablesIndicesSize = lf_IndicesSize(m_MetadataSet.VarsIndices);
        const size_t attributesIndicesSize = lf_IndicesSize(m_MetadataSet.AttributesIndices);

        // first pre-allocate
        const size_t serializedIndicesSize =
            8 * 4 + pgIndicesSize + variablesIndicesSize + attributesIndicesSize;

        m_SerializedIndices.reserve(serializedIndicesSize + 4);
        m_SerializedIndices.resize(serializedIndicesSize + 4);

        const uint32_t rank32 = static_cast<uint32_t>(rank);
        const uint64_t size64 = static_cast<uint64_t>(serializedIndicesSize);
        const uint64_t variablesIndexOffset = static_cast<uint64_t>(pgIndicesSize + 36);
        const uint64_t attributesIndexOffset =
            static_cast<uint64_t>(pgIndicesSize + 36 + variablesIndicesSize);

        size_t position = 0;
        helper::CopyToBuffer(m_SerializedIndices, position, &rank32);
        helper::CopyToBuffer(m_SerializedIndices, position, &size64);
        helper::CopyToBuffer(m_SerializedIndices, position, &variablesIndexOffset);
        helper::CopyToBuffer(m_SerializedIndices, position, &attributesIndexOffset);
        helper::CopyToBuffer(m_SerializedIndices, position, &m_MetadataSet.DataPGCount);

        helper::CopyToBuffer(m_SerializedIndices, position, m_MetadataSet.PGIndex.Buffer.data(),
                             m_MetadataSet.PGIndex.Buffer.size());
        lf_SerializeIndices(m_MetadataSet.VarsIndices, position);
        lf_SerializeIndices(m_MetadataSet.AttributesIndices, position);
    };

    auto lf_DeserializeIndices =
        [&](std::unordered_map<std::string, std::vector<SerialElementIndex>> &deserialized,
            const int rankSource, const std::vector<char> &serialized, const size_t position,
            const size_t endPosition, const bool isRankConstant)

    {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        size_t localPosition = position;
        while (localPosition < endPosition)
        {
            size_t indexPosition = localPosition;
            const ElementIndexHeader header =
                ReadElementIndexHeader(serialized, indexPosition, helper::IsLittleEndian());
            const size_t bufferSize = static_cast<size_t>(header.Length) + 4;

            if (isRankConstant && deserialized.count(header.Name) == 1)
            {
                localPosition += bufferSize;
                continue;
            }
            std::vector<BP3Base::SerialElementIndex> *deserializedIndexes = nullptr;
            auto search = deserialized.find(header.Name);
            if (search == deserialized.end())
            {
                // key (variable name) is not in the map, add it
                // mutex portion
                {
                    std::lock_guard<std::mutex> lock(m_Mutex);
                    deserializedIndexes =
                        &(deserialized
                              .emplace(std::piecewise_construct, std::forward_as_tuple(header.Name),
                                       std::forward_as_tuple(
                                           size, SerialElementIndex(header.MemberID, bufferSize)))
                              .first->second);
                }
            }
            else
            {
                // variable name is already in the map, just get the pointer
                // of the vector of metadata indices
                deserializedIndexes = &(search->second);
            }

            SerialElementIndex &index = deserializedIndexes->at(static_cast<size_t>(rankSource));
            helper::InsertToBuffer(index.Buffer, &serialized[localPosition], bufferSize);
            localPosition += bufferSize;
        }
    };

    auto lf_DeserializeAllIndices = [&](const int rankSource, const std::vector<size_t> headerInfo,
                                        const std::vector<char> &serialized, const size_t position)

    {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        const size_t rankIndicesSize = headerInfo[0];
        const size_t variablesIndexOffset = headerInfo[1] + position;
        const size_t attributesIndexOffset = headerInfo[2] + position;
        pgCount += headerInfo[3];
        size_t localPosition = position + 36;

        const size_t pgIndexLength = variablesIndexOffset - localPosition;
        // first deserialize pg indices
        {
            std::lock_guard<std::mutex> lock(m_Mutex);
            helper::InsertToBuffer(m_PGRankIndices, &serialized[localPosition], pgIndexLength);
        }
        // deserialize variable indices
        localPosition = variablesIndexOffset;
        size_t endPosition = attributesIndexOffset;

        lf_DeserializeIndices(m_VariableRankIndices, rankSource, serialized, localPosition,
                              endPosition, false);

        // deserialize attributes indices
        localPosition = attributesIndexOffset;
        endPosition = rankIndicesSize + 4 + position;
        // attributes are constant and unique across ranks
        lf_DeserializeIndices(m_AttributesRankIndices, rankSource, serialized, localPosition,
                              endPosition, true);
    };

    auto lf_SortMergeIndices =
        [&](const std::unordered_map<std::string, std::vector<SerialElementIndex>>
                &deserializedIndices) {
            PERFSTUBS_SCOPED_TIMER_FUNC();
            auto &position = bufferSTL.m_Position;
            auto &buffer = bufferSTL.m_Buffer;

            size_t countPosition = position;

            const uint32_t totalCountU32 = static_cast<uint32_t>(deserializedIndices.size());
            helper::CopyToBuffer(buffer, countPosition, &totalCountU32);
            position += 12; // skip for length

            MergeSerializeIndices(deserializedIndices, comm, bufferSTL);

            // Write length
            const uint64_t totalLengthU64 = static_cast<uint64_t>(position - countPosition - 8);
            helper::CopyToBuffer(buffer, countPosition, &totalLengthU64);
        };

    // BODY of function starts here
    std::vector<size_t> indexPositions(3);
    // pgIndex
    indexPositions[0] = bufferSTL.m_AbsolutePosition;

    lf_SerializeAllIndices(comm, rank); // Set m_SerializedIndices

    size_t countPosition = bufferSTL.m_Position;

    comm.GathervVectors(m_SerializedIndices, bufferSTL.m_Buffer, bufferSTL.m_Position, 0);

    // use bufferSTL (will resize) to GatherV
    const size_t extraSize = 16 + 12 + 12 + m_MetadataSet.MiniFooterSize;

    try
    {
        bufferSTL.m_Buffer.reserve(bufferSTL.m_Position + extraSize); // to avoid power of 2 growth
        bufferSTL.m_Buffer.resize(bufferSTL.m_Position + extraSize);
    }
    catch (...)
    {
        helper::ThrowNested<std::runtime_error>(
            "Toolkit::Format", "bp::bp3::BP3Serializer", "AggregateCollectiveMetadataIndices",
            "buffer overflow when resizing to " + std::to_string(bufferSTL.m_Position + extraSize) +
                " bytes");
    }

    // deserialize, it's all local inside rank 0
    if (rank == 0)
    {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        const size_t serializedSize = bufferSTL.m_Position;
        const std::vector<char> &serialized = bufferSTL.m_Buffer;
        size_t serializedPosition = 0;
        std::vector<size_t> headerInfo(4);
        const bool isLittleEndian = helper::IsLittleEndian();

        // if (m_Parameters.Threads == 1)
        {
            while (serializedPosition < serializedSize)
            {
                size_t localPosition = serializedPosition;

                const int rankSource = static_cast<int>(
                    helper::ReadValue<uint32_t>(serialized, localPosition, isLittleEndian));

                for (auto i = 0; i < 4; ++i)
                {
                    headerInfo[i] = static_cast<size_t>(
                        helper::ReadValue<uint64_t>(serialized, localPosition, isLittleEndian));
                }

                lf_DeserializeAllIndices(rankSource, headerInfo, serialized, serializedPosition);
                serializedPosition += headerInfo[0] + 4;
            }
        }
    }

    // now merge (and sort variables and attributes) indices
    if (rank == 0)
    {
        PERFSTUBS_SCOPED_TIMER_FUNC();
        auto &position = bufferSTL.m_Position;
        auto &buffer = bufferSTL.m_Buffer;
        position = countPosition; // back to pg count position

        const uint64_t pgCount64 = static_cast<uint64_t>(pgCount);
        const uint64_t pgLength64 = static_cast<uint64_t>(m_PGRankIndices.size());

        helper::CopyToBuffer(buffer, position, &pgCount64);
        helper::CopyToBuffer(buffer, position, &pgLength64);
        helper::CopyToBuffer(buffer, position, m_PGRankIndices.data(), m_PGRankIndices.size());

        indexPositions[1] = bufferSTL.m_AbsolutePosition + position;
        lf_SortMergeIndices(m_VariableRankIndices);
        indexPositions[2] = bufferSTL.m_AbsolutePosition + position;
        lf_SortMergeIndices(m_AttributesRankIndices);
    }

    return indexPositions;
}

#define declare_template_instantiation(T)                                                          \
    void BP3Serializer::DoPutAttributeInData(const core::Attribute<T> &attribute,                  \
                                             Stats<T> &stats) noexcept                             \
    {                                                                                              \
        PutAttributeInDataCommon(attribute, stats);                                                \
    }
ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

#define declare_template_instantiation(T)                                                          \
    template void BP3Serializer::PutVariableMetadata(                                              \
        const core::Variable<T> &, const typename core::Variable<T>::BPInfo &, const bool,         \
        typename core::Variable<T>::Span *) noexcept;                                              \
                                                                                                   \
    template void BP3Serializer::PutVariablePayload(                                               \
        const core::Variable<T> &, const typename core::Variable<T>::BPInfo &, const bool,         \
        typename core::Variable<T>::Span *) noexcept;

ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

#define declare_template_instantiation(T)                                                          \
    template void BP3Serializer::PutSpanMetadata(                                                  \
        const core::Variable<T> &, const typename core::Variable<T>::Span &) noexcept;

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace format
} // end namespace adios2
