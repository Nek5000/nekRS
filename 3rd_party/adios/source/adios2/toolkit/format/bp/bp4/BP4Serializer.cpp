/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Serializer.cpp
 *
 *  Created on: Aug 1, 2018
 *  Author: William F Godoy godoywf@ornl.gov
 *          Lipeng Wan wanl@ornl.gov
 *          Norbert Podhorszki pnb@ornl.gov
 */

#include "BP4Serializer.h"
#include "BP4Serializer.tcc"

#include <chrono>
#include <future>
#include <string>
#include <tuple>
#include <vector>

#include "adios2/helper/adiosFunctions.h" //helper::GetDataType<T>, helper::ReadValue<T>

#ifdef _WIN32
#pragma warning(disable : 4503) // Windows complains about SubFileInfoMap levels
#endif

namespace adios2
{
namespace format
{

BP4Serializer::BP4Serializer(helper::Comm const &comm)
: BPBase(comm), BP4Base(comm), BPSerializer(comm, 4)
{
}

/*generate the header for the metadata index file*/
void BP4Serializer::MakeHeader(BufferSTL &b, const std::string fileType, const bool isActive)
{
    auto lf_CopyVersionChar = [](const std::string version, std::vector<char> &buffer,
                                 size_t &position) {
        helper::CopyToBuffer(buffer, position, version.c_str());
    };

    auto &buffer = b.m_Buffer;
    auto &position = b.m_Position;
    auto &absolutePosition = b.m_AbsolutePosition;
    if (position > 0)
    {
        helper::Throw<std::invalid_argument>("Toolkit", "format::bp::BP4Serializer", "MakeHeader",
                                             "can only be called for an empty "
                                             "buffer. This one for " +
                                                 fileType + " already has content of " +
                                                 std::to_string(position) + " bytes.");
    }

    if (b.GetAvailableSize() < 64)
    {
        b.Resize(position + 64, "BP4Serializer::MakeHeader " + fileType);
    }

    const std::string majorVersion(std::to_string(ADIOS2_VERSION_MAJOR));
    const char minorVersionChar = '0' + ADIOS2_VERSION_MINOR;
    const std::string minorVersion(1, minorVersionChar);
    const std::string patchVersion(std::to_string(ADIOS2_VERSION_PATCH));

    // byte 0-31: Readable tag
    if (position != m_VersionTagPosition)
    {
        helper::Throw<std::runtime_error>("Toolkit", "format::bp::BP4Serializer", "MakeHeader",
                                          "Version Tag "
                                          "position mismatch");
    }
    std::string versionLongTag("ADIOS-BP v" + majorVersion + "." + minorVersion + "." +
                               patchVersion + " ");
    size_t maxTypeLen = m_VersionTagLength - versionLongTag.size();
    const std::string fileTypeStr = fileType.substr(0, maxTypeLen);
    versionLongTag += fileTypeStr;
    const size_t versionLongTagSize = versionLongTag.size();
    if (versionLongTagSize < m_VersionTagLength)
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), versionLongTagSize);
        position += m_VersionTagLength - versionLongTagSize;
    }
    else if (versionLongTagSize > m_VersionTagLength)
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), m_VersionTagLength);
    }
    else
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), m_VersionTagLength);
    }

    // byte 32-35: MAJOR MINOR PATCH Unused

    lf_CopyVersionChar(majorVersion, buffer, position);
    lf_CopyVersionChar(minorVersion, buffer, position);
    lf_CopyVersionChar(patchVersion, buffer, position);
    ++position;

    // Note: Reader does process and use bytes 36-38 in
    // BP4Deserialize.cpp::ParseMetadataIndex().
    // Order and position must match there.

    // byte 36: endianness
    if (position != m_EndianFlagPosition)
    {
        helper::Throw<std::runtime_error>("Toolkit", "format::bp::BP4Serializer", "MakeHeader",
                                          "Endian Flag "
                                          "position mismatch");
    }
    const uint8_t endianness = helper::IsLittleEndian() ? 0 : 1;
    helper::CopyToBuffer(buffer, position, &endianness);

    // byte 37: BP Version 4
    if (position != m_BPVersionPosition)
    {
        helper::Throw<std::runtime_error>("Toolkit", "format::bp::BP4Serializer", "MakeHeader",
                                          "Active Flag "
                                          "position mismatch");
    }
    const uint8_t version = 4;
    helper::CopyToBuffer(buffer, position, &version);

    // byte 38: Active flag (used in Index Table only)
    if (position != m_ActiveFlagPosition)
    {
        helper::Throw<std::runtime_error>("Toolkit", "format::bp::BP4Serializer", "MakeHeader",
                                          "Active Flag "
                                          "position mismatch");
    }
    const uint8_t activeFlag = (isActive ? 1 : 0);
    helper::CopyToBuffer(buffer, position, &activeFlag);

    // byte 39: unused
    position += 1;

    // byte 40-63: unused
    position += 24;
    absolutePosition = position;
}

void BP4Serializer::PutProcessGroupIndex(const std::string &ioName, const std::string hostLanguage,
                                         const std::vector<std::string> &transportsTypes) noexcept
{
    m_Profiler.Start("buffering");
    std::vector<char> &metadataBuffer = m_MetadataSet.PGIndex.Buffer;

    std::vector<char> &dataBuffer = m_Data.m_Buffer;
    size_t &dataPosition = m_Data.m_Position;
    const size_t pgBeginPosition = dataPosition;

    // write a block identifier [PGI
    const char pgi[] = "[PGI"; //  don't write \0!
    helper::CopyToBuffer(dataBuffer, dataPosition, pgi, sizeof(pgi) - 1);

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
    helper::InsertU64(metadataBuffer, m_Data.m_AbsolutePosition + m_PreDataFileLength);

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
    m_Data.m_AbsolutePosition += dataPosition - pgBeginPosition;
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

size_t BP4Serializer::CloseData(core::IO &io)
{
    m_Profiler.Start("buffering");
    size_t dataEndsAt = m_Data.m_Position;
    if (!m_IsClosed)
    {
        if (m_MetadataSet.DataPGIsOpen)
        {
            SerializeDataBuffer(io);
        }
        dataEndsAt = m_Data.m_Position;

        SerializeMetadataInData(false, false);

        if (m_Profiler.m_IsActive)
        {
            m_Profiler.m_Bytes.at("buffering") = m_Data.m_AbsolutePosition;
        }

        m_Aggregator.Close();
        m_IsClosed = true;
    }

    m_Profiler.Stop("buffering");
    return dataEndsAt;
}

size_t BP4Serializer::CloseStream(core::IO &io, const bool addMetadata)
{
    m_Profiler.Start("buffering");
    if (m_MetadataSet.DataPGIsOpen)
    {
        SerializeDataBuffer(io);
    }
    size_t dataEndsAt = m_Data.m_Position;
    SerializeMetadataInData(false, addMetadata);

    if (m_Profiler.m_IsActive)
    {
        m_Profiler.m_Bytes.at("buffering") += m_Data.m_Position;
    }
    m_Profiler.Stop("buffering");

    return dataEndsAt;
}

/* Reset the local metadata indices */
void BP4Serializer::ResetAllIndices()
{
    m_MetadataSet.PGIndex.Buffer.resize(0);
    m_MetadataSet.PGIndex.LastUpdatedPosition = 0;
    m_MetadataSet.DataPGCount = 0;
    m_MetadataSet.DataPGLengthPosition = 0;
    m_MetadataSet.DataPGVarsCount = 0;
    m_MetadataSet.DataPGVarsCountPosition = 0;

    // for (auto &variableIndexPair : m_MetadataSet.VarsIndices)
    // {
    //     const std::string &variableName = variableIndexPair.first;
    //     SerialElementIndex &index = variableIndexPair.second;
    //     const size_t headersize = 15 + 8 + variableName.size();
    //     index.Buffer.resize(headersize);
    //     index.Count = 0;
    //     index.LastUpdatedPosition = headersize;
    //     index.Valid = false; // reset the flag to indicate the variable is
    //     not
    //                          // valid after the "endstep" call
    // }

    // for (auto &attributesIndexPair : m_MetadataSet.AttributesIndices)
    // {
    //     const std::string &attributesName = attributesIndexPair.first;
    //     SerialElementIndex &index = attributesIndexPair.second;
    //     const size_t headersize = 15 + 8 + attributesName.size();
    //     index.Buffer.resize(headersize);
    //     index.Count = 0;
    //     index.Valid = false; // reset the flag to indicate the variable is
    //     not
    //                          // valid after the "endstep" call
    // }
    m_MetadataSet.AttributesIndices.clear();
    m_MetadataSet.VarsIndices.clear();
}

/* Reset the metadata index table*/
void BP4Serializer::ResetMetadataIndexTable() { m_MetadataIndexTable.clear(); }

void BP4Serializer::AggregateCollectiveMetadata(helper::Comm const &comm, BufferSTL &bufferSTL,
                                                const bool inMetadataBuffer)
{
    m_Profiler.Start("buffering");
    m_Profiler.Start("meta_sort_merge");

    AggregateCollectiveMetadataIndices(comm, bufferSTL);

    int rank = comm.Rank();
    if (rank == 0)
    {
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
void BP4Serializer::SerializeDataBuffer(core::IO &io) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    // vars count and Length (only for PG)
    helper::CopyToBuffer(buffer, m_MetadataSet.DataPGVarsCountPosition,
                         &m_MetadataSet.DataPGVarsCount);
    // without record itself and vars count
    // Note: m_MetadataSet.DataPGVarsCount has been incremented by 4
    // in previous CopyToBuffer operation!
    const uint64_t varsLength = position - m_MetadataSet.DataPGVarsCountPosition - 8;
    helper::CopyToBuffer(buffer, m_MetadataSet.DataPGVarsCountPosition, &varsLength);

    // each attribute is only written to output once
    size_t attributesSizeInData = GetAttributesSizeInData(io);
    if (attributesSizeInData)
    {
        attributesSizeInData += 12; // count + length + end ID
        // Take care not to shrink the buffer size in this.
        // Otherwise, growing the buffer again later on is expensive
        // - not because of reallocation, but because of (re)initialization
        // with zeros by BufferSTL::Resize.
        const size_t minSize = position + attributesSizeInData + 4;
        if (m_Data.m_Buffer.size() < minSize)
        {
            m_Data.Resize(minSize, "when writing Attributes in rank=0\n");
        }
        PutAttributes(io);
    }
    else
    {
        const size_t minSize = position + 12 + 4;
        if (m_Data.m_Buffer.size() < minSize)
        {
            m_Data.Resize(minSize, "for empty Attributes\n");
        }
        // Attribute index header for zero attributes: 0, 0LL
        // Resize() already takes care of this
        position += 12;
        absolutePosition += 12;
    }

    // write a block identifier PGI]
    const char pgiend[] = "PGI]"; // no \0
    helper::CopyToBuffer(buffer, position, pgiend, sizeof(pgiend) - 1);
    absolutePosition += sizeof(pgiend) - 1;

    // Finish writing pg group length INCLUDING the record itself and
    // including the closing padding but NOT the opening [PGI
    const uint64_t dataPGLength = position - m_MetadataSet.DataPGLengthPosition;
    helper::CopyToBuffer(buffer, m_MetadataSet.DataPGLengthPosition, &dataPGLength);

    m_MetadataSet.DataPGIsOpen = false;
}

void BP4Serializer::SerializeMetadataInData(const bool updateAbsolutePosition, const bool inData)
{
    auto lf_SetIndexCountLength = [](std::unordered_map<std::string, SerialElementIndex> &indices,
                                     uint32_t &count, uint64_t &length) {
        count = static_cast<uint32_t>(indices.size());
        length = 0;
        for (auto &indexPair : indices) // set each index length
        {
            auto &indexBuffer = indexPair.second.Buffer;
            // const uint32_t indexLength =
            //     static_cast<uint32_t>(indexBuffer.size() - 4);
            // size_t indexLengthPosition = 0;
            // helper::CopyToBuffer(indexBuffer, indexLengthPosition,
            //                      &indexLength);

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

void BP4Serializer::AggregateCollectiveMetadataIndices(helper::Comm const &comm,
                                                       BufferSTL &outBufferSTL)
{
    int rank = comm.Rank();

    BufferSTL inBufferSTL;

    // pre-allocate with rank 0 data
    // size_t pgCount = 0; //< tracks global PG count
    if (rank == 0)
    {
        // assumes that things are more or less balanced
        m_PGIndicesInfo.clear();
        // m_PGRankIndices.reserve(m_MetadataSet.PGIndex.Buffer.size());

        m_VariableIndicesInfo.clear();
        // m_VariableRankIndices.reserve(m_MetadataSet.VarsIndices.size());

        m_AttributesIndicesInfo.clear();
        // m_AttributesRankIndices.reserve(m_MetadataSet.AttributesIndices.size());
    }

    auto lf_IndicesSize =
        [&](const std::unordered_map<std::string, SerialElementIndex> &indices) -> size_t

    {
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
        for (const auto &indexPair : indices)
        {
            const auto &buffer = indexPair.second.Buffer;
            helper::CopyToBuffer(m_SerializedIndices, position, buffer.data(), buffer.size());
        }
    };

    auto lf_SerializeAllIndices = [&](helper::Comm const &comm, const int rank) {
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

    auto lf_LocatePGIndices =
        [&](std::unordered_map<size_t, std::vector<std::tuple<size_t, size_t, size_t>>>
                &pgIndicesInfo,
            const int rankSource, const std::vector<char> &serialized, const size_t position,
            const size_t endPosition) {
            size_t stepStartPosition = position;
            size_t stepBuffersize = 0;
            size_t pgCountPerStep = 0;
            size_t localPosition = position;
            uint32_t currentStep = 0;

            while (localPosition < endPosition)
            {
                size_t indexPosition = localPosition;
                const ProcessGroupIndex header = ReadProcessGroupIndexHeader(
                    serialized, indexPosition, helper::IsLittleEndian());
                if (header.Step == currentStep)
                {
                    stepBuffersize += header.Length + 2;
                    pgCountPerStep++;
                }
                else
                {
                    // found a new step
                    if (currentStep == 0)
                    {
                        // start of going through a new pg buffer
                        stepStartPosition = localPosition;
                        stepBuffersize = header.Length + 2;
                        pgCountPerStep = 1;
                        currentStep = header.Step;
                    }
                    else
                    {
                        // record the pg info of previous step
                        std::tuple<size_t, size_t, size_t> stepPGIndexTuple =
                            std::make_tuple(pgCountPerStep, stepStartPosition, stepBuffersize);
                        auto search = pgIndicesInfo.find(currentStep);
                        if (search == pgIndicesInfo.end())
                        {
                            // the time step hasn't been added to the
                            // unordered_map, add it
                            pgIndicesInfo.emplace(
                                currentStep, std::vector<std::tuple<size_t, size_t, size_t>>());
                            pgIndicesInfo[currentStep].push_back(stepPGIndexTuple);
                        }
                        else
                        {
                            pgIndicesInfo[currentStep].push_back(stepPGIndexTuple);
                        }
                        stepStartPosition = localPosition;
                        stepBuffersize = header.Length + 2;
                        pgCountPerStep = 1;
                        currentStep = header.Step;
                    }
                }
                localPosition += header.Length + 2;
                if (localPosition >= endPosition)
                {
                    // record the pg info of the last step
                    std::tuple<size_t, size_t, size_t> stepPGIndexTuple =
                        std::make_tuple(pgCountPerStep, stepStartPosition, stepBuffersize);
                    auto search = pgIndicesInfo.find(currentStep);
                    if (search == pgIndicesInfo.end())
                    {
                        // the time step hasn't been added to the unordered_map,
                        // add it
                        pgIndicesInfo.emplace(currentStep,
                                              std::vector<std::tuple<size_t, size_t, size_t>>());
                        pgIndicesInfo[currentStep].push_back(stepPGIndexTuple);
                    }
                    else
                    {
                        pgIndicesInfo[currentStep].push_back(stepPGIndexTuple);
                    }
                }
            }
        };

    auto lf_GetCharacteristics = [&](const std::vector<char> &buffer, size_t &position,
                                     const uint8_t dataType, uint8_t &count, uint32_t &length,
                                     uint32_t &timeStep)

    {
        const DataTypes dataTypeEnum = static_cast<DataTypes>(dataType);
        const bool isLittleEndian = helper::IsLittleEndian();

        switch (dataTypeEnum)
        {

#define make_case(T)                                                                               \
    case (TypeTraits<T>::type_enum): {                                                             \
        size_t irrelevant;                                                                         \
        const auto characteristics = ReadElementIndexCharacteristics<T>(                           \
            buffer, position, TypeTraits<T>::type_enum, irrelevant, true, isLittleEndian);         \
        count = characteristics.EntryCount;                                                        \
        length = characteristics.EntryLength;                                                      \
        timeStep = characteristics.Statistics.Step;                                                \
        break;                                                                                     \
    }
            ADIOS2_FOREACH_STDTYPE_1ARG(make_case)
#undef make_case

        case (type_string_array): {
            size_t irrelevant;
            const auto characteristics = ReadElementIndexCharacteristics<std::string>(
                buffer, position, type_string_array, irrelevant, true, isLittleEndian);
            count = characteristics.EntryCount;
            length = characteristics.EntryLength;
            timeStep = characteristics.Statistics.Step;
            break;
        }

        default:
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BP4Serializer", "AggregateCollectiveMetadataIndices",
                "type " + std::to_string(dataType) + " not supported in BP4 Metadata Merge");

        } // end switch
    };

    auto lf_LocateVarIndices =
        [&](std::unordered_map<
                size_t, std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>>
                &indicesInfo,
            const int rankSource, const std::vector<char> &serialized, const size_t position,
            const size_t endPosition)

    {
        // uint32_t currentStep = 0;

        size_t localPosition = position;
        while (localPosition < endPosition)
        {
            // std::cout << "var localPosition: " << localPosition << std::endl;
            size_t indexPosition = localPosition;
            const ElementIndexHeader header =
                ReadElementIndexHeader(serialized, indexPosition, helper::IsLittleEndian());

            // std::cout << "var indexPosition after ReadElementIndexHeader: "
            // << indexPosition << std::endl;

            uint8_t count = 0;
            uint32_t length = 0;
            uint32_t timeStep = 0;

            lf_GetCharacteristics(serialized, indexPosition, header.DataType, count, length,
                                  timeStep);

            size_t varIndexBufferSize = static_cast<size_t>(header.Length) + 4;

            auto stepSearch = indicesInfo.find(timeStep);
            if (stepSearch == indicesInfo.end())
            {

                // the time step hasn't been added to the unordered_map, add it
                std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>
                    varIndexInfo;
                varIndexInfo.emplace(header.Name, std::vector<std::tuple<size_t, size_t>>());
                std::tuple<size_t, size_t> varIndexTuple =
                    std::make_tuple(localPosition, varIndexBufferSize);
                varIndexInfo[header.Name].push_back(varIndexTuple);
                indicesInfo.emplace(timeStep, varIndexInfo);
            }
            else
            {

                auto varSearch = indicesInfo[timeStep].find(header.Name);
                if (varSearch == indicesInfo[timeStep].end())
                {
                    // found a new variable at this step
                    indicesInfo[timeStep].emplace(header.Name,
                                                  std::vector<std::tuple<size_t, size_t>>());
                    std::tuple<size_t, size_t> varIndexTuple =
                        std::make_tuple(localPosition, varIndexBufferSize);
                    indicesInfo[timeStep][header.Name].push_back(varIndexTuple);
                }
                else
                {
                    // variable already exists, insert the location info of this
                    // variable for this rank
                    std::tuple<size_t, size_t> varIndexTuple =
                        std::make_tuple(localPosition, varIndexBufferSize);
                    indicesInfo[timeStep][header.Name].push_back(varIndexTuple);
                }
            }

            localPosition += varIndexBufferSize;
        }
    };

    auto lf_LocateAttrIndices =
        [&](std::unordered_map<
                size_t, std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>>
                &indicesInfo,
            const int rankSource, const std::vector<char> &serialized, const size_t position,
            const size_t endPosition)

    {
        size_t localPosition = position;
        while (localPosition < endPosition)
        {
            size_t indexPosition = localPosition;
            const ElementIndexHeader header =
                ReadElementIndexHeader(serialized, indexPosition, helper::IsLittleEndian());

            uint8_t count = 0;
            uint32_t length = 0;
            uint32_t timeStep = 0;

            lf_GetCharacteristics(serialized, indexPosition, header.DataType, count, length,
                                  timeStep);

            size_t attrIndexBufferSize = static_cast<size_t>(header.Length) + 4;

            if (indicesInfo[timeStep][header.Name].size() == 1)
            {
                localPosition += attrIndexBufferSize;
                continue;
            }

            auto stepSearch = indicesInfo.find(timeStep);
            if (stepSearch == indicesInfo.end())
            {
                // the time step hasn't been added to the unordered_map, add it
                std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>
                    attrIndexInfo;
                attrIndexInfo.emplace(header.Name, std::vector<std::tuple<size_t, size_t>>());
                std::tuple<size_t, size_t> attrIndexTuple =
                    std::make_tuple(localPosition, attrIndexBufferSize);
                attrIndexInfo[header.Name].push_back(attrIndexTuple);
                indicesInfo.emplace(timeStep, attrIndexInfo);
            }
            else
            {
                // found a new step of a new rank
                auto attrSearch = indicesInfo[timeStep].find(header.Name);
                if (attrSearch == indicesInfo[timeStep].end())
                {
                    // found a new attribute at this step
                    indicesInfo[timeStep].emplace(header.Name,
                                                  std::vector<std::tuple<size_t, size_t>>());
                    std::tuple<size_t, size_t> attrIndexTuple =
                        std::make_tuple(localPosition, attrIndexBufferSize);
                    indicesInfo[timeStep][header.Name].push_back(attrIndexTuple);
                }
                else
                {
                    // attribute already exists, insert the location info of
                    // this attribute for this rank
                    std::tuple<size_t, size_t> attrIndexTuple =
                        std::make_tuple(localPosition, attrIndexBufferSize);
                    indicesInfo[timeStep][header.Name].push_back(attrIndexTuple);
                }
            }

            localPosition += attrIndexBufferSize;
        }
    };

    auto lf_LocateAllIndices = [&](const int rankSource, const std::vector<size_t> headerInfo,
                                   const std::vector<char> &serialized, const size_t position)

    {
        const size_t rankIndicesSize = headerInfo[0];
        const size_t variablesIndexOffset = headerInfo[1] + position;
        const size_t attributesIndexOffset = headerInfo[2] + position;
        size_t localPosition = position + 36;

        size_t endPosition = variablesIndexOffset;
        // first deserialize pg indices
        lf_LocatePGIndices(m_PGIndicesInfo, rankSource, serialized, localPosition, endPosition);

        // deserialize variable indices
        localPosition = variablesIndexOffset;
        endPosition = attributesIndexOffset;

        lf_LocateVarIndices(m_VariableIndicesInfo, rankSource, serialized, localPosition,
                            endPosition);

        // deserialize attributes indices
        localPosition = attributesIndexOffset;
        endPosition = rankIndicesSize + 4 + position;
        // attributes are constant and unique across ranks
        lf_LocateAttrIndices(m_AttributesIndicesInfo, rankSource, serialized, localPosition,
                             endPosition);
    };

    auto lf_SortMergeIndices =
        [&](const std::unordered_map<size_t, std::vector<std::tuple<size_t, size_t, size_t>>>
                &pgIndicesInfo,
            const std::unordered_map<
                size_t, std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>>
                &varIndicesInfo,
            const std::unordered_map<
                size_t, std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>>
                &attrIndicesInfo,
            const std::vector<char> &serialized) {
            auto &position = outBufferSTL.m_Position;
            auto &buffer = outBufferSTL.m_Buffer;

            std::vector<size_t> timeSteps;
            timeSteps.reserve(pgIndicesInfo.size());
            for (auto const &pair : pgIndicesInfo)
            {
                timeSteps.push_back(pair.first);
            }
            std::sort(timeSteps.begin(), timeSteps.end());
            m_MetadataIndexTable[rank] = {};
            for (auto t : timeSteps)
            {
                std::vector<uint64_t> ptrs;

                const uint64_t pgIndexStart = position;
                ptrs.push_back(pgIndexStart);
                std::vector<std::tuple<size_t, size_t, size_t>> perStepPGIndicesInfo =
                    pgIndicesInfo.at(t);
                size_t perStepPGCountPosition = position;
                position += 16; // skip the pgcount and pglength
                uint64_t perStepPGCountU64 = 0;

                for (auto &item : perStepPGIndicesInfo)
                {
                    size_t start = std::get<1>(item);
                    size_t length = std::get<2>(item);
                    std::copy(serialized.begin() + start, serialized.begin() + start + length,
                              buffer.begin() + position);
                    position += length;
                    perStepPGCountU64 += std::get<0>(item);
                }
                uint64_t perStepPGLengthU64 = position - perStepPGCountPosition - 16;
                helper::CopyToBuffer(buffer, perStepPGCountPosition, &perStepPGCountU64);
                helper::CopyToBuffer(buffer, perStepPGCountPosition, &perStepPGLengthU64);

                const uint64_t variablesIndexStart = position;
                ptrs.push_back(variablesIndexStart);
                const auto itvars = varIndicesInfo.find(t);
                if (itvars != varIndicesInfo.end())
                {
                    const auto perStepVarIndicesInfo = itvars->second;
                    size_t perStepVarCountPosition = position;
                    const uint32_t perStepVarCountU32 =
                        static_cast<uint32_t>(perStepVarIndicesInfo.size());
                    helper::CopyToBuffer(buffer, perStepVarCountPosition, &perStepVarCountU32);
                    position += 12; // skip for count and length

                    for (auto const &pair : perStepVarIndicesInfo)
                    {
                        const size_t entryLengthPosition = position;
                        size_t headerStartPosition = std::get<0>(pair.second[0]);
                        size_t localPosition = headerStartPosition;
                        ElementIndexHeader header =
                            ReadElementIndexHeader(serialized, localPosition);
                        size_t headerSize = localPosition - headerStartPosition;

                        position += headerSize; // skip the header
                        uint64_t setsCount = 0;
                        for (auto &item : pair.second)
                        {
                            size_t start = std::get<0>(item);
                            size_t length = std::get<1>(item);
                            std::copy(serialized.begin() + start + headerSize,
                                      serialized.begin() + start + length,
                                      buffer.begin() + position);
                            position += length - headerSize;
                            setsCount += header.CharacteristicsSetsCount;
                        }
                        const uint32_t entryLength =
                            static_cast<uint32_t>(position - entryLengthPosition - 4);
                        size_t backPosition = entryLengthPosition;
                        helper::CopyToBuffer(buffer, backPosition, &entryLength);
                        helper::CopyToBuffer(buffer, backPosition,
                                             &serialized[headerStartPosition + 4],
                                             headerSize - 8 - 4);
                        helper::CopyToBuffer(buffer, backPosition, &setsCount);
                    }
                    const uint64_t perStepVarLengthU64 =
                        static_cast<uint64_t>(position - perStepVarCountPosition - 8);
                    helper::CopyToBuffer(buffer, perStepVarCountPosition, &perStepVarLengthU64);
                }
                else
                {
                    size_t perStepVarCountPosition = position;
                    const uint32_t perStepVarCountU32 = static_cast<uint32_t>(0);
                    helper::CopyToBuffer(buffer, perStepVarCountPosition, &perStepVarCountU32);
                    const uint64_t perStepVarLengthU64 = static_cast<uint64_t>(0);
                    helper::CopyToBuffer(buffer, perStepVarCountPosition, &perStepVarLengthU64);
                    position += 12; // skip for count and length
                }

                const uint64_t attributesIndexStart = position;
                ptrs.push_back(attributesIndexStart);

                const auto itattrs = attrIndicesInfo.find(t);
                if (itattrs != attrIndicesInfo.end())
                {
                    const auto perStepAttrIndicesInfo = itattrs->second;
                    size_t perStepAttrCountPosition = position;
                    const uint32_t perStepAttrCountU32 =
                        static_cast<uint32_t>(perStepAttrIndicesInfo.size());
                    helper::CopyToBuffer(buffer, perStepAttrCountPosition, &perStepAttrCountU32);
                    position += 12; // skip for length

                    for (auto const &pair : perStepAttrIndicesInfo)
                    {
                        for (auto &item : pair.second)
                        {
                            size_t start = std::get<0>(item);
                            size_t length = std::get<1>(item);
                            std::copy(serialized.begin() + start,
                                      serialized.begin() + start + length,
                                      buffer.begin() + position);
                            position += length;
                        }
                    }
                    const uint64_t perStepAttrLengthU64 =
                        static_cast<uint64_t>(position - perStepAttrCountPosition - 8);
                    helper::CopyToBuffer(buffer, perStepAttrCountPosition, &perStepAttrLengthU64);
                }
                else
                {
                    size_t perStepAttrCountPosition = position;
                    const uint32_t perStepAttrCountU32 = static_cast<uint32_t>(0);
                    helper::CopyToBuffer(buffer, perStepAttrCountPosition, &perStepAttrCountU32);
                    const uint64_t perStepAttrLengthU64 = static_cast<uint64_t>(0);
                    helper::CopyToBuffer(buffer, perStepAttrCountPosition, &perStepAttrLengthU64);
                    position += 12; // skip for count and length
                }

                const uint64_t currentStepEndPos = position;
                ptrs.push_back(currentStepEndPos);
                m_MetadataIndexTable[rank][t] = ptrs;
            }
        };

    // BODY of function starts here
    lf_SerializeAllIndices(comm, rank); // Set m_SerializedIndices

    comm.GathervVectors(m_SerializedIndices, inBufferSTL.m_Buffer, inBufferSTL.m_Position, 0);

    // deserialize, it's all local inside rank 0
    if (rank == 0)
    {
        const size_t serializedSize = inBufferSTL.m_Position;
        const std::vector<char> &serialized = inBufferSTL.m_Buffer;
        size_t serializedPosition = 0;
        std::vector<size_t> headerInfo(4);
        const bool isLittleEndian = helper::IsLittleEndian();

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

            lf_LocateAllIndices(rankSource, headerInfo, serialized, serializedPosition);
            serializedPosition += headerInfo[0] + 4;
        }
    }

    // now merge (and sort variables and attributes) indices
    if (rank == 0)
    {
        auto &buffer = outBufferSTL.m_Buffer;
        const std::vector<char> &serialized = inBufferSTL.m_Buffer;

        size_t totalStep = m_PGIndicesInfo.size();
        size_t perStepExtraSize = 16 + 12 + 12;
        const size_t totalExtraSize = totalStep * perStepExtraSize;
        buffer.reserve(outBufferSTL.m_Position + inBufferSTL.m_Position + totalExtraSize);
        buffer.resize(outBufferSTL.m_Position + inBufferSTL.m_Position + totalExtraSize);
        lf_SortMergeIndices(m_PGIndicesInfo, m_VariableIndicesInfo, m_AttributesIndicesInfo,
                            serialized);
    }
}

#define declare_template_instantiation(T)                                                          \
    void BP4Serializer::DoPutAttributeInData(const core::Attribute<T> &attribute,                  \
                                             Stats<T> &stats) noexcept                             \
    {                                                                                              \
        PutAttributeInDataCommon(attribute, stats);                                                \
    }
ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

//------------------------------------------------------------------------------
// Explicit instantiation of only public templates

#define declare_template_instantiation(T)                                                          \
    template void BP4Serializer::PutVariablePayload(                                               \
        const core::Variable<T> &, const typename core::Variable<T>::BPInfo &, const bool,         \
        typename core::Variable<T>::Span *) noexcept;                                              \
                                                                                                   \
    template void BP4Serializer::PutVariableMetadata(                                              \
        const core::Variable<T> &, const typename core::Variable<T>::BPInfo &, const bool,         \
        typename core::Variable<T>::Span *) noexcept;

ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

//------------------------------------------------------------------------------

#define declare_template_instantiation(T)                                                          \
    template void BP4Serializer::PutSpanMetadata(                                                  \
        const core::Variable<T> &, const typename core::Variable<T>::BPInfo &,                     \
        const typename core::Variable<T>::Span &) noexcept;

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace format
} // end namespace adios2
