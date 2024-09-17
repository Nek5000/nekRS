/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPSerializer.cpp
 *
 *  Created on: Sep 16, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BPSerializer.h"
#include "BPSerializer.tcc"

#ifdef _WIN32
#pragma warning(disable : 4503) // Windows complains about MergeSerializeIndex
                                // long types
#endif

namespace adios2
{
namespace format
{

std::mutex BPSerializer::m_Mutex;

// PUBLIC
BPSerializer::BPSerializer(const helper::Comm &comm, const uint8_t version)
: BPBase(comm), m_Version(version)
{
}

void BPSerializer::SerializeData(core::IO &io, const bool advanceStep)
{
    m_Profiler.Start("buffering");
    SerializeDataBuffer(io);
    if (advanceStep)
    {
        ++m_MetadataSet.TimeStep;
        ++m_MetadataSet.CurrentStep;
    }
    m_Profiler.Stop("buffering");
}

std::string BPSerializer::GetRankProfilingJSON(
    const std::vector<std::string> &transportsTypes,
    const std::vector<profiling::IOChrono *> &transportsProfilers) noexcept
{
    auto lf_WriterTimer = [](std::string &rankLog, const profiling::Timer &timer) {
        rankLog += ", \"" + timer.m_Process + "_" + timer.GetShortUnits() +
                   "\": " + std::to_string(timer.m_ProcessTime);
    };

    // prepare string dictionary per rank
    std::string rankLog("{ \"rank\": " + std::to_string(m_RankMPI));

    auto &profiler = m_Profiler;

    std::string timeDate(profiler.m_Timers.at("buffering").m_LocalTimeDate);
    timeDate.pop_back();
    // avoid whitespace
    std::replace(timeDate.begin(), timeDate.end(), ' ', '_');

    rankLog += ", \"start\": \"" + timeDate + "\"";
    rankLog += ", \"threads\": " + std::to_string(m_Parameters.Threads);
    rankLog += ", \"bytes\": " + std::to_string(profiler.m_Bytes.at("buffering"));

    for (const auto &timerPair : profiler.m_Timers)
    {
        const profiling::Timer &timer = timerPair.second;
        rankLog += ", \"" + timer.m_Process + "_" + timer.GetShortUnits() +
                   "\": " + std::to_string(timer.m_ProcessTime);
    }

    const size_t transportsSize = transportsTypes.size();

    for (unsigned int t = 0; t < transportsSize; ++t)
    {
        rankLog += ", \"transport_" + std::to_string(t) + "\": { ";
        rankLog += "\"type\": \"" + transportsTypes[t] + "\"";

        for (const auto &transportTimerPair : transportsProfilers[t]->m_Timers)
        {
            lf_WriterTimer(rankLog, transportTimerPair.second);
        }
        rankLog += "}";
    }
    rankLog += " }"; // end rank entry

    return rankLog;
}

std::vector<char> BPSerializer::AggregateProfilingJSON(const std::string &rankLog) const
{
    // Gather sizes
    const size_t rankLogSize = rankLog.size();
    std::vector<size_t> rankLogsSizes = m_Comm.GatherValues(rankLogSize);

    // Gatherv JSON per rank
    std::vector<char> profilingJSON(3);
    const std::string header("[\n");
    const std::string footer("\n]\n");
    size_t gatheredSize = 0;
    size_t position = 0;

    if (m_RankMPI == 0) // pre-allocate in destination
    {
        gatheredSize = std::accumulate(rankLogsSizes.begin(), rankLogsSizes.end(), size_t(0));

        profilingJSON.resize(gatheredSize + header.size() + footer.size() - 2);
        helper::CopyToBuffer(profilingJSON, position, header.c_str(), header.size());
    }

    m_Comm.GathervArrays(rankLog.c_str(), rankLog.size(), rankLogsSizes.data(),
                         rankLogsSizes.size(), &profilingJSON[position]);

    if (m_RankMPI == 0) // add footer to close JSON
    {
        position += gatheredSize - 2;
        helper::CopyToBuffer(profilingJSON, position, footer.c_str(), footer.size());
    }

    return profilingJSON;
}

void BPSerializer::PutNameRecord(const std::string name, std::vector<char> &buffer) noexcept
{
    const uint16_t length = static_cast<uint16_t>(name.size());
    helper::InsertToBuffer(buffer, &length);
    helper::InsertToBuffer(buffer, name.c_str(), name.size());
}

void BPSerializer::PutNameRecord(const std::string name, std::vector<char> &buffer,
                                 size_t &position) noexcept
{
    const uint16_t length = static_cast<uint16_t>(name.length());
    helper::CopyToBuffer(buffer, position, &length);
    helper::CopyToBuffer(buffer, position, name.c_str(), length);
}

void BPSerializer::PutDimensionsRecord(const Dims &localDimensions, const Dims &globalDimensions,
                                       const Dims &offsets, std::vector<char> &buffer) noexcept
{
    if (offsets.empty() && globalDimensions.empty())
    { // local array
        for (const auto localDimension : localDimensions)
        {
            helper::InsertU64(buffer, localDimension);
            buffer.insert(buffer.end(), 2 * sizeof(uint64_t), '\0');
        }
    }
    else if (offsets.empty())
    { // joined array has no offsets but has global dims
        for (unsigned int d = 0; d < localDimensions.size(); ++d)
        {
            helper::InsertU64(buffer, localDimensions[d]);
            helper::InsertU64(buffer, globalDimensions[d]);
            buffer.insert(buffer.end(), sizeof(uint64_t), '\0');
        }
    }
    else
    { // global array
        for (unsigned int d = 0; d < localDimensions.size(); ++d)
        {
            helper::InsertU64(buffer, localDimensions[d]);
            helper::InsertU64(buffer, globalDimensions[d]);
            helper::InsertU64(buffer, offsets[d]);
        }
    }
}

void BPSerializer::PutDimensionsRecord(const Dims &localDimensions, const Dims &globalDimensions,
                                       const Dims &offsets, std::vector<char> &buffer,
                                       size_t &position, const bool isCharacteristic) noexcept
{
    auto lf_CopyDimension = [](std::vector<char> &buffer, size_t &position, const size_t dimension,
                               const bool isCharacteristic) {
        if (!isCharacteristic)
        {
            constexpr char no = 'n';
            helper::CopyToBuffer(buffer, position, &no);
        }

        const uint64_t dimension64 = static_cast<uint64_t>(dimension);

        helper::CopyToBuffer(buffer, position, &dimension64);
    };

    // BODY Starts here
    if (offsets.empty() && globalDimensions.empty())
    {
        unsigned int globalBoundsSkip = 18;
        if (isCharacteristic)
        {
            globalBoundsSkip = 16;
        }

        for (const auto &localDimension : localDimensions)
        {
            lf_CopyDimension(buffer, position, localDimension, isCharacteristic);
            position += globalBoundsSkip;
        }
    }
    else if (offsets.empty())
    {
        // joined array has no offsets but has global dims
        size_t zeroOffset = 0;
        for (unsigned int d = 0; d < localDimensions.size(); ++d)
        {
            lf_CopyDimension(buffer, position, localDimensions[d], isCharacteristic);
            lf_CopyDimension(buffer, position, globalDimensions[d], isCharacteristic);
            lf_CopyDimension(buffer, position, zeroOffset, isCharacteristic);
        }
    }
    else
    {
        for (unsigned int d = 0; d < localDimensions.size(); ++d)
        {
            lf_CopyDimension(buffer, position, localDimensions[d], isCharacteristic);
            lf_CopyDimension(buffer, position, globalDimensions[d], isCharacteristic);
            lf_CopyDimension(buffer, position, offsets[d], isCharacteristic);
        }
    }
}

void BPSerializer::PutMinifooter(const uint64_t pgIndexStart, const uint64_t variablesIndexStart,
                                 const uint64_t attributesIndexStart, std::vector<char> &buffer,
                                 size_t &position, const bool addSubfiles)
{
    auto lf_CopyVersionChar = [](const std::string version, std::vector<char> &buffer,
                                 size_t &position) {
        helper::CopyToBuffer(buffer, position, version.c_str());
    };

    const std::string majorVersion(std::to_string(ADIOS2_VERSION_MAJOR));
    const char minorVersionChar = '0' + ADIOS2_VERSION_MINOR;
    const std::string minorVersion(1, minorVersionChar);
    const std::string patchVersion(std::to_string(ADIOS2_VERSION_PATCH));

    const std::string versionLongTag("ADIOS-BP v" + majorVersion + "." + minorVersion + "." +
                                     patchVersion);
    const size_t versionLongTagSize = versionLongTag.size();
    if (versionLongTagSize < 24)
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), versionLongTagSize);
        position += 24 - versionLongTagSize;
    }
    else
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), 24);
    }

    lf_CopyVersionChar(majorVersion, buffer, position);
    lf_CopyVersionChar(minorVersion, buffer, position);
    lf_CopyVersionChar(patchVersion, buffer, position);
    ++position;

    helper::CopyToBuffer(buffer, position, &pgIndexStart);
    helper::CopyToBuffer(buffer, position, &variablesIndexStart);
    helper::CopyToBuffer(buffer, position, &attributesIndexStart);

    const uint8_t endianness = helper::IsLittleEndian() ? 0 : 1;
    helper::CopyToBuffer(buffer, position, &endianness);

    if (addSubfiles)
    {
        const uint8_t zeros1 = 0;
        helper::CopyToBuffer(buffer, position, &zeros1);
        helper::CopyToBuffer(buffer, position, &m_Version);
    }
    else
    {
        const uint16_t zeros2 = 0;
        helper::CopyToBuffer(buffer, position, &zeros2);
    }
    helper::CopyToBuffer(buffer, position, &m_Version);
}

void BPSerializer::UpdateOffsetsInMetadata()
{
    auto lf_UpdatePGIndexOffsets = [&]() {
        auto &buffer = m_MetadataSet.PGIndex.Buffer;
        size_t &currentPosition = m_MetadataSet.PGIndex.LastUpdatedPosition;
        const bool isLittleEndian = helper::IsLittleEndian();

        while (currentPosition < buffer.size())
        {
            ProcessGroupIndex pgIndex =
                ReadProcessGroupIndexHeader(buffer, currentPosition, isLittleEndian);

            const uint64_t updatedOffset =
                pgIndex.Offset + static_cast<uint64_t>(m_Data.m_AbsolutePosition);
            currentPosition -= sizeof(uint64_t);
            helper::CopyToBuffer(buffer, currentPosition, &updatedOffset);
        }
    };

    auto lf_UpdateIndexOffsets = [&](SerialElementIndex &index) {
        auto &buffer = index.Buffer;

        // First get the type:
        size_t headerPosition = 0;
        ElementIndexHeader header =
            ReadElementIndexHeader(buffer, headerPosition, helper::IsLittleEndian());
        const DataTypes dataTypeEnum = static_cast<DataTypes>(header.DataType);

        size_t &currentPosition = index.LastUpdatedPosition;

        while (currentPosition < buffer.size())
        {
            switch (dataTypeEnum)
            {

            case (type_string): {
                // do nothing, strings are obtained from metadata
                currentPosition = buffer.size();
                break;
            }

#define make_case(T)                                                                               \
    case (TypeTraits<T>::type_enum): {                                                             \
        UpdateIndexOffsetsCharacteristics<T>(currentPosition, TypeTraits<T>::type_enum, buffer);   \
        break;                                                                                     \
    }
                ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(make_case)
#undef make_case

            default:
                // TODO: complex, long double
                helper::Throw<std::invalid_argument>(
                    "Toolkit", "format::bp::BPSerializer", "UpdateOffsetsInMetadata",
                    "type " + std::to_string(header.DataType) +
                        " not supported in updating aggregated offsets");

            } // end switch
        }
    };

    // BODY OF FUNCTION STARTS HERE
    if (m_Aggregator.m_IsAggregator)
    {
        return;
    }

    // PG Indices
    lf_UpdatePGIndexOffsets();

    // Variable Indices
    for (auto &varIndexPair : m_MetadataSet.VarsIndices)
    {
        SerialElementIndex &index = varIndexPair.second;
        lf_UpdateIndexOffsets(index);
    }
}

void BPSerializer::MergeSerializeIndices(
    const std::unordered_map<std::string, std::vector<SerialElementIndex>> &nameRankIndices,
    helper::Comm const &comm, BufferSTL &bufferSTL)
{
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
                "Toolkit", "format::bp::BPSerializer", "MergeSerializeIndices",
                "type " + std::to_string(dataType) + " not supported in BP3 Metadata Merge");

        } // end switch
    };

    auto lf_MergeRankSerial = [&](const std::vector<SerialElementIndex> &indices,
                                  BufferSTL &bufferSTL) {
        auto &bufferOut = bufferSTL.m_Buffer;
        auto &positionOut = bufferSTL.m_Position;

        // extract header
        ElementIndexHeader header;
        // index non-empty buffer
        size_t firstRank = 0;
        // index positions per rank
        std::vector<size_t> positions(indices.size(), 0);
        // merge index length
        size_t headerSize = 0;

        const bool isLittleEndian = helper::IsLittleEndian();

        for (size_t r = 0; r < indices.size(); ++r)
        {
            const auto &buffer = indices[r].Buffer;
            if (buffer.empty())
            {
                continue;
            }
            size_t &position = positions[r];

            header = ReadElementIndexHeader(buffer, position, isLittleEndian);
            firstRank = r;

            headerSize = position;
            break;
        }

        if (header.DataType == std::numeric_limits<uint8_t>::max() - 1)
        {
            helper::Throw<std::runtime_error>(
                "Toolkit", "format::bp::BPSerializer", "MergeSerializeIndices",
                "invalid data type for variable " + header.Name + "when writing metadata index");
        }

        // move all positions to headerSize
        for (size_t r = 0; r < indices.size(); ++r)
        {
            const auto &buffer = indices[r].Buffer;
            if (buffer.empty())
            {
                continue;
            }
            positions[r] = headerSize;
        }

        uint64_t setsCount = 0;
        unsigned int currentTimeStep = 1;
        bool marching = true;

        const size_t entryLengthPosition = positionOut;
        positionOut += headerSize;

        while (marching)
        {
            marching = false;

            for (size_t r = firstRank; r < indices.size(); ++r)
            {
                const auto &buffer = indices[r].Buffer;
                if (buffer.empty())
                {
                    continue;
                }

                auto &position = positions[r];
                if (position < buffer.size())
                {
                    marching = true;
                }
                else
                {
                    continue;
                }

                uint8_t count = 0;
                uint32_t length = 0;
                uint32_t timeStep = static_cast<uint32_t>(currentTimeStep);

                while (timeStep == currentTimeStep)
                {
                    size_t localPosition = position;
                    lf_GetCharacteristics(buffer, localPosition, header.DataType, count, length,
                                          timeStep);

                    if (timeStep != currentTimeStep)
                    {
                        break;
                    }

                    ++setsCount;

                    helper::CopyToBuffer(bufferOut, positionOut, &buffer[position], length + 5);

                    position += length + 5;

                    if (position >= buffer.size())
                    {
                        break;
                    }
                }
            }
            ++currentTimeStep;
        }

        const uint32_t entryLength = static_cast<uint32_t>(positionOut - entryLengthPosition - 4);

        size_t backPosition = entryLengthPosition;
        helper::CopyToBuffer(bufferOut, backPosition, &entryLength);
        helper::CopyToBuffer(bufferOut, backPosition, &indices[firstRank].Buffer[4],
                             headerSize - 8 - 4);
        helper::CopyToBuffer(bufferOut, backPosition, &setsCount);
    };

    auto lf_MergeRank = [&](const std::vector<SerialElementIndex> &indices, BufferSTL &bufferSTL) {
        ElementIndexHeader header;
        size_t firstRank = 0;
        // index positions per rank
        std::vector<size_t> positions(indices.size(), 0);
        // merge index length
        size_t headerSize = 0;

        const bool isLittleEndian = helper::IsLittleEndian();

        for (size_t r = 0; r < indices.size(); ++r)
        {
            const auto &buffer = indices[r].Buffer;
            if (buffer.empty())
            {
                continue;
            }
            size_t &position = positions[r];

            header = ReadElementIndexHeader(buffer, position, isLittleEndian);
            firstRank = r;

            headerSize = position;
            break;
        }

        if (header.DataType == std::numeric_limits<uint8_t>::max() - 1)
        {
            helper::Throw<std::runtime_error>("Toolkit", "format::bp::BPSerializer",
                                              "MergeSerializeIndices",
                                              "invalid data type for variable " + header.Name +
                                                  "when writing collective metadata");
        }

        // move all positions to headerSize
        for (size_t r = 0; r < indices.size(); ++r)
        {
            const auto &buffer = indices[r].Buffer;
            if (buffer.empty())
            {
                continue;
            }
            positions[r] = headerSize;
        }

        uint64_t setsCount = 0;
        unsigned int currentTimeStep = 1;
        bool marching = true;
        std::vector<char> sorted;

        while (marching)
        {
            marching = false;

            for (size_t r = firstRank; r < indices.size(); ++r)
            {
                const auto &buffer = indices[r].Buffer;
                if (buffer.empty())
                {
                    continue;
                }

                auto &position = positions[r];
                if (position < buffer.size())
                {
                    marching = true;
                }
                else
                {
                    continue;
                }

                uint8_t count = 0;
                uint32_t length = 0;
                uint32_t timeStep = static_cast<uint32_t>(currentTimeStep);

                while (timeStep == currentTimeStep)
                {
                    size_t localPosition = position;
                    lf_GetCharacteristics(buffer, localPosition, header.DataType, count, length,
                                          timeStep);

                    if (timeStep != currentTimeStep)
                    {
                        break;
                    }

                    ++setsCount;

                    // here copy to sorted buffer
                    helper::InsertToBuffer(sorted, &buffer[position], length + 5);

                    position += length + 5;

                    if (position >= buffer.size())
                    {
                        break;
                    }
                }
            }
            ++currentTimeStep;
        }

        const uint32_t entryLength = static_cast<uint32_t>(headerSize + sorted.size() - 4);
        // Copy header to metadata buffer, need mutex here
        {
            std::lock_guard<std::mutex> lock(m_Mutex);
            auto &buffer = bufferSTL.m_Buffer;
            auto &position = bufferSTL.m_Position;

            helper::CopyToBuffer(buffer, position, &entryLength);
            helper::CopyToBuffer(buffer, position, &indices[firstRank].Buffer[4],
                                 headerSize - 8 - 4);
            helper::CopyToBuffer(buffer, position, &setsCount);
            helper::CopyToBuffer(buffer, position, sorted.data(), sorted.size());
        }
    };

    auto lf_MergeRankRange =
        [&](const std::unordered_map<std::string, std::vector<SerialElementIndex>> &nameRankIndices,
            const std::vector<std::string> &names, const size_t start, const size_t end,
            BufferSTL &bufferSTL)

    {
        for (auto i = start; i < end; ++i)
        {
            auto itIndex = nameRankIndices.find(names[i]);
            lf_MergeRank(itIndex->second, bufferSTL);
        }
    };

    // BODY OF FUNCTION STARTS HERE
    if (m_Parameters.Threads == 1) // enforcing serial version for now
    {
        for (const auto &rankIndices : nameRankIndices)
        {
            lf_MergeRankSerial(rankIndices.second, bufferSTL);
        }
        return;
    }

    const size_t elements = nameRankIndices.size();
    const size_t stride = elements / m_Parameters.Threads;        // elements per thread
    const size_t last = stride + elements % m_Parameters.Threads; // remainder to last

    std::vector<std::thread> threads;
    threads.reserve(m_Parameters.Threads);

    // copy names in order to use threads
    std::vector<std::string> names;
    names.reserve(nameRankIndices.size());

    for (const auto &nameRankIndexPair : nameRankIndices)
    {
        names.push_back(nameRankIndexPair.first);
    }

    for (unsigned int t = 0; t < m_Parameters.Threads; ++t)
    {
        const size_t start = stride * t;
        size_t end = start + stride;

        if (t == m_Parameters.Threads - 1)
        {
            end = start + last;
        }

        threads.push_back(std::thread(lf_MergeRankRange, std::ref(nameRankIndices), std::ref(names),
                                      start, end, std::ref(bufferSTL)));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }
}

uint32_t BPSerializer::GetFileIndex() const noexcept
{
    if (m_Aggregator.m_IsActive)
    {
        return static_cast<uint32_t>(m_Aggregator.m_SubStreamIndex);
    }

    return static_cast<uint32_t>(m_RankMPI);
}

size_t BPSerializer::GetAttributesSizeInData(core::IO &io) const noexcept
{
    size_t attributesSizeInData = 0;

    auto &attributes = io.GetAttributes();

    for (const auto &attribute : attributes)
    {
        const DataType type = attribute.second->m_Type;

        // each attribute is only written to output once
        // so filter out the ones already written
        auto it = m_SerializedAttributes.find(attribute.first);
        if (it != m_SerializedAttributes.end())
        {
            continue;
        }

        if (type == DataType::Struct)
        {
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        const std::string name = attribute.first;                                                  \
        const core::Attribute<T> &attribute = *io.InquireAttribute<T>(name);                       \
        attributesSizeInData += GetAttributeSizeInData<T>(attribute);                              \
    }
        ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_type)
#undef declare_type
    }

    return attributesSizeInData;
}

void BPSerializer::PutAttributes(core::IO &io)
{
    const auto &attributes = io.GetAttributes();

    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    const size_t attributesCountPosition = position;

    // count is known ahead of time
    const uint32_t attributesCount = static_cast<uint32_t>(attributes.size());
    helper::CopyToBuffer(buffer, position, &attributesCount);

    // will go back
    const size_t attributesLengthPosition = position;
    position += 8; // skip attributes length

    absolutePosition += position - attributesCountPosition;

    uint32_t memberID = 0;

    for (const auto &attributePair : attributes)
    {
        const std::string name(attributePair.first);
        const DataType type(attributePair.second->m_Type);

        // each attribute is only written to output once
        // so filter out the ones already written
        auto it = m_SerializedAttributes.find(name);
        if (it != m_SerializedAttributes.end())
        {
            continue;
        }

        if (type == DataType::None)
        {
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        Stats<T> stats;                                                                            \
        stats.Offset = absolutePosition + m_PreDataFileLength;                                     \
        stats.MemberID = memberID;                                                                 \
        stats.Step = m_MetadataSet.TimeStep;                                                       \
        stats.FileIndex = GetFileIndex();                                                          \
        core::Attribute<T> &attribute = *io.InquireAttribute<T>(name);                             \
        PutAttributeInData(attribute, stats);                                                      \
        PutAttributeInIndex(attribute, stats);                                                     \
    }
        ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_type)
#undef declare_type

        ++memberID;
    }

    // complete attributes length
    const uint64_t attributesLength = static_cast<uint64_t>(position - attributesLengthPosition);

    size_t backPosition = attributesLengthPosition;
    helper::CopyToBuffer(buffer, backPosition, &attributesLength);
}

BPSerializer::SerialElementIndex &
BPSerializer::GetSerialElementIndex(const std::string &name,
                                    std::unordered_map<std::string, SerialElementIndex> &indices,
                                    bool &isNew) const noexcept
{
    auto itName = indices.find(name);
    if (itName == indices.end())
    {
        indices.emplace(name, SerialElementIndex(static_cast<uint32_t>(indices.size())));
        isNew = true;
        return indices.at(name);
    }

    isNew = false;
    return itName->second;
}

#define declare_template_instantiation(T)                                                          \
    template void BPSerializer::PutAttributeCharacteristicValueInIndex(                            \
        uint8_t &, const core::Attribute<T> &, std::vector<char> &) noexcept;                      \
                                                                                                   \
    template size_t BPSerializer::GetAttributeSizeInData(const core::Attribute<T> &)               \
        const noexcept;                                                                            \
                                                                                                   \
    template void BPSerializer::PutAttributeInData(const core::Attribute<T> &,                     \
                                                   Stats<T> &) noexcept;                           \
    template void BPSerializer::PutAttributeInIndex(const core::Attribute<T> &attribute,           \
                                                    const Stats<T> &stats) noexcept;

ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

#define declare_template_instantiation(T)                                                          \
                                                                                                   \
    template void BPSerializer::PutCharacteristicRecord(const uint8_t, uint8_t &, const T &,       \
                                                        std::vector<char> &) noexcept;             \
                                                                                                   \
    template void BPSerializer::PutCharacteristicRecord(const uint8_t, uint8_t &, const T &,       \
                                                        std::vector<char> &, size_t &) noexcept;   \
                                                                                                   \
    template void BPSerializer::PutCharacteristicOperation(                                        \
        const core::Variable<T> &, const typename core::Variable<T>::BPInfo &,                     \
        std::vector<char> &) noexcept;                                                             \
                                                                                                   \
    template void BPSerializer::PutOperationPayloadInBuffer(                                       \
        const core::Variable<T> &, const typename core::Variable<T>::BPInfo &);

ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

#define declare_template_instantiation(T)                                                          \
    template void BPSerializer::PutPayloadInBuffer(const core::Variable<T> &,                      \
                                                   const typename core::Variable<T>::BPInfo &,     \
                                                   const bool) noexcept;

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace format
} // end namespace adios2
