/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP3Serializer.tcc
 *
 *  Created on: Apr 11, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP3_BP3SERIALIZER_TCC_
#define ADIOS2_TOOLKIT_FORMAT_BP3_BP3SERIALIZER_TCC_

#include "BP3Serializer.h"

#include <algorithm> // std::all_of, std::fill_n
#include <array>

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace format
{

template <class T>
void BP3Serializer::PutVariableMetadata(const core::Variable<T> &variable,
                                        const typename core::Variable<T>::BPInfo &blockInfo,
                                        const bool sourceRowMajor,
                                        typename core::Variable<T>::Span *span) noexcept
{
    auto lf_SetOffset = [&](uint64_t &offset) {
        if (m_Aggregator.m_IsActive && !m_Aggregator.m_IsAggregator)
        {
            offset = static_cast<uint64_t>(m_Data.m_Position);
        }
        else
        {
            offset = static_cast<uint64_t>(m_Data.m_AbsolutePosition + m_PreDataFileLength);
        }
    };

    m_Profiler.Start("buffering");

    Stats<T> stats = GetBPStats<T>(variable.m_SingleValue, blockInfo, sourceRowMajor);

    // Get new Index or point to existing index
    bool isNew = true; // flag to check if variable is new
    SerialElementIndex &variableIndex =
        GetSerialElementIndex(variable.m_Name, m_MetadataSet.VarsIndices, isNew);
    stats.MemberID = variableIndex.MemberID;

    lf_SetOffset(stats.Offset);
    PutVariableMetadataInData(variable, blockInfo, stats, span);
    lf_SetOffset(stats.PayloadOffset);
    if (span != nullptr)
    {
        span->m_PayloadPosition = m_Data.m_Position;
    }

    // write to metadata  index
    PutVariableMetadataInIndex(variable, blockInfo, stats, isNew, variableIndex, span);
    ++m_MetadataSet.DataPGVarsCount;

    m_Profiler.Stop("buffering");
}

template <class T>
inline void BP3Serializer::PutVariablePayload(const core::Variable<T> &variable,
                                              const typename core::Variable<T>::BPInfo &blockInfo,
                                              const bool sourceRowMajor,
                                              typename core::Variable<T>::Span *span) noexcept
{
    m_Profiler.Start("buffering");
    if (span != nullptr)
    {
        const size_t blockSize = helper::GetTotalSize(blockInfo.Count);
        if (span->m_Value != T{})
        {
            T *itBegin = reinterpret_cast<T *>(m_Data.m_Buffer.data() + m_Data.m_Position);

            // TODO: does std::fill_n have a bug in gcc or due to optimizations
            // this is impossible due to memory alignment? This seg faults in
            // Release mode only . Even RelWithDebInfo works, replacing with
            // explicit loop below using access operator [] std::fill_n(itBegin,
            // blockSize, span->m_Value);

            for (size_t i = 0; i < blockSize; ++i)
            {
                itBegin[i] = span->m_Value;
            }
        }

        m_Data.m_Position += blockSize * sizeof(T);
        m_Data.m_AbsolutePosition += blockSize * sizeof(T);
        m_Profiler.Stop("buffering");
        return;
    }

    if (blockInfo.Operations.empty())
    {
        PutPayloadInBuffer(variable, blockInfo, sourceRowMajor);
    }
    else
    {
        PutOperationPayloadInBuffer(variable, blockInfo);
    }

    m_Profiler.Stop("buffering");
}

template <class T>
void BP3Serializer::PutSpanMetadata(const core::Variable<T> &variable,
                                    const typename core::Variable<T>::Span &span) noexcept
{
    if (m_Parameters.StatsLevel > 0)
    {
        // Get Min/Max from populated data
        m_Profiler.Start("minmax");
        T min, max;
        helper::GetMinMaxThreads(span.Data(), span.Size(), min, max, m_Parameters.Threads,
                                 variable.m_MemSpace);
        m_Profiler.Stop("minmax");

        // Put min/max in variable index
        SerialElementIndex &variableIndex = m_MetadataSet.VarsIndices.at(variable.m_Name);
        auto &buffer = variableIndex.Buffer;

        const size_t minPosition = span.m_MinMaxMetadataPositions.first;
        const size_t maxPosition = span.m_MinMaxMetadataPositions.second;
        std::copy(&min, &min + 1, reinterpret_cast<T *>(buffer.data() + minPosition));
        std::copy(&max, &max + 1, reinterpret_cast<T *>(buffer.data() + maxPosition));
    }
}

// PRIVATE
template <class T>
size_t BP3Serializer::PutAttributeHeaderInData(const core::Attribute<T> &attribute,
                                               Stats<T> &stats) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;

    // will go back to write length
    const size_t attributeLengthPosition = position;
    position += 4; // skip length

    helper::CopyToBuffer(buffer, position, &stats.MemberID);
    PutNameRecord(attribute.m_Name, buffer, position);
    position += 2; // skip path

    constexpr int8_t no = 'n';
    helper::CopyToBuffer(buffer, position,
                         &no); // not associated with a Variable

    return attributeLengthPosition;
}

template <class T>
void BP3Serializer::PutAttributeLengthInData(const core::Attribute<T> &attribute, Stats<T> &stats,
                                             const size_t attributeLengthPosition) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    // back to attribute length
    size_t backPosition = attributeLengthPosition;
    uint32_t len = static_cast<uint32_t>(position - attributeLengthPosition);
    helper::CopyToBuffer(buffer, backPosition, &len);

    absolutePosition += position - attributeLengthPosition;
}

template <>
inline void BP3Serializer::PutAttributeInDataCommon(const core::Attribute<std::string> &attribute,
                                                    Stats<std::string> &stats) noexcept
{
    const size_t attributeLengthPosition = PutAttributeHeaderInData(attribute, stats);

    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    uint8_t dataType = TypeTraits<std::string>::type_enum;
    if (!attribute.m_IsSingleValue)
    {
        dataType = type_string_array;
    }
    helper::CopyToBuffer(buffer, position, &dataType);

    // here record payload offset
    stats.PayloadOffset = absolutePosition + position - attributeLengthPosition;

    if (dataType == type_string)
    {
        const uint32_t dataSize = static_cast<uint32_t>(attribute.m_DataSingleValue.size());
        helper::CopyToBuffer(buffer, position, &dataSize);
        helper::CopyToBuffer(buffer, position, attribute.m_DataSingleValue.data(),
                             attribute.m_DataSingleValue.size());
    }
    else if (dataType == type_string_array)
    {
        const uint32_t elements = static_cast<uint32_t>(attribute.m_Elements);
        helper::CopyToBuffer(buffer, position, &elements);

        for (size_t s = 0; s < attribute.m_Elements; ++s)
        {
            // include zero terminated
            const std::string element(attribute.m_DataArray[s] + '\0');

            const uint32_t elementSize = static_cast<uint32_t>(element.size());

            helper::CopyToBuffer(buffer, position, &elementSize);
            helper::CopyToBuffer(buffer, position, element.data(), element.size());
        }
    }

    PutAttributeLengthInData(attribute, stats, attributeLengthPosition);
}

template <class T>
void BP3Serializer::PutAttributeInDataCommon(const core::Attribute<T> &attribute,
                                             Stats<T> &stats) noexcept
{
    const size_t attributeLengthPosition = PutAttributeHeaderInData(attribute, stats);

    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    uint8_t dataType = TypeTraits<T>::type_enum;
    helper::CopyToBuffer(buffer, position, &dataType);

    // here record payload offset
    stats.PayloadOffset = absolutePosition + position - attributeLengthPosition;

    const uint32_t dataSize = static_cast<uint32_t>(attribute.m_Elements * sizeof(T));
    helper::CopyToBuffer(buffer, position, &dataSize);

    if (attribute.m_IsSingleValue) // single value
    {
        helper::CopyToBuffer(buffer, position, &attribute.m_DataSingleValue);
    }
    else // array
    {
        helper::CopyToBuffer(buffer, position, attribute.m_DataArray.data(), attribute.m_Elements);
    }

    PutAttributeLengthInData(attribute, stats, attributeLengthPosition);
}

template <>
inline BP3Serializer::Stats<std::string>
BP3Serializer::GetBPStats(const bool /*singleValue*/,
                          const typename core::Variable<std::string>::BPInfo & /*blockInfo*/,
                          const bool /*isRowMajor*/) noexcept
{
    Stats<std::string> stats;
    stats.Step = m_MetadataSet.TimeStep;
    stats.FileIndex = GetFileIndex();
    return stats;
}

template <class T>
BP3Serializer::Stats<T>
BP3Serializer::GetBPStats(const bool singleValue,
                          const typename core::Variable<T>::BPInfo &blockInfo,
                          const bool isRowMajor) noexcept
{
    Stats<T> stats;
    stats.Step = m_MetadataSet.TimeStep;
    stats.FileIndex = GetFileIndex();

    // added to support span
    if (blockInfo.Data == nullptr)
    {
        stats.Min = {};
        stats.Max = {};
        return stats;
    }

    if (singleValue)
    {
        stats.Value = *blockInfo.Data;
        stats.Min = stats.Value;
        stats.Max = stats.Value;
        return stats;
    }

    if (m_Parameters.StatsLevel > 0)
    {
        m_Profiler.Start("minmax");
        if (blockInfo.MemoryStart.empty())
        {
            const std::size_t valuesSize = helper::GetTotalSize(blockInfo.Count);
            helper::GetMinMaxThreads(blockInfo.Data, valuesSize, stats.Min, stats.Max,
                                     m_Parameters.Threads, blockInfo.MemSpace);
        }
        else // non-contiguous memory min/max
        {
            helper::GetMinMaxSelection(blockInfo.Data, blockInfo.MemoryCount, blockInfo.MemoryStart,
                                       blockInfo.Count, isRowMajor, stats.Min, stats.Max,
                                       blockInfo.MemSpace);
        }
        m_Profiler.Stop("minmax");
    }

    return stats;
}

template <class T>
void BP3Serializer::PutVariableMetadataInData(const core::Variable<T> &variable,
                                              const typename core::Variable<T>::BPInfo &blockInfo,
                                              const Stats<T> &stats,
                                              const typename core::Variable<T>::Span *span) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    // for writing length at the end
    const size_t varLengthPosition = position;
    position += 8; // skip var length (8)

    helper::CopyToBuffer(buffer, position, &stats.MemberID);

    PutNameRecord(variable.m_Name, buffer, position);
    position += 2; // skip path

    const uint8_t dataType = TypeTraits<T>::type_enum;
    helper::CopyToBuffer(buffer, position, &dataType);

    constexpr char no = 'n'; // isDimension
    helper::CopyToBuffer(buffer, position, &no);

    const uint8_t dimensions = static_cast<uint8_t>(variable.m_Count.size());
    helper::CopyToBuffer(buffer, position, &dimensions); // count

    // 27 is from 9 bytes for each: var y/n + local, var y/n + global dimension,
    // var y/n + global offset, changed for characteristic
    uint16_t dimensionsLength = 27 * dimensions;
    helper::CopyToBuffer(buffer, position, &dimensionsLength); // length

    PutDimensionsRecord(variable.m_Count, variable.m_Shape, variable.m_Start, buffer, position);

    // CHARACTERISTICS
    PutVariableCharacteristics(variable, blockInfo, stats, buffer, position);

    // here align pointer for span
    if (span != nullptr)
    {
        const size_t padLengthPosition = position;
        constexpr std::array<uint8_t, 5> zeros = {0, 0, 0, 0, 0};
        // skip 1 for paddingLength and 4 for VMD] ending
        helper::CopyToBuffer(buffer, position, zeros.data(), 5);
        // here check for the next aligned pointer
        const size_t extraBytes = m_Data.Align<T>();
        const std::string pad = std::string(extraBytes, '\0') + "VMD]";

        size_t backPosition = padLengthPosition;
        const uint8_t padLength = static_cast<uint8_t>(pad.size());
        helper::CopyToBuffer(buffer, backPosition, &padLength);
        helper::CopyToBuffer(buffer, backPosition, pad.c_str(), pad.size());

        position += extraBytes;
    }

    // Back to varLength including payload size and pad
    // not need to remove its own size (8) from length from bpdump
    const uint64_t varLength = static_cast<uint64_t>(
        position - varLengthPosition + helper::PayloadSize(blockInfo.Data, blockInfo.Count));

    size_t backPosition = varLengthPosition;
    helper::CopyToBuffer(buffer, backPosition, &varLength);

    absolutePosition += position - varLengthPosition;
}

template <>
inline void BP3Serializer::PutVariableMetadataInData(
    const core::Variable<std::string> &variable,
    const typename core::Variable<std::string>::BPInfo &blockInfo, const Stats<std::string> &stats,
    const typename core::Variable<std::string>::Span * /*span*/) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    // for writing length at the end
    const size_t varLengthPosition = position;
    position += 8; // skip var length (8)

    helper::CopyToBuffer(buffer, position, &stats.MemberID);

    PutNameRecord(variable.m_Name, buffer, position);
    position += 2; // skip path

    const uint8_t dataType = TypeTraits<std::string>::type_enum;
    helper::CopyToBuffer(buffer, position, &dataType);

    constexpr char no = 'n'; // is dimension is deprecated
    helper::CopyToBuffer(buffer, position, &no);

    const uint8_t dimensions = static_cast<uint8_t>(blockInfo.Count.size());
    helper::CopyToBuffer(buffer, position, &dimensions); // count

    uint16_t dimensionsLength = 27 * dimensions;
    helper::CopyToBuffer(buffer, position, &dimensionsLength); // length

    PutDimensionsRecord(blockInfo.Count, blockInfo.Shape, blockInfo.Start, buffer, position);

    position += 5; // skipping characteristics

    // Back to varLength including payload size
    // not need to remove its own size (8) from length from bpdump
    const uint64_t varLength = static_cast<uint64_t>(
        position - varLengthPosition + helper::PayloadSize(blockInfo.Data, blockInfo.Count));

    size_t backPosition = varLengthPosition;
    helper::CopyToBuffer(buffer, backPosition, &varLength);

    absolutePosition += position - varLengthPosition;
}

template <class T>
void BP3Serializer::PutVariableMetadataInIndex(const core::Variable<T> &variable,
                                               const typename core::Variable<T>::BPInfo &blockInfo,
                                               const Stats<T> &stats, const bool isNew,
                                               SerialElementIndex &index,
                                               typename core::Variable<T>::Span *span) noexcept
{
    auto &buffer = index.Buffer;

    if (isNew) // write variable header
    {
        buffer.insert(buffer.end(), 4, '\0'); // skip var length (4)
        helper::InsertToBuffer(buffer, &stats.MemberID);
        buffer.insert(buffer.end(), 2, '\0'); // skip group name
        PutNameRecord(variable.m_Name, buffer);
        buffer.insert(buffer.end(), 2, '\0'); // skip path

        const uint8_t dataType = TypeTraits<T>::type_enum;
        helper::InsertToBuffer(buffer, &dataType);

        // Characteristics Sets Count in Metadata
        index.Count = 1;
        helper::InsertToBuffer(buffer, &index.Count);

        // For updating absolute offsets in agreggation
        index.LastUpdatedPosition = buffer.size();
    }
    else // update characteristics sets count
    {
        if (m_Parameters.StatsLevel > 0)
        {
            ++index.Count;
            // fixed since group and path are not printed
            size_t setsCountPosition = 15 + variable.m_Name.size();
            helper::CopyToBuffer(buffer, setsCountPosition, &index.Count);
        }
    }

    PutVariableCharacteristics(variable, blockInfo, stats, buffer, span);
}

template <class T>
void BP3Serializer::PutBoundsRecord(const bool singleValue, const Stats<T> &stats,
                                    uint8_t &characteristicsCounter,
                                    std::vector<char> &buffer) noexcept
{
    if (singleValue)
    {
        PutCharacteristicRecord(characteristic_value, characteristicsCounter, stats.Min, buffer);
    }
    else
    {
        if (m_Parameters.StatsLevel > 0) // default verbose
        {
            PutCharacteristicRecord(characteristic_min, characteristicsCounter, stats.Min, buffer);
            PutCharacteristicRecord(characteristic_max, characteristicsCounter, stats.Max, buffer);
        }
    }
}

template <class T>
void BP3Serializer::PutBoundsRecord(const bool singleValue, const Stats<T> &stats,
                                    uint8_t &characteristicsCounter, std::vector<char> &buffer,
                                    size_t &position) noexcept
{
    if (singleValue)
    {
        const uint8_t id = characteristic_value;
        helper::CopyToBuffer(buffer, position, &id);
        // special case required by bpdump
        const uint16_t length = sizeof(T);
        helper::CopyToBuffer(buffer, position, &length);
        helper::CopyToBuffer(buffer, position, &stats.Min);
        ++characteristicsCounter;
    }
    else
    {
        if (m_Parameters.StatsLevel > 0) // default min and max only
        {
            PutCharacteristicRecord(characteristic_min, characteristicsCounter, stats.Min, buffer,
                                    position);

            PutCharacteristicRecord(characteristic_max, characteristicsCounter, stats.Max, buffer,
                                    position);
        }
    }
}

template <>
inline void BP3Serializer::PutVariableCharacteristics(
    const core::Variable<std::string> &variable,
    const core::Variable<std::string>::BPInfo &blockInfo, const Stats<std::string> &stats,
    std::vector<char> &buffer, typename core::Variable<std::string>::Span * /*span*/) noexcept
{
    const size_t characteristicsCountPosition = buffer.size();
    // skip characteristics count(1) + length (4)
    buffer.insert(buffer.end(), 5, '\0');
    uint8_t characteristicsCounter = 0;

    PutCharacteristicRecord(characteristic_time_index, characteristicsCounter, stats.Step, buffer);

    PutCharacteristicRecord(characteristic_file_index, characteristicsCounter, stats.FileIndex,
                            buffer);

    uint8_t characteristicID = characteristic_value;
    helper::InsertToBuffer(buffer, &characteristicID);
    PutNameRecord(*blockInfo.Data, buffer);
    ++characteristicsCounter;

    characteristicID = characteristic_dimensions;
    helper::InsertToBuffer(buffer, &characteristicID);
    const uint8_t dimensions = static_cast<uint8_t>(blockInfo.Count.size());
    helper::InsertToBuffer(buffer, &dimensions); // count
    const uint16_t dimensionsLength = static_cast<uint16_t>(24 * dimensions);
    helper::InsertToBuffer(buffer, &dimensionsLength); // length
    PutDimensionsRecord(blockInfo.Count, blockInfo.Shape, blockInfo.Start, buffer);
    ++characteristicsCounter;

    PutCharacteristicRecord(characteristic_offset, characteristicsCounter, stats.Offset, buffer);

    PutCharacteristicRecord(characteristic_payload_offset, characteristicsCounter,
                            stats.PayloadOffset, buffer);

    // END OF CHARACTERISTICS

    // Back to characteristics count and length
    size_t backPosition = characteristicsCountPosition;
    helper::CopyToBuffer(buffer, backPosition,
                         &characteristicsCounter); // count (1)

    // remove its own length (4) + characteristic counter (1)
    const uint32_t characteristicsLength =
        static_cast<uint32_t>(buffer.size() - characteristicsCountPosition - 4 - 1);

    helper::CopyToBuffer(buffer, backPosition,
                         &characteristicsLength); // length
}

template <class T>
void BP3Serializer::PutVariableCharacteristics(const core::Variable<T> &variable,
                                               const typename core::Variable<T>::BPInfo &blockInfo,
                                               const Stats<T> &stats, std::vector<char> &buffer,
                                               typename core::Variable<T>::Span *span) noexcept
{
    // going back at the end
    const size_t characteristicsCountPosition = buffer.size();
    // skip characteristics count(1) + length (4)
    buffer.insert(buffer.end(), 5, '\0');
    uint8_t characteristicsCounter = 0;

    // DIMENSIONS
    PutCharacteristicRecord(characteristic_time_index, characteristicsCounter, stats.Step, buffer);

    PutCharacteristicRecord(characteristic_file_index, characteristicsCounter, stats.FileIndex,
                            buffer);

    if (blockInfo.Data != nullptr || span != nullptr)
    {
        if (m_Parameters.StatsLevel > 0 && span != nullptr)
        {
            span->m_MinMaxMetadataPositions.first = buffer.size() + 1;
            span->m_MinMaxMetadataPositions.second = buffer.size() + 2 + sizeof(T);
        }

        PutBoundsRecord(variable.m_SingleValue, stats, characteristicsCounter, buffer);
    }

    uint8_t characteristicID = characteristic_dimensions;
    helper::InsertToBuffer(buffer, &characteristicID);
    const uint8_t dimensions = static_cast<uint8_t>(blockInfo.Count.size());
    helper::InsertToBuffer(buffer, &dimensions); // count
    const uint16_t dimensionsLength = static_cast<uint16_t>(24 * dimensions);
    helper::InsertToBuffer(buffer, &dimensionsLength); // length
    PutDimensionsRecord(blockInfo.Count, blockInfo.Shape, blockInfo.Start, buffer);
    ++characteristicsCounter;

    PutCharacteristicRecord(characteristic_offset, characteristicsCounter, stats.Offset, buffer);

    PutCharacteristicRecord(characteristic_payload_offset, characteristicsCounter,
                            stats.PayloadOffset, buffer);

    if (blockInfo.Operations.size())
    {
        const bool isZeroCount = std::all_of(blockInfo.Count.begin(), blockInfo.Count.end(),
                                             [](const size_t i) { return i == 0; });

        // do not compress if count dimensions are all zero
        if (!isZeroCount)
        {
            characteristicID = characteristic_transform_type;
            helper::InsertToBuffer(buffer, &characteristicID);
            PutCharacteristicOperation(variable, blockInfo, buffer);
            ++characteristicsCounter;
        }
    }

    // END OF CHARACTERISTICS

    // Back to characteristics count and length
    size_t backPosition = characteristicsCountPosition;
    helper::CopyToBuffer(buffer, backPosition,
                         &characteristicsCounter); // count (1)

    // remove its own length (4) + characteristic counter (1)
    const uint32_t characteristicsLength =
        static_cast<uint32_t>(buffer.size() - characteristicsCountPosition - 4 - 1);

    helper::CopyToBuffer(buffer, backPosition,
                         &characteristicsLength); // length
}

template <class T>
void BP3Serializer::PutVariableCharacteristics(const core::Variable<T> &variable,
                                               const typename core::Variable<T>::BPInfo &blockInfo,
                                               const Stats<T> &stats, std::vector<char> &buffer,
                                               size_t &position) noexcept
{
    // going back at the end
    const size_t characteristicsCountPosition = position;
    // skip characteristics count(1) + length (4)
    position += 5;
    uint8_t characteristicsCounter = 0;

    // DIMENSIONS
    uint8_t characteristicID = characteristic_dimensions;
    helper::CopyToBuffer(buffer, position, &characteristicID);

    const uint8_t dimensions = static_cast<uint8_t>(blockInfo.Count.size());
    helper::CopyToBuffer(buffer, position, &dimensions); // count
    const uint16_t dimensionsLength = static_cast<uint16_t>(24 * dimensions);
    helper::CopyToBuffer(buffer, position, &dimensionsLength); // length
    PutDimensionsRecord(blockInfo.Count, blockInfo.Shape, blockInfo.Start, buffer, position, true);
    ++characteristicsCounter;

    // VALUE for SCALAR or STAT min, max for ARRAY
    if (blockInfo.Data != nullptr)
    {
        PutBoundsRecord(variable.m_SingleValue, stats, characteristicsCounter, buffer, position);
    }
    // END OF CHARACTERISTICS

    // Back to characteristics count and length
    size_t backPosition = characteristicsCountPosition;
    helper::CopyToBuffer(buffer, backPosition, &characteristicsCounter);

    // remove its own length (4) + characteristic counter (1)
    const uint32_t characteristicsLength =
        static_cast<uint32_t>(position - characteristicsCountPosition - 4 - 1);
    helper::CopyToBuffer(buffer, backPosition, &characteristicsLength);
}

} // end namespace format
} // end namespace adios2

#endif // ADIOS2_TOOLKIT_FORMAT_BP3_BP3SERIALIZER_TCC_
