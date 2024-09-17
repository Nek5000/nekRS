/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Serializer.tcc
 *
 *  Created on: Aug 1, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 *              Lipeng Wan wanl@ornl.gov
 *              Norbert Podhorszki pnb@ornl.gov
 *
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP4_BP4SERIALIZER_TCC_
#define ADIOS2_TOOLKIT_FORMAT_BP4_BP4SERIALIZER_TCC_

#include "BP4Serializer.h"

#include <algorithm> // std::all_of
#include <array>
#include <iostream>

#include "adios2/helper/adiosFunctions.h"

#include <iostream>
#include <stddef.h>

namespace adios2
{
namespace format
{

template <class T>
inline void BP4Serializer::PutVariableMetadata(const core::Variable<T> &variable,
                                               const typename core::Variable<T>::BPInfo &blockInfo,
                                               const bool sourceRowMajor,
                                               typename core::Variable<T>::Span *span) noexcept
{
    auto lf_SetOffset = [&](uint64_t &offset) {
        if (m_Aggregator.m_IsActive && !m_Aggregator.m_IsAggregator)
        {
            offset = static_cast<uint64_t>(m_Data.m_Position + m_PreDataFileLength);
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
    variableIndex.Valid = true; // flag to indicate this variable is put at current step
    stats.MemberID = variableIndex.MemberID;

    lf_SetOffset(stats.Offset);
    m_LastVarLengthPosInBuffer = PutVariableMetadataInData(variable, blockInfo, stats, span);
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
inline void BP4Serializer::PutVariablePayload(const core::Variable<T> &variable,
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
        const bool isZeroCount = std::all_of(blockInfo.Count.begin(), blockInfo.Count.end(),
                                             [](const size_t i) { return i == 0; });
        if (!isZeroCount)
        {
            PutOperationPayloadInBuffer(variable, blockInfo);
        }
    }

    /* Now we can update the varLength including payload size including the
     * closing padding but NOT the opening [VMD
     */
    const uint64_t varLength =
        static_cast<uint64_t>(m_Data.m_Position - m_LastVarLengthPosInBuffer);
    size_t backPosition = m_LastVarLengthPosInBuffer;
    helper::CopyToBuffer(m_Data.m_Buffer, backPosition, &varLength);

    m_Profiler.Stop("buffering");
}

template <class T>
void BP4Serializer::PutSpanMetadata(const core::Variable<T> &variable,
                                    const typename core::Variable<T>::BPInfo &blockInfo,
                                    const typename core::Variable<T>::Span &span) noexcept
{
    if (m_Parameters.StatsLevel > 0)
    {
        // Get Min/Max from populated data
        m_Profiler.Start("minmax");
        Stats<T> stats;
        stats.Min = {};
        stats.Max = {};
        stats.SubBlockInfo = helper::DivideBlock(blockInfo.Count, m_Parameters.StatsBlockSize,
                                                 helper::BlockDivisionMethod::Contiguous);
        // set stats MinMaxs with the correct size
        helper::GetMinMaxSubblocks(span.Data(), blockInfo.Count, stats.SubBlockInfo, stats.MinMaxs,
                                   stats.Min, stats.Max, m_Parameters.Threads, blockInfo.MemSpace);
        m_Profiler.Stop("minmax");

        // Put min/max blocks in variable index
        SerialElementIndex &variableIndex = m_MetadataSet.VarsIndices.at(variable.m_Name);
        auto &buffer = variableIndex.Buffer;

        size_t minMaxPosition = span.m_MinMaxMetadataPositions.first;
        // here repopulate the metadata buffer
        uint8_t dummyCounter = 0;
        PutBoundsRecord(false, stats, dummyCounter, buffer, minMaxPosition);
    }
}

// PRIVATE
template <class T>
size_t BP4Serializer::PutAttributeHeaderInData(const core::Attribute<T> &attribute, Stats<T> &stats,
                                               const char *headerID,
                                               const size_t headerIDLength) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;

    helper::CopyToBuffer(buffer, position, headerID, headerIDLength);

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
void BP4Serializer::PutAttributeLengthInData(const core::Attribute<T> &attribute, Stats<T> &stats,
                                             const size_t attributeLengthPosition) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;

    // back to attribute length
    size_t backPosition = attributeLengthPosition;
    uint32_t len = static_cast<uint32_t>(position - attributeLengthPosition);
    helper::CopyToBuffer(buffer, backPosition, &len);
}

template <>
inline void BP4Serializer::PutAttributeInDataCommon(const core::Attribute<std::string> &attribute,
                                                    Stats<std::string> &stats) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    const size_t mdBeginPosition = position;

    // write a block identifier [AMD
    const char amd[] = "[AMD"; // no \0
    const size_t attributeLengthPosition =
        PutAttributeHeaderInData(attribute, stats, amd, sizeof(amd) - 1);

    uint8_t dataType = TypeTraits<std::string>::type_enum;
    if (!attribute.m_IsSingleValue)
    {
        dataType = type_string_array;
    }
    helper::CopyToBuffer(buffer, position, &dataType);

    // here record payload offset
    stats.PayloadOffset = absolutePosition + position - mdBeginPosition + m_PreDataFileLength;

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

    // write a block identifier AMD]
    const char amdend[] = "AMD]"; // no \0
    helper::CopyToBuffer(buffer, position, amdend, sizeof(amdend) - 1);

    PutAttributeLengthInData(attribute, stats, attributeLengthPosition);
    absolutePosition += position - mdBeginPosition;
}

template <class T>
void BP4Serializer::PutAttributeInDataCommon(const core::Attribute<T> &attribute,
                                             Stats<T> &stats) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    const size_t mdBeginPosition = position;

    // write a block identifier [AMD
    const char amd[] = "[AMD"; // no \0
    const size_t attributeLengthPosition =
        PutAttributeHeaderInData(attribute, stats, amd, sizeof(amd) - 1);

    uint8_t dataType = TypeTraits<T>::type_enum;
    helper::CopyToBuffer(buffer, position, &dataType);

    // here record payload offset
    stats.PayloadOffset = absolutePosition + position - mdBeginPosition + m_PreDataFileLength;

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

    // write a block identifier AMD]
    const char amdend[] = "AMD]"; // no \0
    helper::CopyToBuffer(buffer, position, amdend, sizeof(amdend) - 1);

    PutAttributeLengthInData(attribute, stats, attributeLengthPosition);
    absolutePosition += position - mdBeginPosition;
}

template <>
inline BP4Serializer::Stats<std::string>
BP4Serializer::GetBPStats(const bool /*singleValue*/,
                          const typename core::Variable<std::string>::BPInfo & /*blockInfo*/,
                          const bool /*isRowMajor*/) noexcept
{
    Stats<std::string> stats;
    stats.Step = m_MetadataSet.TimeStep;
    stats.FileIndex = GetFileIndex();
    return stats;
}

template <class T>
BP4Serializer::Stats<T>
BP4Serializer::GetBPStats(const bool singleValue,
                          const typename core::Variable<T>::BPInfo &blockInfo,
                          const bool isRowMajor) noexcept
{
    Stats<T> stats;
    stats.Step = m_MetadataSet.TimeStep;
    stats.FileIndex = GetFileIndex();

    // support span
    if (blockInfo.Data == nullptr && m_Parameters.StatsLevel > 0)
    {
        stats.Min = {};
        stats.Max = {};
        stats.SubBlockInfo = helper::DivideBlock(blockInfo.Count, m_Parameters.StatsBlockSize,
                                                 helper::BlockDivisionMethod::Contiguous);
        // set stats MinMaxs with the correct size
        helper::GetMinMaxSubblocks(blockInfo.Data, blockInfo.Count, stats.SubBlockInfo,
                                   stats.MinMaxs, stats.Min, stats.Max, m_Parameters.Threads,
                                   blockInfo.MemSpace);
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
            stats.SubBlockInfo = helper::DivideBlock(blockInfo.Count, m_Parameters.StatsBlockSize,
                                                     helper::BlockDivisionMethod::Contiguous);
            helper::GetMinMaxSubblocks(blockInfo.Data, blockInfo.Count, stats.SubBlockInfo,
                                       stats.MinMaxs, stats.Min, stats.Max, m_Parameters.Threads,
                                       blockInfo.MemSpace);
        }
        else
        {
            // non-contiguous memory min/max
            helper::GetMinMaxSelection(blockInfo.Data, blockInfo.MemoryCount, blockInfo.MemoryStart,
                                       blockInfo.Count, isRowMajor, stats.Min, stats.Max,
                                       blockInfo.MemSpace);
        }
        m_Profiler.Stop("minmax");
    }

    return stats;
}

template <class T>
size_t BP4Serializer::PutVariableMetadataInData(
    const core::Variable<T> &variable, const typename core::Variable<T>::BPInfo &blockInfo,
    const Stats<T> &stats, const typename core::Variable<T>::Span *span) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    const size_t mdBeginPosition = position;

    // write a block identifier [VMD
    const char vmd[] = "[VMD"; //  don't write \0!
    helper::CopyToBuffer(buffer, position, vmd, sizeof(vmd) - 1);

    // for writing length at the end
    const size_t varLengthPosition = position;
    position += 8; // skip var length (8)

    helper::CopyToBuffer(buffer, position, &stats.MemberID);

    PutNameRecord(variable.m_Name, buffer, position);

    // Layout can be 'K' = same as the language default, 'C' for row major and
    // 'F' for column major -- these are Numpy-like flags
    const char layout = 'K';
    helper::CopyToBuffer(buffer, position, &layout);

    // unused byte, write a 0 length to skip it
    const uint8_t zero8 = 0;
    helper::CopyToBuffer(buffer, position, &zero8);

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
    PutVariableCharacteristicsInData(variable, blockInfo, stats, buffer, position);

    // here align pointer for span
    const size_t padLengthPosition = position;
    constexpr std::array<uint8_t, 5> zeros = {0, 0, 0, 0, 0};
    // skip 1 for paddingLength and 4 for VMD] ending
    helper::CopyToBuffer(buffer, position, zeros.data(), 5);
    // here check for the next aligned pointer
    const size_t extraBytes = span == nullptr ? 0 : m_Data.Align<T>();
    const std::string pad = span == nullptr ? "VMD]" : std::string(extraBytes, '\0') + "VMD]";

    size_t backPosition = padLengthPosition;
    const uint8_t padLength = static_cast<uint8_t>(pad.size());
    helper::CopyToBuffer(buffer, backPosition, &padLength);
    helper::CopyToBuffer(buffer, backPosition, pad.c_str(), pad.size());

    position += extraBytes;

    absolutePosition += position - mdBeginPosition;
    return varLengthPosition;
}

template <>
inline size_t BP4Serializer::PutVariableMetadataInData(
    const core::Variable<std::string> &variable,
    const typename core::Variable<std::string>::BPInfo &blockInfo, const Stats<std::string> &stats,
    const typename core::Variable<std::string>::Span * /*span*/) noexcept
{
    auto &buffer = m_Data.m_Buffer;
    auto &position = m_Data.m_Position;
    auto &absolutePosition = m_Data.m_AbsolutePosition;

    const size_t mdBeginPosition = position;

    // write a block identifier [VMD
    const char vmd[] = "[VMD"; // no \0
    helper::CopyToBuffer(buffer, position, vmd, sizeof(vmd) - 1);

    // for writing length at the end
    const size_t varLengthPosition = position;
    position += 8; // skip var length (8)

    helper::CopyToBuffer(buffer, position, &stats.MemberID);

    PutNameRecord(variable.m_Name, buffer, position);

    // Layout can be 'K' = same as the language default, 'C' for row major and
    // 'F' for column major -- these are Numpy-like flags
    const char layout = 'K';
    helper::CopyToBuffer(buffer, position, &layout);

    // unused byte, write a 0 length to skip it
    const uint8_t zero8 = 0;
    helper::CopyToBuffer(buffer, position, &zero8);

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

    // write a block identifier VMD]
    // byte 4 and then four characters (as len+str without terminating 0)
    const char vmdend[] = "\4VMD]"; // no \0
    helper::CopyToBuffer(buffer, position, vmdend, sizeof(vmdend) - 1);

    absolutePosition += position - mdBeginPosition;
    return varLengthPosition;
}

template <class T>
void BP4Serializer::PutVariableMetadataInIndex(const core::Variable<T> &variable,
                                               const typename core::Variable<T>::BPInfo &blockInfo,
                                               const Stats<T> &stats, const bool isNew,
                                               SerialElementIndex &index,
                                               typename core::Variable<T>::Span *span) noexcept
{
    auto &buffer = index.Buffer;

    if (index.CurrentStep != stats.Step) // create a new variable header for a new step
    {
        size_t indexLengthPosition = buffer.size();
        index.CurrentHeaderPosition = buffer.size();

        buffer.insert(buffer.end(), 4, '\0'); // skip var length (4)
        helper::InsertToBuffer(buffer, &stats.MemberID);
        buffer.insert(buffer.end(), 2, '\0'); // skip group name
        PutNameRecord(variable.m_Name, buffer);

        // Layout can be 'K' = same as the language default, 'C' for row major
        // and 'F' for column major -- these are Numpy-like flags
        const char layout = 'K';
        buffer.insert(buffer.end(), 1, layout);

        // unused byte, write a 0 length to skip it
        buffer.insert(buffer.end(), 1, '\0'); // unused byte

        const uint8_t dataType = TypeTraits<T>::type_enum;
        helper::InsertToBuffer(buffer, &dataType);

        // Characteristics Sets Count in Metadata
        index.Count = 1;
        helper::InsertToBuffer(buffer, &index.Count);

        // For updating absolute offsets in agreggation
        index.LastUpdatedPosition = buffer.size();

        PutVariableCharacteristics(variable, blockInfo, stats, buffer, span);
        const uint32_t indexLength = static_cast<uint32_t>(buffer.size() - indexLengthPosition - 4);

        helper::CopyToBuffer(buffer, indexLengthPosition, &indexLength);

        index.CurrentStep = stats.Step;
    }
    else // update characteristics sets length and count
    {
        size_t currentIndexStartPosition = buffer.size();
        PutVariableCharacteristics(variable, blockInfo, stats, buffer, span);
        uint32_t currentIndexLength =
            static_cast<uint32_t>(buffer.size() - currentIndexStartPosition);

        size_t localPosition = index.CurrentHeaderPosition;
        uint32_t preIndexLength =
            helper::ReadValue<uint32_t>(buffer, localPosition, helper::IsLittleEndian());

        uint32_t newIndexLength = preIndexLength + currentIndexLength;

        localPosition = index.CurrentHeaderPosition; // back to beginning of the header
        helper::CopyToBuffer(buffer, localPosition, &newIndexLength);

        ++index.Count;
        // fixed since group and path are not printed
        size_t setsCountPosition = index.CurrentHeaderPosition + 15 + variable.m_Name.size();
        helper::CopyToBuffer(buffer, setsCountPosition, &index.Count);
    }
}

template <class T>
void BP4Serializer::PutBoundsRecord(const bool singleValue, const Stats<T> &stats,
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
            // Record entire Min-Max subblock arrays
            const uint8_t id = characteristic_minmax;
            uint16_t M = static_cast<uint16_t>(stats.MinMaxs.size() / 2);
            if (M == 0)
            {
                M = 1;
            }
            helper::InsertToBuffer(buffer, &id);
            helper::InsertToBuffer(buffer, &M);
            helper::InsertToBuffer(buffer, &stats.Min);
            helper::InsertToBuffer(buffer, &stats.Max);

            if (M > 1)
            {
                uint8_t method = static_cast<uint8_t>(stats.SubBlockInfo.DivisionMethod);
                helper::InsertToBuffer(buffer, &method);
                uint64_t subBlockSize = static_cast<uint64_t>(stats.SubBlockInfo.SubBlockSize);
                helper::InsertToBuffer(buffer, &subBlockSize);

                for (auto const d : stats.SubBlockInfo.Div)
                {
                    helper::InsertToBuffer(buffer, &d);
                }
                // insert min+max (alternating) elements (2*M values)
                for (auto const &m : stats.MinMaxs)
                {
                    helper::InsertToBuffer(buffer, &m);
                }
            }
            ++characteristicsCounter;
        }
    }
}

template <class T>
void BP4Serializer::PutBoundsRecord(const bool singleValue, const Stats<T> &stats,
                                    uint8_t &characteristicsCounter, std::vector<char> &buffer,
                                    size_t &position) noexcept
{
    if (singleValue)
    {
        PutCharacteristicRecord(characteristic_value, characteristicsCounter, stats.Min, buffer,
                                position);
    }
    else
    {
        if (m_Parameters.StatsLevel > 0) // default min and max only
        {
            // Record entire Min-Max subblock arrays
            const uint8_t id = characteristic_minmax;
            uint16_t M = static_cast<uint16_t>(stats.MinMaxs.size() / 2);
            if (M == 0)
            {
                M = 1;
            }
            helper::CopyToBuffer(buffer, position, &id);
            helper::CopyToBuffer(buffer, position, &M);
            helper::CopyToBuffer(buffer, position, &stats.Min);
            helper::CopyToBuffer(buffer, position, &stats.Max);

            if (M > 1)
            {
                uint8_t method = static_cast<uint8_t>(stats.SubBlockInfo.DivisionMethod);
                helper::CopyToBuffer(buffer, position, &method);
                uint64_t subBlockSize = static_cast<uint64_t>(stats.SubBlockInfo.SubBlockSize);
                helper::CopyToBuffer(buffer, position, &subBlockSize);

                for (auto const d : stats.SubBlockInfo.Div)
                {
                    helper::CopyToBuffer(buffer, position, &d);
                }
                // insert min+max (alternating) elements (2*M values)
                for (auto const &m : stats.MinMaxs)
                {
                    helper::CopyToBuffer(buffer, position, &m);
                }
            }
            ++characteristicsCounter;
        }
    }
}

template <>
inline void BP4Serializer::PutVariableCharacteristics(
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
void BP4Serializer::PutVariableCharacteristics(const core::Variable<T> &variable,
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

    uint8_t characteristicID = characteristic_dimensions;
    helper::InsertToBuffer(buffer, &characteristicID);
    const uint8_t dimensions = static_cast<uint8_t>(blockInfo.Count.size());
    helper::InsertToBuffer(buffer, &dimensions); // count
    const uint16_t dimensionsLength = static_cast<uint16_t>(24 * dimensions);
    helper::InsertToBuffer(buffer, &dimensionsLength); // length
    PutDimensionsRecord(blockInfo.Count, blockInfo.Shape, blockInfo.Start, buffer);
    ++characteristicsCounter;

    if (blockInfo.Data != nullptr || span != nullptr)
    {
        // minmax array depends on number of dimensions so this must
        // come after dimensions characteristics

        // TODO conflict between BP3 and BP4 with Span
        if (m_Parameters.StatsLevel > 0 && span != nullptr)
        {
            // store the characteristic position to be filled by PutSpanMetadata
            span->m_MinMaxMetadataPositions.first = buffer.size();
            span->m_MinMaxMetadataPositions.second = buffer.size();
        }

        PutBoundsRecord(variable.m_SingleValue, stats, characteristicsCounter, buffer);
    }

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
void BP4Serializer::PutVariableCharacteristicsInData(
    const core::Variable<T> &variable, const typename core::Variable<T>::BPInfo &blockInfo,
    const Stats<T> &stats, std::vector<char> &buffer, size_t &position) noexcept
{
    // going back at the end
    const size_t characteristicsCountPosition = position;
    // skip characteristics count(1) + length (4)
    position += 5;
    uint8_t characteristicsCounter = 0;

    // Write STAT min, max characteristics for an ARRAY
    // VALUE variable data is not written into characteristics
    // in the data file (only in metadata file in other function)
    if (blockInfo.Data != nullptr && !variable.m_SingleValue)
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

#endif // ADIOS2_TOOLKIT_FORMAT_BP4_BP4SERIALIZER_TCC_
