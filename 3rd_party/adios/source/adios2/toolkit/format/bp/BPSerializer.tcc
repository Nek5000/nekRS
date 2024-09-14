/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPSerializer.tcc
 *
 *  Created on: Sep 16, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_TCC_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_TCC_

#include "BPSerializer.h"

namespace adios2
{
namespace format
{

template <class T>
inline void
BPSerializer::PutAttributeCharacteristicValueInIndex(uint8_t &characteristicsCounter,
                                                     const core::Attribute<T> &attribute,
                                                     std::vector<char> &buffer) noexcept
{
    const uint8_t characteristicID = CharacteristicID::characteristic_value;

    helper::InsertToBuffer(buffer, &characteristicID);

    if (attribute.m_IsSingleValue) // single value
    {
        helper::InsertToBuffer(buffer, &attribute.m_DataSingleValue);
    }
    else // array
    {
        helper::InsertToBuffer(buffer, attribute.m_DataArray.data(), attribute.m_Elements);
    }
    ++characteristicsCounter;
}

template <class T>
void BPSerializer::PutCharacteristicRecord(const uint8_t characteristicID,
                                           uint8_t &characteristicsCounter, const T &value,
                                           std::vector<char> &buffer) noexcept
{
    const uint8_t id = characteristicID;
    helper::InsertToBuffer(buffer, &id);
    helper::InsertToBuffer(buffer, &value);
    ++characteristicsCounter;
}

template <class T>
void BPSerializer::PutCharacteristicRecord(const uint8_t characteristicID,
                                           uint8_t &characteristicsCounter, const T &value,
                                           std::vector<char> &buffer, size_t &position) noexcept
{
    const uint8_t id = characteristicID;
    helper::CopyToBuffer(buffer, position, &id);
    helper::CopyToBuffer(buffer, position, &value);
    ++characteristicsCounter;
}

template <class T>
inline void BPSerializer::PutPayloadInBuffer(const core::Variable<T> &variable,
                                             const typename core::Variable<T>::BPInfo &blockInfo,
                                             const bool sourceRowMajor) noexcept
{
    const size_t blockSize = helper::GetTotalSize(blockInfo.Count);
    m_Profiler.Start("memcpy");

#ifdef ADIOS2_HAVE_GPU_SUPPORT
    if (blockInfo.MemSpace == MemorySpace::GPU)
    {
        helper::CopyFromGPUToBuffer(m_Data.m_Buffer, m_Data.m_Position, blockInfo.Data,
                                    blockInfo.MemSpace, blockSize);
        m_Profiler.Stop("memcpy");
        m_Data.m_AbsolutePosition += blockSize * sizeof(T);
        return;
    }
#endif

    if (!blockInfo.MemoryStart.empty())
    {
        helper::CopyMemoryBlock(reinterpret_cast<T *>(m_Data.m_Buffer.data() + m_Data.m_Position),
                                blockInfo.Start, blockInfo.Count, sourceRowMajor, blockInfo.Data,
                                blockInfo.Start, blockInfo.Count, sourceRowMajor, false, Dims(),
                                Dims(), blockInfo.MemoryStart, blockInfo.MemoryCount);
        m_Data.m_Position += blockSize * sizeof(T);
    }
    else
    {
        helper::CopyToBufferThreads(m_Data.m_Buffer, m_Data.m_Position, blockInfo.Data, blockSize,
                                    m_Parameters.Threads);
    }
    m_Profiler.Stop("memcpy");
    m_Data.m_AbsolutePosition += blockSize * sizeof(T); // payload size
}

// PRIVATE
template <class T>
void BPSerializer::UpdateIndexOffsetsCharacteristics(size_t &currentPosition,
                                                     const DataTypes dataType,
                                                     std::vector<char> &buffer)
{
    const bool isLittleEndian = helper::IsLittleEndian();

    helper::ReadValue<uint8_t>(buffer, currentPosition, isLittleEndian);

    const uint32_t characteristicsLength =
        helper::ReadValue<uint32_t>(buffer, currentPosition, isLittleEndian);

    const size_t endPosition = currentPosition + static_cast<size_t>(characteristicsLength);

    size_t dimensionsSize = 0; // get it from dimensions characteristics

    while (currentPosition < endPosition)
    {
        const uint8_t id = helper::ReadValue<uint8_t>(buffer, currentPosition, isLittleEndian);

        switch (id)
        {
        case (characteristic_time_index): {
            currentPosition += sizeof(uint32_t);
            break;
        }

        case (characteristic_file_index): {
            currentPosition += sizeof(uint32_t);
            break;
        }

        case (characteristic_value): {
            if (dataType == type_string)
            {
                // first get the length of the string
                const size_t length = static_cast<size_t>(

                    helper::ReadValue<uint16_t>(buffer, currentPosition, isLittleEndian));

                currentPosition += length;
            }
            // using this function only for variables
            // TODO string array if string arrays are supported in the future
            else
            {
                currentPosition += sizeof(T);
            }

            break;
        }
        case (characteristic_min): {
            currentPosition += sizeof(T);
            break;
        }
        case (characteristic_max): {
            currentPosition += sizeof(T);
            break;
        }
        case (characteristic_minmax): {
            // first get the number of subblocks
            const uint16_t M = helper::ReadValue<uint16_t>(buffer, currentPosition);
            currentPosition += 2 * sizeof(T); // block min/max
            if (M > 1)
            {
                currentPosition += 1 + 8; // method (byte), blockSize (uint64_t)
                currentPosition += dimensionsSize * sizeof(uint16_t); // N-dim division
                currentPosition += 2 * M * sizeof(T);                 // M * min/max
            }
            break;
        }
        case (characteristic_offset): {
            const uint64_t currentOffset =
                helper::ReadValue<uint64_t>(buffer, currentPosition, isLittleEndian);

            const uint64_t updatedOffset =
                currentOffset + static_cast<uint64_t>(m_Data.m_AbsolutePosition);

            currentPosition -= sizeof(uint64_t);
            helper::CopyToBuffer(buffer, currentPosition, &updatedOffset);
            break;
        }
        case (characteristic_payload_offset): {
            const uint64_t currentPayloadOffset =
                helper::ReadValue<uint64_t>(buffer, currentPosition, isLittleEndian);

            const uint64_t updatedPayloadOffset =
                currentPayloadOffset + static_cast<uint64_t>(m_Data.m_AbsolutePosition);

            currentPosition -= sizeof(uint64_t);
            helper::CopyToBuffer(buffer, currentPosition, &updatedPayloadOffset);
            break;
        }
        case (characteristic_dimensions): {
            dimensionsSize = static_cast<size_t>(
                helper::ReadValue<uint8_t>(buffer, currentPosition, isLittleEndian));

            currentPosition += 3 * sizeof(uint64_t) * dimensionsSize + 2; // 2 is for length
            break;
        }
        case (characteristic_transform_type): {
            const size_t typeLength = static_cast<size_t>(
                helper::ReadValue<uint8_t>(buffer, currentPosition, isLittleEndian));
            // skip over operator name (transform type) string
            currentPosition += typeLength;

            // skip over pre-data type (1) and dimensionsSize (1)
            currentPosition += 2;

            const uint16_t dimensionsLength =
                helper::ReadValue<uint16_t>(buffer, currentPosition, isLittleEndian);
            // skip over dimensions
            currentPosition += dimensionsLength;

            const size_t metadataLength = static_cast<size_t>(
                helper::ReadValue<uint16_t>(buffer, currentPosition, isLittleEndian));
            // skip over operator metadata
            currentPosition += metadataLength;

            break;
        }
        default: {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BPSerializer", "UpdateIndexOffsetsCharacteristics",
                "characteristic ID " + std::to_string(id) + " not supported when updating offsets");
        }

        } // end id switch
    }     // end while
}

template <class T>
inline size_t
BPSerializer::GetAttributeSizeInData(const core::Attribute<T> &attribute) const noexcept
{
    size_t size = 14 + attribute.m_Name.size() + 10;
    size += 4 + sizeof(T) * attribute.m_Elements;
    return size;
}

template <class T>
void BPSerializer::PutAttributeInData(const core::Attribute<T> &attribute, Stats<T> &stats) noexcept
{
    DoPutAttributeInData(attribute, stats);
}

template <class T>
void BPSerializer::PutAttributeInIndex(const core::Attribute<T> &attribute,
                                       const Stats<T> &stats) noexcept
{
    SerialElementIndex index(stats.MemberID);
    auto &buffer = index.Buffer;

    // index.Valid = true; // when the attribute is put, set this flag to true
    size_t indexLengthPosition = buffer.size();

    buffer.insert(buffer.end(), 4, '\0'); // skip attribute length (4)
    helper::InsertToBuffer(buffer, &stats.MemberID);
    buffer.insert(buffer.end(), 2, '\0'); // skip group name
    PutNameRecord(attribute.m_Name, buffer);
    buffer.insert(buffer.end(), 2, '\0'); // skip path

    uint8_t dataType = TypeTraits<T>::type_enum; // dataType

    if (dataType == type_string && !attribute.m_IsSingleValue)
    {
        dataType = type_string_array;
    }

    helper::InsertToBuffer(buffer, &dataType);

    // Characteristics Sets Count in Metadata
    index.Count = 1;
    helper::InsertToBuffer(buffer, &index.Count);

    // START OF CHARACTERISTICS
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
    constexpr uint8_t dimensions = 1;
    helper::InsertToBuffer(buffer, &dimensions); // count
    constexpr uint16_t dimensionsLength = 24;
    helper::InsertToBuffer(buffer, &dimensionsLength); // length
    PutDimensionsRecord({attribute.m_Elements}, {}, {}, buffer);
    ++characteristicsCounter;

    // VALUE
    PutAttributeCharacteristicValueInIndex(characteristicsCounter, attribute, buffer);

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

    // Remember this attribute and its serialized piece
    // should not affect BP3 as it's recalculated
    const uint32_t indexLength = static_cast<uint32_t>(buffer.size() - indexLengthPosition - 4);

    helper::CopyToBuffer(buffer, indexLengthPosition, &indexLength);
    m_MetadataSet.AttributesIndices.emplace(attribute.m_Name, index);
    m_SerializedAttributes.emplace(attribute.m_Name);
}

// operations related functions
template <class T>
void BPSerializer::PutCharacteristicOperation(const core::Variable<T> &variable,
                                              const typename core::Variable<T>::BPInfo &blockInfo,
                                              std::vector<char> &buffer) noexcept
{
    const std::string type = blockInfo.Operations[0]->m_TypeString;
    const uint8_t typeLength = static_cast<uint8_t>(type.size());
    helper::InsertToBuffer(buffer, &typeLength);
    helper::InsertToBuffer(buffer, type.c_str(), type.size());

    // pre-transform type
    const uint8_t dataType = TypeTraits<T>::type_enum;
    helper::InsertToBuffer(buffer, &dataType);
    // pre-transform dimensions
    const uint8_t dimensions = static_cast<uint8_t>(blockInfo.Count.size());
    helper::InsertToBuffer(buffer, &dimensions); // count
    const uint16_t dimensionsLength = static_cast<uint16_t>(24 * dimensions);
    helper::InsertToBuffer(buffer, &dimensionsLength); // length
    PutDimensionsRecord(blockInfo.Count, blockInfo.Shape, blockInfo.Start, buffer);

    // here put the metadata info depending on operation
    const uint64_t inputSize =
        static_cast<uint64_t>(helper::GetTotalSize(blockInfo.Count) * sizeof(T));

    // fixed size only stores inputSize 8-bytes and outputSize 8-bytes
    constexpr uint16_t metadataSize = 16;
    helper::InsertToBuffer(buffer, &metadataSize);
    helper::InsertToBuffer(buffer, &inputSize);

    m_OutputSizeMetadataPosition = buffer.size();

    constexpr uint64_t outputSize = 0;
    helper::InsertToBuffer(buffer, &outputSize);
}

template <class T>
void BPSerializer::PutOperationPayloadInBuffer(const core::Variable<T> &variable,
                                               const typename core::Variable<T>::BPInfo &blockInfo)
{
    size_t outputSize = blockInfo.Operations[0]->Operate(
        reinterpret_cast<char *>(blockInfo.Data), blockInfo.Start, blockInfo.Count, variable.m_Type,
        m_Data.m_Buffer.data() + m_Data.m_Position);

    if (outputSize == 0) // the operator was not applied
        outputSize = helper::CopyMemoryWithOpHeader(
            reinterpret_cast<char *>(blockInfo.Data), blockInfo.Count, variable.m_Type,
            m_Data.m_Buffer.data() + m_Data.m_Position, blockInfo.Operations[0]->GetHeaderSize(),
            blockInfo.MemSpace);

    m_Data.m_Position += outputSize;
    m_Data.m_AbsolutePosition += outputSize;

    // update metadata
    bool isFound = false;
    SerialElementIndex &variableIndex =
        GetSerialElementIndex(variable.m_Name, m_MetadataSet.VarsIndices, isFound);

    size_t backPosition = m_OutputSizeMetadataPosition;

    helper::CopyToBuffer(variableIndex.Buffer, backPosition, &outputSize);
}

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_TCC_ */
