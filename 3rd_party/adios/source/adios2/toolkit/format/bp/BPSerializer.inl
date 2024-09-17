/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPSerializer.inl
 *
 *  Created on: Sep 16, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_INL_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_INL_
#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_H_
#error "Inline file should only be included from its header, never on its own"
#endif

namespace adios2
{
namespace format
{

template <>
inline void BPSerializer::PutAttributeCharacteristicValueInIndex(
    uint8_t &characteristicsCounter,
    const core::Attribute<std::string> &attribute,
    std::vector<char> &buffer) noexcept
{
    const uint8_t characteristicID =
        static_cast<uint8_t>(CharacteristicID::characteristic_value);

    helper::InsertToBuffer(buffer, &characteristicID);

    if (attribute.m_IsSingleValue) // Single string
    {
        const uint16_t dataSize =
            static_cast<uint16_t>(attribute.m_DataSingleValue.size());
        helper::InsertToBuffer(buffer, &dataSize);
        helper::InsertToBuffer(buffer, attribute.m_DataSingleValue.data(),
                               attribute.m_DataSingleValue.size());
    }
    else // string array
    {
        for (size_t s = 0; s < attribute.m_Elements; ++s)
        {
            // without zero terminated character
            const std::string element(attribute.m_DataArray[s]);

            const uint16_t elementSize = static_cast<uint16_t>(element.size());

            helper::InsertToBuffer(buffer, &elementSize);
            helper::InsertToBuffer(buffer, element.data(), element.size());
        }
    }
    ++characteristicsCounter;
}

template <>
inline size_t BPSerializer::GetAttributeSizeInData(
    const core::Attribute<std::string> &attribute) const noexcept
{
    // index header
    size_t size = 14 + attribute.m_Name.size() + 10;

    if (attribute.m_IsSingleValue)
    {
        size += 4 + attribute.m_DataSingleValue.size();
    }
    else
    {
        size += 4;
        for (const std::string &dataString : attribute.m_DataArray)
        {
            size += 4 + dataString.size();
        }
    }
    return size;
}

template <>
inline void BPSerializer::PutPayloadInBuffer(
    const core::Variable<std::string> &variable,
    const typename core::Variable<std::string>::BPInfo &blockInfo,
    const bool /* sourceRowMajor*/) noexcept
{
    PutNameRecord(*blockInfo.Data, m_Data.m_Buffer, m_Data.m_Position);
    m_Data.m_AbsolutePosition += blockInfo.Data->size() + 2;
}

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_INL_ */
