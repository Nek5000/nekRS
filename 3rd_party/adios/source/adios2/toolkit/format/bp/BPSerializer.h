/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPSerializer.h
 *
 *  Created on: Sep 16, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_H_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_H_

#include "adios2/toolkit/format/bp/BPBase.h"

#include <mutex>

namespace adios2
{
namespace format
{

class BPSerializer : virtual public BPBase
{
public:
    BPSerializer(helper::Comm const &comm, const uint8_t version);

    virtual ~BPSerializer() = default;

    /**
     * Serializes the data buffer and closes current process group
     * @param io : attributes written in first step
     * @param advanceStep true: advances step, false: doesn't advance
     */
    void SerializeData(core::IO &io, const bool advanceStep);

    /**
     * Get a string with profiling information for this rank
     * @param name stream name
     * @param transportsTypes list of transport types
     * @param transportsProfilers list of references to transport profilers
     */
    std::string
    GetRankProfilingJSON(const std::vector<std::string> &transportsTypes,
                         const std::vector<profiling::IOChrono *> &transportsProfilers) noexcept;

    /**
     * Forms the final profiling.json string aggregating from all ranks
     * @param rankProfilingJSON
     * @return profiling.json
     */
    std::vector<char> AggregateProfilingJSON(const std::string &rankLog) const;

    void UpdateOffsetsInMetadata();

protected:
    /** BP format version */
    const uint8_t m_Version;

    static std::mutex m_Mutex;

    virtual void SerializeDataBuffer(core::IO &io) noexcept = 0;

    template <class T>
    void PutAttributeCharacteristicValueInIndex(uint8_t &characteristicsCounter,
                                                const core::Attribute<T> &attribute,
                                                std::vector<char> &buffer) noexcept;

    template <class T>
    void PutCharacteristicRecord(const uint8_t characteristicID, uint8_t &characteristicsCounter,
                                 const T &value, std::vector<char> &buffer) noexcept;

    template <class T>
    void PutCharacteristicRecord(const uint8_t characteristicID, uint8_t &characteristicsCounter,
                                 const T &value, std::vector<char> &buffer,
                                 size_t &position) noexcept;

    template <class T>
    void PutPayloadInBuffer(const core::Variable<T> &variable,
                            const typename core::Variable<T>::BPInfo &blockInfo,
                            const bool sourceRowMajor) noexcept;

    void PutNameRecord(const std::string name, std::vector<char> &buffer) noexcept;

    void PutNameRecord(const std::string name, std::vector<char> &buffer,
                       size_t &position) noexcept;

    void PutDimensionsRecord(const Dims &localDimensions, const Dims &globalDimensions,
                             const Dims &offsets, std::vector<char> &buffer) noexcept;

    void PutDimensionsRecord(const Dims &localDimensions, const Dims &globalDimensions,
                             const Dims &offsets, std::vector<char> &buffer, size_t &position,
                             const bool isCharacteristic = false) noexcept;

    void PutMinifooter(const uint64_t pgIndexStart, const uint64_t variablesIndexStart,
                       const uint64_t attributesIndexStart, std::vector<char> &buffer,
                       size_t &position, const bool addSubfiles = false);

    void MergeSerializeIndices(
        const std::unordered_map<std::string, std::vector<SerialElementIndex>> &nameRankIndices,
        helper::Comm const &comm, BufferSTL &bufferSTL);

    template <class T>
    void UpdateIndexOffsetsCharacteristics(size_t &currentPosition, const DataTypes dataType,
                                           std::vector<char> &buffer);

    uint32_t GetFileIndex() const noexcept;

    size_t GetAttributesSizeInData(core::IO &io) const noexcept;

    template <class T>
    size_t GetAttributeSizeInData(const core::Attribute<T> &attribute) const noexcept;

    void PutAttributes(core::IO &io);

    SerialElementIndex &
    GetSerialElementIndex(const std::string &name,
                          std::unordered_map<std::string, SerialElementIndex> &indices,
                          bool &isNew) const noexcept;

    template <class T>
    void PutAttributeInData(const core::Attribute<T> &attribute, Stats<T> &stats) noexcept;
    template <class T>
    void PutAttributeInIndex(const core::Attribute<T> &attribute, const Stats<T> &stats) noexcept;

#define declare_template_instantiation(T)                                                          \
    virtual void DoPutAttributeInData(const core::Attribute<T> &attribute,                         \
                                      Stats<T> &stats) noexcept = 0;

    ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

    // Operations related functions
    template <class T>
    void PutCharacteristicOperation(const core::Variable<T> &variable,
                                    const typename core::Variable<T>::BPInfo &blockInfo,
                                    std::vector<char> &buffer) noexcept;

    template <class T>
    void PutOperationPayloadInBuffer(const core::Variable<T> &variable,
                                     const typename core::Variable<T>::BPInfo &blockInfo);

private:
    size_t m_OutputSizeMetadataPosition;
};

} // end namespace format
} // end namespace adios2

#include "BPSerializer.inl"

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BPSERIALIZER_H_ */
