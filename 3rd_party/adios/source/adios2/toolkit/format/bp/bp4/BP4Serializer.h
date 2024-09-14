/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Serializer.h
 *
 *  Created on: Aug 1, 2018
 *  Author: William F Godoy godoywf@ornl.gov
 *          Lipeng Wan wanl@ornl.gov
 *          Norbert Podhorszki pnb@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BP4_BP4SERIALIZER_H_
#define ADIOS2_TOOLKIT_FORMAT_BP_BP4_BP4SERIALIZER_H_

#include "BP4Base.h"

#include "adios2/core/Attribute.h"
#include "adios2/core/IO.h"

#include "adios2/toolkit/format/bp/BPSerializer.h"

#ifdef _WIN32
#pragma warning(disable : 4250)
#endif

namespace adios2
{
namespace format
{

class BP4Serializer : public BP4Base, public BPSerializer
{

public:
    /**
     * Unique constructor
     * @param comm multi-process communicator for BP1 Aggregator
     */
    BP4Serializer(helper::Comm const &comm);

    ~BP4Serializer() = default;

    /** Writes a 64 byte header into the data/metadata buffer.
     *  Must be called only when the buffer is empty.
     *  @param buffer the data or metadata buffer
     *  @param fileType a small string up to 8 characters that is
     *  concatenated to the version string
     */
    void MakeHeader(BufferSTL &b, const std::string fileType, const bool isActive);

    /**
     * Writes a process group index PGIndex and list of methods (from
     * transports). Done at Open or Advance.
     * @param ioName from IO class, identify Process Group with IO name
     * @param hostLanguage from ADIOS class passed to IO
     * @param transportsTypes passed to get list of transport "bp methods"
     */
    void PutProcessGroupIndex(const std::string &ioName, const std::string hostLanguage,
                              const std::vector<std::string> &transportsTypes) noexcept;

    /**
     * Put in buffer metadata for a given variable
     * @param variable
     */
    template <class T>
    void PutVariableMetadata(const core::Variable<T> &variable,
                             const typename core::Variable<T>::BPInfo &blockInfo,
                             const bool sourceRowMajor = true,
                             typename core::Variable<T>::Span *span = nullptr) noexcept;

    /**
     * Put in buffer variable payload. Expensive part.
     * @param variable payload input from m_PutValues
     */
    template <class T>
    void PutVariablePayload(const core::Variable<T> &variable,
                            const typename core::Variable<T>::BPInfo &blockInfo,
                            const bool sourceRowMajor = true,
                            typename core::Variable<T>::Span *span = nullptr) noexcept;

    template <class T>
    void PutSpanMetadata(const core::Variable<T> &variable,
                         const typename core::Variable<T>::BPInfo &blockInfo,
                         const typename core::Variable<T>::Span &span) noexcept;

    /**
     * Serializes the metadata indices appending it into the data buffer inside
     * m_Data
     * @param updateAbsolutePosition true: adds footer size to absolute position
     * used for Close,
     * false: doesn't update absolute, used for partial buffer
     * @param inData true: serialize in data , false: only update metadata
     * indices, do not serialize in data
     */
    void SerializeMetadataInData(const bool updateAbsolutePosition = true,
                                 const bool inData = true);

    /**
     * Finishes bp buffer by serializing data and adding trailing metadata
     * @param io
     * @return the position of the data buffer before writing the metadata
     * footer into it (in case someone wants to write only the data portion)
     */
    size_t CloseData(core::IO &io);

    /**
     * Closes bp buffer for streaming mode...must reset metadata for the next
     * step
     * @param io object containing all attributes
     * @param addMetadata true: process metadata and add to data buffer
     * (minifooter)
     * @return the position of the data buffer before writing the metadata
     * footer into it (in case someone wants to write only the data portion)
     */
    size_t CloseStream(core::IO &io, const bool addMetadata = true);

    /* Reset all metadata indices at the end of each step */
    void ResetAllIndices();

    /* Reset metadata index table*/
    void ResetMetadataIndexTable();

    /**
     * Aggregate collective metadata
     * @param comm input establishing domain (all or per aggregator)
     * @param bufferSTL buffer to put the metadata
     * @param inMetadataBuffer collective metadata from absolute rank = 0, else
     *                         from aggregators
     */
    void AggregateCollectiveMetadata(helper::Comm const &comm, BufferSTL &bufferSTL,
                                     const bool inMetadataBuffer);

private:
    std::vector<char> m_SerializedIndices;
    std::vector<char> m_GatheredSerializedIndices;

    /** aggregate pg rank indices */
    std::unordered_map<size_t, std::vector<std::tuple<size_t, size_t, size_t>>> m_PGIndicesInfo;
    /** deserialized variable indices per rank (vector index) */
    std::unordered_map<size_t,
                       std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>>
        m_VariableIndicesInfo;
    /** deserialized attribute indices per rank (vector index) */
    std::unordered_map<size_t,
                       std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t>>>>
        m_AttributesIndicesInfo;

    /**
     * Put in BP buffer attribute header, called from PutAttributeInData
     * specialized functions
     * @param attribute input
     * @param stats BP4 Stats
     * @param headerID  short string to identify block ("[AMD")
     * @return attribute length position
     */
    template <class T>
    size_t PutAttributeHeaderInData(const core::Attribute<T> &attribute, Stats<T> &stats,
                                    const char *headerID, const size_t headerIDLength) noexcept;

    /**
     * Called from WriteAttributeInData specialized functions
     * @param attribute input
     * @param stats BP4 stats
     * @param attributeLengthPosition
     */
    template <class T>
    void PutAttributeLengthInData(const core::Attribute<T> &attribute, Stats<T> &stats,
                                  const size_t attributeLengthPosition) noexcept;

#define declare_template_instantiation(T)                                                          \
    void DoPutAttributeInData(const core::Attribute<T> &attribute, Stats<T> &stats) noexcept final;
    ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

    template <class T>
    void PutAttributeInDataCommon(const core::Attribute<T> &attribute, Stats<T> &stats) noexcept;

    /**
     * Get variable statistics
     * @param variable
     * @param isRowMajor
     * @return stats BP4 Stats
     */
    template <class T>
    Stats<T> GetBPStats(const bool singleValue, const typename core::Variable<T>::BPInfo &blockInfo,
                        const bool isRowMajor) noexcept;

    /** @return The position that holds the length of the variable entry
     * (metadat+data length). The actual lengths is know after
     * PutVariablePayload()
     */
    template <class T>
    size_t PutVariableMetadataInData(const core::Variable<T> &variable,
                                     const typename core::Variable<T>::BPInfo &blockInfo,
                                     const Stats<T> &stats,
                                     const typename core::Variable<T>::Span *span) noexcept;

    template <class T>
    void PutVariableMetadataInIndex(const core::Variable<T> &variable,
                                    const typename core::Variable<T>::BPInfo &blockInfo,
                                    const Stats<T> &stats, const bool isNew,
                                    SerialElementIndex &index,
                                    typename core::Variable<T>::Span *span) noexcept;

    template <class T>
    void PutVariableCharacteristics(const core::Variable<T> &variable,
                                    const typename core::Variable<T>::BPInfo &blockInfo,
                                    const Stats<T> &stats, std::vector<char> &buffer,
                                    typename core::Variable<T>::Span *span) noexcept;

    template <class T>
    void PutVariableCharacteristicsInData(const core::Variable<T> &variable,
                                          const typename core::Variable<T>::BPInfo &blockInfo,
                                          const Stats<T> &stats, std::vector<char> &buffer,
                                          size_t &position) noexcept;

    /** Writes min max */
    template <class T>
    void PutBoundsRecord(const bool singleValue, const Stats<T> &stats,
                         uint8_t &characteristicsCounter, std::vector<char> &buffer) noexcept;

    /** Overloaded version for data buffer */
    template <class T>
    void PutBoundsRecord(const bool singleValue, const Stats<T> &stats,
                         uint8_t &characteristicsCounter, std::vector<char> &buffer,
                         size_t &position) noexcept;

    /**
     * Wraps up the data buffer serialization in m_HeapBuffer and fills the pg
     * length, vars count, vars
     * length and attributes count and attributes length
     * @param io object containing all attributes
     */
    void SerializeDataBuffer(core::IO &io) noexcept final;

    void AggregateCollectiveMetadataIndices(helper::Comm const &comm, BufferSTL &bufferSTL);
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_UTILITIES_FORMAT_BP4_BP4Serializer_H_ */
