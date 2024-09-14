/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Deserializer.h
 *
 *  Created on: Aug 1, 2018
 *  Author: William F Godoy godoywf@ornl.gov
 *          Lipeng Wan wanl@ornl.gov
 *          Norbert Podhorszki pnb@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BP4_BP4DESERIALIZER_H_
#define ADIOS2_TOOLKIT_FORMAT_BP_BP4_BP4DESERIALIZER_H_

#include "BP4Base.h"

#include <mutex>
#include <set>
#include <utility> //std::pair
#include <vector>

#include "adios2/core/IO.h"
#include "adios2/core/Variable.h"
#include "adios2/helper/adiosFunctions.h" //VariablesSubFileInfo, BlockOperation

namespace adios2
{
namespace format
{

class BP4Deserializer : virtual public BP4Base
{

public:
    /** BP Minifooter fields */
    Minifooter m_Minifooter;

    bool m_WriterIsActive = false;

    /**
     * Unique constructor
     * @param comm multi-process communicator
     */
    BP4Deserializer(helper::Comm const &comm);

    ~BP4Deserializer() = default;

    void ParseMetadataIndex(BufferSTL &bufferSTL, const size_t absoluteStartPos,
                            const bool hasHeader, const bool oneStepOnly);

    /* Return the position in the buffer where processing ends. The processing
     * is controlled by the number of records in the Index, which may be less
     * than the actual entries in the metadata in a streaming situation (where
     * writer has just written metadata for step K+1,...,K+L while the index
     * contains K steps when the reader looks at it).
     */
    size_t ParseMetadata(const BufferSTL &bufferSTL, core::Engine &engine,
                         const bool firstStep = true);

    /**
     * Used to get the variable payload data for the current selection (dims and
     * steps), used in single buffer for streaming
     * @param variable
     * @param bufferSTL bp buffer input that contains metadata and data
     */
    template <class T>
    void GetSyncVariableDataFromStream(core::Variable<T> &variable, BufferSTL &bufferSTL) const;

    /**
     * Initializes a block inside variable.m_BlocksInfo
     * @param variable input
     * @param data user data pointer
     * @return a reference inside variable.m_BlocksInfo (invalidated if called
     * twice)
     */
    template <class T>
    typename core::Variable<T>::BPInfo &InitVariableBlockInfo(core::Variable<T> &variable,
                                                              T *data) const;

    /**
     * Sets read block information from the available metadata information
     * @param variable
     * @param blockInfo
     */
    template <class T>
    void SetVariableBlockInfo(core::Variable<T> &variable,
                              typename core::Variable<T>::BPInfo &blockInfo) const;

    /**
     * Prepares the information to get raw data from the transport manager for a
     * required substream box (block)
     * @param variable input Variable
     * @param blockInfo input blockInfo with information about Get request
     * @param subStreamBoxInfo contains information (e.g. bounds, operation,
     * etc.) about the available box (block) to be accessed by the Transport
     * Manager.
     * @param buffer output to be passed to Transport Manager for current box
     * @param payloadSize output to be passed to Transport Manager for current
     * box
     * @param payloadStart output to be passed to Transport Manager for current
     * box
     * @param threadID assign different thread ID to have independent raw memory
     * spaces per thread, default = 0
     */
    template <class T>
    void PreDataRead(core::Variable<T> &variable, typename core::Variable<T>::BPInfo &blockInfo,
                     const helper::SubStreamBoxInfo &subStreamBoxInfo, char *&buffer,
                     size_t &payloadSize, size_t &payloadOffset, const size_t threadID = 0);

    template <class T>
    void PostDataRead(core::Variable<T> &variable, typename core::Variable<T>::BPInfo &blockInfo,
                      const helper::SubStreamBoxInfo &subStreamBoxInfo,
                      const bool isRowMajorDestination, const size_t threadID = 0);

    void BackCompatDecompress(const helper::SubStreamBoxInfo &subStreamBoxInfo,
                              const size_t threadID = 0);

    /**
     * Clips and assigns memory to blockInfo.Data from a contiguous memory
     * input
     * @param blockInfo
     * @param contiguousMemory
     * @param blockBox
     * @param intersectionBox
     */
    template <class T>
    void ClipContiguousMemory(typename core::Variable<T>::BPInfo &blockInfo,
                              const std::vector<char> &contiguousMemory, const Box<Dims> &blockBox,
                              const Box<Dims> &intersectionBox) const;

    /**
     * Gets a value directly from metadata (if Variable is single value)
     * @param variable
     * @param data
     */
    template <class T>
    void GetValueFromMetadata(core::Variable<T> &variable, T *data) const;

    template <class T>
    std::map<size_t, std::vector<typename core::Variable<T>::BPInfo>>
    AllStepsBlocksInfo(const core::Variable<T> &variable) const;

    template <class T>
    std::vector<std::vector<typename core::Variable<T>::BPInfo>>
    AllRelativeStepsBlocksInfo(const core::Variable<T> &variable) const;

    template <class T>
    std::vector<typename core::Variable<T>::BPInfo> BlocksInfo(const core::Variable<T> &variable,
                                                               const size_t step) const;

    // TODO : Will deprecate all function below
    std::map<std::string, helper::SubFileInfoMap> PerformGetsVariablesSubFileInfo(core::IO &io);

    // TODO : will deprecate
    template <class T>
    std::map<std::string, helper::SubFileInfoMap>
    GetSyncVariableSubFileInfo(const core::Variable<T> &variable) const;

    // TODO : will deprecate
    template <class T>
    void GetDeferredVariable(core::Variable<T> &variable, T *data);

    // TODO : will deprecate
    template <class T>
    helper::SubFileInfoMap GetSubFileInfo(const core::Variable<T> &variable) const;

    // TODO : will deprecate
    void ClipMemory(const std::string &variableName, core::IO &io,
                    const std::vector<char> &contiguousMemory, const Box<Dims> &blockBox,
                    const Box<Dims> &intersectionBox) const;

    bool ReadActiveFlag(std::vector<char> &buffer);

    // TODO: will deprecate
    bool m_PerformedGets = false;

private:
    std::map<std::string, helper::SubFileInfoMap> m_DeferredVariablesMap;

    static std::mutex m_Mutex;

    void ParseMinifooter(const BufferSTL &bufferSTL);

    // void ParsePGIndex(const BufferSTL &bufferSTL, const core::IO &io);
    void ParsePGIndexPerStep(const BufferSTL &bufferSTL, const std::string hostLanguage,
                             size_t submetadatafileId, size_t step);

    // void ParseVariablesIndex(const BufferSTL &bufferSTL, core::IO &io);
    void ParseVariablesIndexPerStep(const BufferSTL &bufferSTL, core::Engine &engine,
                                    size_t submetadatafileId, size_t step);

    // void ParseAttributesIndex(const BufferSTL &bufferSTL, core::IO &io);
    void ParseAttributesIndexPerStep(const BufferSTL &bufferSTL, core::Engine &engine,
                                     size_t submetadatafileId, size_t step);

    /**
     * Reads a variable index element (serialized) and calls IO.DefineVariable
     * to deserialize the Variable metadata
     * @param header serialize
     * @param io
     * @param buffer
     * @param position
     */

    template <class T>
    void DefineVariableInEngineIOPerStep(const ElementIndexHeader &header, core::Engine &engine,
                                         const std::vector<char> &buffer, size_t position,
                                         size_t step) const;

    template <class T>
    void DefineAttributeInEngineIO(const ElementIndexHeader &header, core::Engine &engine,
                                   const std::vector<char> &buffer, size_t position) const;

    template <class T>
    std::vector<typename core::Variable<T>::BPInfo>
    BlocksInfoCommon(const core::Variable<T> &variable,
                     const std::vector<size_t> &blocksIndexOffsets) const;

    const helper::BlockOperationInfo &InitPostOperatorBlockData(
        const std::vector<helper::BlockOperationInfo> &blockOperationsInfo) const;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BP4_BP4DESERIALIZER_H_ */
