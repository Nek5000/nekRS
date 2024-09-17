/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Deserializer.tcc
 *
 *  Created on: Aug 1, 2018
 *  Author: William F Godoy godoywf@ornl.gov
 *          Lipeng Wan wanl@ornl.gov
 *          Norbert Podhorszki pnb@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP1_BP4DESERIALIZER_TCC_
#define ADIOS2_TOOLKIT_FORMAT_BP1_BP4DESERIALIZER_TCC_

#include "BP4Deserializer.h"

#include <algorithm> //std::reverse
#include <iostream>
#include <unordered_set>

#include "adios2/helper/adiosFunctions.h"
#include "adios2/operator/OperatorFactory.h"

namespace adios2
{
namespace format
{

template <class T>
void BP4Deserializer::GetSyncVariableDataFromStream(core::Variable<T> &variable,
                                                    BufferSTL &bufferSTL) const
{
    auto itStep = variable.m_AvailableStepBlockIndexOffsets.find(variable.m_StepsStart + 1);

    if (itStep == variable.m_AvailableStepBlockIndexOffsets.end())
    {
        variable.m_Data = nullptr;
        return;
    }

    auto &buffer = bufferSTL.m_Buffer;
    size_t position = itStep->second.front();

    size_t irrelevant;
    const Characteristics<T> characteristics = ReadElementIndexCharacteristics<T>(
        buffer, position, TypeTraits<T>::type_enum, irrelevant, false, m_Minifooter.IsLittleEndian);

    const size_t payloadOffset = characteristics.Statistics.PayloadOffset;
    variable.m_Data = reinterpret_cast<T *>(&buffer[payloadOffset]);
}

template <class T>
typename core::Variable<T>::BPInfo &
BP4Deserializer::InitVariableBlockInfo(core::Variable<T> &variable, T *data) const
{
    const size_t stepsStart = variable.m_StepsStart;
    const size_t stepsCount = variable.m_StepsCount;

    const auto &indices = variable.m_AvailableStepBlockIndexOffsets;
    const size_t maxStep = indices.rbegin()->first;
    if (stepsStart + 1 > maxStep)
    {
        helper::Throw<std::invalid_argument>(
            "Toolkit", "format::bp::BP4Deserializer", "InitVariableBlockInfo",
            "steps start " + std::to_string(stepsStart) +
                " from SetStepsSelection or BeginStep is larger than "
                "the maximum available step " +
                std::to_string(maxStep - 1) + " for variable " + variable.m_Name +
                ", in call to Get");
    }

    auto itStep = std::next(indices.begin(), stepsStart);
    // BlocksInfo() expects absolute step, stepsStart is relative
    const size_t absStep = itStep->first;

    // Check that we have enough steps from stepsStart
    for (size_t i = 0; i < stepsCount; ++i)
    {
        if (itStep == indices.end())
        {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BP4Deserializer", "InitVariableBlockInfo",
                "offset " + std::to_string(i) + " from steps start " + std::to_string(stepsStart) +
                    " in variable " + variable.m_Name +
                    " is beyond the largest available step = " + std::to_string(maxStep - 1) +
                    ", check Variable SetStepSelection argument stepsCount "
                    "(random access), or "
                    "number of BeginStep calls (streaming), in call to Get");
        }
        ++itStep;
    }

    if (variable.m_SelectionType == SelectionType::WriteBlock)
    {
        // BlocksInfo() expects absolute step, stepsStart is relative
        // BlocksInfo() adds +1 to match the step starting from 1
        // but absStep already is the actual step in the map
        const std::vector<typename core::Variable<T>::BPInfo> blocksInfo =
            BlocksInfo(variable, absStep - 1);

        if (variable.m_BlockID >= blocksInfo.size())
        {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BP4Deserializer", "InitVariableBlockInfo",
                "invalid blockID " + std::to_string(variable.m_BlockID) + " from steps start " +
                    std::to_string(stepsStart) + " in variable " + variable.m_Name +
                    ", check argument to Variable<T>::SetBlockID, in call "
                    "to Get");
        }

        // switch to bounding box for global array
        if (variable.m_ShapeID == ShapeID::GlobalArray)
        {
            const Dims &start = blocksInfo[variable.m_BlockID].Start;
            const Dims &count = blocksInfo[variable.m_BlockID].Count;
            // TODO check if we need to reverse dimensions
            variable.SetSelection({start, count});
        }
        else if (variable.m_ShapeID == ShapeID::LocalArray)
        {
            // TODO keep Count for block updated
            variable.m_Count = blocksInfo[variable.m_BlockID].Count;
        }
    }

    // create block info
    return variable.SetBlockInfo(data, stepsStart, stepsCount);
}

template <class T>
void BP4Deserializer::SetVariableBlockInfo(core::Variable<T> &variable,
                                           typename core::Variable<T>::BPInfo &blockInfo) const
{
    auto lf_SetSubStreamInfoOperations = [&](const BPOpInfo &bpOpInfo, const size_t payloadOffset,
                                             helper::SubStreamBoxInfo &subStreamInfo,
                                             const bool isRowMajor)

    {
        helper::BlockOperationInfo blockOperation;
        blockOperation.PayloadOffset = payloadOffset;
        blockOperation.PreShape = bpOpInfo.PreShape;
        blockOperation.PreStart = bpOpInfo.PreStart;
        blockOperation.PreCount = bpOpInfo.PreCount;
        blockOperation.Info["PreDataType"] = ToString(helper::GetDataType<T>());
        // TODO: need to verify it's a match with PreDataType
        // std::to_string(static_cast<size_t>(bp4OpInfo.PreDataType));
        blockOperation.Info["Type"] = bpOpInfo.Type;
        blockOperation.PreSizeOf = sizeof(T);

        // read metadata from supported type and populate Info
        if (m_Minifooter.ADIOSVersion >= 2008000)
        {
            std::memcpy(&blockOperation.PayloadSize, bpOpInfo.Metadata.data() + 8, 8);
        }
        else
        {
            // Files made before 2.8 have incompatible compression operators
            // Add backward compatible fixes here to parse compression metadata
            // directly from the BP4 metadata format
            std::shared_ptr<BPBackCompatOperation> bpOp = SetBPBackCompatOperation(bpOpInfo.Type);
            if (bpOp)
            {
                bpOp->GetMetadata(bpOpInfo.Metadata, blockOperation.Info);
                blockOperation.PayloadSize =
                    static_cast<size_t>(std::stoull(blockOperation.Info.at("OutputSize")));
            }
            else
            {
                std::memcpy(&blockOperation.PayloadSize, bpOpInfo.Metadata.data() + 8, 8);
            }
        }

        subStreamInfo.OperationsInfo.push_back(std::move(blockOperation));
    };

    auto lf_SetSubStreamInfoLocalArray =
        [&](const std::string &variableName, const Box<Dims> &selectionBox,
            typename core::Variable<T>::BPInfo &blockInfo, const size_t step,
            const size_t blockIndexOffset, const BufferSTL &bufferSTL, const bool isRowMajor)

    {
        const std::vector<char> &buffer = bufferSTL.m_Buffer;

        size_t position = blockIndexOffset;
        size_t irrelevant;

        const Characteristics<T> blockCharacteristics =
            ReadElementIndexCharacteristics<T>(buffer, position, TypeTraits<T>::type_enum,
                                               irrelevant, false, m_Minifooter.IsLittleEndian);
        // check if they intersect
        helper::SubStreamBoxInfo subStreamInfo;

        if (helper::GetTotalSize(blockCharacteristics.Count) == 0)
        {
            subStreamInfo.ZeroBlock = true;
        }

        // if selection box start is not empty = local selection
        subStreamInfo.BlockBox = helper::StartEndBox(Dims(blockCharacteristics.Count.size(), 0),
                                                     blockCharacteristics.Count);

        if (!selectionBox.first.empty()) // selection start count was defined
        {
            subStreamInfo.IntersectionBox =
                helper::IntersectionBox(selectionBox, subStreamInfo.BlockBox);
        }
        else // read the entire block
        {
            subStreamInfo.IntersectionBox = subStreamInfo.BlockBox;
        }

        if (subStreamInfo.IntersectionBox.first.empty() ||
            subStreamInfo.IntersectionBox.second.empty())
        {
            return;
        }

        const size_t dimensions = blockCharacteristics.Count.size();
        if (dimensions != blockInfo.Count.size())
        {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BP4Deserializer", "SetVariableBlockInfo",
                "block Count (available) and "
                "selection Count (requested) number of dimensions, do not "
                "match "
                "when reading local array variable " +
                    variableName + ", in call to Get");
        }

        const Dims readInCount = m_ReverseDimensions ? Dims(blockCharacteristics.Count.rbegin(),
                                                            blockCharacteristics.Count.rend())
                                                     : blockCharacteristics.Count;

        const Dims blockInfoStart =
            blockInfo.Start.empty() ? Dims(blockInfo.Count.size(), 0) : blockInfo.Start;

        for (size_t i = 0; i < dimensions; ++i)
        {
            if (blockInfoStart[i] + blockInfo.Count[i] > readInCount[i])
            {
                helper::Throw<std::invalid_argument>(
                    "Toolkit", "format::bp::BP4Deserializer", "SetVariableBlockInfo",
                    "selection Start " + helper::DimsToString(blockInfoStart) + " and Count " +
                        helper::DimsToString(blockInfo.Count) +
                        " (requested) is out of bounds of (available) local"
                        " Count " +
                        helper::DimsToString(readInCount) +
                        " , when reading local array variable " + variableName +
                        ", in call to Get");
            }
        }

        subStreamInfo.Seeks.first =
            sizeof(T) * helper::LinearIndex(subStreamInfo.BlockBox,
                                            subStreamInfo.IntersectionBox.first, isRowMajor);

        subStreamInfo.Seeks.second =
            sizeof(T) * (helper::LinearIndex(subStreamInfo.BlockBox,
                                             subStreamInfo.IntersectionBox.second, isRowMajor) +
                         1);

        const size_t payloadOffset = blockCharacteristics.Statistics.PayloadOffset;

        const BPOpInfo &bpOp = blockCharacteristics.Statistics.Op;
        // if they intersect get info Seeks (first: start, second:
        // count) depending on operation info
        if (bpOp.IsActive)
        {
            lf_SetSubStreamInfoOperations(bpOp, payloadOffset, subStreamInfo, m_IsRowMajor);
        }
        else
        {
            // make it absolute if no operations
            subStreamInfo.Seeks.first += payloadOffset;
            subStreamInfo.Seeks.second += payloadOffset;
        }
        subStreamInfo.SubStreamID = static_cast<size_t>(blockCharacteristics.Statistics.FileIndex);

        blockInfo.StepBlockSubStreamsInfo[step].push_back(std::move(subStreamInfo));
    };

    auto lf_SetSubStreamInfoGlobalArray =
        [&](const std::string &variableName, const Box<Dims> &selectionBox,
            typename core::Variable<T>::BPInfo &blockInfo, const size_t step,
            const size_t blockIndexOffset, const BufferSTL &bufferSTL, const bool isRowMajor)

    {
        const std::vector<char> &buffer = bufferSTL.m_Buffer;

        size_t position = blockIndexOffset;
        size_t irrelevant;

        const Characteristics<T> blockCharacteristics =
            ReadElementIndexCharacteristics<T>(buffer, position, TypeTraits<T>::type_enum,
                                               irrelevant, false, m_Minifooter.IsLittleEndian);

        // check if they intersect
        helper::SubStreamBoxInfo subStreamInfo;

        if (helper::GetTotalSize(blockCharacteristics.Count) == 0)
        {
            subStreamInfo.ZeroBlock = true;
        }

        subStreamInfo.BlockBox =
            helper::StartEndBox(blockCharacteristics.Start, blockCharacteristics.Count);
        subStreamInfo.IntersectionBox =
            helper::IntersectionBox(selectionBox, subStreamInfo.BlockBox);

        if (subStreamInfo.IntersectionBox.first.empty() ||
            subStreamInfo.IntersectionBox.second.empty())
        {
            return;
        }

        // relative position
        subStreamInfo.Seeks.first =
            sizeof(T) * helper::LinearIndex(subStreamInfo.BlockBox,
                                            subStreamInfo.IntersectionBox.first, isRowMajor);

        subStreamInfo.Seeks.second =
            sizeof(T) * (helper::LinearIndex(subStreamInfo.BlockBox,
                                             subStreamInfo.IntersectionBox.second, isRowMajor) +
                         1);

        const size_t payloadOffset = blockCharacteristics.Statistics.PayloadOffset;
        const auto &bp4Op = blockCharacteristics.Statistics.Op;
        // if they intersect get info Seeks (first: start, second:
        // count) depending on operation info
        if (bp4Op.IsActive)
        {
            lf_SetSubStreamInfoOperations(bp4Op, payloadOffset, subStreamInfo, m_IsRowMajor);
        }
        else
        {
            // make it absolute if no operations
            subStreamInfo.Seeks.first += payloadOffset;
            subStreamInfo.Seeks.second += payloadOffset;
        }
        subStreamInfo.SubStreamID = static_cast<size_t>(blockCharacteristics.Statistics.FileIndex);

        blockInfo.StepBlockSubStreamsInfo[step].push_back(std::move(subStreamInfo));
    };

    // BODY OF FUNCTIONS STARTS HERE
    const std::map<size_t, std::vector<size_t>> &indices =
        variable.m_AvailableStepBlockIndexOffsets;

    const Box<Dims> selectionBox =
        helper::StartEndBox(blockInfo.Start, blockInfo.Count, m_ReverseDimensions);

    auto itStep = std::next(indices.begin(), blockInfo.StepsStart);

    for (size_t i = 0; i < blockInfo.StepsCount; ++i)
    {
        const size_t step = itStep->first;
        const std::vector<size_t> &blockOffsets = itStep->second;

        if (variable.m_ShapeID == ShapeID::GlobalArray)
        {
            // Check that dimensions and boundaries are okay
            Dims readInShape = variable.m_AvailableShapes[step];
            const size_t dimensions = readInShape.size();
            if (dimensions != blockInfo.Shape.size())
            {
                helper::Throw<std::invalid_argument>(
                    "Toolkit", "format::bp::BP4Deserializer", "SetVariableBlockInfo",
                    "variable Shape (available) and "
                    "selection Shape (requested) number of dimensions, do not "
                    "match in step " +
                        std::to_string(step) + " when reading global array variable " +
                        variable.m_Name + ", in call to Get");
            }

            for (size_t i = 0; i < dimensions; ++i)
            {
                if (blockInfo.Start[i] + blockInfo.Count[i] > readInShape[i])
                {
                    helper::Throw<std::invalid_argument>(
                        "Toolkit", "format::bp::BP4Deserializer", "SetVariableBlockInfo",
                        "selection Start " + helper::DimsToString(blockInfo.Start) + " and Count " +
                            helper::DimsToString(blockInfo.Count) +
                            " (requested) is out of bounds of (available) "
                            "Shape " +
                            helper::DimsToString(readInShape) +
                            " , when reading global array variable " + variable.m_Name +
                            " in step " + std::to_string(step) + ", in call to Get");
                }
            }
            // Get intersections with each block
            for (const size_t blockOffset : blockOffsets)
            {
                lf_SetSubStreamInfoGlobalArray(variable.m_Name, selectionBox, blockInfo, step,
                                               blockOffset, m_Metadata, m_IsRowMajor);
            }
        }
        else if (variable.m_ShapeID == ShapeID::LocalArray)
        {
            lf_SetSubStreamInfoLocalArray(variable.m_Name, selectionBox, blockInfo, step,
                                          blockOffsets[blockInfo.BlockID], m_Metadata,
                                          m_IsRowMajor);
        }

        ++itStep;
    }
}

template <class T>
void BP4Deserializer::GetValueFromMetadata(core::Variable<T> &variable, T *data) const
{
    const auto &buffer = m_Metadata.m_Buffer;

    const typename core::Variable<T>::BPInfo &blockInfo = InitVariableBlockInfo(variable, data);

    const size_t stepsStart = blockInfo.StepsStart;
    const size_t stepsCount = blockInfo.StepsCount;
    const std::map<size_t, std::vector<size_t>> &indices =
        variable.m_AvailableStepBlockIndexOffsets;
    auto itStep = std::next(indices.begin(), stepsStart);

    size_t dataCounter = 0;
    for (size_t s = 0; s < stepsCount; ++s)
    {
        const std::vector<size_t> &positions = itStep->second;

        // global values only read one block per step
        const size_t blocksStart =
            (variable.m_ShapeID == ShapeID::GlobalArray) ? blockInfo.Start.front() : 0;

        const size_t blocksCount =
            (variable.m_ShapeID == ShapeID::GlobalArray) ? variable.m_Count.front() : 1;

        if (blocksStart + blocksCount > positions.size())
        {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::bp::BP4Deserializer", "GetValueFromMetadata",
                "selection Start {" + std::to_string(blocksStart) + "} and Count {" +
                    std::to_string(blocksCount) +
                    "} (requested) is out of bounds of (available) Shape {" +
                    std::to_string(positions.size()) + "} for relative step " + std::to_string(s) +
                    " , when reading 1D global array variable " + variable.m_Name +
                    ", in call to Get");
        }

        size_t irrelevant;
        for (size_t b = blocksStart; b < blocksStart + blocksCount; ++b)
        {
            size_t localPosition = positions[b];
            const Characteristics<T> characteristics = ReadElementIndexCharacteristics<T>(
                buffer, localPosition, type_string, irrelevant, false, m_Minifooter.IsLittleEndian);

            data[dataCounter] = characteristics.Statistics.Value;
            ++dataCounter;
        }

        ++itStep;
    }

    variable.m_Value = data[0];
}

template <class T>
void BP4Deserializer::PreDataRead(core::Variable<T> &variable,
                                  typename core::Variable<T>::BPInfo &blockInfo,
                                  const helper::SubStreamBoxInfo &subStreamBoxInfo, char *&buffer,
                                  size_t &payloadSize, size_t &payloadOffset, const size_t threadID)
{
    if (subStreamBoxInfo.OperationsInfo.size() > 0)
    {
        const helper::BlockOperationInfo &blockOperationInfo =
            InitPostOperatorBlockData(subStreamBoxInfo.OperationsInfo);

        m_ThreadBuffers[threadID][1].resize(blockOperationInfo.PayloadSize, '\0');

        buffer = m_ThreadBuffers[threadID][1].data();

        payloadSize = blockOperationInfo.PayloadSize;
        payloadOffset = blockOperationInfo.PayloadOffset;
    }
    else
    {
        payloadOffset = subStreamBoxInfo.Seeks.first;
        payloadSize = subStreamBoxInfo.Seeks.second - payloadOffset;
        m_ThreadBuffers[threadID][0].resize(payloadSize);

        buffer = m_ThreadBuffers[threadID][0].data();
    }
}

template <class T>
void BP4Deserializer::PostDataRead(core::Variable<T> &variable,
                                   typename core::Variable<T>::BPInfo &blockInfo,
                                   const helper::SubStreamBoxInfo &subStreamBoxInfo,
                                   const bool isRowMajorDestination, const size_t threadID)
{
    if (subStreamBoxInfo.OperationsInfo.size() > 0)
    {
        if (m_Minifooter.ADIOSVersion >= 2008000)
        {
            const helper::BlockOperationInfo &blockOperationInfo =
                InitPostOperatorBlockData(subStreamBoxInfo.OperationsInfo);

            const size_t preOpPayloadSize =
                helper::GetTotalSize(blockOperationInfo.PreCount) * blockOperationInfo.PreSizeOf;
            m_ThreadBuffers[threadID][0].resize(preOpPayloadSize);

            // get original block back
            char *preOpData = m_ThreadBuffers[threadID][0].data();
            const char *postOpData = m_ThreadBuffers[threadID][1].data();

            std::shared_ptr<core::Operator> op = nullptr;
            for (auto &o : blockInfo.Operations)
            {
                if (o->m_Category == "compress" || o->m_Category == "plugin")
                {
                    op = o;
                    break;
                }
            }
            core::Decompress(postOpData, blockOperationInfo.PayloadSize, preOpData,
                             blockInfo.MemSpace, op);

            // clip block to match selection
            helper::ClipVector(m_ThreadBuffers[threadID][0], subStreamBoxInfo.Seeks.first,
                               subStreamBoxInfo.Seeks.second);
        }
        else
        {
            // Files made before 2.8 have incompatible compression operators
            // Add backward compatible fixes in the function below
            BackCompatDecompress(subStreamBoxInfo, threadID);
        }
    }

#ifdef ADIOS2_HAVE_ENDIAN_REVERSE
    const bool endianReverse =
        (helper::IsLittleEndian() != m_Minifooter.IsLittleEndian) ? true : false;
#else
    constexpr bool endianReverse = false;
#endif

    const Dims blockInfoStart =
        (variable.m_ShapeID == ShapeID::LocalArray && blockInfo.Start.empty())
            ? Dims(blockInfo.Count.size(), 0)
            : blockInfo.Start;

    if (!blockInfo.MemoryStart.empty())
    {
        if (endianReverse)
        {
            helper::Throw<std::invalid_argument>("Toolkit", "format::bp::BP4Deserializer",
                                                 "PostDataRead",
                                                 "endianReverse "
                                                 "not supported with MemorySelection");
        }

        if (m_ReverseDimensions)
        {
            helper::Throw<std::invalid_argument>("Toolkit", "format::bp::BP4Deserializer",
                                                 "PostDataRead",
                                                 "ReverseDimensions not supported with "
                                                 "MemorySelection");
        }

        helper::DimsArray intersectStart(subStreamBoxInfo.IntersectionBox.first);
        helper::DimsArray intersectCount(subStreamBoxInfo.IntersectionBox.second);
        helper::DimsArray blockStart(subStreamBoxInfo.BlockBox.first);
        helper::DimsArray blockCount(subStreamBoxInfo.BlockBox.second);
        helper::DimsArray memoryStart(blockInfoStart);
        for (size_t d = 0; d < intersectStart.size(); d++)
        {
            // change {intersect,block}Count from [start, end] to {start, count}
            intersectCount[d] -= (intersectStart[d] - 1);
            blockCount[d] -= (blockStart[d] - 1);
            // shift everything by MemoryStart
            intersectStart[d] += blockInfo.MemoryStart[d];
            blockStart[d] += blockInfo.MemoryStart[d];
        }
        helper::NdCopy(m_ThreadBuffers[threadID][0].data(), intersectStart, intersectCount, true,
                       true, reinterpret_cast<char *>(blockInfo.Data), intersectStart,
                       intersectCount, true, true, sizeof(T), intersectStart, blockCount,
                       memoryStart, helper::DimsArray(blockInfo.MemoryCount), false);
    }
    else
    {
        helper::ClipContiguousMemory(blockInfo.Data, blockInfoStart, blockInfo.Count,
                                     m_ThreadBuffers[threadID][0].data(), subStreamBoxInfo.BlockBox,
                                     subStreamBoxInfo.IntersectionBox, m_IsRowMajor,
                                     m_ReverseDimensions, endianReverse, blockInfo.MemSpace);
    }
}

void BP4Deserializer::BackCompatDecompress(const helper::SubStreamBoxInfo &subStreamBoxInfo,
                                           const size_t threadID)
{
    // Files made before 2.8 have incompatible compression operators
    // Add backward compatible fixes here
    const helper::BlockOperationInfo &blockOperationInfo =
        InitPostOperatorBlockData(subStreamBoxInfo.OperationsInfo);

    const size_t preOpPayloadSize =
        helper::GetTotalSize(blockOperationInfo.PreCount) * blockOperationInfo.PreSizeOf;
    m_ThreadBuffers[threadID][0].resize(preOpPayloadSize);

    std::string opType = blockOperationInfo.Info.at("Type");

    // get original block back
    char *preOpData = m_ThreadBuffers[threadID][0].data();
    const char *postOpData = m_ThreadBuffers[threadID][1].data();
    // get the right bp4Op
    std::shared_ptr<BPBackCompatOperation> bp4Op = SetBPBackCompatOperation(opType);

    if (bp4Op)
    {
        bp4Op->GetData(postOpData, blockOperationInfo, preOpData);
        // clip block to match selection
        helper::ClipVector(m_ThreadBuffers[threadID][0], subStreamBoxInfo.Seeks.first,
                           subStreamBoxInfo.Seeks.second);
    }
    else
    {
        helper::Throw<std::runtime_error>(
            "Toolkit", "format::bp::BP4Deserializer", "PostDataRead",
            "This file was created by pre-ADIOS 2.8.0 using "
            "compression type " +
                opType +
                ", for which there is no backward compatible reader in this "
                "ADIOS version");
    }
}

template <class T>
std::map<size_t, std::vector<typename core::Variable<T>::BPInfo>>
BP4Deserializer::AllStepsBlocksInfo(const core::Variable<T> &variable) const
{
    std::map<size_t, std::vector<typename core::Variable<T>::BPInfo>> allStepsBlocksInfo;

    for (const auto &pair : variable.m_AvailableStepBlockIndexOffsets)
    {
        const size_t step = pair.first;
        const std::vector<size_t> &blockPositions = pair.second;
        // bp4 index starts at 1
        allStepsBlocksInfo[step - 1] = BlocksInfoCommon(variable, blockPositions);
    }
    return allStepsBlocksInfo;
}

template <class T>
std::vector<std::vector<typename core::Variable<T>::BPInfo>>
BP4Deserializer::AllRelativeStepsBlocksInfo(const core::Variable<T> &variable) const
{
    std::vector<std::vector<typename core::Variable<T>::BPInfo>> allRelativeStepsBlocksInfo(
        variable.m_AvailableStepBlockIndexOffsets.size());

    size_t relativeStep = 0;
    for (const auto &pair : variable.m_AvailableStepBlockIndexOffsets)
    {
        const std::vector<size_t> &blockPositions = pair.second;
        allRelativeStepsBlocksInfo[relativeStep] = BlocksInfoCommon(variable, blockPositions);
        ++relativeStep;
    }
    return allRelativeStepsBlocksInfo;
}

template <class T>
std::vector<typename core::Variable<T>::BPInfo>
BP4Deserializer::BlocksInfo(const core::Variable<T> &variable, const size_t step) const
{
    // bp4 format starts at 1
    auto itStep = variable.m_AvailableStepBlockIndexOffsets.find(step + 1);
    if (itStep == variable.m_AvailableStepBlockIndexOffsets.end())
    {
        return std::vector<typename core::Variable<T>::BPInfo>();
    }
    return BlocksInfoCommon(variable, itStep->second);
}

template <class T>
void BP4Deserializer::ClipContiguousMemory(typename core::Variable<T>::BPInfo &blockInfo,
                                           const std::vector<char> &contiguousMemory,
                                           const Box<Dims> &blockBox,
                                           const Box<Dims> &intersectionBox) const
{
    helper::ClipContiguousMemory(blockInfo.Data, blockInfo.Start, blockInfo.Count, contiguousMemory,
                                 blockBox, intersectionBox, m_IsRowMajor, m_ReverseDimensions);
}

// PRIVATE

template <>
inline void BP4Deserializer::DefineVariableInEngineIOPerStep<std::string>(
    const ElementIndexHeader &header, core::Engine &engine, const std::vector<char> &buffer,
    size_t position, size_t step) const
{
    const size_t initialPosition = position;
    size_t irrelevant;

    const Characteristics<std::string> characteristics =
        ReadElementIndexCharacteristics<std::string>(
            buffer, position, static_cast<DataTypes>(header.DataType), irrelevant, false,
            m_Minifooter.IsLittleEndian);

    const std::string variableName =
        header.Path.empty() ? header.Name : header.Path + PathSeparator + header.Name;

    core::Variable<std::string> *variable = nullptr;
    variable = engine.m_IO.InquireVariable<std::string>(variableName);
    if (variable)
    {
        size_t endPositionCurrentStep =
            initialPosition -
            (header.Name.size() + header.GroupName.size() + header.Path.size() + 23) +
            static_cast<size_t>(header.Length) + 4;
        position = initialPosition;
        // variable->m_AvailableStepsCount = step;
        ++variable->m_AvailableStepsCount;
        // std::cout << variable->m_Name << ", " <<
        // variable->m_AvailableStepsCount << std::endl;
        while (position < endPositionCurrentStep)
        {
            const size_t subsetPosition = position;
            size_t irrelevant;

            // read until step is found
            const Characteristics<std::string> subsetCharacteristics =
                ReadElementIndexCharacteristics<std::string>(
                    buffer, position, static_cast<DataTypes>(header.DataType), irrelevant, false,
                    m_Minifooter.IsLittleEndian);

            if (subsetCharacteristics.EntryShapeID == ShapeID::LocalValue)
            {
                if (subsetPosition == initialPosition)
                {
                    // reset shape and count
                    variable->m_Shape[0] = 1;
                    variable->m_Count[0] = 1;
                }
                else
                {
                    ++variable->m_Shape[0];
                    ++variable->m_Count[0];
                }
            }

            variable->m_AvailableStepBlockIndexOffsets[step].push_back(subsetPosition);
            position = subsetPosition + subsetCharacteristics.EntryLength + 5;
        }
        return;
    }

    if (characteristics.Statistics.IsValue)
    {
        std::lock_guard<std::mutex> lock(m_Mutex);
        variable = &engine.m_IO.DefineVariable<std::string>(variableName);
        engine.RegisterCreatedVariable(variable);
        variable->m_Value = characteristics.Statistics.Value; // assigning first step

        if (characteristics.EntryShapeID == ShapeID::LocalValue)
        {
            variable->m_Shape = {1};
            variable->m_Start = {0};
            variable->m_Count = {1};
            variable->m_ShapeID = ShapeID::LocalValue;
        }
    }
    else
    {
        helper::Throw<std::runtime_error>("Toolkit", "format::bp::BP4Deserializer",
                                          "DefineVariableInEngineIOPerStep",
                                          "variable " + variableName +
                                              " of type string can't be an array, when "
                                              "parsing metadata in call to Open");
    }

    // going back to get variable index position
    variable->m_IndexStart =
        initialPosition - (header.Name.size() + header.GroupName.size() + header.Path.size() + 23);

    const size_t endPosition = variable->m_IndexStart + static_cast<size_t>(header.Length) + 4;

    position = initialPosition;

    size_t currentStep = 0; // Starts at 1 in bp file
    std::set<uint32_t> stepsFound;
    variable->m_AvailableStepsCount = 0;
    while (position < endPosition)
    {
        const size_t subsetPosition = position;
        size_t irrelevant;

        // read until step is found
        const Characteristics<std::string> subsetCharacteristics =
            ReadElementIndexCharacteristics<std::string>(
                buffer, position, static_cast<DataTypes>(header.DataType), irrelevant, false,
                m_Minifooter.IsLittleEndian);

        const bool isNextStep = stepsFound.insert(subsetCharacteristics.Statistics.Step).second;

        // if new step is inserted
        if (isNextStep)
        {
            currentStep = subsetCharacteristics.Statistics.Step;
            ++variable->m_AvailableStepsCount;
            if (subsetCharacteristics.EntryShapeID == ShapeID::LocalValue)
            {
                // reset shape and count
                variable->m_Shape[0] = 1;
                variable->m_Count[0] = 1;
            }
        }
        else
        {
            if (subsetCharacteristics.EntryShapeID == ShapeID::LocalValue)
            {
                ++variable->m_Shape[0];
                ++variable->m_Count[0];
            }
        }

        variable->m_AvailableStepBlockIndexOffsets[currentStep].push_back(subsetPosition);
        position = subsetPosition + subsetCharacteristics.EntryLength + 5;
    }

    if (variable->m_ShapeID == ShapeID::LocalValue)
    {
        variable->m_ShapeID = ShapeID::GlobalArray;
        variable->m_SingleValue = true;
    }
    /* Update variable's starting step, which is always 0 */
    variable->m_StepsStart = 0;

    // update variable Engine for read streaming functions
    variable->m_Engine = &engine;
}

/* Define the variable of each step when parsing the metadata */
template <class T>
void BP4Deserializer::DefineVariableInEngineIOPerStep(const ElementIndexHeader &header,
                                                      core::Engine &engine,
                                                      const std::vector<char> &buffer,
                                                      size_t position, size_t step) const
{
    const size_t initialPosition = position;
    size_t irrelevant;

    const Characteristics<T> characteristics = ReadElementIndexCharacteristics<T>(
        buffer, position, static_cast<DataTypes>(header.DataType), irrelevant, false,
        m_Minifooter.IsLittleEndian);

    const std::string variableName =
        header.Path.empty() ? header.Name : header.Path + PathSeparator + header.Name;

    core::Variable<T> *variable = nullptr;
    {
        // to prevent conflict with DefineVariable
        std::lock_guard<std::mutex> lock(m_Mutex);
        variable = engine.m_IO.InquireVariable<T>(variableName);
    }

    if (variable)
    {
        size_t endPositionCurrentStep =
            initialPosition -
            (header.Name.size() + header.GroupName.size() + header.Path.size() + 23) +
            static_cast<size_t>(header.Length) + 4;
        position = initialPosition;
        // variable->m_AvailableStepsCount = step;
        ++variable->m_AvailableStepsCount;
        while (position < endPositionCurrentStep)
        {
            const size_t subsetPosition = position;
            size_t joinedArrayStartValuePos;

            // read until step is found
            const Characteristics<T> subsetCharacteristics = ReadElementIndexCharacteristics<T>(
                buffer, position, static_cast<DataTypes>(header.DataType), joinedArrayStartValuePos,
                false, m_Minifooter.IsLittleEndian);

            const T blockMin = characteristics.Statistics.IsValue
                                   ? subsetCharacteristics.Statistics.Value
                                   : subsetCharacteristics.Statistics.Min;
            const T blockMax = characteristics.Statistics.IsValue
                                   ? subsetCharacteristics.Statistics.Value
                                   : subsetCharacteristics.Statistics.Max;

            if (helper::LessThan(blockMin, variable->m_Min))
            {
                variable->m_Min = blockMin;
            }

            if (helper::GreaterThan(blockMax, variable->m_Max))
            {
                variable->m_Max = blockMax;
            }

            if (subsetCharacteristics.EntryShapeID == ShapeID::LocalValue)
            {
                if (subsetPosition == initialPosition)
                {
                    // reset shape and count
                    variable->m_Shape[0] = 1;
                    variable->m_Count[0] = 1;
                }
                else
                {
                    ++variable->m_Shape[0];
                    ++variable->m_Count[0];
                }
            }
            else if (subsetCharacteristics.EntryShapeID == ShapeID::JoinedArray)
            {
                Dims shape = m_ReverseDimensions ? Dims(subsetCharacteristics.Shape.rbegin(),
                                                        subsetCharacteristics.Shape.rend())
                                                 : subsetCharacteristics.Shape;
                const Dims count = m_ReverseDimensions ? Dims(subsetCharacteristics.Count.rbegin(),
                                                              subsetCharacteristics.Count.rend())
                                                       : subsetCharacteristics.Count;
                uint64_t newStart;
                if (variable->m_AvailableShapes[step].empty())
                {
                    shape[variable->m_JoinedDimPos] = count[variable->m_JoinedDimPos];
                    newStart = 0;
                }
                else
                {
                    newStart = static_cast<uint64_t>(
                        variable->m_AvailableShapes[step][variable->m_JoinedDimPos]);

                    // shape increases with each block
                    shape[variable->m_JoinedDimPos] =
                        variable->m_AvailableShapes[step][variable->m_JoinedDimPos] +
                        count[variable->m_JoinedDimPos];
                }
                // big hack: modify the metada in place, update the Start[i]
                // of the block
                char *src = reinterpret_cast<char *>(&newStart);
                char *dst = const_cast<char *>(&buffer[joinedArrayStartValuePos]);
                memcpy(dst, src, sizeof(uint64_t));

                variable->m_Shape = shape;
                variable->m_AvailableShapes[step] = shape;
            }
            else if (subsetCharacteristics.EntryShapeID == ShapeID::GlobalArray)
            {
                // Shape definition is by the last block now, not the first
                // block
                const Dims shape = m_ReverseDimensions ? Dims(subsetCharacteristics.Shape.rbegin(),
                                                              subsetCharacteristics.Shape.rend())
                                                       : subsetCharacteristics.Shape;
                variable->m_Shape = shape;
                variable->m_AvailableShapes[step] = shape;
            }

            variable->m_AvailableStepBlockIndexOffsets[step].push_back(subsetPosition);
            position = subsetPosition + subsetCharacteristics.EntryLength + 5;
        }
        return;
    }

    {
        std::lock_guard<std::mutex> lock(m_Mutex);
        switch (characteristics.EntryShapeID)
        {
        case (ShapeID::GlobalValue): {
            variable = &engine.m_IO.DefineVariable<T>(variableName);
            break;
        }
        case (ShapeID::GlobalArray): {
            const Dims shape = m_ReverseDimensions ? Dims(characteristics.Shape.rbegin(),
                                                          characteristics.Shape.rend())
                                                   : characteristics.Shape;

            variable =
                &engine.m_IO.DefineVariable<T>(variableName, shape, Dims(shape.size(), 0), shape);
            variable->m_AvailableShapes[characteristics.Statistics.Step] = variable->m_Shape;
            break;
        }
        case (ShapeID::JoinedArray): {
            Dims shape = m_ReverseDimensions
                             ? Dims(characteristics.Shape.rbegin(), characteristics.Shape.rend())
                             : characteristics.Shape;
            const Dims count = m_ReverseDimensions ? Dims(characteristics.Count.rbegin(),
                                                          characteristics.Count.rend())
                                                   : characteristics.Count;
            size_t joinedDimPos = 0;
            for (size_t i = 0; i < shape.size(); ++i)
            {
                if (shape[i] == JoinedDim)
                {
                    joinedDimPos = i;
                    shape[i] = 0; // will increase with each block
                }
            }

            variable =
                &engine.m_IO.DefineVariable<T>(variableName, shape, Dims(shape.size(), 0), shape);
            variable->m_AvailableShapes[characteristics.Statistics.Step] = variable->m_Shape;
            variable->m_JoinedDimPos = joinedDimPos;
            break;
        }
        case (ShapeID::LocalValue): {
            variable = &engine.m_IO.DefineVariable<T>(variableName, {1}, {0}, {1});
            variable->m_ShapeID = ShapeID::LocalValue;
            break;
        }
        case (ShapeID::LocalArray): {
            const Dims count = m_ReverseDimensions ? Dims(characteristics.Count.rbegin(),
                                                          characteristics.Count.rend())
                                                   : characteristics.Count;
            variable = &engine.m_IO.DefineVariable<T>(variableName, {}, {}, count);
            break;
        }
        default:
            helper::Throw<std::runtime_error>("Toolkit", "format::bp::BP4Deserializer",
                                              "DefineVariableInEngineIOPerStep",
                                              "invalid ShapeID or not yet supported for variable " +
                                                  variableName + ", in call to Open");
        } // end switch

        engine.RegisterCreatedVariable(variable);
        if (characteristics.Statistics.IsValue)
        {
            variable->m_Value = characteristics.Statistics.Value;
            variable->m_Min = characteristics.Statistics.Value;
            variable->m_Max = characteristics.Statistics.Value;
        }
        else
        {
            variable->m_Min = characteristics.Statistics.Min;
            variable->m_Max = characteristics.Statistics.Max;
        }
    } // end mutex lock

    // going back to get variable index position
    variable->m_IndexStart =
        initialPosition - (header.Name.size() + header.GroupName.size() + header.Path.size() + 23);

    const size_t endPosition = variable->m_IndexStart + static_cast<size_t>(header.Length) + 4;

    position = initialPosition;

    size_t currentStep = 0; // Starts at 1 in bp file
    std::set<uint32_t> stepsFound;
    variable->m_AvailableStepsCount = 0;
    while (position < endPosition)
    {
        const size_t subsetPosition = position;
        size_t joinedArrayStartValuePos;

        // read until step is found
        const Characteristics<T> subsetCharacteristics = ReadElementIndexCharacteristics<T>(
            buffer, position, static_cast<DataTypes>(header.DataType), joinedArrayStartValuePos,
            false, m_Minifooter.IsLittleEndian);

        const T blockMin = characteristics.Statistics.IsValue
                               ? subsetCharacteristics.Statistics.Value
                               : subsetCharacteristics.Statistics.Min;
        const T blockMax = characteristics.Statistics.IsValue
                               ? subsetCharacteristics.Statistics.Value
                               : subsetCharacteristics.Statistics.Max;

        const bool isNextStep = stepsFound.insert(subsetCharacteristics.Statistics.Step).second;

        if (isNextStep)
        {
            currentStep = subsetCharacteristics.Statistics.Step;
            ++variable->m_AvailableStepsCount;
            if (subsetCharacteristics.EntryShapeID == ShapeID::LocalValue)
            {
                // reset shape and count
                variable->m_Shape[0] = 1;
                variable->m_Count[0] = 1;
            }
        }
        else
        {
            if (subsetCharacteristics.EntryShapeID == ShapeID::LocalValue)
            {
                ++variable->m_Shape[0];
                ++variable->m_Count[0];
            }
        }

        if (subsetCharacteristics.EntryShapeID == ShapeID::JoinedArray)
        {
            Dims shape = m_ReverseDimensions ? Dims(subsetCharacteristics.Shape.rbegin(),
                                                    subsetCharacteristics.Shape.rend())
                                             : subsetCharacteristics.Shape;
            const Dims count = m_ReverseDimensions ? Dims(subsetCharacteristics.Count.rbegin(),
                                                          subsetCharacteristics.Count.rend())
                                                   : subsetCharacteristics.Count;

            uint64_t newStart;
            if (isNextStep)
            {
                shape[variable->m_JoinedDimPos] = count[variable->m_JoinedDimPos];
                newStart = 0;
            }
            else
            {
                newStart = static_cast<uint64_t>(
                    variable->m_AvailableShapes[currentStep][variable->m_JoinedDimPos]);
                // shape increase with each block
                shape[variable->m_JoinedDimPos] =
                    variable->m_AvailableShapes[currentStep][variable->m_JoinedDimPos] +
                    count[variable->m_JoinedDimPos];
            }
            // big hack: modify the metada in place, update the Start[i]
            // of the block.
            char *src = reinterpret_cast<char *>(&newStart);
            char *dst = const_cast<char *>(&buffer[joinedArrayStartValuePos]);
            memcpy(dst, src, sizeof(uint64_t));

            variable->m_Shape = shape;
            variable->m_AvailableShapes[currentStep] = shape;
        }

        // Shape definition is by the last block now, not the first block
        if (subsetCharacteristics.EntryShapeID == ShapeID::GlobalArray)
        {
            const Dims shape = m_ReverseDimensions ? Dims(subsetCharacteristics.Shape.rbegin(),
                                                          subsetCharacteristics.Shape.rend())
                                                   : subsetCharacteristics.Shape;
            variable->m_Shape = shape;
            variable->m_AvailableShapes[currentStep] = shape;
        }

        // update min max for global values only if new step is found
        if ((isNextStep && subsetCharacteristics.EntryShapeID == ShapeID::GlobalValue) ||
            (subsetCharacteristics.EntryShapeID != ShapeID::GlobalValue))
        {
            if (helper::LessThan(blockMin, variable->m_Min))
            {
                variable->m_Min = blockMin;
            }

            if (helper::GreaterThan(blockMax, variable->m_Max))
            {
                variable->m_Max = blockMax;
            }
        }

        variable->m_AvailableStepBlockIndexOffsets[currentStep].push_back(subsetPosition);
        position = subsetPosition + subsetCharacteristics.EntryLength + 5;
    }

    if (variable->m_ShapeID == ShapeID::LocalValue)
    {
        variable->m_ShapeID = ShapeID::GlobalArray;
        // variable will be presented as a 1D global array, but contents exist
        // in metadata
        variable->m_SingleValue = true;
    }
    /* Update variable's starting step, which is always 0 */
    variable->m_StepsStart = 0;

    // update variable Engine for read streaming functions
    variable->m_Engine = &engine;
}

template <class T>
void BP4Deserializer::DefineAttributeInEngineIO(const ElementIndexHeader &header,
                                                core::Engine &engine,
                                                const std::vector<char> &buffer,
                                                size_t position) const
{
    size_t irrelevant;
    const Characteristics<T> characteristics = ReadElementIndexCharacteristics<T>(
        buffer, position, static_cast<DataTypes>(header.DataType), irrelevant, false,
        m_Minifooter.IsLittleEndian);

    std::string attributeName(header.Name);
    if (!header.Path.empty())
    {
        attributeName = header.Path + PathSeparator + header.Name;
    }

    if (characteristics.Statistics.IsValue)
    {
        engine.m_IO.DefineAttribute<T>(attributeName, characteristics.Statistics.Value, "", "",
                                       true);
    }
    else
    {
        engine.m_IO.DefineAttribute<T>(attributeName, characteristics.Statistics.Values.data(),
                                       characteristics.Statistics.Values.size(), "", "", true);
    }
}

template <class T>
std::map<std::string, helper::SubFileInfoMap>
BP4Deserializer::GetSyncVariableSubFileInfo(const core::Variable<T> &variable) const
{
    std::map<std::string, helper::SubFileInfoMap> variableSubFileInfo;
    variableSubFileInfo[variable.m_Name] = GetSubFileInfo(variable);
    return variableSubFileInfo;
}

template <class T>
void BP4Deserializer::GetDeferredVariable(core::Variable<T> &variable, T *data)
{
    variable.m_Data = data;
    m_DeferredVariablesMap[variable.m_Name] = helper::SubFileInfoMap();
}

template <class T>
helper::SubFileInfoMap BP4Deserializer::GetSubFileInfo(const core::Variable<T> &variable) const
{
    helper::SubFileInfoMap infoMap;

    const auto &buffer = m_Metadata.m_Buffer;

    const size_t stepStart = variable.m_StepsStart + 1;
    const size_t stepEnd = stepStart + variable.m_StepsCount; // exclusive

    const Box<Dims> selectionBox =
        helper::StartEndBox(variable.m_Start, variable.m_Count, m_ReverseDimensions);

    for (size_t step = stepStart; step < stepEnd; ++step)
    {
        auto itBlockStarts = variable.m_AvailableStepBlockIndexOffsets.find(step);
        if (itBlockStarts == variable.m_AvailableStepBlockIndexOffsets.end())
        {
            continue;
        }

        const std::vector<size_t> &blockStarts = itBlockStarts->second;

        // blockPosition gets updated by Read, can't be const
        for (size_t blockPosition : blockStarts)
        {
            size_t irrelevant;
            const Characteristics<T> blockCharacteristics =
                ReadElementIndexCharacteristics<T>(buffer, blockPosition, TypeTraits<T>::type_enum,
                                                   irrelevant, false, m_Minifooter.IsLittleEndian);

            // check if they intersect
            helper::SubFileInfo info;
            info.BlockBox =
                helper::StartEndBox(blockCharacteristics.Start, blockCharacteristics.Count);
            info.IntersectionBox = helper::IntersectionBox(selectionBox, info.BlockBox);

            if (info.IntersectionBox.first.empty() || info.IntersectionBox.second.empty())
            {
                continue;
            }
            // if they intersect get info Seeks (first: start, second:
            // count)
            info.Seeks.first =
                blockCharacteristics.Statistics.PayloadOffset +
                helper::LinearIndex(info.BlockBox, info.IntersectionBox.first, m_IsRowMajor) *
                    sizeof(T);

            info.Seeks.second =
                blockCharacteristics.Statistics.PayloadOffset +
                (helper::LinearIndex(info.BlockBox, info.IntersectionBox.second, m_IsRowMajor) +
                 1) *
                    sizeof(T);

            const size_t fileIndex = static_cast<size_t>(blockCharacteristics.Statistics.FileIndex);

            infoMap[fileIndex][step].push_back(std::move(info));
        }
    }

    return infoMap;
}

// PRIVATE
template <class T>
std::vector<typename core::Variable<T>::BPInfo>
BP4Deserializer::BlocksInfoCommon(const core::Variable<T> &variable,
                                  const std::vector<size_t> &blocksIndexOffsets) const
{
    std::vector<typename core::Variable<T>::BPInfo> blocksInfo;
    blocksInfo.reserve(blocksIndexOffsets.size());

    size_t n = 0;
    for (const size_t blockIndexOffset : blocksIndexOffsets)
    {
        size_t position = blockIndexOffset;
        size_t irrelevant;

        const Characteristics<T> blockCharacteristics = ReadElementIndexCharacteristics<T>(
            m_Metadata.m_Buffer, position, TypeTraits<T>::type_enum, irrelevant, false,
            m_Minifooter.IsLittleEndian);

        typename core::Variable<T>::BPInfo blockInfo;
        blockInfo.Shape = blockCharacteristics.Shape;
        blockInfo.Start = blockCharacteristics.Start;
        blockInfo.Count = blockCharacteristics.Count;
        blockInfo.WriterID = blockCharacteristics.Statistics.FileIndex;
        blockInfo.IsReverseDims = m_ReverseDimensions;

        if (m_ReverseDimensions)
        {
            std::reverse(blockInfo.Shape.begin(), blockInfo.Shape.end());
            std::reverse(blockInfo.Start.begin(), blockInfo.Start.end());
            std::reverse(blockInfo.Count.begin(), blockInfo.Count.end());
        }

        if (blockCharacteristics.Statistics.IsValue) // value
        {
            blockInfo.IsValue = true;
            blockInfo.Value = blockCharacteristics.Statistics.Value;
        }
        else // array
        {
            blockInfo.IsValue = false;
            blockInfo.Min = blockCharacteristics.Statistics.Min;
            blockInfo.Max = blockCharacteristics.Statistics.Max;
            blockInfo.MinMaxs = blockCharacteristics.Statistics.MinMaxs;
            blockInfo.SubBlockInfo = blockCharacteristics.Statistics.SubBlockInfo;
        }
        if (blockInfo.Shape.size() == 1 && blockInfo.Shape.front() == LocalValueDim)
        {
            blockInfo.Shape = Dims{blocksIndexOffsets.size()};
            blockInfo.Count = Dims{1};
            blockInfo.Start = Dims{n};
            blockInfo.Min = blockCharacteristics.Statistics.Value;
            blockInfo.Max = blockCharacteristics.Statistics.Value;
        }
        // bp index starts at 1
        blockInfo.Step = static_cast<size_t>(blockCharacteristics.Statistics.Step - 1);
        blockInfo.BlockID = n;
        blocksInfo.push_back(blockInfo);
        ++n;
    }
    return blocksInfo;
}

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP4_BP4DESERIALIZER_TCC_ */
