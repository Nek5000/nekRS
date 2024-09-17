/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderGeneric.cpp
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#include "SscReaderGeneric.tcc"

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

SscReaderGeneric::SscReaderGeneric(IO &io, const std::string &name, const Mode mode, MPI_Comm comm)
: SscReaderBase(io, name, mode, comm)
{
}

void SscReaderGeneric::BeginStepConsequentFixed()
{
    MPI_Waitall(static_cast<int>(m_MpiRequests.size()), m_MpiRequests.data(), MPI_STATUS_IGNORE);
    m_MpiRequests.clear();
}

void SscReaderGeneric::BeginStepFlexible(StepStatus &status)
{
    m_AllReceivingWriterRanks.clear();
    m_Buffer.resize(1);
    m_Buffer[0] = 0;
    m_GlobalWritePattern.clear();
    m_GlobalWritePattern.resize(m_StreamSize);
    m_LocalReadPattern.clear();
    m_GlobalWritePatternBuffer.clear();
    bool finalStep = SyncWritePattern();
    if (finalStep)
    {
        status = StepStatus::EndOfStream;
        return;
    }
    MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, m_StreamComm, &m_MpiWin);
}

StepStatus SscReaderGeneric::BeginStep(const StepMode stepMode, const float timeoutSeconds,
                                       const bool readerLocked)
{

    m_ReaderSelectionsLocked = readerLocked;

    ++m_CurrentStep;

    m_StepBegun = true;

    if (m_CurrentStep == 0 || m_WriterDefinitionsLocked == false ||
        m_ReaderSelectionsLocked == false)
    {
        if (m_Threading && m_EndStepThread.joinable())
        {
            m_EndStepThread.join();
        }
        else
        {
            BeginStepFlexible(m_StepStatus);
        }
        if (m_StepStatus == StepStatus::EndOfStream)
        {
            return StepStatus::EndOfStream;
        }
    }
    else
    {
        BeginStepConsequentFixed();
    }

    for (const auto &r : m_GlobalWritePattern)
    {
        for (auto &v : r)
        {
            if (v.shapeId == ShapeID::GlobalValue || v.shapeId == ShapeID::LocalValue)
            {
                std::vector<char> value(v.bufferCount);
                if (m_CurrentStep == 0 || m_WriterDefinitionsLocked == false ||
                    m_ReaderSelectionsLocked == false)
                {
                    std::memcpy(value.data(), v.value.data(), v.value.size());
                }
                else
                {
                    std::memcpy(value.data(), m_Buffer.data() + v.bufferStart, v.bufferCount);
                }
                if (v.type == DataType::String)
                {
                    auto variable = m_IO.InquireVariable<std::string>(v.name);
                    if (variable)
                    {
                        variable->m_Value = std::string(value.begin(), value.end());
                        variable->m_Min = std::string(value.begin(), value.end());
                        variable->m_Max = std::string(value.begin(), value.end());
                    }
                }
#define declare_type(T)                                                                            \
    else if (v.type == helper::GetDataType<T>())                                                   \
    {                                                                                              \
        auto variable = m_IO.InquireVariable<T>(v.name);                                           \
        if (variable)                                                                              \
        {                                                                                          \
            std::memcpy(reinterpret_cast<char *>(&variable->m_Min), value.data(), value.size());   \
            std::memcpy(reinterpret_cast<char *>(&variable->m_Max), value.data(), value.size());   \
            std::memcpy(reinterpret_cast<char *>(&variable->m_Value), value.data(), value.size()); \
        }                                                                                          \
    }
                ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
                else
                {
                    helper::Log("Engine", "SscReaderGeneric", "BeginStep", "unknown data type",
                                m_ReaderRank, m_ReaderRank, 0, m_Verbosity, helper::FATALERROR);
                }
            }
        }
    }

    if (m_Buffer[0] == 1)
    {
        return StepStatus::EndOfStream;
    }

    return StepStatus::OK;
}

size_t SscReaderGeneric::CurrentStep() { return m_CurrentStep; }

void SscReaderGeneric::EndStepFixed()
{
    if (m_CurrentStep == 0)
    {
        MPI_Win_free(&m_MpiWin);
        SyncReadPattern();
    }
    for (const auto &i : m_AllReceivingWriterRanks)
    {
        m_MpiRequests.emplace_back();
        MPI_Irecv(m_Buffer.data() + i.second.first, static_cast<int>(i.second.second), MPI_CHAR,
                  i.first, 0, m_StreamComm, &m_MpiRequests.back());
    }
}

void SscReaderGeneric::EndStepFirstFlexible()
{
    MPI_Win_free(&m_MpiWin);
    SyncReadPattern();
    BeginStepFlexible(m_StepStatus);
}

void SscReaderGeneric::EndStepConsequentFlexible()
{
    MPI_Win_free(&m_MpiWin);
    BeginStepFlexible(m_StepStatus);
}

void SscReaderGeneric::EndStep(const bool readerLocked)
{
    m_ReaderSelectionsLocked = readerLocked;
    PerformGets();

    if (m_WriterDefinitionsLocked && m_ReaderSelectionsLocked)
    {
        EndStepFixed();
    }
    else
    {
        if (m_CurrentStep == 0)
        {
            if (m_Threading)
            {
                m_EndStepThread = std::thread(&SscReaderGeneric::EndStepFirstFlexible, this);
            }
            else
            {
                MPI_Win_free(&m_MpiWin);
                SyncReadPattern();
            }
        }
        else
        {
            if (m_Threading)
            {
                m_EndStepThread = std::thread(&SscReaderGeneric::EndStepConsequentFlexible, this);
            }
            else
            {
                MPI_Win_free(&m_MpiWin);
            }
        }
    }

    m_StepBegun = false;
}

void SscReaderGeneric::PerformGets()
{

    if (m_CurrentStep == 0 || m_WriterDefinitionsLocked == false ||
        m_ReaderSelectionsLocked == false)
    {
        ssc::Deserialize(m_GlobalWritePatternBuffer, m_GlobalWritePattern, m_IO, false, false,
                         false, m_StructDefinitions);
        size_t oldSize = m_AllReceivingWriterRanks.size();
        m_AllReceivingWriterRanks = ssc::CalculateOverlap(m_GlobalWritePattern, m_LocalReadPattern);
        CalculatePosition(m_GlobalWritePattern, m_AllReceivingWriterRanks);
        size_t newSize = m_AllReceivingWriterRanks.size();
        if (oldSize != newSize)
        {
            size_t totalDataSize = 0;
            for (auto i : m_AllReceivingWriterRanks)
            {
                totalDataSize += i.second.second;
            }
            m_Buffer.resize(totalDataSize);
            for (const auto &i : m_AllReceivingWriterRanks)
            {
                MPI_Win_lock(MPI_LOCK_SHARED, i.first, 0, m_MpiWin);
                MPI_Get(m_Buffer.data() + i.second.first, static_cast<int>(i.second.second),
                        MPI_CHAR, i.first, 0, static_cast<int>(i.second.second), MPI_CHAR,
                        m_MpiWin);
                MPI_Win_unlock(i.first, m_MpiWin);
            }
        }

        for (auto &br : m_LocalReadPattern)
        {
            if (br.performed)
            {
                continue;
            }
            for (const auto &i : m_AllReceivingWriterRanks)
            {
                const auto &v = m_GlobalWritePattern[i.first];
                for (const auto &b : v)
                {
                    if (b.name == br.name)
                    {
                        if (b.type == DataType::None)
                        {
                            helper::Log("Engine", "SscReaderGeneric", "PerformGets",
                                        "unknown data type", m_ReaderRank, m_ReaderRank, 0,
                                        m_Verbosity, helper::FATALERROR);
                        }
                        else if (b.type == DataType::String)
                        {
                            *reinterpret_cast<std::string *>(br.data) =
                                std::string(b.value.begin(), b.value.end());
                        }
                        else
                        {
                            if (b.shapeId == ShapeID::GlobalArray ||
                                b.shapeId == ShapeID::LocalArray)
                            {
                                bool empty = false;
                                for (const auto c : b.count)
                                {
                                    if (c == 0)
                                    {
                                        empty = true;
                                    }
                                }
                                if (empty)
                                {
                                    continue;
                                }
                                helper::NdCopy(
                                    m_Buffer.data<char>() + b.bufferStart,
                                    helper::CoreDims(b.start), helper::CoreDims(b.count), true,
                                    true, reinterpret_cast<char *>(br.data),
                                    helper::CoreDims(br.start), helper::CoreDims(br.count), true,
                                    true, static_cast<int>(b.elementSize),
                                    helper::CoreDims(b.start), helper::CoreDims(b.count),
                                    helper::CoreDims(br.memStart), helper::CoreDims(br.memCount));
                            }
                            else if (b.shapeId == ShapeID::GlobalValue ||
                                     b.shapeId == ShapeID::LocalValue)
                            {
                                std::memcpy(br.data, m_Buffer.data() + b.bufferStart,
                                            b.bufferCount);
                            }
                        }
                    }
                }
            }
            br.performed = true;
        }
    }
}

bool SscReaderGeneric::SyncWritePattern()
{

    ssc::BroadcastMetadata(m_GlobalWritePatternBuffer, m_WriterMasterStreamRank, m_StreamComm);

    if (m_GlobalWritePatternBuffer[0] == 1)
    {
        return true;
    }

    m_WriterDefinitionsLocked = m_GlobalWritePatternBuffer[1];

    ssc::Deserialize(m_GlobalWritePatternBuffer, m_GlobalWritePattern, m_IO, true, true, true,
                     m_StructDefinitions);

    if (m_Verbosity >= 10 && m_ReaderRank == 0)
    {
        ssc::PrintBlockVecVec(m_GlobalWritePattern, "Global Write Pattern");
    }
    return false;
}

void SscReaderGeneric::SyncReadPattern()
{

    ssc::Buffer localBuffer(8);
    localBuffer.value<uint64_t>() = 8;

    ssc::SerializeVariables(m_LocalReadPattern, localBuffer, m_StreamRank);

    ssc::Buffer globalBuffer;

    ssc::AggregateMetadata(localBuffer, globalBuffer, m_ReaderComm, false,
                           m_ReaderSelectionsLocked);

    ssc::BroadcastMetadata(globalBuffer, m_ReaderMasterStreamRank, m_StreamComm);

    ssc::Deserialize(m_GlobalWritePatternBuffer, m_GlobalWritePattern, m_IO, true, true, true,
                     m_StructDefinitions);

    m_AllReceivingWriterRanks = ssc::CalculateOverlap(m_GlobalWritePattern, m_LocalReadPattern);
    CalculatePosition(m_GlobalWritePattern, m_AllReceivingWriterRanks);

    size_t totalDataSize = 0;
    for (auto i : m_AllReceivingWriterRanks)
    {
        totalDataSize += i.second.second;
    }
    m_Buffer.resize(totalDataSize);

    if (m_Verbosity >= 20)
    {
        for (int i = 0; i < m_ReaderSize; ++i)
        {
            MPI_Barrier(m_ReaderComm);
            if (i == m_ReaderRank)
            {
                ssc::PrintBlockVec(m_LocalReadPattern, "\n\nGlobal Read Pattern on Rank " +
                                                           std::to_string(m_ReaderRank));
            }
        }
        MPI_Barrier(m_ReaderComm);
    }
}

void SscReaderGeneric::CalculatePosition(ssc::BlockVecVec &bvv, ssc::RankPosMap &allRanks)
{

    size_t bufferPosition = 0;

    for (int rank = 0; rank < static_cast<int>(bvv.size()); ++rank)
    {
        bool hasOverlap = false;
        for (const auto &r : allRanks)
        {
            if (r.first == rank)
            {
                hasOverlap = true;
                break;
            }
        }
        if (hasOverlap)
        {
            allRanks[rank].first = bufferPosition;
            auto &bv = bvv[rank];
            for (auto &b : bv)
            {
                b.bufferStart += bufferPosition;
            }
            size_t currentRankTotalSize = ssc::TotalDataSize(bv);
            allRanks[rank].second = currentRankTotalSize + 1;
            bufferPosition += currentRankTotalSize + 1;
        }
    }
}

void SscReaderGeneric::Close(const int transportIndex)
{
    if (!m_StepBegun)
    {
        BeginStep(StepMode::Read, -1.0, m_ReaderSelectionsLocked);
    }
}

#define declare_type(T)                                                                            \
    std::vector<typename Variable<T>::BPInfo> SscReaderGeneric::BlocksInfo(                        \
        const Variable<T> &variable, const size_t step) const                                      \
    {                                                                                              \
        return BlocksInfoCommon(variable, step);                                                   \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void SscReaderGeneric::GetDeferredDeltaCommon(VariableBase &variable, void *data)
{

    Dims vStart = variable.m_Start;
    Dims vCount = variable.m_Count;
    Dims vShape = variable.m_Shape;
    Dims vMemStart = variable.m_MemoryStart;
    Dims vMemCount = variable.m_MemoryCount;

    if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
    {
        std::reverse(vStart.begin(), vStart.end());
        std::reverse(vCount.begin(), vCount.end());
        std::reverse(vShape.begin(), vShape.end());
        std::reverse(vMemStart.begin(), vMemStart.end());
        std::reverse(vMemCount.begin(), vMemCount.end());
    }

    m_LocalReadPattern.emplace_back();
    auto &b = m_LocalReadPattern.back();
    b.name = variable.m_Name;
    b.type = variable.m_Type;
    b.elementSize = variable.m_ElementSize;
    b.shapeId = variable.m_ShapeID;
    b.start = vStart;
    b.count = vCount;
    b.shape = vShape;
    b.memStart = vMemStart;
    b.memCount = vMemCount;
    b.bufferStart = 0;
    b.bufferCount = 0;
    b.data = data;
    b.performed = false;

    for (const auto &d : b.count)
    {
        if (d == 0)
        {
            helper::Throw<std::invalid_argument>("Engine", "SscReader", "GetDeferredDeltaCommon",
                                                 "SetSelection count dimensions cannot be 0");
        }
    }
}

void SscReaderGeneric::GetDeferred(VariableBase &variable, void *data)
{

    if (variable.m_Type == DataType::String)
    {
        auto *dataString = reinterpret_cast<std::string *>(data);
        if (m_CurrentStep == 0 || m_WriterDefinitionsLocked == false ||
            m_ReaderSelectionsLocked == false)
        {
            GetDeferredDeltaCommon(variable, data);
        }
        else
        {
            for (const auto &i : m_AllReceivingWriterRanks)
            {
                const auto &v = m_GlobalWritePattern[i.first];
                for (const auto &b : v)
                {
                    if (b.name == variable.m_Name)
                    {
                        *dataString = std::string(b.value.begin(), b.value.end());
                    }
                }
            }
        }
        return;
    }

    helper::DimsArray vStart(variable.m_Start);
    helper::DimsArray vCount(variable.m_Count);
    helper::DimsArray vShape(variable.m_Shape);
    helper::DimsArray vMemStart(variable.m_MemoryStart);
    helper::DimsArray vMemCount(variable.m_MemoryCount);

    if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
    {
        std::reverse(vStart.begin(), vStart.end());
        std::reverse(vCount.begin(), vCount.end());
        std::reverse(vShape.begin(), vShape.end());
        std::reverse(vMemStart.begin(), vMemStart.end());
        std::reverse(vMemCount.begin(), vMemCount.end());
    }

    if (m_CurrentStep == 0 || m_WriterDefinitionsLocked == false ||
        m_ReaderSelectionsLocked == false)
    {
        GetDeferredDeltaCommon(variable, data);
    }
    else
    {

        for (const auto &i : m_AllReceivingWriterRanks)
        {
            const auto &v = m_GlobalWritePattern[i.first];
            for (const auto &b : v)
            {
                if (b.name == variable.m_Name)
                {
                    bool empty = false;
                    for (const auto c : b.count)
                    {
                        if (c == 0)
                        {
                            empty = true;
                        }
                    }
                    if (empty)
                    {
                        continue;
                    }

                    if (b.shapeId == ShapeID::GlobalArray || b.shapeId == ShapeID::LocalArray)
                    {
                        helper::NdCopy(m_Buffer.data<char>() + b.bufferStart,
                                       helper::CoreDims(b.start), helper::CoreDims(b.count), true,
                                       true, reinterpret_cast<char *>(data), vStart, vCount, true,
                                       true, static_cast<int>(variable.m_ElementSize),
                                       helper::CoreDims(b.start), helper::CoreDims(b.count),
                                       vMemStart, vMemCount);
                    }
                    else if (b.shapeId == ShapeID::GlobalValue || b.shapeId == ShapeID::LocalValue)
                    {
                        std::memcpy(data, m_Buffer.data() + b.bufferStart, b.bufferCount);
                    }
                    else
                    {
                        helper::Log("Engine", "SscReaderGeneric", "GetDeferredCommon",
                                    "unknown ShapeID", m_ReaderRank, m_ReaderRank, 0, m_Verbosity,
                                    helper::LogMode::FATALERROR);
                    }
                }
            }
        }
    }
}

std::vector<VariableStruct::BPInfo> SscReaderGeneric::BlocksInfo(const VariableStruct &variable,
                                                                 const size_t step) const
{
    std::vector<VariableStruct::BPInfo> ret;
    size_t blockID = 0;
    for (int i = 0; i < static_cast<int>(m_GlobalWritePattern.size()); ++i)
    {
        for (auto &v : m_GlobalWritePattern[i])
        {
            if (v.name == variable.m_Name)
            {
                ret.emplace_back();
                auto &b = ret.back();
                b.Start = v.start;
                b.Count = v.count;
                b.Shape = v.shape;
                b.Step = m_CurrentStep;
                b.StepsStart = m_CurrentStep;
                b.StepsCount = 1;
                b.WriterID = i;
                b.BlockID = blockID;
                if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
                {
                    b.IsReverseDims = true;
                }
                ++blockID;
            }
        }
    }
    return ret;
}
}
}
}
}
