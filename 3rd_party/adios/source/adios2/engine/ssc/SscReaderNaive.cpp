/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderNaive.cpp
 *
 *  Created on: Mar 7, 2022
 *      Author: Jason Wang
 */

#include "SscReaderNaive.tcc"

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

SscReaderNaive::SscReaderNaive(IO &io, const std::string &name, const Mode mode, MPI_Comm comm)
: SscReaderBase(io, name, mode, comm)
{
}

StepStatus SscReaderNaive::BeginStep(const StepMode stepMode, const float timeoutSeconds,
                                     const bool readerLocked)
{

    m_Buffer.clear();
    m_BlockMap.clear();

    ++m_CurrentStep;

    int globalSize;

    if (m_ReaderRank == 0)
    {
        MPI_Recv(&globalSize, 1, MPI_INT, m_WriterMasterStreamRank, 0, m_StreamComm,
                 MPI_STATUS_IGNORE);
        m_Buffer.resize(globalSize);
        //        MPI_Recv(m_Buffer.data(), globalSize, MPI_CHAR,
        //        m_WriterMasterStreamRank, 0, m_StreamComm, MPI_STATUS_IGNORE);
        // TODO: revert when the crusher MPI bug is fixed
        ssc::Buffer tmp(globalSize);
        MPI_Recv(tmp.data(), globalSize, MPI_CHAR, m_WriterMasterStreamRank, 0, m_StreamComm,
                 MPI_STATUS_IGNORE);
        std::memcpy(m_Buffer.data(), tmp.data(), globalSize);
    }

    MPI_Bcast(&globalSize, 1, MPI_INT, 0, m_ReaderComm);
    if (m_ReaderRank != 0)
    {
        m_Buffer.resize(globalSize);
    }
    if (globalSize == 1)
    {
        return StepStatus::EndOfStream;
    }
    MPI_Bcast(m_Buffer.data(), globalSize, MPI_CHAR, 0, m_ReaderComm);

    uint64_t pos = 0;
    while (pos < m_Buffer.size())
    {
        uint64_t start = pos;
        uint64_t end = pos + m_Buffer.value<uint64_t>(pos);
        pos += m_Buffer.value<uint64_t>(pos + 8);

        while (pos < end)
        {
            uint8_t shapeId = m_Buffer[pos];
            ++pos;

            if (shapeId == 65)
            {
                DeserializeStructDefinitions(m_Buffer, pos, m_IO, true, m_StructDefinitions);
            }
            else if (shapeId == 66)
            {
                DeserializeAttribute(m_Buffer, pos, m_IO, true);
            }
            else
            {
                pos += 4;
                ssc::BlockInfo b;
                DeserializeVariable(m_Buffer, static_cast<ShapeID>(shapeId), pos, b, m_IO, true,
                                    m_StructDefinitions);
                b.bufferStart += start;
                m_BlockMap[b.name].push_back(b);
                if (b.shapeId == ShapeID::GlobalValue || b.shapeId == ShapeID::LocalValue)
                {
                    std::vector<char> value(b.bufferCount);
                    std::memcpy(value.data(), b.value.data(), b.value.size());

                    if (b.type == DataType::String)
                    {
                        auto variable = m_IO.InquireVariable<std::string>(b.name);
                        if (variable)
                        {
                            variable->m_Value = std::string(value.begin(), value.end());
                            variable->m_Min = std::string(value.begin(), value.end());
                            variable->m_Max = std::string(value.begin(), value.end());
                        }
                    }
#define declare_type(T)                                                                            \
    else if (b.type == helper::GetDataType<T>())                                                   \
    {                                                                                              \
        auto variable = m_IO.InquireVariable<T>(b.name);                                           \
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
                        helper::Log("Engine", "SscReaderNaive", "BeginStep", "unknown data type",
                                    m_ReaderRank, m_ReaderRank, 0, m_Verbosity, helper::FATALERROR);
                    }
                }
            }
        }
    }

    return StepStatus::OK;
}

size_t SscReaderNaive::CurrentStep() { return m_CurrentStep; }

void SscReaderNaive::EndStep(const bool readerLocked) {}

void SscReaderNaive::PerformGets() {}

void SscReaderNaive::Close(const int transportIndex) {}

#define declare_type(T)                                                                            \
    std::vector<typename Variable<T>::BPInfo> SscReaderNaive::BlocksInfo(                          \
        const Variable<T> &variable, const size_t step) const                                      \
    {                                                                                              \
        return BlocksInfoCommon(variable, step);                                                   \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void SscReaderNaive::GetDeferred(VariableBase &variable, void *data)
{

    if (variable.m_Type == DataType::String)
    {
        auto *variableString = dynamic_cast<Variable<std::string> *>(&variable);
        auto *dataString = reinterpret_cast<std::string *>(data);
        for (const auto &b : m_BlockMap[variable.m_Name])
        {
            if (b.name == variable.m_Name)
            {
                *dataString = std::string(m_Buffer.data() + b.bufferStart,
                                          m_Buffer.data() + b.bufferStart + b.bufferCount);
                variableString->m_Value = *dataString;
                variableString->m_Min = *dataString;
                variableString->m_Max = *dataString;
            }
        }
        return;
    }

    helper::DimsArray vStart(variable.m_Start);
    helper::DimsArray vCount(variable.m_Count);
    helper::DimsArray vMemStart(variable.m_MemoryStart);
    helper::DimsArray vMemCount(variable.m_MemoryCount);

    if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
    {
        std::reverse(vStart.begin(), vStart.end());
        std::reverse(vCount.begin(), vCount.end());
        std::reverse(vMemStart.begin(), vMemStart.end());
        std::reverse(vMemCount.begin(), vMemCount.end());
    }

    for (const auto &b : m_BlockMap[variable.m_Name])
    {
        if (b.shapeId == ShapeID::GlobalArray || b.shapeId == ShapeID::LocalArray)
        {
            helper::NdCopy(m_Buffer.data<char>() + b.bufferStart, helper::CoreDims(b.start),
                           helper::CoreDims(b.count), true, true, reinterpret_cast<char *>(data),
                           vStart, vCount, true, true, static_cast<int>(variable.m_ElementSize),
                           helper::CoreDims(b.start), helper::CoreDims(b.count), vMemStart,
                           vMemCount);
        }
        else if (b.shapeId == ShapeID::GlobalValue || b.shapeId == ShapeID::LocalValue)
        {
            std::memcpy(data, m_Buffer.data() + b.bufferStart, b.bufferCount);
        }
    }
}

std::vector<VariableStruct::BPInfo> SscReaderNaive::BlocksInfo(const VariableStruct &variable,
                                                               const size_t step) const
{
    std::vector<VariableStruct::BPInfo> ret;
    size_t blockID = 0;
    auto it = m_BlockMap.find(variable.m_Name);
    if (it != m_BlockMap.end())
    {
        for (const auto &v : it->second)
        {
            ret.emplace_back();
            auto &b = ret.back();
            b.Start = v.start;
            b.Count = v.count;
            b.Shape = v.shape;
            b.Step = m_CurrentStep;
            b.StepsStart = m_CurrentStep;
            b.StepsCount = 1;
            b.BlockID = blockID;
            if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
            {
                b.IsReverseDims = true;
            }
            ++blockID;
        }
    }
    return ret;
}

}
}
}
}
