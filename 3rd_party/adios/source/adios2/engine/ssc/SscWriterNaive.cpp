/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscWriterNaive.cpp
 *
 *  Created on: Mar 7, 2022
 *      Author: Jason Wang
 */

#include "SscWriterNaive.h"

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

SscWriterNaive::SscWriterNaive(IO &io, const std::string &name, const Mode mode, MPI_Comm comm)
: SscWriterBase(io, name, mode, comm)
{
}

StepStatus SscWriterNaive::BeginStep(const StepMode mode, const float timeoutSeconds,
                                     const bool writerLocked)
{
    ++m_CurrentStep;

    m_Buffer.clear();
    m_Buffer.resize(16);
    m_Metadata.clear();

    return StepStatus::OK;
}

size_t SscWriterNaive::CurrentStep() { return m_CurrentStep; }

void SscWriterNaive::PerformPuts() {}

void SscWriterNaive::EndStep(const bool writerLocked)
{
    m_Buffer.value<uint64_t>() = m_Buffer.size();
    m_Buffer.value<uint64_t>(8) = m_Buffer.size();

    if (m_WriterRank == 0)
    {
        ssc::SerializeStructDefinitions(m_IO.m_ADIOS.m_StructDefinitions, m_Buffer);
    }

    if (m_WriterRank == m_WriterSize - 1)
    {
        ssc::SerializeAttributes(m_IO, m_Buffer);
    }

    ssc::SerializeVariables(m_Metadata, m_Buffer, m_WriterRank);

    int localSize = static_cast<int>(m_Buffer.value<uint64_t>());
    std::vector<int> localSizes(m_WriterSize);
    MPI_Gather(&localSize, 1, MPI_INT, localSizes.data(), 1, MPI_INT, 0, m_WriterComm);
    int globalSize = std::accumulate(localSizes.begin(), localSizes.end(), 0);
    ssc::Buffer globalBuffer(globalSize);

    std::vector<int> displs(m_WriterSize);
    for (size_t i = 1; i < static_cast<size_t>(m_WriterSize); ++i)
    {
        displs[i] = displs[i - 1] + localSizes[i - 1];
    }

    MPI_Gatherv(m_Buffer.data(), localSize, MPI_CHAR, globalBuffer.data(), localSizes.data(),
                displs.data(), MPI_CHAR, 0, m_WriterComm);

    if (m_WriterRank == 0)
    {
        MPI_Send(&globalSize, 1, MPI_INT, m_ReaderMasterStreamRank, 0, m_StreamComm);
        MPI_Send(globalBuffer.data(), globalSize, MPI_CHAR, m_ReaderMasterStreamRank, 0,
                 m_StreamComm);
    }
}

void SscWriterNaive::Close(const int transportIndex)
{

    int globalSize = 1;
    ssc::Buffer globalBuffer(globalSize);
    if (m_WriterRank == 0)
    {
        MPI_Send(&globalSize, 1, MPI_INT, m_ReaderMasterStreamRank, 0, m_StreamComm);
        MPI_Send(globalBuffer.data(), globalSize, MPI_CHAR, m_ReaderMasterStreamRank, 0,
                 m_StreamComm);
    }
}

void SscWriterNaive::PutDeferred(VariableBase &variable, const void *data)
{

    if (variable.m_Type == DataType::String)
    {
        const auto dataString = reinterpret_cast<const std::string *>(data);
        m_Metadata.emplace_back();
        auto &b = m_Metadata.back();
        b.name = variable.m_Name;
        b.type = DataType::String;
        b.shapeId = variable.m_ShapeID;
        b.shape = variable.m_Shape;
        b.start = variable.m_Start;
        b.count = variable.m_Count;
        b.bufferStart = m_Buffer.size();
        b.bufferCount = dataString->size();
        m_Buffer.resize(b.bufferStart + b.bufferCount);
        std::memcpy(m_Buffer.data() + b.bufferStart, dataString->data(), dataString->size());
        b.value.resize(dataString->size());
        std::memcpy(b.value.data(), dataString->data(), dataString->size());
        return;
    }

    if ((variable.m_ShapeID == ShapeID::GlobalValue || variable.m_ShapeID == ShapeID::LocalValue ||
         variable.m_Type == DataType::String) &&
        m_WriterRank != 0)
    {
        return;
    }

    Dims vStart = variable.m_Start;
    Dims vCount = variable.m_Count;
    Dims vShape = variable.m_Shape;

    if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
    {
        std::reverse(vStart.begin(), vStart.end());
        std::reverse(vCount.begin(), vCount.end());
        std::reverse(vShape.begin(), vShape.end());
    }

    bool found = false;
    for (const auto &b : m_Metadata)
    {
        if (b.name == variable.m_Name && ssc::AreSameDims(vStart, b.start) &&
            ssc::AreSameDims(vCount, b.count) && ssc::AreSameDims(vShape, b.shape))
        {
            std::memcpy(m_Buffer.data() + b.bufferStart, data, b.bufferCount);
            found = true;
        }
    }

    if (!found)
    {
        m_Metadata.emplace_back();
        auto &b = m_Metadata.back();
        b.name = variable.m_Name;
        b.type = variable.m_Type;
        b.shapeId = variable.m_ShapeID;
        b.shape = vShape;
        b.start = vStart;
        b.count = vCount;
        b.elementSize = variable.m_ElementSize;
        b.bufferStart = m_Buffer.size();
        b.bufferCount = ssc::TotalDataSize(b.count, b.elementSize, b.shapeId);
        m_Buffer.resize(b.bufferStart + b.bufferCount);
        std::memcpy(m_Buffer.data() + b.bufferStart, data, b.bufferCount);
        if (b.shapeId == ShapeID::GlobalValue || b.shapeId == ShapeID::LocalValue)
        {
            b.value.resize(variable.m_ElementSize);
            std::memcpy(b.value.data(), data, b.bufferCount);
        }
        if (variable.m_Type == DataType::Struct)
        {
            b.structDef = reinterpret_cast<VariableStruct *>(&variable)
                              ->m_WriteStructDefinition->StructName();
        }
    }
}

}
}
}
}
