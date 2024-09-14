/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * InlineWriter.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: Aron Helser aron.helser@kitware.com
 */

#include "InlineWriter.h"
#include "InlineReader.h"
#include "InlineWriter.tcc"

#include "adios2/helper/adiosFunctions.h"
#include <adios2-perfstubs-interface.h>

#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

InlineWriter::InlineWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("InlineWriter", io, name, mode, std::move(comm))
{
    PERFSTUBS_SCOPED_TIMER("InlineWriter::Open");
    m_WriterRank = m_Comm.Rank();
    Init();
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << " Open(" << m_Name << ")." << std::endl;
    }
    m_IsOpen = true;
}

InlineWriter::~InlineWriter()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

const InlineReader *InlineWriter::GetReader() const
{
    const auto &engine_map = m_IO.GetEngines();
    if (engine_map.size() == 1)
    {
        // it should be fine for a writer to be created and start running,
        // without the reader having been created. This is necessary to run
        // correctly with ParaView Catalyst Live.
        return nullptr;
    }
    else if (engine_map.size() > 2)
    {
        helper::Throw<std::runtime_error>("Engine", "InlineWriter", "GetReader",
                                          "There must be only one inline writer and at most "
                                          "one inline reader.");
    }

    std::shared_ptr<Engine> e = engine_map.begin()->second;
    if (e->OpenMode() == adios2::Mode::Write)
    {
        e = engine_map.rbegin()->second;
    }

    const auto reader = dynamic_cast<InlineReader *>(e.get());
    if (!reader)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "InlineWriter", "GetReader",
            "dynamic_cast<InlineReader*> failed; this is very likely a bug.");
    }
    return reader;
}

StepStatus InlineWriter::BeginStep(StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER("InlineWriter::BeginStep");
    if (m_InsideStep)
    {
        helper::Throw<std::runtime_error>("Engine", "InlineWriter", "BeginStep",
                                          "InlineWriter::BeginStep was called but the "
                                          "writer is already inside a step");
    }

    auto reader = GetReader();
    if (reader && reader->IsInsideStep())
    {
        m_InsideStep = false;
        return StepStatus::NotReady;
    }
    m_InsideStep = true;
    if (m_CurrentStep == static_cast<size_t>(-1))
    {
        m_CurrentStep = 0; // 0 is the first step
    }
    else
    {
        ++m_CurrentStep;
    }
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << "   BeginStep() new step " << m_CurrentStep
                  << "\n";
    }

    // m_BlocksInfo for all variables should be cleared at this point,
    // whether they were read in the last step or not.
    ResetVariables();

    return StepStatus::OK;
}

void InlineWriter::ResetVariables()
{
    auto availVars = m_IO.GetAvailableVariables();
    for (auto &varPair : availVars)
    {
        const auto &name = varPair.first;
        const DataType type = m_IO.InquireVariableType(name);

        if (type == DataType::Struct)
        {
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        Variable<T> &variable = FindVariable<T>(name, "in call to BeginStep");                     \
        variable.m_BlocksInfo.clear();                                                             \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    }
    m_ResetVariables = false;
}

size_t InlineWriter::CurrentStep() const { return m_CurrentStep; }

void InlineWriter::PerformPuts()
{
    PERFSTUBS_SCOPED_TIMER("InlineWriter::PerformPuts");
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << "     PerformPuts()\n";
    }
    m_ResetVariables = true;
}

void InlineWriter::EndStep()
{
    PERFSTUBS_SCOPED_TIMER("InlineWriter::EndStep");
    if (!m_InsideStep)
    {
        helper::Throw<std::runtime_error>("Engine", "InlineWriter", "EndStep",
                                          "InlineWriter::EndStep() cannot be called "
                                          "without a call to BeginStep() first");
    }
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << " EndStep() Step " << m_CurrentStep
                  << std::endl;
    }
    m_InsideStep = false;
}

void InlineWriter::Flush(const int)
{
    PERFSTUBS_SCOPED_TIMER("InlineWriter::Flush");
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << "   Flush()\n";
    }
}

bool InlineWriter::IsInsideStep() const { return m_InsideStep; }

// PRIVATE

#define declare_type(T)                                                                            \
    void InlineWriter::DoPutSync(Variable<T> &variable, const T *data)                             \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineWriter::DoPutSync");                                         \
        PutSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void InlineWriter::DoPutDeferred(Variable<T> &variable, const T *data)                         \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineWriter::DoPutDeferred");                                     \
        PutDeferredCommon(variable, data);                                                         \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void InlineWriter::Init()
{
    InitParameters();
    InitTransports();
}

void InlineWriter::InitParameters()
{
    for (const auto &pair : m_IO.m_Parameters)
    {
        std::string key(pair.first);
        std::transform(key.begin(), key.end(), key.begin(), ::tolower);

        std::string value(pair.second);

        if (key == "verbose")
        {
            m_Verbosity = std::stoi(value);
            if (m_Verbosity < 0 || m_Verbosity > 5)
                helper::Throw<std::invalid_argument>("Engine", "InlineWriter", "InitParameters",
                                                     "Method verbose argument must be an "
                                                     "integer in the range [0,5], in call to "
                                                     "Open or Engine constructor");
        }
    }
}

void InlineWriter::InitTransports()
{
    // Nothing to process from m_IO.m_TransportsParameters
}

void InlineWriter::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("InlineWriter::DoClose");
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << " Close(" << m_Name << ")\n";
    }
    // end of stream
    m_CurrentStep = static_cast<size_t>(-1);
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
