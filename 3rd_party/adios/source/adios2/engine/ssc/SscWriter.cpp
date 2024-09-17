/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscWriter.cpp
 *
 *  Created on: Nov 1, 2018
 *      Author: Jason Wang
 */

#include "SscWriter.h"
#include "SscWriterGeneric.h"
#include "SscWriterNaive.h"
#include "adios2/helper/adiosCommMPI.h"
#include "adios2/helper/adiosString.h"
#include <adios2-perfstubs-interface.h>

namespace adios2
{
namespace core
{
namespace engine
{

SscWriter::SscWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("SscWriter", io, name, mode, std::move(comm))
{
    PERFSTUBS_SCOPED_TIMER_FUNC();

    helper::GetParameter(m_IO.m_Parameters, "EngineMode", m_EngineMode);
    helper::GetParameter(m_IO.m_Parameters, "Verbose", m_Verbosity);

    helper::Log("Engine", "SscWriter", "SscWriter", m_EngineMode,
                m_Verbosity >= 10 ? m_Comm.Rank() : 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    if (m_EngineMode == "generic")
    {
        m_EngineInstance =
            std::make_shared<ssc::SscWriterGeneric>(io, name, mode, CommAsMPI(m_Comm));
    }
    else if (m_EngineMode == "naive")
    {
        m_EngineInstance = std::make_shared<ssc::SscWriterNaive>(io, name, mode, CommAsMPI(m_Comm));
    }
    m_IsOpen = true;
}

SscWriter::~SscWriter()
{
    if (adios2::core::Engine::m_IsOpen)
    {
        DestructorClose(adios2::core::Engine::m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus SscWriter::BeginStep(StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();

    auto ret = m_EngineInstance->BeginStep(mode, timeoutSeconds, m_WriterDefinitionsLocked);

    helper::Log("Engine", "SscWriter", "BeginStep", std::to_string(CurrentStep()),
                m_Verbosity >= 10 ? m_Comm.Rank() : 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    return ret;
}

size_t SscWriter::CurrentStep() const { return m_EngineInstance->CurrentStep(); }

void SscWriter::PerformPuts()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    m_EngineInstance->PerformPuts();
}

void SscWriter::EndStep()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();

    helper::Log("Engine", "SscWriter", "EndStep", std::to_string(CurrentStep()),
                m_Verbosity >= 10 ? m_Comm.Rank() : 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    m_EngineInstance->EndStep(m_WriterDefinitionsLocked);
}

void SscWriter::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();

    helper::Log("Engine", "SscWriter", "Close", m_Name, m_Verbosity >= 10 ? m_Comm.Rank() : 0,
                m_Comm.Rank(), 5, m_Verbosity, helper::LogMode::INFO);

    m_EngineInstance->Close(transportIndex);
}

#define declare_type(T)                                                                            \
    void SscWriter::DoPutSync(Variable<T> &variable, const T *data)                                \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER_FUNC();                                                             \
        helper::Log("Engine", "SscWriter", "DoPutSync", variable.m_Name,                           \
                    m_Verbosity >= 10 ? m_Comm.Rank() : 0, m_Comm.Rank(), 5, m_Verbosity,          \
                    helper::LogMode::INFO);                                                        \
        m_EngineInstance->PutDeferred(variable, data);                                             \
        m_EngineInstance->PerformPuts();                                                           \
    }                                                                                              \
    void SscWriter::DoPutDeferred(Variable<T> &variable, const T *data)                            \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER_FUNC();                                                             \
        helper::Log("Engine", "SscWriter", "DoPutDeferred", variable.m_Name,                       \
                    m_Verbosity >= 10 ? m_Comm.Rank() : 0, m_Comm.Rank(), 5, m_Verbosity,          \
                    helper::LogMode::INFO);                                                        \
        m_EngineInstance->PutDeferred(variable, data);                                             \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void SscWriter::DoPutStructSync(VariableStruct &variable, const void *data)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    helper::Log("Engine", "SscWriter", "DoPutSync", variable.m_Name,
                m_Verbosity >= 10 ? m_Comm.Rank() : 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);
    m_EngineInstance->PutDeferred(variable, data);
    m_EngineInstance->PerformPuts();
}
void SscWriter::DoPutStructDeferred(VariableStruct &variable, const void *data)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    helper::Log("Engine", "SscWriter", "DoPutDeferred", variable.m_Name,
                m_Verbosity >= 10 ? m_Comm.Rank() : 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);
    m_EngineInstance->PutDeferred(variable, data);
}

void SscWriter::Flush(const int transportIndex) {}

} // end namespace engine
} // end namespace core
} // end namespace adios2
