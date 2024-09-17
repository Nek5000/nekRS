/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * InlineReader.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: Aron Helser aron.helser@kitware.com
 */

#include "InlineReader.h"
#include "InlineReader.tcc"

#include "adios2/helper/adiosFunctions.h" // CSVToVector
#include <adios2-perfstubs-interface.h>

#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

InlineReader::InlineReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("InlineReader", io, name, mode, std::move(comm))
{
    PERFSTUBS_SCOPED_TIMER("InlineReader::Open");
    m_ReaderRank = m_Comm.Rank();
    Init();
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Reader " << m_ReaderRank << " Open(" << m_Name << ") in constructor"
                  << std::endl;
    }
    m_IsOpen = true;
}

InlineReader::~InlineReader()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

const InlineWriter *InlineReader::GetWriter() const
{
    const auto &engine_map = m_IO.GetEngines();
    if (engine_map.size() != 2)
    {
        helper::Throw<std::runtime_error>("Engine", "InlineReader", "GetWriter",
                                          "There must be exactly one reader and one "
                                          "writer for the inline engine.");
    }

    std::shared_ptr<Engine> e = engine_map.begin()->second;
    if (e->OpenMode() == adios2::Mode::Read)
    {
        e = engine_map.rbegin()->second;
    }

    const auto writer = dynamic_cast<InlineWriter *>(e.get());
    if (!writer)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "InlineReader", "GetWriter",
            "dynamic_cast<InlineWriter*> failed; this is very likely a bug.");
    }
    return writer;
}

StepStatus InlineReader::BeginStep(const StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER("InlineReader::BeginStep");
    if (m_InsideStep)
    {
        helper::Throw<std::runtime_error>("Engine", "InlineReader", "BeginStep",
                                          "InlineReader::BeginStep was called but the "
                                          "reader is already inside a step");
    }
    // Reader should be on step that writer just completed
    auto writer = GetWriter();
    if (writer->IsInsideStep())
    {
        m_InsideStep = false;
        return StepStatus::NotReady;
    }
    m_CurrentStep = writer->CurrentStep();
    if (m_CurrentStep == static_cast<size_t>(-1))
    {
        m_InsideStep = false;
        return StepStatus::EndOfStream;
    }
    m_InsideStep = true;

    if (m_Verbosity == 5)
    {
        std::cout << "Inline Reader " << m_ReaderRank << "   BeginStep() new step " << m_CurrentStep
                  << "\n";
    }

    return StepStatus::OK;
}

void InlineReader::PerformGets()
{
    PERFSTUBS_SCOPED_TIMER("InlineReader::PerformGets");
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Reader " << m_ReaderRank << "     PerformGets()\n";
    }
    SetDeferredVariablePointers();
}

size_t InlineReader::CurrentStep() const
{
    // Reader should be on same step as writer
    // added here since it's not really necessary to use beginstep/endstep for
    // this engine's reader so this ensures we do report the correct step
    const auto writer = GetWriter();
    return writer->CurrentStep();
}

void InlineReader::EndStep()
{
    PERFSTUBS_SCOPED_TIMER("InlineReader::EndStep");
    if (!m_InsideStep)
    {
        helper::Throw<std::runtime_error>("Engine", "InlineReader", "EndStep",
                                          "InlineReader::EndStep() cannot be called "
                                          "without a call to BeginStep() first");
    }
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Reader " << m_ReaderRank << " EndStep() Step " << m_CurrentStep
                  << std::endl;
    }
    if (!m_DeferredVariables.empty())
    {
        SetDeferredVariablePointers();
    }
    m_InsideStep = false;
}

bool InlineReader::IsInsideStep() const { return m_InsideStep; }

// PRIVATE

#define declare_type(T)                                                                            \
    void InlineReader::DoGetSync(Variable<T> &variable, T *data)                                   \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineReader::DoGetSync");                                         \
        GetSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void InlineReader::DoGetDeferred(Variable<T> &variable, T *data)                               \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineReader::DoGetDeferred");                                     \
        GetDeferredCommon(variable, data);                                                         \
    }                                                                                              \
    typename Variable<T>::BPInfo *InlineReader::DoGetBlockSync(Variable<T> &variable)              \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineReader::DoGetBlockSync");                                    \
        return GetBlockSyncCommon(variable);                                                       \
    }                                                                                              \
    typename Variable<T>::BPInfo *InlineReader::DoGetBlockDeferred(Variable<T> &variable)          \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineReader::DoGetBlockDeferred");                                \
        return GetBlockDeferredCommon(variable);                                                   \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

// Design note: Returns a copy. Instead, could return a reference, then
// Engine::Get() would not need an Info parameter passed in - binding could
// retrieve the current Core Info object at a later time.
// See note on binding Engine::BlocksInfo
#define declare_type(T)                                                                            \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>>                                    \
    InlineReader::DoAllStepsBlocksInfo(const Variable<T> &variable) const                          \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineReader::AllStepsBlockInfo");                                 \
        return std::map<size_t, std::vector<typename Variable<T>::BPInfo>>();                      \
    }                                                                                              \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> InlineReader::DoBlocksInfo(                          \
        const Variable<T> &variable, const size_t step) const                                      \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("InlineReader::DoBlocksInfo");                                      \
        return variable.m_BlocksInfo;                                                              \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void InlineReader::Init()
{
    InitParameters();
    InitTransports();
}

void InlineReader::InitParameters()
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
                helper::Throw<std::invalid_argument>("Engine", "InlineReader", "InitParameters",
                                                     "Method verbose argument must be an "
                                                     "integer in the range [0,5], in call to "
                                                     "Open or Engine constructor");
        }
    }
}

void InlineReader::InitTransports()
{
    // Nothing to process from m_IO.m_TransportsParameters
}

void InlineReader::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("InlineReader::DoClose");
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Reader " << m_ReaderRank << " Close(" << m_Name << ")\n";
    }
}

void InlineReader::SetDeferredVariablePointers()
{
    // need to set core::Variable::BPInfo::BufferP for each deferred variable
    // to the ptr stored in core::Variable::BPInfo::Data
    // this will make Variable::BPInfo::Data() work correctly for the user
    for (const auto &varName : m_DeferredVariables)
    {
        const DataType type = m_IO.InquireVariableType(varName);
        if (type == DataType::Struct)
        {
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        Variable<T> &variable = FindVariable<T>(varName, "in call to EndStep");                    \
        for (auto &info : variable.m_BlocksInfo)                                                   \
        {                                                                                          \
            info.BufferP = info.Data;                                                              \
        }                                                                                          \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    }
    m_DeferredVariables.clear();
}

#define declare_type(T) template void InlineReader::Get<T>(Variable<T> &, T **) const;
ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

} // end namespace engine
} // end namespace core
} // end namespace adios2
