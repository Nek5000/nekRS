/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Engine.cpp
 *
 *  Created on: Dec 19, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "Engine.h"
#include "Engine.tcc"

#include <stdexcept>
#include <utility>

#include "adios2/core/IO.h"

namespace adios2
{
namespace core
{

Engine::Engine(const std::string engineType, IO &io, const std::string &name, const Mode openMode,
               helper::Comm comm)
: m_EngineType(engineType), m_IO(io), m_Name(name), m_OpenMode(openMode), m_Comm(std::move(comm)),
  m_UserOptions(io.m_ADIOS.GetUserOptions())
{
    m_FailVerbose = (m_Comm.Rank() == 0);
}

Engine::~Engine()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
}

Engine::operator bool() const noexcept { return !m_IsClosed; }

IO &Engine::GetIO() noexcept { return m_IO; }

Mode Engine::OpenMode() const noexcept { return m_OpenMode; }

StepStatus Engine::BeginStep()
{
    if (m_OpenMode == Mode::Read)
    {
        return BeginStep(StepMode::Read, -1.0);
    }
    else
    {
        return BeginStep(StepMode::Append, -1.0);
    }
}

StepStatus Engine::BeginStep(StepMode mode, const float timeoutSeconds)
{
    ThrowUp("BeginStep");
    return StepStatus::OtherError;
}

size_t Engine::CurrentStep() const
{
    ThrowUp("CurrentStep");
    return 0;
}

void Engine::EndStep() { ThrowUp("EndStep"); }
void Engine::PerformPuts() { ThrowUp("PerformPuts"); }
void Engine::PerformGets() { ThrowUp("PerformGets"); }
void Engine::PerformDataWrite() { return; }

void Engine::Close(const int transportIndex)
{
    DoClose(transportIndex);

    m_IsOpen = false;

    if (transportIndex == -1)
    {
        m_Comm.Free("freeing comm in Engine " + m_Name + ", in call to Close");
        m_IsClosed = true;
    }
}

/**
 * Called if destructor is called on an open engine.  Should warn or take any
 * non-complex measure that might help recover.
 */
void Engine::DestructorClose(bool Verbose) noexcept
{
    if (Verbose)
    {
        std::cerr << "Engine \"" << m_Name << "\" destroyed without a prior Close()." << std::endl;
        std::cerr << "This may have negative consequences." << std::endl;
    }
};

void Engine::Flush(const int /*transportIndex*/) { ThrowUp("Flush"); }

size_t Engine::Steps() const { return DoSteps(); }

void Engine::LockWriterDefinitions() noexcept { m_WriterDefinitionsLocked = true; }

bool Engine::BetweenStepPairs() const { return m_BetweenStepPairs; }

void Engine::LockReaderSelections() noexcept { m_ReaderSelectionsLocked = true; }

size_t Engine::DebugGetDataBufferSize() const
{
    ThrowUp("DebugGetDataBufferSize");
    return 0;
}

void Engine::Put(VariableStruct &variable, const void *data, const Mode launch)
{
    CommonChecks(variable, data, {Mode::Write, Mode::Append}, "in call to Put");

    switch (launch)
    {
    case Mode::Deferred:
        DoPutStructDeferred(variable, data);
        break;
    case Mode::Sync:
        DoPutStructSync(variable, data);
        break;
    default:
        helper::Throw<std::invalid_argument>("Core", "Engine", "Put",
                                             "invalid launch Mode for variable " + variable.m_Name +
                                                 ", only Mode::Deferred and Mode::Sync are valid");
    }
}

void Engine::Get(VariableStruct &variable, void *data, const Mode launch)
{
    CommonChecks(variable, data, {Mode::Read, Mode::ReadRandomAccess}, "in call to Get");

    switch (launch)
    {
    case Mode::Deferred:
        DoGetStructDeferred(variable, data);
        break;
    case Mode::Sync:
        DoGetStructSync(variable, data);
        break;
    default:
        helper::Throw<std::invalid_argument>("Core", "Engine", "Get",
                                             "invalid launch Mode for variable " + variable.m_Name +
                                                 ", only Mode::Deferred and Mode::Sync are valid");
    }
}

void Engine::EnterComputationBlock() noexcept {}
void Engine::ExitComputationBlock() noexcept {}

// PROTECTED
void Engine::Init() {}
void Engine::InitParameters() {}
void Engine::InitTransports() {}

void Engine::NotifyEngineAttribute(std::string name, DataType type) noexcept {}

// if not overriden, default to name/type version
void Engine::NotifyEngineAttribute(std::string name, AttributeBase *attr, void *Data) noexcept
{
    NotifyEngineAttribute(name, attr->m_Type);
}

void Engine::NotifyEngineNoVarsQuery() {}

// DoPut*
#define declare_type(T)                                                                            \
    void Engine::DoPut(Variable<T> &, typename Variable<T>::Span &, const bool, const T &)         \
    {                                                                                              \
        ThrowUp("DoPut");                                                                          \
    }
ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

#define declare_type(T)                                                                            \
    void Engine::DoPutSync(Variable<T> &, const T *) { ThrowUp("DoPutSync"); }                     \
    void Engine::DoPutDeferred(Variable<T> &, const T *) { ThrowUp("DoPutDeferred"); }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void Engine::DoPutStructSync(VariableStruct &, const void *) { ThrowUp("DoPutStructSync"); }
void Engine::DoPutStructDeferred(VariableStruct &, const void *) { ThrowUp("DoPutStructDeferred"); }

// DoGet*
#define declare_type(T)                                                                            \
    void Engine::DoGetSync(Variable<T> &, T *) { ThrowUp("DoGetSync"); }                           \
    void Engine::DoGetDeferred(Variable<T> &, T *) { ThrowUp("DoGetDeferred"); }                   \
    typename Variable<T>::BPInfo *Engine::DoGetBlockSync(Variable<T> &v) { return nullptr; }       \
    typename Variable<T>::BPInfo *Engine::DoGetBlockDeferred(Variable<T> &v) { return nullptr; }

void Engine::RegisterCreatedVariable(const VariableBase *var) { m_CreatedVars.insert(var); }

void Engine::RemoveCreatedVars()
{
    for (auto &VarRec : m_CreatedVars)
    {
        m_IO.RemoveVariable(VarRec->m_Name);
    }
    m_CreatedVars.clear();
}

void Engine::DoGetAbsoluteSteps(const VariableBase &variable, std::vector<size_t> &keys) const
{
    ThrowUp("DoGetAbsoluteSteps");
    return;
}

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void Engine::DoGetStructSync(VariableStruct &, void *) { ThrowUp("DoGetSync for Struct Variable"); }
void Engine::DoGetStructDeferred(VariableStruct &, void *)
{
    ThrowUp("DoGetDeferred for Struct Variable");
}

#define declare_type(T)                                                                            \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> Engine::DoAllStepsBlocksInfo(      \
        const Variable<T> &variable) const                                                         \
    {                                                                                              \
        ThrowUp("DoAllStepsBlocksInfo");                                                           \
        return std::map<size_t, std::vector<typename Variable<T>::BPInfo>>();                      \
    }                                                                                              \
                                                                                                   \
    std::vector<std::vector<typename Variable<T>::BPInfo>> Engine::DoAllRelativeStepsBlocksInfo(   \
        const Variable<T> &variable) const                                                         \
    {                                                                                              \
        ThrowUp("DoAllRelativeStepsBlocksInfo");                                                   \
        return std::vector<std::vector<typename Variable<T>::BPInfo>>();                           \
    }                                                                                              \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> Engine::DoBlocksInfo(const Variable<T> &variable,    \
                                                                   const size_t step) const        \
    {                                                                                              \
        ThrowUp("DoBlocksInfo");                                                                   \
        return std::vector<typename Variable<T>::BPInfo>();                                        \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

std::map<size_t, std::vector<VariableStruct::BPInfo>>
Engine::DoAllStepsBlocksInfoStruct(const VariableStruct &variable) const
{
    ThrowUp("DoAllStepsBlocksInfo");
    return std::map<size_t, std::vector<VariableStruct::BPInfo>>();
}

std::vector<std::vector<VariableStruct::BPInfo>>
Engine::DoAllRelativeStepsBlocksInfoStruct(const VariableStruct &variable) const
{
    ThrowUp("DoAllRelativeStepsBlocksInfo");
    return std::vector<std::vector<VariableStruct::BPInfo>>();
}

std::vector<VariableStruct::BPInfo> Engine::DoBlocksInfoStruct(const VariableStruct &variable,
                                                               const size_t step) const
{
    ThrowUp("DoBlocksInfo");
    return std::vector<VariableStruct::BPInfo>();
}

#define declare_type(T, L)                                                                         \
    T *Engine::DoBufferData_##L(const int bufferIdx, const size_t payloadPosition,                 \
                                const size_t bufferID) noexcept                                    \
    {                                                                                              \
        T *data = nullptr;                                                                         \
        return data;                                                                               \
    }

ADIOS2_FOREACH_PRIMITVE_STDTYPE_2ARGS(declare_type)
#undef declare_type

size_t Engine::DoSteps() const
{
    ThrowUp("DoSteps");
    return MaxSizeT;
}

// PRIVATE
void Engine::ThrowUp(const std::string function) const
{
    helper::Throw<std::invalid_argument>(
        "Core", "Engine", "ThrowUp", "Engine " + m_EngineType + " does not support " + function);
}

void Engine::CheckOpenModes(const std::set<Mode> &modes, const std::string hint) const
{
    if (modes.count(m_OpenMode) == 0)
    {
        helper::Throw<std::invalid_argument>("Core", "Engine", "CheckOpenModes",
                                             "Engine open mode not valid for " + hint);
    }
}

void Engine::CommonChecks(VariableBase &variable, const void *data, const std::set<Mode> &modes,
                          const std::string hint) const
{
    variable.CheckDimensions(hint);
    CheckOpenModes(modes, " for variable " + variable.m_Name + ", " + hint);

    // If no dimension has a zero count then there must be data to write.
    if (std::find(variable.m_Count.begin(), variable.m_Count.end(), 0) == variable.m_Count.end())
    {
        helper::CheckForNullptr(data, "for data argument in non-zero count block, " + hint);
    }
}

std::map<size_t, std::vector<VariableStruct::BPInfo>>
Engine::AllStepsBlocksInfoStruct(const VariableStruct &variable) const
{
    return DoAllStepsBlocksInfoStruct(variable);
}

std::vector<std::vector<VariableStruct::BPInfo>>
Engine::AllRelativeStepsBlocksInfoStruct(const VariableStruct &variable) const
{
    return DoAllRelativeStepsBlocksInfoStruct(variable);
}

std::vector<VariableStruct::BPInfo> Engine::BlocksInfoStruct(const VariableStruct &variable,
                                                             const size_t step) const
{
    return DoBlocksInfoStruct(variable, step);
}

// PUBLIC TEMPLATE FUNCTIONS EXPANSION WITH SCOPED TYPES
#define declare_template_instantiation(T)                                                          \
                                                                                                   \
    template void Engine::Put<T>(Variable<T> &, const T *, const Mode);                            \
    template void Engine::Put<T>(const std::string &, const T *, const Mode);                      \
                                                                                                   \
    template void Engine::Put<T>(Variable<T> &, const T &, const Mode);                            \
    template void Engine::Put<T>(const std::string &, const T &, const Mode);                      \
                                                                                                   \
    template void Engine::Get<T>(Variable<T> &, T *, const Mode);                                  \
    template void Engine::Get<T>(const std::string &, T *, const Mode);                            \
                                                                                                   \
    template void Engine::Get<T>(Variable<T> &, T &, const Mode);                                  \
    template void Engine::Get<T>(const std::string &, T &, const Mode);                            \
                                                                                                   \
    template void Engine::Get<T>(Variable<T> &, std::vector<T> &, const Mode);                     \
    template void Engine::Get<T>(const std::string &, std::vector<T> &, const Mode);               \
                                                                                                   \
    template typename Variable<T>::BPInfo *Engine::Get<T>(Variable<T> &, const Mode);              \
    template typename Variable<T>::BPInfo *Engine::Get<T>(const std::string &, const Mode);        \
                                                                                                   \
    template Variable<T> &Engine::FindVariable(const std::string &variableName,                    \
                                               const std::string hint);                            \
                                                                                                   \
    template std::map<size_t, std::vector<typename Variable<T>::BPInfo>>                           \
    Engine::AllStepsBlocksInfo(const Variable<T> &) const;                                         \
                                                                                                   \
    template std::vector<std::vector<typename Variable<T>::BPInfo>>                                \
    Engine::AllRelativeStepsBlocksInfo(const Variable<T> &) const;                                 \
                                                                                                   \
    template std::vector<typename Variable<T>::BPInfo> Engine::BlocksInfo(const Variable<T> &,     \
                                                                          const size_t) const;     \
    template std::vector<size_t> Engine::GetAbsoluteSteps(const Variable<T> &) const;

ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

#define declare_template_instantiation(T)                                                          \
    template typename Variable<T>::Span &Engine::Put(Variable<T> &, const bool, const T &);        \
    template void Engine::Get<T>(core::Variable<T> &, T **) const;

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace core
} // end namespace adios2
