/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP3Reader.cpp
 *
 *  Created on: Feb 27, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BP3Reader.h"
#include "BP3Reader.tcc"

#include "adios2/helper/adiosComm.h"
#include <adios2-perfstubs-interface.h>

namespace adios2
{
namespace core
{
namespace engine
{

BP3Reader::BP3Reader(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("BP3", io, name, mode, std::move(comm)), m_BP3Deserializer(m_Comm),
  m_FileManager(io, m_Comm), m_SubFileManager(io, m_Comm)
{
    PERFSTUBS_SCOPED_TIMER("BP3Reader::Open");
    Init();
    m_IsOpen = true;
}

BP3Reader::~BP3Reader()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus BP3Reader::BeginStep(StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER("BP3Reader::BeginStep");
    if (mode != StepMode::Read)
    {
        helper::Throw<std::invalid_argument>("Engine", "BP3Reader", "BeginStep",
                                             "mode is not supported yet, "
                                             "only Read is valid for "
                                             "engine BP3 with adios2::Mode::Read, in call to "
                                             "BeginStep");
    }

    if (!m_BP3Deserializer.m_DeferredVariables.empty())
    {
        helper::Throw<std::invalid_argument>("Engine", "BP3Reader", "BeginStep",
                                             "existing variables subscribed with "
                                             "GetDeferred, did you forget to call "
                                             "PerformGets() or EndStep()?, in call to BeginStep");
    }

    if (m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "BP3Reader", "BeginStep",
                                        "BeginStep() is called a second time "
                                        "without an intervening EndStep()");
    }

    m_BetweenStepPairs = true;
    if (m_FirstStep)
    {
        m_FirstStep = false;
    }
    else
    {
        ++m_CurrentStep;
    }

    // used to inquire for variables in streaming mode
    m_IO.m_ReadStreaming = true;
    m_IO.m_EngineStep = m_CurrentStep;

    if (m_CurrentStep >= m_BP3Deserializer.m_MetadataSet.StepsCount)
    {
        m_IO.m_ReadStreaming = false;
        return StepStatus::EndOfStream;
    }

    m_IO.ResetVariablesStepSelection(false, "in call to BP3 Reader BeginStep");

    return StepStatus::OK;
}

size_t BP3Reader::CurrentStep() const { return m_CurrentStep; }

void BP3Reader::EndStep()
{
    if (!m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "BP3Reader", "EndStep",
                                        "EndStep() is called without a successful BeginStep()");
    }
    m_BetweenStepPairs = false;
    PERFSTUBS_SCOPED_TIMER("BP3Reader::EndStep");
    PerformGets();
}

void BP3Reader::PerformGets()
{
    PERFSTUBS_SCOPED_TIMER("BP3Reader::PerformGets");
    if (m_BP3Deserializer.m_DeferredVariables.empty())
    {
        return;
    }

    for (const std::string &name : m_BP3Deserializer.m_DeferredVariables)
    {
        const DataType type = m_IO.InquireVariableType(name);

        if (type == DataType::Struct)
        {
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        Variable<T> &variable = FindVariable<T>(name, "in call to PerformGets, EndStep or Close"); \
        for (auto &blockInfo : variable.m_BlocksInfo)                                              \
        {                                                                                          \
            m_BP3Deserializer.SetVariableBlockInfo(variable, blockInfo);                           \
        }                                                                                          \
        ReadVariableBlocks(variable);                                                              \
        variable.m_BlocksInfo.clear();                                                             \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    }

    m_BP3Deserializer.m_DeferredVariables.clear();
}

// PRIVATE
void BP3Reader::Init()
{
    if (m_OpenMode != Mode::Read)
    {
        helper::Throw<std::invalid_argument>("Engine", "BP3Reader", "Init",
                                             "BPFileReader only supports OpenMode::Read from" +
                                                 m_Name);
    }

    // if IO was involved in reading before this flag may be true now
    m_IO.m_ReadStreaming = false;

    helper::RaiseLimitNoFile();
    InitTransports();
    InitBuffer();
}

void BP3Reader::InitTransports()
{
    if (m_IO.m_TransportsParameters.empty())
    {
        Params defaultTransportParameters;
        defaultTransportParameters["transport"] = "File";
        m_IO.m_TransportsParameters.push_back(defaultTransportParameters);
    }
    // TODO Set Parameters

    if (m_BP3Deserializer.m_RankMPI == 0)
    {
        const bool profile = m_BP3Deserializer.m_Profiler.m_IsActive;
        try
        {
            m_FileManager.OpenFiles({m_Name}, adios2::Mode::Read, m_IO.m_TransportsParameters,
                                    profile);
        }
        catch (...)
        {
            const std::string bpName = helper::AddExtension(m_Name, ".bp");
            m_FileManager.OpenFiles({bpName}, adios2::Mode::Read, m_IO.m_TransportsParameters,
                                    profile);
        }
    }
}

void BP3Reader::InitBuffer()
{
    if (m_BP3Deserializer.m_RankMPI == 0)
    {
        const size_t fileSize = m_FileManager.GetFileSize();
        // handle single bp files from ADIOS 1.x by getting onl the metadata in
        // buffer

        // Load/Read Minifooter
        const size_t miniFooterSize = m_BP3Deserializer.m_MetadataSet.MiniFooterSize;
        if (fileSize < miniFooterSize)
        {
            std::string err = "The size of the input file " + m_Name + "(" +
                              std::to_string(fileSize) +
                              " bytes) is less than the minimum BP3 header "
                              "size, which is " +
                              std::to_string(miniFooterSize) + " bytes." +
                              " It is unlikely that this is a .bp file.";
            helper::Throw<std::logic_error>("Engine", "BP3Reader", "Init", err);
        }
        const size_t miniFooterStart = helper::GetDistance(
            fileSize, miniFooterSize, " fileSize < miniFooterSize, in call to Open");

        m_BP3Deserializer.m_Metadata.Resize(
            miniFooterSize, "allocating metadata buffer to inspect bp minifooter, in call to "
                            "Open");

        m_FileManager.ReadFile(m_BP3Deserializer.m_Metadata.m_Buffer.data(), miniFooterSize,
                               miniFooterStart);

        // Load/Read Metadata
        const size_t metadataStart = m_BP3Deserializer.MetadataStart(m_BP3Deserializer.m_Metadata);
        const size_t metadataSize = helper::GetDistance(
            fileSize, metadataStart, " fileSize < miniFooterSize, in call to Open");

        m_BP3Deserializer.m_Metadata.Resize(metadataSize,
                                            "allocating metadata buffer, in call to Open");

        m_FileManager.ReadFile(m_BP3Deserializer.m_Metadata.m_Buffer.data(), metadataSize,
                               metadataStart);
    }

    // broadcast metadata buffer to all ranks from zero
    m_Comm.BroadcastVector(m_BP3Deserializer.m_Metadata.m_Buffer);

    // fills IO with available Variables and Attributes
    m_BP3Deserializer.ParseMetadata(m_BP3Deserializer.m_Metadata, *this);
    // caches attributes associated with variables
    m_IO.SetPrefixedNames(false);
}

#define declare_type(T)                                                                            \
    void BP3Reader::DoGetSync(Variable<T> &variable, T *data)                                      \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Reader::Get");                                                  \
        GetSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void BP3Reader::DoGetDeferred(Variable<T> &variable, T *data)                                  \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Reader::Get");                                                  \
        GetDeferredCommon(variable, data);                                                         \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void BP3Reader::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP3Reader::Close");
    PerformGets();
    RemoveCreatedVars();
    m_SubFileManager.CloseFiles();
    m_FileManager.CloseFiles();
}

#define declare_type(T)                                                                            \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> BP3Reader::DoAllStepsBlocksInfo(   \
        const Variable<T> &variable) const                                                         \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Reader::AllStepsBlocksInfo");                                   \
        return m_BP3Deserializer.AllStepsBlocksInfo(variable);                                     \
    }                                                                                              \
                                                                                                   \
    std::vector<std::vector<typename Variable<T>::BPInfo>>                                         \
    BP3Reader::DoAllRelativeStepsBlocksInfo(const Variable<T> &variable) const                     \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Reader::AllRelativeStepsBlocksInfo");                           \
        return m_BP3Deserializer.AllRelativeStepsBlocksInfo(variable);                             \
    }                                                                                              \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> BP3Reader::DoBlocksInfo(const Variable<T> &variable, \
                                                                      const size_t step) const     \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Reader::BlocksInfo");                                           \
        return m_BP3Deserializer.BlocksInfo(variable, step);                                       \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

size_t BP3Reader::DoSteps() const { return m_BP3Deserializer.m_MetadataSet.StepsCount; }

} // end namespace engine
} // end namespace core
} // end namespace adios2
