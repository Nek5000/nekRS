/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP3Writer.cpp
 *
 *  Created on: Dec 19, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BP3Writer.h"
#include "BP3Writer.tcc"

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosFunctions.h" //CheckIndexRange
#include "adios2/toolkit/transport/file/FileFStream.h"

#include <adios2-perfstubs-interface.h>

namespace adios2
{
namespace core
{
namespace engine
{

BP3Writer::BP3Writer(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("BP3", io, name, mode, std::move(comm)), m_BP3Serializer(m_Comm),
  m_FileDataManager(io, m_Comm), m_FileMetadataManager(io, m_Comm)
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::Open");
    m_IO.m_ReadStreaming = false;
    Init();
    m_IsOpen = true;
}

BP3Writer::~BP3Writer()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus BP3Writer::BeginStep(StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::BeginStep");
    m_BP3Serializer.m_DeferredVariables.clear();
    m_BP3Serializer.m_DeferredVariablesDataSize = 0;
    m_IO.m_ReadStreaming = false;
    m_DidBeginStep = true;
    return StepStatus::OK;
}

size_t BP3Writer::CurrentStep() const { return m_BP3Serializer.m_MetadataSet.CurrentStep; }

void BP3Writer::PerformPuts()
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::PerformPuts");
    if (m_BP3Serializer.m_DeferredVariables.empty())
    {
        return;
    }

    m_BP3Serializer.ResizeBuffer(m_BP3Serializer.m_DeferredVariablesDataSize,
                                 "in call to PerformPuts");

    for (const std::string &variableName : m_BP3Serializer.m_DeferredVariables)
    {
        const DataType type = m_IO.InquireVariableType(variableName);
        if (type == DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        Variable<T> &variable =                                                                    \
            FindVariable<T>(variableName, "in call to PerformPuts, EndStep or Close");             \
        PerformPutCommon(variable);                                                                \
    }

        ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }
    m_BP3Serializer.m_DeferredVariables.clear();
    m_BP3Serializer.m_DeferredVariablesDataSize = 0;
}

void BP3Writer::EndStep()
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::EndStep");
    if (m_BP3Serializer.m_DeferredVariables.size() > 0)
    {
        PerformPuts();
    }

    // true: advances step
    m_BP3Serializer.SerializeData(m_IO, true);

    const size_t currentStep = CurrentStep();
    const size_t flushStepsCount = m_BP3Serializer.m_Parameters.FlushStepsCount;

    if (currentStep % flushStepsCount == 0)
    {
        Flush();
    }

    if (m_BP3Serializer.m_RankMPI == 0)
    {
        m_IO.m_ADIOS.RecordOutputStep(m_Name, UnknownStep, UnknownTime);
    }
}

void BP3Writer::Flush(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::Flush");
    DoFlush(false, transportIndex);
    m_BP3Serializer.ResetBuffer(m_BP3Serializer.m_Data);

    if (m_BP3Serializer.m_Parameters.CollectiveMetadata)
    {
        WriteCollectiveMetadataFile();
    }
}

// PRIVATE
void BP3Writer::Init()
{
    InitParameters();
    if (m_BP3Serializer.m_Parameters.NumAggregators <
        static_cast<unsigned int>(m_BP3Serializer.m_SizeMPI))
    {
        m_BP3Serializer.m_Aggregator.Init(m_BP3Serializer.m_Parameters.NumAggregators,
                                          m_BP3Serializer.m_Parameters.NumAggregators, m_Comm);
    }
    InitTransports();
    InitBPBuffer();
}

#define declare_type(T)                                                                            \
    void BP3Writer::DoPut(Variable<T> &variable, typename Variable<T>::Span &span,                 \
                          const bool initialize, const T &value)                                   \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Writer::Put");                                                  \
        PutCommon(variable, span, 0, value);                                                       \
    }

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

#define declare_type(T)                                                                            \
    void BP3Writer::DoPutSync(Variable<T> &variable, const T *data)                                \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Writer::Put");                                                  \
        PutSyncCommon(variable, variable.SetBlockInfo(data, CurrentStep()));                       \
        variable.m_BlocksInfo.pop_back();                                                          \
    }                                                                                              \
    void BP3Writer::DoPutDeferred(Variable<T> &variable, const T *data)                            \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP3Writer::Put");                                                  \
        PutDeferredCommon(variable, data);                                                         \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void BP3Writer::InitParameters()
{
    m_BP3Serializer.Init(m_IO.m_Parameters, "in call to BP3::Open for writing");
    m_BP3Serializer.ResizeBuffer(m_BP3Serializer.m_Parameters.InitialBufferSize,
                                 "in call to BP3::Open to write");
}

void BP3Writer::InitTransports()
{
    // TODO need to add support for aggregators here later
    if (m_IO.m_TransportsParameters.empty())
    {
        Params defaultTransportParameters;
        defaultTransportParameters["transport"] = "File";
        m_IO.m_TransportsParameters.push_back(defaultTransportParameters);
    }

    // only consumers will interact with transport managers
    std::vector<std::string> bpSubStreamNames;

    if (m_BP3Serializer.m_Aggregator.m_IsAggregator)
    {
        // Names passed to IO AddTransport option with key "Name"
        const std::vector<std::string> transportsNames =
            m_FileDataManager.GetFilesBaseNames(m_Name, m_IO.m_TransportsParameters);

        // /path/name.bp.dir/name.bp.rank
        bpSubStreamNames = m_BP3Serializer.GetBPSubStreamNames(transportsNames);
    }

    m_BP3Serializer.m_Profiler.Start("mkdir");
    m_FileDataManager.MkDirsBarrier(bpSubStreamNames, m_IO.m_TransportsParameters,
                                    m_BP3Serializer.m_Parameters.NodeLocal);
    m_BP3Serializer.m_Profiler.Stop("mkdir");

    if (m_BP3Serializer.m_Aggregator.m_IsAggregator)
    {
        if (m_BP3Serializer.m_Parameters.AsyncOpen)
        {
            for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
            {
                m_IO.m_TransportsParameters[i]["asyncopen"] = "true";
            }
        }
        m_FileDataManager.OpenFiles(bpSubStreamNames, m_OpenMode, m_IO.m_TransportsParameters,
                                    m_BP3Serializer.m_Profiler.m_IsActive);
    }
}

void BP3Writer::InitBPBuffer()
{
    if (m_OpenMode == Mode::Append)
    {
        helper::Throw<std::invalid_argument>("Engine", "BP3Writer", "InitBPBuffer",
                                             "Mode::Append is only available in "
                                             "BP4; it is not implemented "
                                             "for BP3 files.");
        // TODO: Get last pg timestep and update timestep counter in
    }
    else
    {
        m_BP3Serializer.PutProcessGroupIndex(
            m_IO.m_Name, (m_IO.m_ArrayOrder == ArrayOrdering::RowMajor) ? "C++" : "Fortran",
            m_FileDataManager.GetTransportsTypes());
    }
}

void BP3Writer::DoFlush(const bool isFinal, const int transportIndex)
{
    if (m_BP3Serializer.m_Aggregator.m_IsActive)
    {
        AggregateWriteData(isFinal, transportIndex);
    }
    else
    {
        WriteData(isFinal, transportIndex);
    }
}

void BP3Writer::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::Close");
    if (m_BP3Serializer.m_DeferredVariables.size() > 0)
    {
        PerformPuts();
    }

    DoFlush(true, transportIndex);

    if (m_BP3Serializer.m_Aggregator.m_IsAggregator)
    {
        m_FileDataManager.CloseFiles(transportIndex);
    }

    if (m_BP3Serializer.m_Parameters.CollectiveMetadata && m_FileDataManager.AllTransportsClosed())
    {
        WriteCollectiveMetadataFile(true);
    }

    if (m_BP3Serializer.m_Profiler.m_IsActive && m_FileDataManager.AllTransportsClosed())
    {
        WriteProfilingJSONFile();
    }

    m_BP3Serializer.DeleteBuffers();

    if (!m_DidBeginStep && m_BP3Serializer.m_RankMPI == 0)
    {
        m_IO.m_ADIOS.RecordOutputStep(m_Name, UnknownStep, UnknownTime);
    }
}

void BP3Writer::WriteProfilingJSONFile()
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::WriteProfilingJSONFile");
    auto transportTypes = m_FileDataManager.GetTransportsTypes();

    // find first File type output, where we can write the profile
    int fileTransportIdx = -1;
    for (size_t i = 0; i < transportTypes.size(); ++i)
    {
        if (transportTypes[i].compare(0, 4, "File") == 0)
        {
            fileTransportIdx = static_cast<int>(i);
        }
    }

    auto transportProfilers = m_FileDataManager.GetTransportsProfilers();

    auto transportTypesMD = m_FileMetadataManager.GetTransportsTypes();
    auto transportProfilersMD = m_FileMetadataManager.GetTransportsProfilers();

    transportTypes.insert(transportTypes.end(), transportTypesMD.begin(), transportTypesMD.end());

    transportProfilers.insert(transportProfilers.end(), transportProfilersMD.begin(),
                              transportProfilersMD.end());

    const std::string lineJSON(
        m_BP3Serializer.GetRankProfilingJSON(transportTypes, transportProfilers) + ",\n");

    const std::vector<char> profilingJSON(m_BP3Serializer.AggregateProfilingJSON(lineJSON));

    if (m_BP3Serializer.m_RankMPI == 0)
    {
        transport::FileFStream profilingJSONStream(m_Comm);
        std::string profileFileName;
        if (fileTransportIdx > -1)
        {
            // write profile to <filename.bp>.dir/profiling.json
            auto bpBaseNames = m_BP3Serializer.GetBPBaseNames({m_Name});
            profileFileName = bpBaseNames[fileTransportIdx] + "/profiling.json";
        }
        else
        {
            // write profile to <filename.bp>_profiling.json
            auto transportsNames =
                m_FileMetadataManager.GetFilesBaseNames(m_Name, m_IO.m_TransportsParameters);

            auto bpMetadataFileNames = m_BP3Serializer.GetBPMetadataFileNames(transportsNames);
            profileFileName = bpMetadataFileNames[0] + "_profiling.json";
        }
        profilingJSONStream.Open(profileFileName, Mode::Write);
        profilingJSONStream.Write(profilingJSON.data(), profilingJSON.size());
        profilingJSONStream.Close();
    }
}

void BP3Writer::WriteCollectiveMetadataFile(const bool isFinal)
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::WriteCollectiveMetadataFile");
    m_BP3Serializer.AggregateCollectiveMetadata(m_Comm, m_BP3Serializer.m_Metadata, true);

    if (m_BP3Serializer.m_RankMPI == 0)
    {
        // first init metadata files
        const std::vector<std::string> transportsNames =
            m_FileMetadataManager.GetFilesBaseNames(m_Name, m_IO.m_TransportsParameters);

        const std::vector<std::string> bpMetadataFileNames =
            m_BP3Serializer.GetBPMetadataFileNames(transportsNames);

        m_FileMetadataManager.OpenFiles(bpMetadataFileNames, m_OpenMode,
                                        m_IO.m_TransportsParameters,
                                        m_BP3Serializer.m_Profiler.m_IsActive);

        m_FileMetadataManager.WriteFiles(m_BP3Serializer.m_Metadata.m_Buffer.data(),
                                         m_BP3Serializer.m_Metadata.m_Position);
        m_FileMetadataManager.CloseFiles();

        if (!isFinal)
        {
            m_BP3Serializer.ResetBuffer(m_BP3Serializer.m_Metadata, true);
            m_FileMetadataManager.m_Transports.clear();
        }
    }
}

void BP3Writer::WriteData(const bool isFinal, const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::WriteData");
    size_t dataSize = m_BP3Serializer.m_Data.m_Position;

    if (isFinal)
    {
        m_BP3Serializer.CloseData(m_IO);
        dataSize = m_BP3Serializer.m_Data.m_Position;
    }
    else
    {
        m_BP3Serializer.CloseStream(m_IO);
    }

    m_FileDataManager.WriteFiles(m_BP3Serializer.m_Data.m_Buffer.data(), dataSize, transportIndex);

    m_FileDataManager.FlushFiles(transportIndex);
}

void BP3Writer::AggregateWriteData(const bool isFinal, const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP3Writer::AggregateWriteData");
    m_BP3Serializer.CloseStream(m_IO, false);

    // async?
    for (int r = 0; r < m_BP3Serializer.m_Aggregator.m_Size; ++r)
    {
        aggregator::MPIChain::ExchangeRequests dataRequests =
            m_BP3Serializer.m_Aggregator.IExchange(m_BP3Serializer.m_Data, r);

        aggregator::MPIChain::ExchangeAbsolutePositionRequests absolutePositionRequests =
            m_BP3Serializer.m_Aggregator.IExchangeAbsolutePosition(m_BP3Serializer.m_Data, r);

        if (m_BP3Serializer.m_Aggregator.m_IsAggregator)
        {
            const format::Buffer &bufferSTL =
                m_BP3Serializer.m_Aggregator.GetConsumerBuffer(m_BP3Serializer.m_Data);

            m_FileDataManager.WriteFiles(bufferSTL.Data(), bufferSTL.m_Position, transportIndex);

            m_FileDataManager.FlushFiles(transportIndex);
        }

        m_BP3Serializer.m_Aggregator.WaitAbsolutePosition(absolutePositionRequests, r);

        m_BP3Serializer.m_Aggregator.Wait(dataRequests, r);
        m_BP3Serializer.m_Aggregator.SwapBuffers(r);
    }

    m_BP3Serializer.UpdateOffsetsInMetadata();

    if (isFinal) // Write metadata footer
    {
        format::BufferSTL &bufferSTL = m_BP3Serializer.m_Data;
        m_BP3Serializer.ResetBuffer(bufferSTL, false, false);

        m_BP3Serializer.AggregateCollectiveMetadata(m_BP3Serializer.m_Aggregator.m_Comm, bufferSTL,
                                                    false);

        if (m_BP3Serializer.m_Aggregator.m_IsAggregator)
        {
            m_FileDataManager.WriteFiles(bufferSTL.m_Buffer.data(), bufferSTL.m_Position,
                                         transportIndex);

            m_FileDataManager.FlushFiles(transportIndex);
        }

        m_BP3Serializer.m_Aggregator.Close();
    }

    m_BP3Serializer.m_Aggregator.ResetBuffers();
}

#define declare_type(T, L)                                                                         \
    T *BP3Writer::DoBufferData_##L(const int bufferIdx, const size_t payloadPosition,              \
                                   const size_t bufferID) noexcept                                 \
    {                                                                                              \
        return BufferDataCommon<T>(bufferIdx, payloadPosition, bufferID);                          \
    }

ADIOS2_FOREACH_PRIMITVE_STDTYPE_2ARGS(declare_type)
#undef declare_type

size_t BP3Writer::DebugGetDataBufferSize() const
{
    return m_BP3Serializer.DebugGetDataBufferSize();
}
} // end namespace engine
} // end namespace core
} // end namespace adios2
