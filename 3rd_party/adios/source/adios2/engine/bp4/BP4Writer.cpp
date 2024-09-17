/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Writer.cpp
 *
 *  Created on: Aug 1, 2018
 *      Author: Lipeng Wan wanl@ornl.gov
 */

#include "BP4Writer.h"
#include "BP4Writer.tcc"

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosFunctions.h" //CheckIndexRange
#include "adios2/toolkit/transport/file/FileFStream.h"
#include <adios2-perfstubs-interface.h>

#include <ctime>
#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

BP4Writer::BP4Writer(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("BP4Writer", io, name, mode, std::move(comm)), m_BP4Serializer(m_Comm),
  m_FileDataManager(io, m_Comm), m_FileMetadataManager(io, m_Comm),
  m_FileMetadataIndexManager(io, m_Comm), m_FileDrainer()
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::Open");
    helper::GetParameter(m_IO.m_Parameters, "Verbose", m_Verbosity);
    helper::Log("Engine", "BP4Writer", "Open", m_Name, 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    m_IO.m_ReadStreaming = false;

    Init();

    m_IsOpen = true;
}

BP4Writer::~BP4Writer()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus BP4Writer::BeginStep(StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::BeginStep");
    helper::Log("Engine", "BP4Writer", "BeginStep", std::to_string(CurrentStep()), 0, m_Comm.Rank(),
                5, m_Verbosity, helper::LogMode::INFO);

    m_BP4Serializer.m_DeferredVariables.clear();
    m_BP4Serializer.m_DeferredVariablesDataSize = 0;
    m_IO.m_ReadStreaming = false;
    m_DidBeginStep = true;
    return StepStatus::OK;
}

size_t BP4Writer::CurrentStep() const { return m_BP4Serializer.m_MetadataSet.CurrentStep; }

void BP4Writer::PerformPuts()
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::PerformPuts");
    helper::Log("Engine", "BP4Writer", "PerformPuts", "", 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    if (m_BP4Serializer.m_DeferredVariables.empty())
    {
        return;
    }

    m_BP4Serializer.ResizeBuffer(m_BP4Serializer.m_DeferredVariablesDataSize,
                                 "in call to PerformPuts");

    for (const std::string &variableName : m_BP4Serializer.m_DeferredVariables)
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
    m_BP4Serializer.m_DeferredVariables.clear();
    m_BP4Serializer.m_DeferredVariablesDataSize = 0;
}

void BP4Writer::EndStep()
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::EndStep");
    helper::Log("Engine", "BP4Writer", "EndStep", std::to_string(CurrentStep()), 0, m_Comm.Rank(),
                5, m_Verbosity, helper::LogMode::INFO);

    if (m_BP4Serializer.m_DeferredVariables.size() > 0)
    {
        PerformPuts();
    }

    // true: advances step
    m_BP4Serializer.SerializeData(m_IO, true);

    const size_t currentStep = CurrentStep();
    const size_t flushStepsCount = m_BP4Serializer.m_Parameters.FlushStepsCount;

    if (currentStep % flushStepsCount == 0)
    {
        Flush();
    }

    if (m_BP4Serializer.m_RankMPI == 0)
    {
        m_IO.m_ADIOS.RecordOutputStep(m_Name, UnknownStep, UnknownTime);
    }
}

void BP4Writer::Flush(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::Flush");
    DoFlush(false, transportIndex);
    m_BP4Serializer.ResetBuffer(m_BP4Serializer.m_Data, false, false);

    if (m_BP4Serializer.m_Parameters.CollectiveMetadata)
    {
        WriteCollectiveMetadataFile();
    }
}

// PRIVATE
void BP4Writer::Init()
{
    InitParameters();
    if (m_BP4Serializer.m_Parameters.NumAggregators <
        static_cast<unsigned int>(m_BP4Serializer.m_SizeMPI))
    {
        m_BP4Serializer.m_Aggregator.Init(m_BP4Serializer.m_Parameters.NumAggregators,
                                          m_BP4Serializer.m_Parameters.NumAggregators, m_Comm);
    }
    InitTransports();
    InitBPBuffer();
}

#define declare_type(T)                                                                            \
    void BP4Writer::DoPut(Variable<T> &variable, typename Variable<T>::Span &span,                 \
                          const bool initialize, const T &value)                                   \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Writer::Put");                                                  \
        helper::Log("Engine", "BP4Writer", "Put", variable.m_Name, 0, m_Comm.Rank(), 5,            \
                    m_Verbosity, helper::LogMode::INFO);                                           \
        PutCommon(variable, span, 0, value);                                                       \
    }

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

#define declare_type(T)                                                                            \
    void BP4Writer::DoPutSync(Variable<T> &variable, const T *data)                                \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Writer::Put");                                                  \
        helper::Log("Engine", "BP4Writer", "PutSync", variable.m_Name, 0, m_Comm.Rank(), 5,        \
                    m_Verbosity, helper::LogMode::INFO);                                           \
        PutSyncCommon(variable, variable.SetBlockInfo(data, CurrentStep()));                       \
        variable.m_BlocksInfo.pop_back();                                                          \
    }                                                                                              \
    void BP4Writer::DoPutDeferred(Variable<T> &variable, const T *data)                            \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Writer::Put");                                                  \
        helper::Log("Engine", "BP4Writer", "PutDeferred", variable.m_Name, 0, m_Comm.Rank(), 5,    \
                    m_Verbosity, helper::LogMode::INFO);                                           \
        PutDeferredCommon(variable, data);                                                         \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void BP4Writer::InitParameters()
{
    m_BP4Serializer.Init(m_IO.m_Parameters, "in call to BP4::Open to write");
    m_BP4Serializer.ResizeBuffer(m_BP4Serializer.m_Parameters.InitialBufferSize,
                                 "in call to BP4::Open to write");
    m_WriteToBB = !(m_BP4Serializer.m_Parameters.BurstBufferPath.empty());
    m_DrainBB = m_WriteToBB && m_BP4Serializer.m_Parameters.BurstBufferDrain;
}

void BP4Writer::InitTransports()
{
    // TODO need to add support for aggregators here later
    if (m_IO.m_TransportsParameters.empty())
    {
        Params defaultTransportParameters;
        defaultTransportParameters["transport"] = "File";
        m_IO.m_TransportsParameters.push_back(defaultTransportParameters);
    }

    // only consumers will interact with transport managers
    m_BBName = m_Name;
    if (m_WriteToBB)
    {
        m_BBName = m_BP4Serializer.m_Parameters.BurstBufferPath + PathSeparator + m_Name;
    }

    // Names passed to IO AddTransport option with key "Name"
    const std::vector<std::string> transportsNames =
        m_FileDataManager.GetFilesBaseNames(m_BBName, m_IO.m_TransportsParameters);

    if (m_BP4Serializer.m_Aggregator.m_IsAggregator)
    {
        // /path/name.bp.dir/name.bp.rank
        m_SubStreamNames = m_BP4Serializer.GetBPSubStreamNames(transportsNames);
        if (m_DrainBB)
        {
            const std::vector<std::string> drainTransportNames =
                m_FileDataManager.GetFilesBaseNames(m_Name, m_IO.m_TransportsParameters);
            m_DrainSubStreamNames = m_BP4Serializer.GetBPSubStreamNames(drainTransportNames);
            /* start up BB thread */
            m_FileDrainer.SetVerbose(m_BP4Serializer.m_Parameters.BurstBufferVerbose,
                                     m_BP4Serializer.m_RankMPI);
            m_FileDrainer.Start();
        }
    }

    /* Create the directories either on target or burst buffer if used */
    m_BP4Serializer.m_Profiler.Start("mkdir");
    m_FileDataManager.MkDirsBarrier(m_SubStreamNames, m_IO.m_TransportsParameters,
                                    m_BP4Serializer.m_Parameters.NodeLocal || m_WriteToBB);
    if (m_DrainBB)
    {
        /* Create the directories on target anyway by main thread */
        m_FileDataManager.MkDirsBarrier(m_DrainSubStreamNames, m_IO.m_TransportsParameters,
                                        m_BP4Serializer.m_Parameters.NodeLocal);
    }
    m_BP4Serializer.m_Profiler.Stop("mkdir");

    if (m_BP4Serializer.m_Aggregator.m_IsAggregator)
    {
        if (m_BP4Serializer.m_Parameters.AsyncOpen)
        {
            for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
            {
                m_IO.m_TransportsParameters[i]["asyncopen"] = "true";
            }
        }
        for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
        {
            m_IO.m_TransportsParameters[i].insert({"SingleProcess", "true"});
        }

        m_FileDataManager.OpenFiles(m_SubStreamNames, m_OpenMode, m_IO.m_TransportsParameters,
                                    m_BP4Serializer.m_Profiler.m_IsActive);

        if (m_DrainBB)
        {
            for (const auto &name : m_DrainSubStreamNames)
            {
                m_FileDrainer.AddOperationOpen(name, m_OpenMode);
            }
        }
    }

    if (m_BP4Serializer.m_RankMPI == 0)
    {
        // if (m_BP4Serializer.m_Parameters.CollectiveMetadata)
        //{
        const std::vector<std::string> transportsNames =
            m_FileMetadataManager.GetFilesBaseNames(m_BBName, m_IO.m_TransportsParameters);

        m_MetadataFileNames = m_BP4Serializer.GetBPMetadataFileNames(transportsNames);

        for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
        {
            m_IO.m_TransportsParameters[i].insert({"SingleProcess", "true"});
        }
        m_FileMetadataManager.OpenFiles(m_MetadataFileNames, m_OpenMode,
                                        m_IO.m_TransportsParameters,
                                        m_BP4Serializer.m_Profiler.m_IsActive);

        m_MetadataIndexFileNames = m_BP4Serializer.GetBPMetadataIndexFileNames(transportsNames);

        m_FileMetadataIndexManager.OpenFiles(m_MetadataIndexFileNames, m_OpenMode,
                                             m_IO.m_TransportsParameters,
                                             m_BP4Serializer.m_Profiler.m_IsActive);

        if (m_DrainBB)
        {
            const std::vector<std::string> drainTransportNames =
                m_FileDataManager.GetFilesBaseNames(m_Name, m_IO.m_TransportsParameters);
            m_DrainMetadataFileNames = m_BP4Serializer.GetBPMetadataFileNames(drainTransportNames);
            m_DrainMetadataIndexFileNames =
                m_BP4Serializer.GetBPMetadataIndexFileNames(drainTransportNames);

            for (const auto &name : m_DrainMetadataFileNames)
            {
                m_FileDrainer.AddOperationOpen(name, m_OpenMode);
            }
            for (const auto &name : m_DrainMetadataIndexFileNames)
            {
                m_FileDrainer.AddOperationOpen(name, m_OpenMode);
            }
        }
        //}
    }
}

void BP4Writer::InitBPBuffer()
{
    if (m_OpenMode == Mode::Append)
    {
        // TODO: Get last pg timestep and update timestep counter in
        format::BufferSTL preMetadataIndex;
        size_t preMetadataIndexFileSize;

        if (m_BP4Serializer.m_RankMPI == 0)
        {
            preMetadataIndexFileSize = m_FileMetadataIndexManager.GetFileSize(0);
            preMetadataIndex.m_Buffer.resize(preMetadataIndexFileSize);
            preMetadataIndex.m_Buffer.assign(preMetadataIndex.m_Buffer.size(), '\0');
            preMetadataIndex.m_Position = 0;
            m_FileMetadataIndexManager.ReadFile(preMetadataIndex.m_Buffer.data(),
                                                preMetadataIndexFileSize);
        }
        m_Comm.BroadcastVector(preMetadataIndex.m_Buffer);
        preMetadataIndexFileSize = preMetadataIndex.m_Buffer.size();
        if (preMetadataIndexFileSize > 0)
        {
            size_t position = 0;
            position += 28;
            const uint8_t endianness =
                helper::ReadValue<uint8_t>(preMetadataIndex.m_Buffer, position);
            bool IsLittleEndian = true;
            IsLittleEndian = (endianness == 0) ? true : false;
            if (helper::IsLittleEndian() != IsLittleEndian)
            {
                helper::Throw<std::runtime_error>(
                    "Engine", "BP4Writer", "InitBPBuffer",
                    "previous run generated BigEndian bp file, "
                    "this version of ADIOS2 wasn't compiled "
                    "with the cmake flag -DADIOS2_USE_Endian_Reverse=ON "
                    "explicitly, in call to Open");
            }
            const size_t pos_last_step = preMetadataIndexFileSize - 64;
            position = pos_last_step;
            const uint64_t lastStep =
                helper::ReadValue<uint64_t>(preMetadataIndex.m_Buffer, position, IsLittleEndian);
            m_BP4Serializer.m_MetadataSet.TimeStep += static_cast<uint32_t>(lastStep);
            m_BP4Serializer.m_MetadataSet.CurrentStep += lastStep;

            if (m_BP4Serializer.m_Aggregator.m_IsAggregator)
            {
                m_BP4Serializer.m_PreDataFileLength = m_FileDataManager.GetFileSize(0);
            }
            if (m_BP4Serializer.m_Aggregator.m_IsActive)
            {
                // only use aggregator comm if it is active (created)
                m_BP4Serializer.m_PreDataFileLength =
                    m_BP4Serializer.m_Aggregator.m_Comm.BroadcastValue(
                        m_BP4Serializer.m_PreDataFileLength);
            }

            if (m_BP4Serializer.m_RankMPI == 0)
            {
                // Get the size of existing metadata file
                m_BP4Serializer.m_PreMetadataFileLength = m_FileMetadataManager.GetFileSize(0);
            }
        }
    }

    if (m_BP4Serializer.m_PreDataFileLength == 0)
    {
        /* This is a new file.
         * Make headers in data buffer and metadata buffer (but do not write
         * them yet so that Open() can stay free of writing to disk)
         */
        if (m_BP4Serializer.m_RankMPI == 0)
        {
            m_BP4Serializer.MakeHeader(m_BP4Serializer.m_Metadata, "Metadata", false);
            m_BP4Serializer.MakeHeader(m_BP4Serializer.m_MetadataIndex, "Index Table", true);
        }
        if (m_BP4Serializer.m_Aggregator.m_IsAggregator)
        {
            m_BP4Serializer.MakeHeader(m_BP4Serializer.m_Data, "Data", false);
        }
    }
    else
    {
        if (m_BP4Serializer.m_RankMPI == 0)
        {
            // Set the flag in the header of metadata index table to 1 again
            // to indicate a new run begins
            UpdateActiveFlag(true);
        }
    }

    m_BP4Serializer.PutProcessGroupIndex(
        m_IO.m_Name, (m_IO.m_ArrayOrder == ArrayOrdering::RowMajor) ? "C++" : "Fortran",
        m_FileDataManager.GetTransportsTypes());
}

void BP4Writer::DoFlush(const bool isFinal, const int transportIndex)
{
    if (m_BP4Serializer.m_Aggregator.m_IsActive)
    {
        AggregateWriteData(isFinal, transportIndex);
    }
    else
    {
        WriteData(isFinal, transportIndex);
    }
}

void BP4Writer::DestructorClose(bool Verbose) noexcept
{
    if (Verbose)
    {
        std::cerr << "BP4 Writer \"" << m_Name << "\" Destroyed without a prior Close()."
                  << std::endl;
        std::cerr << "This may result in corrupt output." << std::endl;
    }
    // at least close metadata index file
    UpdateActiveFlag(false);
    m_IsOpen = false;
}

void BP4Writer::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::Close");
    helper::Log("Engine", "BP4Writer", "Close", m_Name, 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    if (m_BP4Serializer.m_DeferredVariables.size() > 0)
    {
        PerformPuts();
    }

    DoFlush(true, transportIndex);

    if (m_BP4Serializer.m_Aggregator.m_IsAggregator)
    {
        m_FileDataManager.CloseFiles(transportIndex);
        // Delete files from temporary storage if draining was on
        if (m_DrainBB)
        {
            for (const auto &name : m_SubStreamNames)
            {
                m_FileDrainer.AddOperationDelete(name);
            }
        }
    }

    if (m_BP4Serializer.m_Parameters.CollectiveMetadata && m_FileDataManager.AllTransportsClosed())
    {
        WriteCollectiveMetadataFile(true);
    }

    if (m_BP4Serializer.m_Profiler.m_IsActive && m_FileDataManager.AllTransportsClosed())
    {
        WriteProfilingJSONFile();
    }
    if (m_BP4Serializer.m_Aggregator.m_IsActive)
    {
        m_BP4Serializer.m_Aggregator.Close();
    }

    if (m_BP4Serializer.m_RankMPI == 0)
    {
        // Update the active flag in index to indicate current run is over.
        UpdateActiveFlag(false);

        // close metadata file
        m_FileMetadataManager.CloseFiles();

        // close metadata index file
        m_FileMetadataIndexManager.CloseFiles();

        // Delete metadata files from temporary storage if draining was on
        if (m_DrainBB)
        {
            for (const auto &name : m_MetadataFileNames)
            {
                m_FileDrainer.AddOperationDelete(name);
            }
            for (const auto &name : m_MetadataIndexFileNames)
            {
                m_FileDrainer.AddOperationDelete(name);
            }
            const std::vector<std::string> transportsNames =
                m_FileDataManager.GetFilesBaseNames(m_BBName, m_IO.m_TransportsParameters);
            for (const auto &name : transportsNames)
            {
                m_FileDrainer.AddOperationDelete(name);
            }
        }
    }

    if (m_BP4Serializer.m_Aggregator.m_IsAggregator && m_DrainBB)
    {
        /* Signal the BB thread that no more work is coming */
        m_FileDrainer.Finish();
    }

    if (!m_DidBeginStep && m_BP4Serializer.m_RankMPI == 0)
    {
        m_IO.m_ADIOS.RecordOutputStep(m_Name, UnknownStep, UnknownTime);
    }
    // m_BP4Serializer.DeleteBuffers();
}

void BP4Writer::WriteProfilingJSONFile()
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::WriteProfilingJSONFile");
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
        m_BP4Serializer.GetRankProfilingJSON(transportTypes, transportProfilers) + ",\n");

    const std::vector<char> profilingJSON(m_BP4Serializer.AggregateProfilingJSON(lineJSON));

    if (m_BP4Serializer.m_RankMPI == 0)
    {
        std::string profileFileName;
        if (m_DrainBB)
        {
            auto bpTargetNames = m_BP4Serializer.GetBPBaseNames({m_Name});
            if (fileTransportIdx > -1)
            {
                profileFileName = bpTargetNames[fileTransportIdx] + "/profiling.json";
            }
            else
            {
                profileFileName = bpTargetNames[0] + "_profiling.json";
            }
            m_FileDrainer.AddOperationWrite(profileFileName, profilingJSON.size(),
                                            profilingJSON.data());
        }
        else
        {
            transport::FileFStream profilingJSONStream(m_Comm);
            auto bpBaseNames = m_BP4Serializer.GetBPBaseNames({m_BBName});
            if (fileTransportIdx > -1)
            {
                profileFileName = bpBaseNames[fileTransportIdx] + "/profiling.json";
            }
            else
            {
                profileFileName = bpBaseNames[0] + "_profiling.json";
            }
            profilingJSONStream.Open(profileFileName, Mode::Write);
            profilingJSONStream.Write(profilingJSON.data(), profilingJSON.size());
            profilingJSONStream.Close();
        }
    }
}

/*write the content of metadata index file*/
void BP4Writer::PopulateMetadataIndexFileContent(format::BufferSTL &b, const uint64_t currentStep,
                                                 const uint64_t mpirank,
                                                 const uint64_t pgIndexStart,
                                                 const uint64_t variablesIndexStart,
                                                 const uint64_t attributesIndexStart,
                                                 const uint64_t currentStepEndPos,
                                                 const uint64_t currentTimeStamp)
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::PopulateMetadataIndexFileContent");
    auto &buffer = b.m_Buffer;
    auto &position = b.m_Position;
    helper::CopyToBuffer(buffer, position, &currentStep);
    helper::CopyToBuffer(buffer, position, &mpirank);
    helper::CopyToBuffer(buffer, position, &pgIndexStart);
    helper::CopyToBuffer(buffer, position, &variablesIndexStart);
    helper::CopyToBuffer(buffer, position, &attributesIndexStart);
    helper::CopyToBuffer(buffer, position, &currentStepEndPos);
    helper::CopyToBuffer(buffer, position, &currentTimeStamp);
    position += 8;
}

void BP4Writer::UpdateActiveFlag(const bool active)
{
    const char activeChar = (active ? '\1' : '\0');
    m_FileMetadataIndexManager.WriteFileAt(&activeChar, 1, m_BP4Serializer.m_ActiveFlagPosition);
    m_FileMetadataIndexManager.FlushFiles();
    m_FileMetadataIndexManager.SeekToFileEnd();
    if (m_DrainBB)
    {
        for (size_t i = 0; i < m_MetadataIndexFileNames.size(); ++i)
        {
            m_FileDrainer.AddOperationWriteAt(m_DrainMetadataIndexFileNames[i],
                                              m_BP4Serializer.m_ActiveFlagPosition, 1, &activeChar);
            m_FileDrainer.AddOperationSeekEnd(m_DrainMetadataIndexFileNames[i]);
        }
    }
}

void BP4Writer::WriteCollectiveMetadataFile(const bool isFinal)
{

    PERFSTUBS_SCOPED_TIMER("BP4Writer::WriteCollectiveMetadataFile");

    if (isFinal && m_BP4Serializer.m_MetadataSet.DataPGCount == 0)
    {
        // If data pg count is zero, it means all metadata
        // has already been written, don't need to write it again.
        return;
    }
    m_BP4Serializer.AggregateCollectiveMetadata(m_Comm, m_BP4Serializer.m_Metadata, true);

    if (m_BP4Serializer.m_RankMPI == 0)
    {

        m_FileMetadataManager.WriteFiles(m_BP4Serializer.m_Metadata.m_Buffer.data(),
                                         m_BP4Serializer.m_Metadata.m_Position);
        m_FileMetadataManager.FlushFiles();

        if (m_DrainBB)
        {
            for (size_t i = 0; i < m_MetadataFileNames.size(); ++i)
            {
                m_FileDrainer.AddOperationCopy(m_MetadataFileNames[i], m_DrainMetadataFileNames[i],
                                               m_BP4Serializer.m_Metadata.m_Position);
            }
        }

        std::time_t currentTimeStamp = std::time(nullptr);

        std::vector<size_t> timeSteps;
        timeSteps.reserve(m_BP4Serializer.m_MetadataIndexTable[m_BP4Serializer.m_RankMPI].size());
        for (auto const &pair : m_BP4Serializer.m_MetadataIndexTable[m_BP4Serializer.m_RankMPI])
        {
            timeSteps.push_back(pair.first);
        }
        std::sort(timeSteps.begin(), timeSteps.end());

        size_t rowsInMetadataIndexTable = timeSteps.size() + 1;
        m_BP4Serializer.m_MetadataIndex.Resize(rowsInMetadataIndexTable * 64, "BP4 Index Table");
        for (auto const &t : timeSteps)
        {
            const uint64_t pgIndexStartMetadataFile =
                m_BP4Serializer.m_MetadataIndexTable[m_BP4Serializer.m_RankMPI][t][0] +
                m_BP4Serializer.m_MetadataSet.MetadataFileLength +
                m_BP4Serializer.m_PreMetadataFileLength;
            const uint64_t varIndexStartMetadataFile =
                m_BP4Serializer.m_MetadataIndexTable[m_BP4Serializer.m_RankMPI][t][1] +
                m_BP4Serializer.m_MetadataSet.MetadataFileLength +
                m_BP4Serializer.m_PreMetadataFileLength;
            const uint64_t attrIndexStartMetadataFile =
                m_BP4Serializer.m_MetadataIndexTable[m_BP4Serializer.m_RankMPI][t][2] +
                m_BP4Serializer.m_MetadataSet.MetadataFileLength +
                m_BP4Serializer.m_PreMetadataFileLength;
            const uint64_t currentStepEndPosMetadataFile =
                m_BP4Serializer.m_MetadataIndexTable[m_BP4Serializer.m_RankMPI][t][3] +
                m_BP4Serializer.m_MetadataSet.MetadataFileLength +
                m_BP4Serializer.m_PreMetadataFileLength;
            PopulateMetadataIndexFileContent(m_BP4Serializer.m_MetadataIndex, t,
                                             m_BP4Serializer.m_RankMPI, pgIndexStartMetadataFile,
                                             varIndexStartMetadataFile, attrIndexStartMetadataFile,
                                             currentStepEndPosMetadataFile, currentTimeStamp);
        }

        m_FileMetadataIndexManager.WriteFiles(m_BP4Serializer.m_MetadataIndex.m_Buffer.data(),
                                              m_BP4Serializer.m_MetadataIndex.m_Position);
        m_FileMetadataIndexManager.FlushFiles();

        m_BP4Serializer.m_MetadataSet.MetadataFileLength += m_BP4Serializer.m_Metadata.m_Position;

        if (m_DrainBB)
        {
            for (size_t i = 0; i < m_MetadataIndexFileNames.size(); ++i)
            {
                m_FileDrainer.AddOperationWrite(m_DrainMetadataIndexFileNames[i],
                                                m_BP4Serializer.m_MetadataIndex.m_Position,
                                                m_BP4Serializer.m_MetadataIndex.m_Buffer.data());
            }
        }
    }
    /*Clear the local indices buffer at the end of each step*/
    m_BP4Serializer.ResetBuffer(m_BP4Serializer.m_Metadata, true, true);

    /* clear the metadata index buffer*/
    m_BP4Serializer.ResetBuffer(m_BP4Serializer.m_MetadataIndex, true, true);

    /* reset the metadata index table*/
    m_BP4Serializer.ResetMetadataIndexTable();
    m_BP4Serializer.ResetAllIndices();
}

void BP4Writer::WriteData(const bool isFinal, const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::WriteData");
    size_t dataSize;

    // write data without footer
    if (isFinal)
    {
        dataSize = m_BP4Serializer.CloseData(m_IO);
    }
    else
    {
        dataSize = m_BP4Serializer.CloseStream(m_IO, false);
    }

    m_FileDataManager.WriteFiles(m_BP4Serializer.m_Data.m_Buffer.data(), dataSize, transportIndex);

    m_FileDataManager.FlushFiles(transportIndex);
    if (m_DrainBB)
    {
        for (size_t i = 0; i < m_SubStreamNames.size(); ++i)
        {
            m_FileDrainer.AddOperationCopy(m_SubStreamNames[i], m_DrainSubStreamNames[i], dataSize);
        }
    }
}

void BP4Writer::AggregateWriteData(const bool isFinal, const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP4Writer::AggregateWriteData");
    m_BP4Serializer.CloseStream(m_IO, false);
    size_t totalBytesWritten = 0;
    const size_t dataBufferSize = m_BP4Serializer.m_Data.m_Position;

    // async?
    for (int r = 0; r < m_BP4Serializer.m_Aggregator.m_Size; ++r)
    {
        aggregator::MPIChain::ExchangeRequests dataRequests =
            m_BP4Serializer.m_Aggregator.IExchange(m_BP4Serializer.m_Data, r);

        aggregator::MPIChain::ExchangeAbsolutePositionRequests absolutePositionRequests =
            m_BP4Serializer.m_Aggregator.IExchangeAbsolutePosition(m_BP4Serializer.m_Data, r);

        if (m_BP4Serializer.m_Aggregator.m_IsAggregator)
        {
            const format::Buffer &bufferSTL =
                m_BP4Serializer.m_Aggregator.GetConsumerBuffer(m_BP4Serializer.m_Data);
            if (bufferSTL.m_Position > 0)
            {
                m_FileDataManager.WriteFiles(bufferSTL.Data(), bufferSTL.m_Position,
                                             transportIndex);

                m_FileDataManager.FlushFiles(transportIndex);

                totalBytesWritten += bufferSTL.m_Position;
            }
        }

        m_BP4Serializer.m_Aggregator.WaitAbsolutePosition(absolutePositionRequests, r);

        m_BP4Serializer.m_Aggregator.Wait(dataRequests, r);
        m_BP4Serializer.m_Aggregator.SwapBuffers(r);
    }

    if (m_DrainBB)
    {
        for (size_t i = 0; i < m_SubStreamNames.size(); ++i)
        {
            m_FileDrainer.AddOperationCopy(m_SubStreamNames[i], m_DrainSubStreamNames[i],
                                           totalBytesWritten);
        }
    }

    m_BP4Serializer.UpdateOffsetsInMetadata();

    if (isFinal) // Write metadata footer
    {
        m_BP4Serializer.m_Aggregator.Close();
    }

    m_BP4Serializer.m_Aggregator.ResetBuffers();

    // Reset Data buffer to its final size on this process at EndStep
    // The aggregation routine has resized it to some incoming process' data
    // size
    m_BP4Serializer.m_Data.Resize(dataBufferSize, "Reset buffersize to final size" +
                                                      std::to_string(dataBufferSize));
}

#define declare_type(T, L)                                                                         \
    T *BP4Writer::DoBufferData_##L(const int bufferIdx, const size_t payloadPosition,              \
                                   const size_t bufferID) noexcept                                 \
    {                                                                                              \
        return BufferDataCommon<T>(bufferIdx, payloadPosition, bufferID);                          \
    }

ADIOS2_FOREACH_PRIMITVE_STDTYPE_2ARGS(declare_type)
#undef declare_type

size_t BP4Writer::DebugGetDataBufferSize() const
{
    return m_BP4Serializer.DebugGetDataBufferSize();
}

void BP4Writer::NotifyEngineAttribute(std::string name, DataType type) noexcept
{
    m_BP4Serializer.m_SerializedAttributes.erase(name);
}

void BP4Writer::NotifyEngineAttribute(std::string name, AttributeBase *attr, void *Data) noexcept
{
    NotifyEngineAttribute(name, attr->m_Type);
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
