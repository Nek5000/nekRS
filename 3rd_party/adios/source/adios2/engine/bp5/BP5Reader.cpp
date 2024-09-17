/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP5Reader.cpp
 *
 */

#include "BP5Reader.h"
#include "BP5Reader.tcc"

#include "adios2/helper/adiosMath.h" // SetWithinLimit
#include "adios2/toolkit/transport/file/FileFStream.h"
#include <adios2-perfstubs-interface.h>

#include <chrono>
#include <cstdio>
#include <errno.h>
#include <iostream>
#include <mutex>
#include <thread>

using TP = std::chrono::high_resolution_clock::time_point;
#define NOW() std::chrono::high_resolution_clock::now();
#define DURATION(T1, T2) static_cast<double>((T2 - T1).count()) / 1000000000.0;

namespace adios2
{
namespace core
{
namespace engine
{

BP5Reader::BP5Reader(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("BP5Reader", io, name, mode, std::move(comm)), m_MDFileManager(io, m_Comm),
  m_DataFileManager(io, m_Comm), m_MDIndexFileManager(io, m_Comm),
  m_FileMetaMetadataManager(io, m_Comm), m_ActiveFlagFileManager(io, m_Comm), m_Remote(),
  m_JSONProfiler(m_Comm)
{
    PERFSTUBS_SCOPED_TIMER("BP5Reader::Open");
    Init();
    m_IsOpen = true;
}

BP5Reader::~BP5Reader()
{
    if (m_BP5Deserializer)
        delete m_BP5Deserializer;
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

void BP5Reader::DestructorClose(bool Verbose) noexcept
{
    // Nothing special needs to be done to "close" a BP5 reader during shutdown
    // if it hasn't already been Closed
    m_IsOpen = false;
}

void BP5Reader::InstallMetadataForTimestep(size_t Step)
{
    size_t pgstart = m_MetadataIndexTable[Step][0];
    size_t Position = pgstart + sizeof(uint64_t); // skip total data size
    const uint64_t WriterCount = m_WriterMap[m_WriterMapIndex[Step]].WriterCount;
    size_t MDPosition = Position + 2 * sizeof(uint64_t) * WriterCount;
    for (size_t WriterRank = 0; WriterRank < WriterCount; WriterRank++)
    {
        // variable metadata for timestep
        size_t ThisMDSize =
            helper::ReadValue<uint64_t>(m_Metadata.Data(), Position, m_Minifooter.IsLittleEndian);
        char *ThisMD = m_Metadata.Data() + MDPosition;
        if ((m_OpenMode == Mode::ReadRandomAccess) || (m_FlattenSteps))
        {
            m_BP5Deserializer->InstallMetaData(ThisMD, ThisMDSize, WriterRank, Step);
        }
        else
        {
            m_BP5Deserializer->InstallMetaData(ThisMD, ThisMDSize, WriterRank);
        }
        MDPosition += ThisMDSize;
    }
    for (size_t WriterRank = 0; WriterRank < WriterCount; WriterRank++)
    {
        // attribute metadata for timestep
        size_t ThisADSize =
            helper::ReadValue<uint64_t>(m_Metadata.Data(), Position, m_Minifooter.IsLittleEndian);
        char *ThisAD = m_Metadata.Data() + MDPosition;
        if (ThisADSize > 0)
            m_BP5Deserializer->InstallAttributeData(ThisAD, ThisADSize);
        MDPosition += ThisADSize;
    }
}

StepStatus BP5Reader::BeginStep(StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER("BP5Reader::BeginStep");

    if (m_OpenMode != Mode::Read)
    {
        helper::Throw<std::logic_error>("Engine", "BP5Reader", "BeginStep",
                                        "BeginStep called in random access mode");
    }
    if (m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "BP5Reader", "BeginStep",
                                        "BeginStep() is called a second time "
                                        "without an intervening EndStep()");
    }

    if (mode != StepMode::Read)
    {
        helper::Throw<std::invalid_argument>(
            "Engine", "BP5Reader", "BeginStep",
            "mode is not supported yet, only Read is valid for engine "
            "BP5Reader, in call to BeginStep");
    }

    StepStatus status = StepStatus::OK;
    if (m_FirstStep)
    {
        if (!m_StepsCount)
        {
            // not steps was found in Open/Init, check for new steps now
            status = CheckForNewSteps(Seconds(timeoutSeconds));
        }
    }
    else
    {
        if (m_CurrentStep + 1 >= m_StepsCount)
        {
            // we processed steps in memory, check for new steps now
            status = CheckForNewSteps(Seconds(timeoutSeconds));
        }
    }

    if (status == StepStatus::OK)
    {
        m_BetweenStepPairs = true;
        if (m_FirstStep)
        {
            m_FirstStep = false;
        }
        else
        {
            ++m_CurrentStep;
        }

        m_IO.m_EngineStep = m_CurrentStep;
        //        SstBlock AttributeBlockList =
        //            SstGetAttributeData(m_Input, SstCurrentStep(m_Input));
        //        i = 0;
        //        while (AttributeBlockList && AttributeBlockList[i].BlockData)
        //        {
        //            m_IO.RemoveAllAttributes();
        //            m_BP5Deserializer->InstallAttributeData(
        //                AttributeBlockList[i].BlockData,
        //                AttributeBlockList[i].BlockSize);
        //            i++;
        //        }

        m_BP5Deserializer->SetupForStep(m_CurrentStep,
                                        m_WriterMap[m_WriterMapIndex[m_CurrentStep]].WriterCount);

        /* Remove all existing variables from previous steps
           It seems easier than trying to update them */
        // m_IO.RemoveAllVariables();

        InstallMetadataForTimestep(m_CurrentStep);
        m_IO.ResetVariablesStepSelection(false, "in call to BP5 Reader BeginStep");

        // caches attributes for each step
        // if a variable name is a prefix
        // e.g. var  prefix = {var/v1, var/v2, var/v3}
        m_IO.SetPrefixedNames(true);
    }

    return status;
}

size_t BP5Reader::CurrentStep() const { return m_CurrentStep; }

void BP5Reader::EndStep()
{
    if (m_OpenMode != Mode::Read)
    {
        helper::Throw<std::logic_error>("Engine", "BP5Reader", "EndStep",
                                        "EndStep called in random access mode");
    }
    if (!m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "BP5Reader", "EndStep",
                                        "EndStep() is called without a successful BeginStep()");
    }
    m_BetweenStepPairs = false;
    PERFSTUBS_SCOPED_TIMER("BP5Reader::EndStep");
    PerformGets();
}

std::pair<double, double> BP5Reader::ReadData(adios2::transportman::TransportMan &FileManager,
                                              const size_t maxOpenFiles, const size_t WriterRank,
                                              const size_t Timestep, const size_t StartOffset,
                                              const size_t Length, char *Destination)
{
    /*
     * Warning: this function is called by multiple threads
     */
    size_t FlushCount = m_MetadataIndexTable[Timestep][2];
    size_t DataPosPos = m_MetadataIndexTable[Timestep][3];
    size_t SubfileNum =
        static_cast<size_t>(m_WriterMap[m_WriterMapIndex[Timestep]].RankToSubfile[WriterRank]);

    // check if subfile is already opened
    TP startSubfile = NOW();
    if (FileManager.m_Transports.count(SubfileNum) == 0)
    {
        const std::string subFileName =
            GetBPSubStreamName(m_Name, SubfileNum, m_Minifooter.HasSubFiles, true);
        if (FileManager.m_Transports.size() >= maxOpenFiles)
        {
            auto m = FileManager.m_Transports.begin();
            FileManager.CloseFiles((int)m->first);
        }
        FileManager.OpenFileID(subFileName, SubfileNum, Mode::Read, m_IO.m_TransportsParameters[0],
                               /*{{"transport", "File"}},*/ false);
        if (!m_WriterIsActive)
        {
            Params transportParameters;
            transportParameters["FailOnEOF"] = "true";
            FileManager.SetParameters(transportParameters, -1);
        }
    }
    TP endSubfile = NOW();
    double timeSubfile = DURATION(startSubfile, endSubfile);

    /* Each block is in exactly one flush. The StartOffset was calculated
       as if all the flushes were in a single contiguous block in file.
    */
    TP startRead = NOW();
    size_t InfoStartPos = DataPosPos + (WriterRank * (2 * FlushCount + 1) * sizeof(uint64_t));
    size_t SumDataSize = 0; // count in contiguous space
    for (size_t flush = 0; flush < FlushCount; flush++)
    {
        size_t ThisDataPos = helper::ReadValue<uint64_t>(m_MetadataIndex.m_Buffer, InfoStartPos,
                                                         m_Minifooter.IsLittleEndian);
        size_t ThisDataSize = helper::ReadValue<uint64_t>(m_MetadataIndex.m_Buffer, InfoStartPos,
                                                          m_Minifooter.IsLittleEndian);

        if (StartOffset < SumDataSize + ThisDataSize)
        {
            // discount offsets of skipped flushes
            size_t Offset = StartOffset - SumDataSize;
            FileManager.ReadFile(Destination, Length, ThisDataPos + Offset, SubfileNum);
            TP endRead = NOW();
            double timeRead = DURATION(startRead, endRead);
            return std::make_pair(timeSubfile, timeRead);
        }
        SumDataSize += ThisDataSize;
    }

    size_t ThisDataPos = helper::ReadValue<uint64_t>(m_MetadataIndex.m_Buffer, InfoStartPos,
                                                     m_Minifooter.IsLittleEndian);
    size_t Offset = StartOffset - SumDataSize;
    FileManager.ReadFile(Destination, Length, ThisDataPos + Offset, SubfileNum);

    TP endRead = NOW();
    double timeRead = DURATION(startRead, endRead);
    return std::make_pair(timeSubfile, timeRead);
}

void BP5Reader::PerformGets()
{
    // if dataIsRemote is true and m_Remote is not true, this is our first time through
    // PerformGets() Either we don't need a remote open (m_dataIsRemote=false), or we need to Open
    // remote file (or die trying)
    if (m_dataIsRemote && !m_Remote)
    {
        bool RowMajorOrdering = (m_IO.m_ArrayOrder == ArrayOrdering::RowMajor);

        // If nothing is pending, don't open
        if (m_BP5Deserializer->PendingGetRequests.size() == 0)
            return;

        if (!m_Parameters.RemoteDataPath.empty())
        {
            m_Remote.Open("localhost", RemoteCommon::ServerPort, m_Parameters.RemoteDataPath,
                          m_OpenMode, RowMajorOrdering);
        }
        else if (getenv("DoRemote"))
        {
            m_Remote.Open("localhost", RemoteCommon::ServerPort, m_Name, m_OpenMode,
                          RowMajorOrdering);
        }
        if (!m_Remote)
        {
            helper::Throw<std::ios_base::failure>(
                "Engine", "BP5Reader", "OpenFiles",
                "Remote file " + m_Name +
                    " cannot be opened. Possible server or file specification error.");
        }
    }

    if (m_Remote)
    {
        PerformRemoteGets();
    }
    else
    {
        PerformLocalGets();
    }

    // clear pending requests inside deserializer
    {
        std::vector<adios2::format::BP5Deserializer::ReadRequest> empty;
        m_BP5Deserializer->FinalizeGets(empty);
    }
}

void BP5Reader::PerformRemoteGets()
{
    // TP startGenerate = NOW();
    auto GetRequests = m_BP5Deserializer->PendingGetRequests;
    for (auto &Req : GetRequests)
    {
        m_Remote.Get(Req.VarName, Req.RelStep, Req.BlockID, Req.Count, Req.Start, Req.Data);
    }
}

void BP5Reader::PerformLocalGets()
{
    auto lf_CompareReqSubfile =
        [&](const adios2::format::BP5Deserializer::ReadRequest &r1,
            const adios2::format::BP5Deserializer::ReadRequest &r2) -> bool {
        return (m_WriterMap[m_WriterMapIndex[r1.Timestep]].RankToSubfile[r1.WriterRank] <
                m_WriterMap[m_WriterMapIndex[r2.Timestep]].RankToSubfile[r2.WriterRank]);
    };

    if (!m_InitialWriterActiveCheckDone)
    {
        CheckWriterActive();
        m_InitialWriterActiveCheckDone = true;
    }
    // TP start = NOW();
    PERFSTUBS_SCOPED_TIMER("BP5Reader::PerformGets");
    m_JSONProfiler.Start("DataRead");
    size_t maxReadSize;

    // TP startGenerate = NOW();
    auto ReadRequests = m_BP5Deserializer->GenerateReadRequests(false, &maxReadSize);
    size_t nRequest = ReadRequests.size();
    // TP endGenerate = NOW();
    // double generateTime = DURATION(startGenerate, endGenerate);

    size_t nextRequest = 0;
    std::mutex mutexReadRequests;

    auto lf_GetNextRequest = [&]() -> size_t {
        std::lock_guard<std::mutex> lockGuard(mutexReadRequests);
        size_t reqidx = MaxSizeT;
        if (nextRequest < nRequest)
        {
            reqidx = nextRequest;
            ++nextRequest;
        }
        if (reqidx <= nRequest)
        {
            m_JSONProfiler.AddBytes("dataread", ReadRequests[reqidx].ReadLength);
        }
        return reqidx;
    };

    auto lf_Reader = [&](const int FileManagerID,
                         const size_t maxOpenFiles) -> std::tuple<double, double, double, size_t> {
        double copyTotal = 0.0;
        double readTotal = 0.0;
        double subfileTotal = 0.0;
        size_t nReads = 0;
        std::vector<char> buf(maxReadSize);

        while (true)
        {
            const auto reqidx = lf_GetNextRequest();
            if (reqidx > nRequest)
            {
                break;
            }
            auto &Req = ReadRequests[reqidx];
            if (!Req.DestinationAddr)
            {
                Req.DestinationAddr = buf.data();
            }
            std::pair<double, double> t =
                ReadData(fileManagers[FileManagerID], maxOpenFiles, Req.WriterRank, Req.Timestep,
                         Req.StartOffset, Req.ReadLength, Req.DestinationAddr);

            TP startCopy = NOW();
            m_BP5Deserializer->FinalizeGet(Req, false);
            TP endCopy = NOW();
            subfileTotal += t.first;
            readTotal += t.second;
            copyTotal += DURATION(startCopy, endCopy);
            ++nReads;
        }
        return std::make_tuple(subfileTotal, readTotal, copyTotal, nReads);
    };

    // TP startRead = NOW();
    // double sortTime = 0.0;
    if (m_Threads > 1 && nRequest > 1)
    {
        // TP startSort = NOW();
        std::sort(ReadRequests.begin(), ReadRequests.end(), lf_CompareReqSubfile);
        // TP endSort = NOW();
        // sortTime = DURATION(startSort, endSort);
        size_t nThreads = (m_Threads < nRequest ? m_Threads : nRequest);

        size_t maxOpenFiles = helper::SetWithinLimit(
            (size_t)m_Parameters.MaxOpenFilesAtOnce / nThreads, (size_t)1, MaxSizeT);

        std::vector<std::future<std::tuple<double, double, double, size_t>>> futures(nThreads - 1);

        // launch Threads-1 threads to process subsets of requests,
        // then main thread process the last subset
        for (size_t tid = 0; tid < nThreads - 1; ++tid)
        {
            futures[tid] = std::async(std::launch::async, lf_Reader, (int)(tid + 1), maxOpenFiles);
        }
        // main thread runs last subset of reads
        /*auto tMain = */ lf_Reader(0, maxOpenFiles);
        /*{
            double tSubfile = std::get<0>(tMain);
            double tRead = std::get<1>(tMain);
            double tCopy = std::get<2>(tMain);
            size_t nReads = std::get<3>(tMain);
            std::cout << " -> PerformGets() thread MAIN total = "
                      << tSubfile + tRead + tCopy << "s, subfile = " << tSubfile
                      << "s, read = " << tRead << "s, copy = " << tCopy
                      << ", nReads = " << nReads << std::endl;
        }*/

        // wait for all async threads
        for (auto &f : futures)
        {
            f.get();
        }
    }
    else
    {
        size_t maxOpenFiles =
            helper::SetWithinLimit((size_t)m_Parameters.MaxOpenFilesAtOnce, (size_t)1, MaxSizeT);
        std::vector<char> buf(maxReadSize);
        for (auto &Req : ReadRequests)
        {
            if (!Req.DestinationAddr)
            {
                Req.DestinationAddr = buf.data();
            }
            m_JSONProfiler.AddBytes("dataread", Req.ReadLength);
            ReadData(m_DataFileManager, maxOpenFiles, Req.WriterRank, Req.Timestep, Req.StartOffset,
                     Req.ReadLength, Req.DestinationAddr);
            m_BP5Deserializer->FinalizeGet(Req, false);
        }
    }
    m_JSONProfiler.Stop("DataRead");
    /*TP end = NOW();
    double t1 = DURATION(start, end);
    double t2 = DURATION(startRead, end);
    std::cout << " -> PerformGets() total = " << t1 << "s, Read loop = " << t2
              << "s, sort = " << sortTime << "s, generate = " << generateTime
              << ", nRequests = " << nRequest << std::endl;*/
}

// PRIVATE
void BP5Reader::Init()
{
    if ((m_OpenMode != Mode::Read) && (m_OpenMode != Mode::ReadRandomAccess))
    {
        helper::Throw<std::invalid_argument>("Engine", "BP5Reader", "Init",
                                             "BPFileReader only supports OpenMode::Read or "
                                             "OpenMode::ReadRandomAccess from" +
                                                 m_Name);
    }

    // if IO was involved in reading before this flag may be true now
    m_IO.m_ReadStreaming = false;
    m_ReaderIsRowMajor = (m_IO.m_ArrayOrder == ArrayOrdering::RowMajor);
    InitParameters();
    InitTransports();
    if (!m_Parameters.SelectSteps.empty())
    {
        m_SelectedSteps.ParseSelection(m_Parameters.SelectSteps);
    }

    /* Do a collective wait for the file(s) to appear within timeout.
       Make sure every process comes to the same conclusion */
    const Seconds timeoutSeconds = Seconds(m_Parameters.OpenTimeoutSecs);

    Seconds pollSeconds = Seconds(m_Parameters.BeginStepPollingFrequencySecs);
    if (pollSeconds > timeoutSeconds)
    {
        pollSeconds = timeoutSeconds;
    }

    TimePoint timeoutInstant = Now() + timeoutSeconds;
    OpenFiles(timeoutInstant, pollSeconds, timeoutSeconds);
    UpdateBuffer(timeoutInstant, pollSeconds / 10, timeoutSeconds);

    // Don't try to open the remote file when we open local metadata.  Do that on demand.
    if (!m_Parameters.RemoteDataPath.empty())
        m_dataIsRemote = true;
    if (getenv("DoRemote"))
        m_dataIsRemote = true;
}

void BP5Reader::InitParameters()
{
    ParseParams(m_IO, m_Parameters);
    if (m_Parameters.OpenTimeoutSecs < 0.0f)
    {
        if (m_OpenMode == Mode::ReadRandomAccess)
        {
            m_Parameters.OpenTimeoutSecs = 0.0f;
        }
        else
        {
            m_Parameters.OpenTimeoutSecs = 3600.0f;
        }
    }

    m_Threads = m_Parameters.Threads;
    if (m_Threads == 0)
    {
        helper::Comm m_NodeComm = m_Comm.GroupByShm("creating per-node comm at BP5 Open(read)");
        unsigned int NodeSize = static_cast<unsigned int>(m_NodeComm.Size());
        unsigned int NodeThreadSize = helper::NumHardwareThreadsPerNode();
        if (NodeThreadSize > 0)
        {
            m_Threads = helper::SetWithinLimit(NodeThreadSize / NodeSize, 1U, 16U);
        }
        else
        {
            m_Threads = helper::SetWithinLimit(8U / NodeSize, 1U, 8U);
        }
    }

    // Create m_Threads-1  extra file managers to be used by threads
    // The main thread uses the DataFileManager pushed here to vector[0]
    fileManagers.push_back(m_DataFileManager);
    for (unsigned int i = 0; i < m_Threads - 1; ++i)
    {
        fileManagers.push_back(
            transportman::TransportMan(transportman::TransportMan(m_IO, singleComm)));
    }

    size_t limit = helper::RaiseLimitNoFile();
    if (m_Parameters.MaxOpenFilesAtOnce > (unsigned int)limit - 8)
    {
        m_Parameters.MaxOpenFilesAtOnce = (unsigned int)limit - 8;
    }
}

bool BP5Reader::SleepOrQuit(const TimePoint &timeoutInstant, const Seconds &pollSeconds)
{
    auto now = Now();
    if (now >= timeoutInstant)
    {
        return false;
    }
    auto remainderTime = timeoutInstant - now;
    auto sleepTime = pollSeconds;
    if (remainderTime < sleepTime)
    {
        sleepTime = remainderTime;
    }
    std::this_thread::sleep_for(sleepTime);
    return true;
}

size_t BP5Reader::OpenWithTimeout(transportman::TransportMan &tm,
                                  const std::vector<std::string> &fileNames,
                                  const TimePoint &timeoutInstant, const Seconds &pollSeconds,
                                  std::string &lasterrmsg /*INOUT*/)
{
    size_t flag = 1; // 0 = OK, opened file, 1 = timeout, 2 = error
    do
    {
        try
        {
            errno = 0;
            const bool profile = false; // m_BP4Deserializer.m_Profiler.m_IsActive;
            tm.OpenFiles(fileNames, adios2::Mode::Read, m_IO.m_TransportsParameters, profile);
            flag = 0; // found file
            break;
        }
        catch (std::ios_base::failure &e)
        {
            lasterrmsg = std::string("errno=" + std::to_string(errno) + ": " + e.what());
            if (errno == ENOENT)
            {
                flag = 1; // timeout
            }
            else
            {
                flag = 2; // fatal error
                break;
            }
        }
    } while (SleepOrQuit(timeoutInstant, pollSeconds));
    return flag;
}

void BP5Reader::OpenFiles(TimePoint &timeoutInstant, const Seconds &pollSeconds,
                          const Seconds &timeoutSeconds)
{
    /* Poll */
    size_t flag = 1; // 0 = OK, opened file, 1 = timeout, 2 = error
    std::string lasterrmsg;
    if (m_Comm.Rank() == 0)
    {
        /* Open the metadata index table */
        const std::string metadataIndexFile(GetBPMetadataIndexFileName(m_Name));

        flag = OpenWithTimeout(m_MDIndexFileManager, {metadataIndexFile}, timeoutInstant,
                               pollSeconds, lasterrmsg);
        if (flag == 0)
        {
            /* Open the metadata file */
            const std::string metadataFile(GetBPMetadataFileName(m_Name));

            /* We found md.idx. If we don't find md.0 immediately  we should
             * wait a little bit hoping for the file system to catch up.
             * This slows down finding the error in file reading mode but
             * it will be more robust in streaming mode
             */
            if (timeoutSeconds == Seconds(0.0))
            {
                timeoutInstant += Seconds(5.0);
            }

            flag = OpenWithTimeout(m_MDFileManager, {metadataFile}, timeoutInstant, pollSeconds,
                                   lasterrmsg);
            if (flag != 0)
            {
                /* Close the metadata index table */
                m_MDIndexFileManager.CloseFiles();
            }
            else
            {
                /* Open the metametadata file */
                const std::string metametadataFile(GetBPMetaMetadataFileName(m_Name));

                /* We found md.idx. If we don't find md.0 immediately  we should
                 * wait a little bit hoping for the file system to catch up.
                 * This slows down finding the error in file reading mode but
                 * it will be more robust in streaming mode
                 */
                if (timeoutSeconds == Seconds(0.0))
                {
                    timeoutInstant += Seconds(5.0);
                }

                flag = OpenWithTimeout(m_FileMetaMetadataManager, {metametadataFile},
                                       timeoutInstant, pollSeconds, lasterrmsg);
                if (flag != 0)
                {
                    /* Close the metametadata index table */
                    m_MDIndexFileManager.CloseFiles();
                    m_MDFileManager.CloseFiles();
                }
            }
        }
    }

    flag = m_Comm.BroadcastValue(flag, 0);
    if (flag == 2)
    {
        if (m_Comm.Rank() == 0 && !lasterrmsg.empty())
        {
            helper::Throw<std::ios_base::failure>("Engine", "BP5Reader", "OpenFiles",
                                                  "File " + m_Name +
                                                      " cannot be opened: " + lasterrmsg);
        }
        else
        {
            helper::Throw<std::ios_base::failure>("Engine", "BP5Reader", "OpenFiles",
                                                  "File " + m_Name + " cannot be opened");
        }
    }
    else if (flag == 1)
    {
        if (m_Comm.Rank() == 0)
        {
            helper::Throw<std::ios_base::failure>(
                "Engine", "BP5Reader", "OpenFiles",
                "File " + m_Name + " could not be found within the " +
                    std::to_string(timeoutSeconds.count()) + "s timeout: " + lasterrmsg);
        }
        else
        {
            helper::Throw<std::ios_base::failure>(
                "Engine", "BP5Reader", "OpenFiles",
                "File " + m_Name + " could not be found within the " +
                    std::to_string(timeoutSeconds.count()) + "s timeout");
        }
    }

    /* At this point we may have an empty index table.
     * The writer has created the file but no content may have been stored yet.
     */
}

MinVarInfo *BP5Reader::MinBlocksInfo(const VariableBase &Var, const size_t Step) const
{
    return m_BP5Deserializer->MinBlocksInfo(Var, Step);
}

bool BP5Reader::VarShape(const VariableBase &Var, const size_t Step, Dims &Shape) const
{
    return m_BP5Deserializer->VarShape(Var, Step, Shape);
}

bool BP5Reader::VariableMinMax(const VariableBase &Var, const size_t Step, MinMaxStruct &MinMax)
{
    return m_BP5Deserializer->VariableMinMax(Var, Step, MinMax);
}

void BP5Reader::InitTransports()
{
    if (m_IO.m_TransportsParameters.empty())
    {
        Params defaultTransportParameters;
        defaultTransportParameters["transport"] = "File";
        m_IO.m_TransportsParameters.push_back(defaultTransportParameters);
    }
}

void BP5Reader::InstallMetaMetaData(format::BufferSTL buffer)
{
    size_t Position = m_MetaMetaDataFileAlreadyProcessedSize;
    while (Position < buffer.m_Buffer.size())
    {
        format::BP5Base::MetaMetaInfoBlock MMI;

        MMI.MetaMetaIDLen =
            helper::ReadValue<uint64_t>(buffer.m_Buffer, Position, m_Minifooter.IsLittleEndian);
        MMI.MetaMetaInfoLen =
            helper::ReadValue<uint64_t>(buffer.m_Buffer, Position, m_Minifooter.IsLittleEndian);
        MMI.MetaMetaID = buffer.Data() + Position;
        MMI.MetaMetaInfo = buffer.Data() + Position + MMI.MetaMetaIDLen;
        m_BP5Deserializer->InstallMetaMetaData(MMI);
        Position += MMI.MetaMetaIDLen + MMI.MetaMetaInfoLen;
    }
    m_MetaMetaDataFileAlreadyProcessedSize = Position;
}

void BP5Reader::UpdateBuffer(const TimePoint &timeoutInstant, const Seconds &pollSeconds,
                             const Seconds &timeoutSeconds)
{
    size_t newIdxSize = 0;
    m_MetadataIndex.Reset(true, false);
    if (m_Comm.Rank() == 0)
    {
        /* Read metadata index table into memory */
        const size_t metadataIndexFileSize = m_MDIndexFileManager.GetFileSize(0);
        newIdxSize = metadataIndexFileSize - m_MDIndexFileAlreadyReadSize;
        if (metadataIndexFileSize > m_MDIndexFileAlreadyReadSize)
        {
            m_MetadataIndex.m_Buffer.resize(newIdxSize);
            m_MDIndexFileManager.ReadFile(m_MetadataIndex.m_Buffer.data(), newIdxSize,
                                          m_MDIndexFileAlreadyReadSize);
        }
        else
        {
            m_MetadataIndex.m_Buffer.resize(0);
        }
    }

    // broadcast metadata index buffer to all ranks from zero
    m_Comm.BroadcastVector(m_MetadataIndex.m_Buffer);
    newIdxSize = m_MetadataIndex.m_Buffer.size();

    size_t parsedIdxSize = 0;
    const auto stepsBefore = m_StepsCount;
    if (newIdxSize > 0)
    {
        /* Parse metadata index table */
        const bool hasHeader = (!m_MDIndexFileAlreadyReadSize);
        parsedIdxSize = ParseMetadataIndex(m_MetadataIndex, 0, hasHeader);
        // now we are sure the index header has been parsed,
        // first step parsing done
        // m_FilteredMetadataInfo is created

        // cut down the index buffer by throwing away the read but unprocessed
        // steps
        m_MetadataIndex.m_Buffer.resize(parsedIdxSize);
        // next time read index file from this position
        m_MDIndexFileAlreadyReadSize += parsedIdxSize;

        // At this point first in time we learned the writer's major and we can
        // create the serializer object
        if (!m_BP5Deserializer)
        {
            m_BP5Deserializer =
                new format::BP5Deserializer(m_WriterIsRowMajor, m_ReaderIsRowMajor,
                                            (m_OpenMode != Mode::Read), (m_FlattenSteps));
            m_BP5Deserializer->m_Engine = this;
        }
    }

    if (m_StepsCount > stepsBefore)
    {
        m_Metadata.Reset(true, false);
        m_MetaMetadata.Reset(true, false);
        if (m_Comm.Rank() == 0)
        {
            // How much metadata do we need to read?
            size_t fileFilteredSize = 0;
            for (auto p : m_FilteredMetadataInfo)
            {
                fileFilteredSize += p.second;
            }

            /* Read metadata file into memory but first make sure
             * it has the content that the index table refers to */
            auto p = m_FilteredMetadataInfo.back();
            uint64_t expectedMinFileSize = p.first + p.second;
            size_t actualFileSize = 0;
            do
            {
                actualFileSize = m_MDFileManager.GetFileSize(0);
                if (actualFileSize >= expectedMinFileSize)
                {
                    break;
                }
            } while (SleepOrQuit(timeoutInstant, pollSeconds));

            if (actualFileSize >= expectedMinFileSize)
            {
                m_JSONProfiler.Start("MetaDataRead");
                m_Metadata.Resize(fileFilteredSize, "allocating metadata buffer, "
                                                    "in call to BP5Reader Open");
                size_t mempos = 0;
                for (auto p : m_FilteredMetadataInfo)
                {
                    m_JSONProfiler.AddBytes("metadataread", p.second);
                    m_MDFileManager.ReadFile(m_Metadata.Data() + mempos, p.second, p.first);
                    mempos += p.second;
                }
                m_MDFileAlreadyReadSize = expectedMinFileSize;
                m_JSONProfiler.Stop("MetaDataRead");
            }
            else
            {
                helper::Throw<std::ios_base::failure>(
                    "Engine", "BP5Reader", "UpdateBuffer",
                    "File " + m_Name +
                        " was found with an index file but md.0 "
                        "has not contained enough data within "
                        "the specified timeout of " +
                        std::to_string(timeoutSeconds.count()) +
                        " seconds. index size = " + std::to_string(newIdxSize) +
                        " metadata size = " + std::to_string(actualFileSize) +
                        " expected size = " + std::to_string(expectedMinFileSize) +
                        ". One reason could be if the reader finds old "
                        "data "
                        "while "
                        "the writer is creating the new files.");
            }

            /* Read new meta-meta-data into memory and append to existing one in
             * memory */
            const size_t metametadataFileSize = m_FileMetaMetadataManager.GetFileSize(0);
            if (metametadataFileSize > m_MetaMetaDataFileAlreadyReadSize)
            {
                const size_t newMMDSize = metametadataFileSize - m_MetaMetaDataFileAlreadyReadSize;
                m_JSONProfiler.Start("MetaMetaDataRead");
                m_JSONProfiler.AddBytes("metametadataread", newMMDSize);
                m_MetaMetadata.Resize(metametadataFileSize, "(re)allocating meta-meta-data buffer, "
                                                            "in call to BP5Reader Open");
                m_FileMetaMetadataManager.ReadFile(m_MetaMetadata.m_Buffer.data() +
                                                       m_MetaMetaDataFileAlreadyReadSize,
                                                   newMMDSize, m_MetaMetaDataFileAlreadyReadSize);
                m_MetaMetaDataFileAlreadyReadSize += newMMDSize;
                m_JSONProfiler.Stop("MetaMetaDataRead");
            }
        }

        // broadcast metadata index buffer to all ranks from zero
        m_Comm.BroadcastVector(m_MetaMetadata.m_Buffer);

        InstallMetaMetaData(m_MetaMetadata);

        size_t inputSize = m_Comm.BroadcastValue(m_Metadata.Size(), 0);

        if (m_Comm.Rank() != 0)
        {
            m_Metadata.Resize(inputSize, "metadata broadcast");
        }

        m_Comm.Bcast(m_Metadata.Data(), inputSize, 0);

        if ((m_OpenMode == Mode::ReadRandomAccess) || m_FlattenSteps)
        {
            for (size_t Step = 0; Step < m_MetadataIndexTable.size(); Step++)
            {
                m_BP5Deserializer->SetupForStep(Step,
                                                m_WriterMap[m_WriterMapIndex[Step]].WriterCount);
                InstallMetadataForTimestep(Step);
            }
        }
    }
}

size_t BP5Reader::ParseMetadataIndex(format::BufferSTL &bufferSTL, const size_t absoluteStartPos,
                                     const bool hasHeader)
{
    const auto &buffer = bufferSTL.m_Buffer;
    size_t &position = bufferSTL.m_Position;

    if (hasHeader)
    {
        // Read header (64 bytes)
        // long version string
        position = m_VersionTagPosition;
        m_Minifooter.VersionTag.assign(&buffer[position], m_VersionTagLength);

        position = m_EndianFlagPosition;
        const uint8_t endianness = helper::ReadValue<uint8_t>(buffer, position);
        m_Minifooter.IsLittleEndian = (endianness == 0) ? true : false;
#ifndef ADIOS2_HAVE_ENDIAN_REVERSE
        if (helper::IsLittleEndian() != m_Minifooter.IsLittleEndian)
        {
            helper::Throw<std::runtime_error>("Engine", "BP5Reader", "ParseMetadataIndex",
                                              "reader found BigEndian bp file, "
                                              "this version of ADIOS2 wasn't compiled "
                                              "with the cmake flag -DADIOS2_USE_Endian_Reverse=ON "
                                              "explicitly, in call to Open");
        }
#endif

        // This has no flag in BP5 header. Always true
        m_Minifooter.HasSubFiles = true;

        // BP version
        position = m_BPVersionPosition;
        m_Minifooter.Version =
            helper::ReadValue<uint8_t>(buffer, position, m_Minifooter.IsLittleEndian);
        if (m_Minifooter.Version != 5)
        {
            helper::Throw<std::runtime_error>("Engine", "BP5Reader", "ParseMetadataIndex",
                                              "ADIOS2 BP5 Engine only supports bp format "
                                              "version 5, found " +
                                                  std::to_string(m_Minifooter.Version) +
                                                  " version");
        }

        // BP minor version, unused
        position = m_BPMinorVersionPosition;
        const uint8_t minorversion =
            helper::ReadValue<uint8_t>(buffer, position, m_Minifooter.IsLittleEndian);
        if (minorversion != m_BP5MinorVersion)
        {
            helper::Throw<std::runtime_error>("Engine", "BP5Reader", "ParseMetadataIndex",
                                              "Current ADIOS2 BP5 Engine only supports version 5." +
                                                  std::to_string(m_BP5MinorVersion) + ", found 5." +
                                                  std::to_string(minorversion) + " version");
        }

        // Writer active flag
        position = m_ActiveFlagPosition;
        const char activeChar =
            helper::ReadValue<uint8_t>(buffer, position, m_Minifooter.IsLittleEndian);
        m_WriterIsActive = (activeChar == '\1' ? true : false);

        position = m_ColumnMajorFlagPosition;
        const uint8_t val =
            helper::ReadValue<uint8_t>(buffer, position, m_Minifooter.IsLittleEndian);
        m_WriterIsRowMajor = val == 'n';

        position = m_FlattenStepsPosition;
        const uint8_t flatten_val =
            helper::ReadValue<uint8_t>(buffer, position, m_Minifooter.IsLittleEndian);
        m_FlattenSteps = (flatten_val != 0);

        if (m_Parameters.IgnoreFlattenSteps)
            m_FlattenSteps = false;

        // move position to first row
        position = m_IndexHeaderSize;
    }

    // set a limit for metadata size in streaming mode
    size_t maxMetadataSizeInMemory = adios2::MaxSizeT;
    if (m_OpenMode == Mode::Read)
    {
        maxMetadataSizeInMemory = 16777216; // 16MB
    }
    size_t metadataSizeToRead = 0;

    // Read each record now
    uint64_t MetadataPosTotalSkip = 0;
    m_MetadataIndexTable.clear();
    m_FilteredMetadataInfo.clear();
    uint64_t minfo_pos = 0;
    uint64_t minfo_size = 0;
    int n = 0;    // a loop counter for current run4
    int nrec = 0; // number of records in current run

    while (position < buffer.size() && metadataSizeToRead < maxMetadataSizeInMemory)
    {

        const unsigned char recordID =
            helper::ReadValue<unsigned char>(buffer, position, m_Minifooter.IsLittleEndian);
        const uint64_t recordLength =
            helper::ReadValue<uint64_t>(buffer, position, m_Minifooter.IsLittleEndian);
        const size_t dbgRecordStartPosition = position;

        switch (recordID)
        {
        case IndexRecord::WriterMapRecord: {
            auto p = m_WriterMap.emplace(m_StepsCount, WriterMapStruct());
            auto &s = p.first->second;
            s.WriterCount = (uint32_t)helper::ReadValue<uint64_t>(buffer, position,
                                                                  m_Minifooter.IsLittleEndian);
            s.AggregatorCount = (uint32_t)helper::ReadValue<uint64_t>(buffer, position,
                                                                      m_Minifooter.IsLittleEndian);
            s.SubfileCount = (uint32_t)helper::ReadValue<uint64_t>(buffer, position,
                                                                   m_Minifooter.IsLittleEndian);
            // Get the process -> subfile map
            s.RankToSubfile.reserve(s.WriterCount);
            for (uint64_t i = 0; i < s.WriterCount; i++)
            {
                const uint64_t subfileIdx =
                    helper::ReadValue<uint64_t>(buffer, position, m_Minifooter.IsLittleEndian);
                s.RankToSubfile.push_back(subfileIdx);
            }
            m_LastMapStep = m_StepsCount;
            m_LastWriterCount = s.WriterCount;
            break;
        }
        case IndexRecord::StepRecord: {
            std::vector<uint64_t> ptrs;
            const uint64_t MetadataPos =
                helper::ReadValue<uint64_t>(buffer, position, m_Minifooter.IsLittleEndian);
            const uint64_t MetadataSize =
                helper::ReadValue<uint64_t>(buffer, position, m_Minifooter.IsLittleEndian);
            const uint64_t FlushCount =
                helper::ReadValue<uint64_t>(buffer, position, m_Minifooter.IsLittleEndian);

            if (!n)
            {
                minfo_pos = MetadataPos; // initialize minfo_pos properly
                MetadataPosTotalSkip = MetadataPos;
            }

            if (m_SelectedSteps.IsSelected(m_AbsStepsInFile))
            {
                m_WriterMapIndex.push_back(m_LastMapStep);

                // pos in metadata in memory
                ptrs.push_back(MetadataPos - MetadataPosTotalSkip);
                ptrs.push_back(MetadataSize);
                ptrs.push_back(FlushCount);
                ptrs.push_back(position);
                // absolute pos in file before read
                ptrs.push_back(MetadataPos);
                m_MetadataIndexTable[m_StepsCount] = ptrs;
#ifdef DUMPDATALOCINFO
                for (uint64_t i = 0; i < m_WriterCount; i++)
                {
                    size_t DataPosPos = ptrs[3];
                    std::cout << "Writer " << i << " data at ";
                    for (uint64_t j = 0; j < FlushCount; j++)
                    {
                        const uint64_t DataPos = helper::ReadValue<uint64_t>(
                            buffer, DataPosPos, m_Minifooter.IsLittleEndian);
                        const uint64_t DataSize = helper::ReadValue<uint64_t>(
                            buffer, DataPosPos, m_Minifooter.IsLittleEndian);
                        std::cout << "loc:" << DataPos << " siz:" << DataSize << "; ";
                    }
                    const uint64_t DataPos = helper::ReadValue<uint64_t>(
                        buffer, DataPosPos, m_Minifooter.IsLittleEndian);
                    std::cout << "loc:" << DataPos << std::endl;
                }
#endif
                minfo_size += MetadataSize;
                metadataSizeToRead += MetadataSize;
                m_StepsCount++;
            }
            else
            {
                MetadataPosTotalSkip += MetadataSize;
                if (minfo_size > 0)
                {
                    m_FilteredMetadataInfo.push_back(std::make_pair(minfo_pos, minfo_size));
                }
                minfo_pos = MetadataPos + MetadataSize;
                minfo_size = 0;
            }

            // skip over the writer -> data file offset records
            position += sizeof(uint64_t) * m_LastWriterCount * ((2 * FlushCount) + 1);
            ++m_AbsStepsInFile;
            ++n;
            break;
        }
        }
        // dbg
        if ((position - dbgRecordStartPosition) != (size_t)recordLength)
        {
            helper::Throw<std::runtime_error>(
                "Engine", "BP5Reader", "ParseMetadataIndex",
                "Record " + std::to_string(nrec) + " (id = " + std::to_string(recordID) +
                    ") has invalid length " + std::to_string(recordLength) + ". We parsed " +
                    std::to_string(position - dbgRecordStartPosition) + " bytes for this record"

            );
        }
        ++nrec;
    }
    if (minfo_size > 0)
    {
        m_FilteredMetadataInfo.push_back(std::make_pair(minfo_pos, minfo_size));
    }

    return position;
}

bool BP5Reader::ReadActiveFlag(std::vector<char> &buffer)
{
    if (buffer.size() < m_ActiveFlagPosition)
    {
        helper::Throw<std::runtime_error>("Engine", "BP5Reader", "ReadActiveFlag",
                                          "called with a buffer smaller than required");
    }
    // Writer active flag
    size_t position = m_ActiveFlagPosition;
    const char activeChar =
        helper::ReadValue<uint8_t>(buffer, position, m_Minifooter.IsLittleEndian);
    m_WriterIsActive = (activeChar == '\1' ? true : false);
    return m_WriterIsActive;
}

bool BP5Reader::CheckWriterActive()
{
    size_t flag = 1;
    if (m_Comm.Rank() == 0)
    {
        auto fsize = m_MDIndexFileManager.GetFileSize(0);
        if (fsize >= m_IndexHeaderSize)
        {
            std::vector<char> header(m_IndexHeaderSize, '\0');
            m_MDIndexFileManager.ReadFile(header.data(), m_IndexHeaderSize, 0, 0);
            bool active = ReadActiveFlag(header);
            flag = (active ? 1 : 0);
        }
    }
    flag = m_Comm.BroadcastValue(flag, 0);
    m_WriterIsActive = (flag > 0);
    return m_WriterIsActive;
}

StepStatus BP5Reader::CheckForNewSteps(Seconds timeoutSeconds)
{
    /* Do a collective wait for a step within timeout.
       Make sure every reader comes to the same conclusion */
    StepStatus retval = StepStatus::OK;

    if (timeoutSeconds < Seconds::zero())
    {
        timeoutSeconds = Seconds(999999999); // max 1 billion seconds wait
    }
    const TimePoint timeoutInstant = Now() + timeoutSeconds;

    auto pollSeconds = Seconds(m_Parameters.BeginStepPollingFrequencySecs);
    if (pollSeconds > timeoutSeconds)
    {
        pollSeconds = timeoutSeconds;
    }

    /* Poll */
    const auto stepsBefore = m_StepsCount;
    do
    {
        UpdateBuffer(timeoutInstant, pollSeconds / 10, timeoutSeconds);
        if (m_StepsCount > stepsBefore)
        {
            break;
        }
        if (!CheckWriterActive())
        {
            /* Race condition: When checking data in UpdateBuffer, new
             * step(s) may have not arrived yet. When checking active flag,
             * the writer may have completed write and terminated. So we may
             * have missed a step or two. */
            UpdateBuffer(timeoutInstant, pollSeconds / 10, timeoutSeconds);
            break;
        }
    } while (SleepOrQuit(timeoutInstant, pollSeconds));

    if (m_StepsCount > stepsBefore)
    {
        /* we have got new steps and new metadata in memory */
        retval = StepStatus::OK;
    }
    else
    {
        m_IO.m_ReadStreaming = false;
        if (m_WriterIsActive)
        {
            retval = StepStatus::NotReady;
        }
        else
        {
            retval = StepStatus::EndOfStream;
        }
    }
    return retval;
}

void BP5Reader::DoGetAbsoluteSteps(const VariableBase &variable, std::vector<size_t> &keys) const
{
    m_BP5Deserializer->GetAbsoluteSteps(variable, keys);
    return;
}

#define declare_type(T)                                                                            \
    void BP5Reader::DoGetSync(Variable<T> &variable, T *data)                                      \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP5Reader::Get");                                                  \
        GetSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void BP5Reader::DoGetDeferred(Variable<T> &variable, T *data)                                  \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP5Reader::Get");                                                  \
        GetDeferredCommon(variable, data);                                                         \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void BP5Reader::DoGetStructSync(VariableStruct &variable, void *data)
{
    PERFSTUBS_SCOPED_TIMER("BP5Reader::Get");
    GetSyncCommon(variable, data);
}

void BP5Reader::DoGetStructDeferred(VariableStruct &variable, void *data)
{
    PERFSTUBS_SCOPED_TIMER("BP5Reader::Get");
    GetDeferredCommon(variable, data);
}

void BP5Reader::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP5Reader::Close");
    if (m_OpenMode == Mode::ReadRandomAccess)
    {
        PerformGets();
    }
    else if (m_BetweenStepPairs)
    {
        EndStep();
    }
    FlushProfiler();
    m_DataFileManager.CloseFiles();
    m_MDFileManager.CloseFiles();
    m_MDIndexFileManager.CloseFiles();
    m_FileMetaMetadataManager.CloseFiles();
    for (unsigned int i = 1; i < m_Threads; ++i)
    {
        fileManagers[i].CloseFiles();
    }
}

#if defined(_WIN32)
#include <windows.h>
#define getpid() GetCurrentProcessId();
#elif defined(__linux) || defined(__APPLE__) || defined(__OpenBSD__) || defined(__FreeBSD__) ||    \
    defined(__NetBSD__) || defined(__DragonFly__) || defined(__CYGWIN__)
#include <sys/types.h>
#include <unistd.h>
#else
#define getpid() (long long)-1;
#endif

void BP5Reader::FlushProfiler()
{

    auto LineJSON = m_JSONProfiler.GetRankProfilingJSON({}, {});
    const std::vector<char> profilingJSON(m_JSONProfiler.AggregateProfilingJSON(LineJSON));

    if (m_RankMPI == 0)
    {
        std::string profileFileName;
        transport::FileFStream profilingJSONStream(m_Comm);
        std::string bpBaseName = m_Name;

        auto PID = getpid();
        std::stringstream PIDstr;
        PIDstr << std::hex << PID;
        // write profile json in /tmp
        profileFileName = "/tmp/" + bpBaseName + "_" + PIDstr.str() + "_profiling.json";

        try
        {
            (void)remove(profileFileName.c_str());
            profilingJSONStream.Open(profileFileName, Mode::Write);
            profilingJSONStream.Write(profilingJSON.data(), profilingJSON.size());
            profilingJSONStream.Close();
        }
        catch (...)
        { // do nothing
        }
    }
}

size_t BP5Reader::DoSteps() const
{
    if (m_FlattenSteps)
        return 1;
    else
        return m_StepsCount;
}

void BP5Reader::NotifyEngineNoVarsQuery()
{
    if (!m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>(
            "Engine", "BP5Reader", "NotifyEngineNoVarsQuery",
            "You've called InquireVariable() when the IO is empty and "
            "outside a BeginStep/EndStep pair.  If this is code that is "
            "newly "
            "transititioning to the BP5 file engine, you may be relying "
            "upon "
            "deprecated behaviour.  If you intend to use ADIOS using the "
            "Begin/EndStep interface, move all InquireVariable calls "
            "inside "
            "the BeginStep/EndStep pair.  If intending to use "
            "random-access "
            "file mode, change your Open() mode parameter to "
            "Mode::ReadRandomAccess.");
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
