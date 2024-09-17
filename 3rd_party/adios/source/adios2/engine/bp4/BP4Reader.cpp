/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Reader.cpp
 *
 *  Created on: Aug 1, 2018
 *      Author: Lipeng Wan wanl@ornl.gov
 */

#include "BP4Reader.h"
#include "BP4Reader.tcc"

#include <adios2-perfstubs-interface.h>

#include <chrono>
#include <errno.h>

namespace adios2
{
namespace core
{
namespace engine
{

BP4Reader::BP4Reader(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("BP4Reader", io, name, mode, std::move(comm)), m_BP4Deserializer(m_Comm),
  m_MDFileManager(io, m_Comm), m_DataFileManager(io, m_Comm), m_MDIndexFileManager(io, m_Comm),
  m_ActiveFlagFileManager(io, m_Comm)
{
    PERFSTUBS_SCOPED_TIMER("BP4Reader::Open");
    helper::GetParameter(m_IO.m_Parameters, "Verbose", m_Verbosity);
    helper::Log("Engine", "BP4Reader", "Open", m_Name, 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);
    Init();
    m_IsOpen = true;
}

BP4Reader::~BP4Reader()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus BP4Reader::BeginStep(StepMode mode, const float timeoutSeconds)
{
    PERFSTUBS_SCOPED_TIMER("BP4Reader::BeginStep");
    helper::Log("Engine", "BP4Reader", "BeginStep", std::to_string(CurrentStep()), 0, m_Comm.Rank(),
                5, m_Verbosity, helper::LogMode::INFO);

    if (mode != StepMode::Read)
    {
        helper::Throw<std::invalid_argument>("Engine", "BP4Reader", "BeginStep",
                                             "mode is not supported yet, "
                                             "only Read is valid for "
                                             "engine BP4Reader, in call to "
                                             "BeginStep");
    }

    if (m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "BP4Reader", "BeginStep",
                                        "BeginStep() is called a second time "
                                        "without an intervening EndStep()");
    }

    if (!m_BP4Deserializer.m_DeferredVariables.empty())
    {
        helper::Throw<std::invalid_argument>("Engine", "BP4Reader", "BeginStep",
                                             "existing variables subscribed with "
                                             "GetDeferred, did you forget to call "
                                             "PerformGets() or EndStep()?, in call to BeginStep");
    }

    // used to inquire for variables in streaming mode
    m_IO.m_ReadStreaming = true;
    StepStatus status = StepStatus::OK;

    if (m_FirstStep)
    {
        if (m_BP4Deserializer.m_MetadataSet.StepsCount == 0)
        {
            status = CheckForNewSteps(Seconds(timeoutSeconds));
        }
    }
    else
    {
        if (m_CurrentStep + 1 >= m_BP4Deserializer.m_MetadataSet.StepsCount)
        {
            status = CheckForNewSteps(Seconds(timeoutSeconds));
        }
    }

    // This should be after getting new steps

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
        m_IO.ResetVariablesStepSelection(false, "in call to BP4 Reader BeginStep");

        // caches attributes for each step
        // if a variable name is a prefix
        // e.g. var  prefix = {var/v1, var/v2, var/v3}
        m_IO.SetPrefixedNames(true);
    }

    return status;
}

size_t BP4Reader::CurrentStep() const { return m_CurrentStep; }

void BP4Reader::EndStep()
{
    helper::Log("Engine", "BP4Reader", "EndStep", std::to_string(CurrentStep()), 0, m_Comm.Rank(),
                5, m_Verbosity, helper::LogMode::INFO);
    if (!m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "BP4Reader", "EndStep",
                                        "EndStep() is called without a successful BeginStep()");
    }
    m_BetweenStepPairs = false;
    PERFSTUBS_SCOPED_TIMER("BP4Reader::EndStep");
    PerformGets();
}

void BP4Reader::PerformGets()
{
    PERFSTUBS_SCOPED_TIMER("BP4Reader::PerformGets");
    helper::Log("Engine", "BP4Reader", "PerformGets", "", 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);
    if (m_BP4Deserializer.m_DeferredVariables.empty())
    {
        return;
    }

    for (const std::string &name : m_BP4Deserializer.m_DeferredVariables)
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
            m_BP4Deserializer.SetVariableBlockInfo(variable, blockInfo);                           \
        }                                                                                          \
        ReadVariableBlocks(variable);                                                              \
        variable.m_BlocksInfo.clear();                                                             \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    }

    m_BP4Deserializer.m_DeferredVariables.clear();
}

// PRIVATE
void BP4Reader::Init()
{
    if (m_OpenMode != Mode::Read)
    {
        helper::Throw<std::invalid_argument>("Engine", "BP4Reader", "Init",
                                             "BPFileReader only "
                                             "supports OpenMode::Read from" +
                                                 m_Name);
    }
    // if IO was involved in reading before this flag may be true now
    m_IO.m_ReadStreaming = false;

    m_BP4Deserializer.Init(m_IO.m_Parameters, "in call to BP4::Open to write");
    InitTransports();

    helper::RaiseLimitNoFile();

    /* Do a collective wait for the file(s) to appear within timeout.
       Make sure every process comes to the same conclusion */
    const Seconds timeoutSeconds(m_BP4Deserializer.m_Parameters.OpenTimeoutSecs);

    Seconds pollSeconds(m_BP4Deserializer.m_Parameters.BeginStepPollingFrequencySecs);
    if (pollSeconds > timeoutSeconds)
    {
        pollSeconds = timeoutSeconds;
    }

    TimePoint timeoutInstant = Now() + timeoutSeconds;

    OpenFiles(timeoutInstant, pollSeconds, timeoutSeconds);
    if (!m_BP4Deserializer.m_Parameters.StreamReader)
    {
        /* non-stream reader gets as much steps as available now */
        InitBuffer(timeoutInstant, pollSeconds / 10, timeoutSeconds);
    }
}

bool BP4Reader::SleepOrQuit(const TimePoint &timeoutInstant, const Seconds &pollSeconds)
{
    auto now = Now();
    if (now + pollSeconds >= timeoutInstant)
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

size_t BP4Reader::OpenWithTimeout(transportman::TransportMan &tm,
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
            const bool profile = m_BP4Deserializer.m_Profiler.m_IsActive;

            for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
            {
                m_IO.m_TransportsParameters[i].insert({"SingleProcess", "true"});
            }
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

void BP4Reader::OpenFiles(TimePoint &timeoutInstant, const Seconds &pollSeconds,
                          const Seconds &timeoutSeconds)
{
    /* Poll */
    size_t flag = 1; // 0 = OK, opened file, 1 = timeout, 2 = error
    std::string lasterrmsg;
    if (m_BP4Deserializer.m_RankMPI == 0)
    {
        /* Open the metadata index table */
        const std::string metadataIndexFile(m_BP4Deserializer.GetBPMetadataIndexFileName(m_Name));
        flag = OpenWithTimeout(m_MDIndexFileManager, {metadataIndexFile}, timeoutInstant,
                               pollSeconds, lasterrmsg);
        if (flag == 0)
        {
            /* Open the metadata file */
            const std::string metadataFile(m_BP4Deserializer.GetBPMetadataFileName(m_Name));

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
        }
    }
    m_Comm.Barrier("wait for rank 0 to open...");
    flag = m_Comm.BroadcastValue(flag, 0);
    if (flag == 2)
    {
        if (m_BP4Deserializer.m_RankMPI == 0 && !lasterrmsg.empty())
        {
            helper::Throw<std::ios_base::failure>("Engine", "BP4Reader", "OpenFiles",
                                                  "File " + m_Name +
                                                      " cannot be opened: " + lasterrmsg);
        }
        else
        {
            helper::Throw<std::ios_base::failure>("Engine", "BP4Reader", "OpenFiles",
                                                  "File " + m_Name + " cannot be opened");
        }
    }
    else if (flag == 1)
    {
        if (m_BP4Deserializer.m_RankMPI == 0)
        {
            helper::Throw<std::ios_base::failure>(
                "Engine", "BP4Reader", "OpenFiles",
                "File " + m_Name + " could not be found within the " +
                    std::to_string(timeoutSeconds.count()) + "s timeout: " + lasterrmsg);
        }
        else
        {
            helper::Throw<std::ios_base::failure>(
                "Engine", "BP4Reader", "OpenFiles",
                "File " + m_Name + " could not be found within the " +
                    std::to_string(timeoutSeconds.count()) + "s timeout");
        }
    }

    /* At this point we may have an empty index table.
     * The writer has created the file but no content may have been stored yet.
     */
}

void BP4Reader::InitTransports()
{
    if (m_IO.m_TransportsParameters.empty())
    {
        Params defaultTransportParameters;
        defaultTransportParameters["transport"] = "File";
        m_IO.m_TransportsParameters.push_back(defaultTransportParameters);
    }
}

/* Count index records to minimum 1 and maximum of N records so that
 * expected metadata size is less then a predetermined constant
 */
void MetadataCalculateMinFileSize(const format::BP4Deserializer &m_BP4Deserializer,
                                  const std::string &IdxFileName, char *buf, size_t idxsize,
                                  bool hasHeader, const size_t mdStartPos, size_t &newIdxSize,
                                  size_t &expectedMinFileSize)
{
    newIdxSize = 0;
    expectedMinFileSize = 0;

    if (hasHeader && idxsize < m_BP4Deserializer.m_IndexRecordSize)
    {
        return;
    }

    /* eliminate header for now for only calculating with records */
    if (hasHeader)
    {
        buf += m_BP4Deserializer.m_IndexRecordSize;
        idxsize -= m_BP4Deserializer.m_IndexRecordSize;
    }

    if (idxsize % m_BP4Deserializer.m_IndexRecordSize != 0)
    {
        helper::Throw<std::runtime_error>("Engine", "BP4Reader", "MetadataCalculateMinFileSize",
                                          "ADIOS Index file " + IdxFileName +
                                              " is assumed to always contain n*" +
                                              std::to_string(m_BP4Deserializer.m_IndexRecordSize) +
                                              " byte-length records. "
                                              "Right now the length of index buffer is " +
                                              std::to_string(idxsize) + " bytes.");
    }

    const size_t nTotalRecords = idxsize / m_BP4Deserializer.m_IndexRecordSize;
    if (nTotalRecords == 0)
    {
        // no (new) step entry in the index, so no metadata is expected
        newIdxSize = 0;
        expectedMinFileSize = 0;
        return;
    }

    size_t nRecords = 1;
    expectedMinFileSize = *(uint64_t *)&(buf[nRecords * m_BP4Deserializer.m_IndexRecordSize - 24]);
    while (nRecords < nTotalRecords)
    {
        const size_t n = nRecords + 1;
        const uint64_t mdEndPos = *(uint64_t *)&(buf[n * m_BP4Deserializer.m_IndexRecordSize - 24]);
        if (mdEndPos - mdStartPos > 16777216)
        {
            break;
        }
        expectedMinFileSize = mdEndPos;
        ++nRecords;
    }
    newIdxSize = nRecords * m_BP4Deserializer.m_IndexRecordSize;
    if (hasHeader)
    {
        newIdxSize += m_BP4Deserializer.m_IndexRecordSize;
    }
}

uint64_t MetadataExpectedMinFileSize(const format::BP4Deserializer &m_BP4Deserializer,
                                     const std::string &IdxFileName, bool hasHeader)
{
    size_t idxsize = m_BP4Deserializer.m_MetadataIndex.m_Buffer.size();
    if (idxsize % 64 != 0)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "BP4Reader", "MetadataExpectedMinFileSize",
            "ADIOS Index file " + IdxFileName +
                " is assumed to always contain n*64 byte-length records. "
                "The file size now is " +
                std::to_string(idxsize) + " bytes.");
    }
    if ((hasHeader &&
         idxsize < m_BP4Deserializer.m_IndexHeaderSize + m_BP4Deserializer.m_IndexRecordSize) ||
        idxsize < m_BP4Deserializer.m_IndexRecordSize)
    {
        // no (new) step entry in the index, so no metadata is expected
        return 0;
    }
    uint64_t lastpos = *(uint64_t *)&(m_BP4Deserializer.m_MetadataIndex.m_Buffer[idxsize - 24]);
    return lastpos;
}

void BP4Reader::InitBuffer(const TimePoint &timeoutInstant, const Seconds &pollSeconds,
                           const Seconds &timeoutSeconds)
{
    size_t newIdxSize = 0;
    // Put all metadata in buffer
    if (m_BP4Deserializer.m_RankMPI == 0)
    {
        /* Read metadata index table into memory */
        const size_t metadataIndexFileSize = m_MDIndexFileManager.GetFileSize(0);
        if (metadataIndexFileSize > 0)
        {
            m_BP4Deserializer.m_MetadataIndex.Resize(metadataIndexFileSize,
                                                     "allocating metadata index buffer, "
                                                     "in call to BPFileReader Open");
            m_MDIndexFileManager.ReadFile(m_BP4Deserializer.m_MetadataIndex.m_Buffer.data(),
                                          metadataIndexFileSize);

            /* Read metadata file into memory but first make sure
             * it has the content that the index table refers to */
            uint64_t expectedMinFileSize =
                MetadataExpectedMinFileSize(m_BP4Deserializer, m_Name, true);
            size_t fileSize = 0;
            do
            {
                fileSize = m_MDFileManager.GetFileSize(0);
                if (fileSize >= expectedMinFileSize)
                {
                    break;
                }
            } while (SleepOrQuit(timeoutInstant, pollSeconds));

            if (fileSize >= expectedMinFileSize)
            {
                m_BP4Deserializer.m_Metadata.Resize(
                    expectedMinFileSize, "allocating metadata buffer, in call to BP4Reader Open");

                m_MDFileManager.ReadFile(m_BP4Deserializer.m_Metadata.m_Buffer.data(),
                                         expectedMinFileSize);
                m_MDFileAlreadyReadSize = expectedMinFileSize;
                m_MDIndexFileAlreadyReadSize = metadataIndexFileSize;
                newIdxSize = metadataIndexFileSize;
            }
            else
            {
                helper::Throw<std::ios_base::failure>(
                    "Engine", "BP4Reader", "InitBuffer",
                    "File " + m_Name +
                        " was found with an index file but md.0 "
                        "has not contained enough data within "
                        "the specified timeout of " +
                        std::to_string(timeoutSeconds.count()) +
                        " seconds. index size = " + std::to_string(metadataIndexFileSize) +
                        " metadata size = " + std::to_string(fileSize) +
                        " expected size = " + std::to_string(expectedMinFileSize) +
                        ". One reason could be if the reader finds old data "
                        "while "
                        "the writer is creating the new files.");
            }
        }
    }

    newIdxSize = m_Comm.BroadcastValue(newIdxSize, 0);

    if (newIdxSize > 0)
    {
        // broadcast buffer to all ranks from zero
        m_Comm.BroadcastVector(m_BP4Deserializer.m_Metadata.m_Buffer);

        // broadcast metadata index buffer to all ranks from zero
        m_Comm.BroadcastVector(m_BP4Deserializer.m_MetadataIndex.m_Buffer);

        /* Parse metadata index table */
        m_BP4Deserializer.ParseMetadataIndex(m_BP4Deserializer.m_MetadataIndex, 0, true, false);
        // now we are sure the index header has been parsed, first step parsing
        // done
        m_IdxHeaderParsed = true;

        // fills IO with Variables and Attributes
        m_MDFileProcessedSize =
            m_BP4Deserializer.ParseMetadata(m_BP4Deserializer.m_Metadata, *this, true);

        /* m_MDFileProcessedSize is the position in the buffer where processing
         * ends. The processing is controlled by the number of records in the
         * Index, which may be less than the actual entries in the metadata in a
         * streaming situation (where writer has just written metadata for step
         * K+1,...,K+L while the index contains K steps when the reader looks at
         * it).
         *
         * In ProcessMetadataForNewSteps(), we will re-read the metadata which
         * is in the buffer but has not been processed yet.
         */
    }
}

size_t BP4Reader::UpdateBuffer(const TimePoint &timeoutInstant, const Seconds &pollSeconds)
{
    std::vector<size_t> sizes(3, 0);
    if (m_BP4Deserializer.m_RankMPI == 0)
    {
        const size_t idxFileSize = m_MDIndexFileManager.GetFileSize(0);
        if (idxFileSize > m_MDIndexFileAlreadyReadSize)
        {
            const size_t maxIdxSize = idxFileSize - m_MDIndexFileAlreadyReadSize;
            std::vector<char> idxbuf(maxIdxSize);
            m_MDIndexFileManager.ReadFile(idxbuf.data(), maxIdxSize, m_MDIndexFileAlreadyReadSize);
            size_t newIdxSize;
            size_t expectedMinFileSize;
            char *buf = idxbuf.data();

            MetadataCalculateMinFileSize(m_BP4Deserializer, m_Name, buf, maxIdxSize,
                                         !m_IdxHeaderParsed, m_MDFileAlreadyReadSize, newIdxSize,
                                         expectedMinFileSize);

            // const uint64_t expectedMinFileSize = MetadataExpectedMinFileSize(
            //    m_BP4Deserializer, m_Name, !m_IdxHeaderParsed);

            if (m_BP4Deserializer.m_MetadataIndex.m_Buffer.size() < newIdxSize)
            {
                m_BP4Deserializer.m_MetadataIndex.Resize(
                    newIdxSize, "re-allocating metadata index buffer, in "
                                "call to BP4Reader::BeginStep/UpdateBuffer");
            }
            m_BP4Deserializer.m_MetadataIndex.Reset(true, false);
            std::copy(idxbuf.begin(), idxbuf.begin() + newIdxSize,
                      m_BP4Deserializer.m_MetadataIndex.m_Buffer.begin());

            /* Wait until as much metadata arrives in the file as much
             * is indicated by the existing index entries
             */

            size_t fileSize = 0;
            do
            {
                fileSize = m_MDFileManager.GetFileSize(0);
                if (fileSize >= expectedMinFileSize)
                {
                    break;
                }
            } while (SleepOrQuit(timeoutInstant, pollSeconds));

            if (fileSize >= expectedMinFileSize)
            {
                /* Read corresponding new metadata (throwing away the old)
                 * There may be unprocessed entries in the metadata if the index
                 * had less steps than the metadata file at the last read.
                 * Those steps are read again here, starting in the beginning of
                 * the buffer now.
                 */
                const size_t newMDSize = expectedMinFileSize - m_MDFileAlreadyReadSize;
                if (m_BP4Deserializer.m_Metadata.m_Buffer.size() < newMDSize)
                {
                    m_BP4Deserializer.m_Metadata.Resize(newMDSize,
                                                        "allocating metadata buffer, in call to "
                                                        "BP4Reader Open");
                }
                m_BP4Deserializer.m_Metadata.Reset(true, false);
                m_MDFileManager.ReadFile(m_BP4Deserializer.m_Metadata.m_Buffer.data(), newMDSize,
                                         m_MDFileAlreadyReadSize);

                m_MDFileAbsolutePos = m_MDFileAlreadyReadSize;
                m_MDFileAlreadyReadSize = expectedMinFileSize;

                m_MDIndexFileAlreadyReadSize += newIdxSize;

                sizes[0] = newIdxSize;
                sizes[1] = m_MDFileAlreadyReadSize;
                sizes[2] = m_MDFileAbsolutePos;
            }
        }
    }

    m_Comm.BroadcastVector(sizes, 0);
    size_t newIdxSize = sizes[0];

    if (newIdxSize > 0)
    {
        if (m_BP4Deserializer.m_RankMPI != 0)
        {
            m_MDFileAlreadyReadSize = sizes[1];
            m_MDFileAbsolutePos = sizes[2];
            m_BP4Deserializer.m_MetadataIndex.Reset(true, false);
            m_BP4Deserializer.m_Metadata.Reset(true, false);
            // we need this pointer in Metadata buffer on all processes
            // for parsing it correctly in ProcessMetadataForNewSteps()
        }

        // broadcast buffer to all ranks from zero
        m_Comm.BroadcastVector(m_BP4Deserializer.m_Metadata.m_Buffer);

        // broadcast metadata index buffer to all ranks from zero
        m_Comm.BroadcastVector(m_BP4Deserializer.m_MetadataIndex.m_Buffer);
    }
    return newIdxSize;
}
void BP4Reader::ProcessMetadataForNewSteps(const size_t newIdxSize)
{
    /* Remove all variables we created in the last step */
    RemoveCreatedVars();

    /* Parse metadata index table (without header) */
    /* We need to skew the index table pointers with the
       size of the already-processed metadata because the memory buffer of
       new metadata starts from 0 */
    m_BP4Deserializer.ParseMetadataIndex(m_BP4Deserializer.m_MetadataIndex, m_MDFileAbsolutePos,
                                         !m_IdxHeaderParsed, true);
    m_IdxHeaderParsed = true;

    // fills IO with Variables and Attributes
    const size_t newProcessedMDSize =
        m_BP4Deserializer.ParseMetadata(m_BP4Deserializer.m_Metadata, *this, false);

    // remember current end position in metadata and index table for next round
    m_MDFileProcessedSize = m_MDFileAbsolutePos + newProcessedMDSize;
    // if (m_BP4Deserializer.m_RankMPI == 0)
    //{
    //    m_MDIndexFileAlreadyReadSize += newIdxSize;
    //}
}

bool BP4Reader::CheckWriterActive()
{
    size_t flag = 0;
    if (m_BP4Deserializer.m_RankMPI == 0)
    {
        std::vector<char> header(m_BP4Deserializer.m_IndexHeaderSize, '\0');
        m_MDIndexFileManager.ReadFile(header.data(), m_BP4Deserializer.m_IndexHeaderSize, 0, 0);
        bool active = m_BP4Deserializer.ReadActiveFlag(header);
        flag = (active ? 1 : 0);
    }
    flag = m_BP4Deserializer.m_Comm.BroadcastValue(flag, 0);
    m_WriterIsActive = (flag > 0);
    return m_WriterIsActive;
}

bool BP4Reader::ProcessNextStepInMemory()
{
    if (m_MDFileAlreadyReadSize > m_MDFileProcessedSize)
    {
        // Hack: processing metadata for multiple new steps only works
        // when pretending not to be in streaming mode
        const bool saveReadStreaming = m_IO.m_ReadStreaming;
        m_IO.m_ReadStreaming = false;
        ProcessMetadataForNewSteps(0);
        m_IO.m_ReadStreaming = saveReadStreaming;
        return true;
    }
    return false;
}

StepStatus BP4Reader::CheckForNewSteps(Seconds timeoutSeconds)
{
    /* Do a collective wait for a step within timeout.
       Make sure every reader comes to the same conclusion */
    StepStatus retval = StepStatus::OK;

    if (ProcessNextStepInMemory())
    {
        return retval;
    }

    if (timeoutSeconds < Seconds::zero())
    {
        timeoutSeconds = Seconds(999999999); // max 1 billion seconds wait
    }
    const TimePoint timeoutInstant = Now() + timeoutSeconds;

    auto pollSeconds = Seconds(m_BP4Deserializer.m_Parameters.BeginStepPollingFrequencySecs);
    if (pollSeconds > timeoutSeconds)
    {
        pollSeconds = timeoutSeconds;
    }

    /* Poll */

    // Hack: processing metadata for multiple new steps only works
    // when pretending not to be in streaming mode
    const bool saveReadStreaming = m_IO.m_ReadStreaming;
    m_IO.m_ReadStreaming = false;
    size_t newIdxSize = 0;

    do
    {
        newIdxSize = UpdateBuffer(timeoutInstant, pollSeconds / 10);
        if (newIdxSize > 0)
        {
            break;
        }
        if (!CheckWriterActive())
        {
            /* Race condition: When checking data in UpdateBuffer, new step(s)
             * may have not arrived yet. When checking active flag, the writer
             * may have completed write and terminated. So we may have missed a
             * step or two. */
            newIdxSize = UpdateBuffer(timeoutInstant, pollSeconds / 10);
            break;
        }
    } while (SleepOrQuit(timeoutInstant, pollSeconds));

    if (newIdxSize > 0)
    {
        /* we have new metadata in memory. Need to parse it now */
        ProcessMetadataForNewSteps(newIdxSize);
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

    m_IO.m_ReadStreaming = saveReadStreaming;

    return retval;
}

#define declare_type(T)                                                                            \
    void BP4Reader::DoGetSync(Variable<T> &variable, T *data)                                      \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Reader::Get");                                                  \
        helper::Log("Engine", "BP4Reader", "GetSync", variable.m_Name, 0, m_Comm.Rank(), 5,        \
                    m_Verbosity, helper::LogMode::INFO);                                           \
        GetSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void BP4Reader::DoGetDeferred(Variable<T> &variable, T *data)                                  \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Reader::Get");                                                  \
        helper::Log("Engine", "BP4Reader", "GetDeferred", variable.m_Name, 0, m_Comm.Rank(), 5,    \
                    m_Verbosity, helper::LogMode::INFO);                                           \
        GetDeferredCommon(variable, data);                                                         \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void BP4Reader::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("BP4Reader::Close");
    helper::Log("Engine", "BP4Reader", "Close", m_Name, 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);
    PerformGets();
    /* Remove all variables we created in the last step */
    RemoveCreatedVars();

    m_DataFileManager.CloseFiles();
    m_MDFileManager.CloseFiles();
    m_MDIndexFileManager.CloseFiles();
}

#define declare_type(T)                                                                            \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> BP4Reader::DoAllStepsBlocksInfo(   \
        const Variable<T> &variable) const                                                         \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Reader::AllStepsBlocksInfo");                                   \
        return m_BP4Deserializer.AllStepsBlocksInfo(variable);                                     \
    }                                                                                              \
                                                                                                   \
    std::vector<std::vector<typename Variable<T>::BPInfo>>                                         \
    BP4Reader::DoAllRelativeStepsBlocksInfo(const Variable<T> &variable) const                     \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Reader::AllRelativeStepsBlocksInfo");                           \
        return m_BP4Deserializer.AllRelativeStepsBlocksInfo(variable);                             \
    }                                                                                              \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> BP4Reader::DoBlocksInfo(const Variable<T> &variable, \
                                                                      const size_t step) const     \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("BP4Reader::BlocksInfo");                                           \
        return m_BP4Deserializer.BlocksInfo(variable, step);                                       \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

size_t BP4Reader::DoSteps() const { return m_BP4Deserializer.m_MetadataSet.StepsCount; }

} // end namespace engine
} // end namespace core
} // end namespace adios2
