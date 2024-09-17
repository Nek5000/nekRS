/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DaosWriter.cpp
 *
 */

#include "DaosWriter.h"
#include "DaosWriter.tcc"

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosFunctions.h" //CheckIndexRange
#include "adios2/helper/adiosMath.h"      // SetWithinLimit
#include "adios2/helper/adiosMemory.h"    // NdCopy
#include "adios2/toolkit/format/buffer/chunk/ChunkV.h"
#include "adios2/toolkit/format/buffer/malloc/MallocV.h"
#include "adios2/toolkit/transport/file/FileFStream.h"
#include <adios2-perfstubs-interface.h>

#include <ctime>
#include <iomanip> // setw
#include <iostream>
#include <memory> // make_shared

namespace adios2
{
namespace core
{
namespace engine
{

using namespace adios2::format;

DaosWriter::DaosWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("DaosWriter", io, name, mode, std::move(comm)), m_BP5Serializer(),
  m_FileDataManager(io, m_Comm), m_FileMetadataManager(io, m_Comm),
  m_FileMetadataIndexManager(io, m_Comm), m_FileMetaMetadataManager(io, m_Comm), m_Profiler(m_Comm)
{
    m_EngineStart = Now();
    PERFSTUBS_SCOPED_TIMER("DaosWriter::Open");
    m_IO.m_ReadStreaming = false;

    Init();
    m_IsOpen = true;
}

StepStatus DaosWriter::BeginStep(StepMode mode, const float timeoutSeconds)
{
    if (m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "DaosWriter", "BeginStep",
                                        "BeginStep() is called a second time "
                                        "without an intervening EndStep()");
    }

    Seconds ts = Now() - m_EngineStart;
    // std::cout << "BEGIN STEP starts at: " << ts.count() << std::endl;
    m_BetweenStepPairs = true;

    if (m_WriterStep > 0)
    {
        m_LastTimeBetweenSteps = Now() - m_EndStepEnd;
        m_TotalTimeBetweenSteps += m_LastTimeBetweenSteps;
        m_AvgTimeBetweenSteps = m_TotalTimeBetweenSteps / m_WriterStep;
        m_ExpectedTimeBetweenSteps = m_LastTimeBetweenSteps;
        if (m_ExpectedTimeBetweenSteps > m_AvgTimeBetweenSteps)
        {
            m_ExpectedTimeBetweenSteps = m_AvgTimeBetweenSteps;
        }
    }

    if ((m_WriterStep == 0) && m_Parameters.UseOneTimeAttributes)
    {
        const auto &attributes = m_IO.GetAttributes();

        for (const auto &attributePair : attributes)
        {
            m_BP5Serializer.OnetimeMarshalAttribute(*(attributePair.second));
        }
    }

    if (m_Parameters.AsyncWrite)
    {
        m_AsyncWriteLock.lock();
        m_flagRush = true;
        m_AsyncWriteLock.unlock();
        TimePoint wait_start = Now();
        if (m_WriteFuture.valid())
        {
            m_Profiler.Start("WaitOnAsync");
            m_WriteFuture.get();
            m_Comm.Barrier();
            AsyncWriteDataCleanup();
            Seconds wait = Now() - wait_start;
            if (m_Comm.Rank() == 0)
            {
                WriteMetadataFileIndex(m_LatestMetaDataPos, m_LatestMetaDataSize);
                if (m_Parameters.verbose > 0)
                {
                    std::cout << "BeginStep, wait on async write was = " << wait.count()
                              << " time since EndStep was = " << m_LastTimeBetweenSteps.count()
                              << " expect next one to be = " << m_ExpectedTimeBetweenSteps.count()
                              << std::endl;
                }
            }
            m_Profiler.Stop("WaitOnAsync");
        }
    }

    if (m_Parameters.BufferVType == (int)BufferVType::MallocVType)
    {
        m_BP5Serializer.InitStep(new MallocV(
            "DaosWriter", false, m_BP5Serializer.m_BufferAlign, m_BP5Serializer.m_BufferBlockSize,
            m_Parameters.InitialBufferSize, m_Parameters.GrowthFactor));
    }
    else
    {
        m_BP5Serializer.InitStep(new ChunkV("DaosWriter", false, m_BP5Serializer.m_BufferAlign,
                                            m_BP5Serializer.m_BufferBlockSize,
                                            m_Parameters.BufferChunkSize));
    }
    m_ThisTimestepDataSize = 0;

    ts = Now() - m_EngineStart;
    // std::cout << "BEGIN STEP ended at: " << ts.count() << std::endl;
    return StepStatus::OK;
}

size_t DaosWriter::CurrentStep() const { return m_WriterStep; }

void DaosWriter::PerformPuts()
{
    PERFSTUBS_SCOPED_TIMER("DaosWriter::PerformPuts");
    m_Profiler.Start("PP");
    m_BP5Serializer.PerformPuts(m_Parameters.AsyncWrite || m_Parameters.DirectIO);
    m_Profiler.Stop("PP");
    return;
}

void DaosWriter::WriteMetaMetadata(
    const std::vector<format::BP5Base::MetaMetaInfoBlock> MetaMetaBlocks)
{
    for (auto &b : MetaMetaBlocks)
    {
        m_FileMetaMetadataManager.WriteFiles((char *)&b.MetaMetaIDLen, sizeof(size_t));
        m_FileMetaMetadataManager.WriteFiles((char *)&b.MetaMetaInfoLen, sizeof(size_t));
        m_FileMetaMetadataManager.WriteFiles((char *)b.MetaMetaID, b.MetaMetaIDLen);
        m_FileMetaMetadataManager.WriteFiles((char *)b.MetaMetaInfo, b.MetaMetaInfoLen);
    }
}

uint64_t DaosWriter::WriteMetadata(const std::vector<core::iovec> &MetaDataBlocks,
                                   const std::vector<core::iovec> &AttributeBlocks)
{
    uint64_t MDataTotalSize = 0;
    uint64_t MetaDataSize = 0;
    std::vector<uint64_t> SizeVector;
    std::vector<uint64_t> AttrSizeVector;
    SizeVector.reserve(MetaDataBlocks.size());
    for (auto &b : MetaDataBlocks)
    {
        MDataTotalSize += sizeof(uint64_t) + b.iov_len;
        SizeVector.push_back(b.iov_len);
    }
    for (auto &b : AttributeBlocks)
    {
        MDataTotalSize += sizeof(uint64_t) + b.iov_len;
        AttrSizeVector.push_back(b.iov_len);
    }
    MetaDataSize = 0;
    m_FileMetadataManager.WriteFiles((char *)&MDataTotalSize, sizeof(uint64_t));
    MetaDataSize += sizeof(uint64_t);
    m_FileMetadataManager.WriteFiles((char *)SizeVector.data(),
                                     sizeof(uint64_t) * SizeVector.size());
    MetaDataSize += sizeof(uint64_t) * AttrSizeVector.size();
    m_FileMetadataManager.WriteFiles((char *)AttrSizeVector.data(),
                                     sizeof(uint64_t) * AttrSizeVector.size());
    MetaDataSize += sizeof(uint64_t) * AttrSizeVector.size();
    for (auto &b : MetaDataBlocks)
    {
        if (!b.iov_base)
            continue;
        m_FileMetadataManager.WriteFiles((char *)b.iov_base, b.iov_len);
        MetaDataSize += b.iov_len;
    }

    for (auto &b : AttributeBlocks)
    {
        if (!b.iov_base)
            continue;
        m_FileMetadataManager.WriteFiles((char *)b.iov_base, b.iov_len);
        MetaDataSize += b.iov_len;
    }

    m_MetaDataPos += MetaDataSize;
    return MetaDataSize;
}

void DaosWriter::AsyncWriteDataCleanup()
{
    if (m_Parameters.AsyncWrite)
    {
        switch (m_Parameters.AggregationType)
        {
        case (int)AggregationType::EveryoneWrites:
        case (int)AggregationType::EveryoneWritesSerial:
            AsyncWriteDataCleanup_EveryoneWrites();
            break;
        case (int)AggregationType::TwoLevelShm:
            AsyncWriteDataCleanup_TwoLevelShm();
            break;
        default:
            break;
        }
    }
}

void DaosWriter::WriteData(format::BufferV *Data)
{
    if (m_Parameters.AsyncWrite)
    {
        switch (m_Parameters.AggregationType)
        {
        case (int)AggregationType::EveryoneWrites:
            WriteData_EveryoneWrites_Async(Data, false);
            break;
        case (int)AggregationType::EveryoneWritesSerial:
            WriteData_EveryoneWrites_Async(Data, true);
            break;
        case (int)AggregationType::TwoLevelShm:
            WriteData_TwoLevelShm_Async(Data);
            break;
        default:
            helper::Throw<std::invalid_argument>("Engine", "DaosWriter", "WriteData",
                                                 "Aggregation method " +
                                                     std::to_string(m_Parameters.AggregationType) +
                                                     "is not supported in BP5");
        }
    }
    else
    {
        switch (m_Parameters.AggregationType)
        {
        case (int)AggregationType::EveryoneWrites:
            WriteData_EveryoneWrites(Data, false);
            break;
        case (int)AggregationType::EveryoneWritesSerial:
            WriteData_EveryoneWrites(Data, true);
            break;
        case (int)AggregationType::TwoLevelShm:
            WriteData_TwoLevelShm(Data);
            break;
        default:
            helper::Throw<std::invalid_argument>("Engine", "DaosWriter", "WriteData",
                                                 "Aggregation method " +
                                                     std::to_string(m_Parameters.AggregationType) +
                                                     "is not supported in BP5");
        }
        delete Data;
    }
}

void DaosWriter::WriteData_EveryoneWrites(format::BufferV *Data, bool SerializedWriters)
{
    const aggregator::MPIChain *a = dynamic_cast<aggregator::MPIChain *>(m_Aggregator);

    // new step writing starts at offset m_DataPos on aggregator
    // others will wait for the position to arrive from the rank below

    if (a->m_Comm.Rank() > 0)
    {
        a->m_Comm.Recv(&m_DataPos, 1, a->m_Comm.Rank() - 1, 0,
                       "Chain token in DaosWriter::WriteData");
    }

    // align to PAGE_SIZE
    m_DataPos += helper::PaddingToAlignOffset(m_DataPos, m_Parameters.StripeSize);
    m_StartDataPos = m_DataPos;

    if (!SerializedWriters && a->m_Comm.Rank() < a->m_Comm.Size() - 1)
    {
        /* Send the token before writing so everyone can start writing asap */
        uint64_t nextWriterPos = m_DataPos + Data->Size();
        a->m_Comm.Isend(&nextWriterPos, 1, a->m_Comm.Rank() + 1, 0,
                        "Chain token in DaosWriter::WriteData");
    }

    m_DataPos += Data->Size();
    std::vector<core::iovec> DataVec = Data->DataVec();
    m_FileDataManager.WriteFileAt(DataVec.data(), DataVec.size(), m_StartDataPos);

    if (SerializedWriters && a->m_Comm.Rank() < a->m_Comm.Size() - 1)
    {
        /* send token now, effectively serializing the writers in the chain */
        uint64_t nextWriterPos = m_DataPos;
        a->m_Comm.Isend(&nextWriterPos, 1, a->m_Comm.Rank() + 1, 0,
                        "Chain token in DaosWriter::WriteData");
    }

    if (a->m_Comm.Size() > 1)
    {
        // at the end, last rank sends back the final data pos to first rank
        // so it can update its data pos
        if (a->m_Comm.Rank() == a->m_Comm.Size() - 1)
        {
            a->m_Comm.Isend(&m_DataPos, 1, 0, 0, "Final chain token in DaosWriter::WriteData");
        }
        if (a->m_Comm.Rank() == 0)
        {
            a->m_Comm.Recv(&m_DataPos, 1, a->m_Comm.Size() - 1, 0,
                           "Chain token in DaosWriter::WriteData");
        }
    }
}

void DaosWriter::WriteMetadataFileIndex(uint64_t MetaDataPos, uint64_t MetaDataSize)
{
    m_FileMetadataManager.FlushFiles();

    // bufsize: Step record
    size_t bufsize =
        1 + (4 + ((FlushPosSizeInfo.size() * 2) + 1) * m_Comm.Size()) * sizeof(uint64_t);
    if (MetaDataPos == 0)
    {
        //  First time, write the headers
        bufsize += m_IndexHeaderSize;
    }
    if (!m_WriterSubfileMap.empty())
    {
        // WriterMap record
        bufsize += 1 + (4 + m_Comm.Size()) * sizeof(uint64_t);
    }

    std::vector<char> buf(bufsize);
    size_t pos = 0;
    uint64_t d;
    unsigned char record;

    if (MetaDataPos == 0)
    {
        //  First time, write the headers
        MakeHeader(buf, pos, "Index Table", true);
    }

    // WriterMap record
    if (!m_WriterSubfileMap.empty())
    {
        record = WriterMapRecord;
        helper::CopyToBuffer(buf, pos, &record, 1); // record type
        d = (3 + m_Comm.Size()) * sizeof(uint64_t);
        helper::CopyToBuffer(buf, pos, &d, 1); // record length
        d = static_cast<uint64_t>(m_Comm.Size());
        helper::CopyToBuffer(buf, pos, &d, 1);
        d = static_cast<uint64_t>(m_Aggregator->m_NumAggregators);
        helper::CopyToBuffer(buf, pos, &d, 1);
        d = static_cast<uint64_t>(m_Aggregator->m_SubStreams);
        helper::CopyToBuffer(buf, pos, &d, 1);
        helper::CopyToBuffer(buf, pos, m_WriterSubfileMap.data(), m_Comm.Size());
        m_WriterSubfileMap.clear();
    }

    // Step record
    record = StepRecord;
    helper::CopyToBuffer(buf, pos, &record, 1); // record type
    d = (3 + ((FlushPosSizeInfo.size() * 2) + 1) * m_Comm.Size()) * sizeof(uint64_t);
    helper::CopyToBuffer(buf, pos, &d, 1); // record length
    helper::CopyToBuffer(buf, pos, &MetaDataPos, 1);
    helper::CopyToBuffer(buf, pos, &MetaDataSize, 1);
    d = static_cast<uint64_t>(FlushPosSizeInfo.size());
    helper::CopyToBuffer(buf, pos, &d, 1);

    for (int writer = 0; writer < m_Comm.Size(); writer++)
    {
        for (size_t flushNum = 0; flushNum < FlushPosSizeInfo.size(); flushNum++)
        {
            // add two numbers here
            helper::CopyToBuffer(buf, pos, &FlushPosSizeInfo[flushNum][2 * writer], 2);
        }
        helper::CopyToBuffer(buf, pos, &m_WriterDataPos[writer], 1);
    }

    m_FileMetadataIndexManager.WriteFiles((char *)buf.data(), buf.size());

#ifdef DUMPDATALOCINFO
    std::cout << "Flush count is :" << FlushPosSizeInfo.size() << std::endl;
    std::cout << "Write Index positions = {" << std::endl;

    for (size_t i = 0; i < m_Comm.Size(); ++i)
    {
        std::cout << "Writer " << i << " has data at: " << std::endl;
        uint64_t eachWriterSize = FlushPosSizeInfo.size() * 2 + 1;
        for (size_t j = 0; j < FlushPosSizeInfo.size(); ++j)
        {
            std::cout << "loc:" << buf[3 + eachWriterSize * i + j * 2]
                      << " siz:" << buf[3 + eachWriterSize * i + j * 2 + 1] << std::endl;
        }
        std::cout << "loc:" << buf[3 + eachWriterSize * (i + 1) - 1] << std::endl;
    }
    std::cout << "}" << std::endl;
#endif
    /* reset for next timestep */
    FlushPosSizeInfo.clear();
}

void DaosWriter::NotifyEngineAttribute(std::string name, DataType type) noexcept
{
    helper::Throw<std::invalid_argument>("DaosWriter", "Engine", "ThrowUp",
                                         "Engine does not support NotifyEngineAttribute");
}

void DaosWriter::NotifyEngineAttribute(std::string name, AttributeBase *Attr, void *data) noexcept
{
    if (!m_Parameters.UseOneTimeAttributes)
    {
        m_MarshalAttributesNecessary = true;
        return;
    }

    m_BP5Serializer.OnetimeMarshalAttribute(*Attr);
    m_MarshalAttributesNecessary = false;
}

void DaosWriter::MarshalAttributes()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    const auto &attributes = m_IO.GetAttributes();

    // if there are no new attributes, nothing to do
    if (!m_MarshalAttributesNecessary)
    {
        return;
    }
    m_MarshalAttributesNecessary = false;

    for (const auto &attributePair : attributes)
    {
        const std::string name(attributePair.first);
        auto baseAttr = &attributePair.second;
        const DataType type((*baseAttr)->m_Type);
        int element_count = -1;

        if (!attributePair.second->m_IsSingleValue)
        {
            element_count = (*baseAttr)->m_Elements;
        }

        if (type == DataType::None)
        {
        }
        else if (type == helper::GetDataType<std::string>())
        {
            core::Attribute<std::string> &attribute = *m_IO.InquireAttribute<std::string>(name);
            void *data_addr;
            if (attribute.m_IsSingleValue)
            {
                data_addr = (void *)attribute.m_DataSingleValue.c_str();
            }
            else
            {
                const char **tmp = (const char **)malloc(sizeof(char *) * element_count);
                for (int i = 0; i < element_count; i++)
                {
                    auto str = &attribute.m_DataArray[i];
                    tmp[i] = str->c_str();
                }
                // tmp will be free'd after final attribute marshalling
                data_addr = (void *)tmp;
            }

            m_BP5Serializer.MarshalAttribute(name.c_str(), type, sizeof(char *), element_count,
                                             data_addr);
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        core::Attribute<T> &attribute = *m_IO.InquireAttribute<T>(name);                           \
        int element_count = -1;                                                                    \
        void *data_addr = &attribute.m_DataSingleValue;                                            \
        if (!attribute.m_IsSingleValue)                                                            \
        {                                                                                          \
            element_count = attribute.m_Elements;                                                  \
            data_addr = attribute.m_DataArray.data();                                              \
        }                                                                                          \
        m_BP5Serializer.MarshalAttribute(attribute.m_Name.c_str(), type, sizeof(T), element_count, \
                                         data_addr);                                               \
    }

        ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type
    }
}

void DaosWriter::EndStep()
{
    /* Seconds ts = Now() - m_EngineStart;
      std::cout << "END STEP starts at: " << ts.count() << std::endl; */
    m_BetweenStepPairs = false;
    PERFSTUBS_SCOPED_TIMER("DaosWriter::EndStep");
    m_Profiler.Start("endstep");

    m_Profiler.Start("close_ts");
    MarshalAttributes();

    // true: advances step
    auto TSInfo = m_BP5Serializer.CloseTimestep(m_WriterStep,
                                                m_Parameters.AsyncWrite || m_Parameters.DirectIO);

    /* TSInfo includes NewMetaMetaBlocks, the MetaEncodeBuffer, the
     * AttributeEncodeBuffer and the data encode Vector */

    m_ThisTimestepDataSize += TSInfo.DataBuffer->Size();
    m_Profiler.Stop("close_ts");

    m_Profiler.Start("AWD");
    // TSInfo destructor would delete the DataBuffer so we need to save it
    // for async IO and let the writer free it up when not needed anymore
    adios2::format::BufferV *databuf = TSInfo.DataBuffer;
    TSInfo.DataBuffer = NULL;
    m_AsyncWriteLock.lock();
    m_flagRush = false;
    m_AsyncWriteLock.unlock();
    WriteData(databuf);
    m_Profiler.Stop("AWD");

    /*
     * Two-step metadata aggregation
     */
    m_Profiler.Start("meta_lvl1");
    std::vector<char> MetaBuffer;
    // core::iovec m{TSInfo.MetaEncodeBuffer->Data(),
    // TSInfo.MetaEncodeBuffer->m_FixedSize};
    core::iovec m{nullptr, 0};
    core::iovec a{nullptr, 0};
    if (TSInfo.AttributeEncodeBuffer)
    {
        a = {TSInfo.AttributeEncodeBuffer->Data(), TSInfo.AttributeEncodeBuffer->m_FixedSize};
    }
    MetaBuffer = m_BP5Serializer.CopyMetadataToContiguous(
        TSInfo.NewMetaMetaBlocks, {m}, {a}, {m_ThisTimestepDataSize}, {m_StartDataPos});

    if (m_Aggregator->m_Comm.Size() > 1)
    { // level 1
        m_Profiler.Start("meta_gather1");
        size_t LocalSize = MetaBuffer.size();
        std::vector<size_t> RecvCounts = m_Aggregator->m_Comm.GatherValues(LocalSize, 0);
        std::vector<char> RecvBuffer;
        if (m_Aggregator->m_Comm.Rank() == 0)
        {
            uint64_t TotalSize = 0;
            for (auto &n : RecvCounts)
                TotalSize += n;
            RecvBuffer.resize(TotalSize);
            /*std::cout << "MD Lvl-1: rank " << m_Comm.Rank() << " gather "
                      << TotalSize << " bytes from aggregator group"
                      << std::endl;*/
        }
        m_Aggregator->m_Comm.GathervArrays(MetaBuffer.data(), LocalSize, RecvCounts.data(),
                                           RecvCounts.size(), RecvBuffer.data(), 0);
        m_Profiler.Stop("meta_gather1");
        if (m_Aggregator->m_Comm.Rank() == 0)
        {
            std::vector<format::BP5Base::MetaMetaInfoBlock> UniqueMetaMetaBlocks;
            std::vector<uint64_t> DataSizes;
            std::vector<uint64_t> WriterDataPositions;
            std::vector<core::iovec> AttributeBlocks;
            auto Metadata = m_BP5Serializer.BreakoutContiguousMetadata(
                RecvBuffer, RecvCounts, UniqueMetaMetaBlocks, AttributeBlocks, DataSizes,
                WriterDataPositions);

            MetaBuffer.clear();
            MetaBuffer = m_BP5Serializer.CopyMetadataToContiguous(
                UniqueMetaMetaBlocks, Metadata, AttributeBlocks, DataSizes, WriterDataPositions);
        }
    } // level 1
    m_Profiler.Stop("meta_lvl1");
    m_Profiler.Start("meta_lvl2");
    // level 2
    if (m_Aggregator->m_Comm.Rank() == 0)
    {
        std::vector<char> RecvBuffer;
        std::vector<char> *buf;
        std::vector<size_t> RecvCounts;
        size_t LocalSize = MetaBuffer.size();
        if (m_CommAggregators.Size() > 1)
        {
            m_Profiler.Start("meta_gather2");
            RecvCounts = m_CommAggregators.GatherValues(LocalSize, 0);
            if (m_CommAggregators.Rank() == 0)
            {
                uint64_t TotalSize = 0;
                for (auto &n : RecvCounts)
                    TotalSize += n;
                RecvBuffer.resize(TotalSize);
                /*std::cout << "MD Lvl-2: rank " << m_Comm.Rank() << " gather "
                          << TotalSize << " bytes from aggregator group"
                          << std::endl;*/
            }

            m_CommAggregators.GathervArrays(MetaBuffer.data(), LocalSize, RecvCounts.data(),
                                            RecvCounts.size(), RecvBuffer.data(), 0);
            buf = &RecvBuffer;
            m_Profiler.Stop("meta_gather2");
        }
        else
        {
            buf = &MetaBuffer;
            RecvCounts.push_back(LocalSize);
        }

        if (m_CommAggregators.Rank() == 0)
        {
            std::vector<format::BP5Base::MetaMetaInfoBlock> UniqueMetaMetaBlocks;
            std::vector<uint64_t> DataSizes;
            std::vector<core::iovec> AttributeBlocks;
            m_WriterDataPos.resize(0);
            auto Metadata = m_BP5Serializer.BreakoutContiguousMetadata(
                *buf, RecvCounts, UniqueMetaMetaBlocks, AttributeBlocks, DataSizes,
                m_WriterDataPos);
            assert(m_WriterDataPos.size() == static_cast<size_t>(m_Comm.Size()));
            WriteMetaMetadata(UniqueMetaMetaBlocks);
            m_LatestMetaDataPos = m_MetaDataPos;
            m_LatestMetaDataSize = WriteMetadata(Metadata, AttributeBlocks);
            // m_LatestMetaDataPos = 0;
            // m_LatestMetaDataSize = 0;
            if (!m_Parameters.AsyncWrite)
            {
                WriteMetadataFileIndex(m_LatestMetaDataPos, m_LatestMetaDataSize);
            }
        }
    } // level 2
    m_Profiler.Stop("meta_lvl2");

    char key[1000];
    int rc;

    sprintf(key, "step%d-rank%d", m_WriterStep, m_Comm.Rank());
    std::cout << __FILE__ << "::" << __func__ << "(), step: " << m_WriterStep << std::endl;
    std::cout << "Rank = " << m_Comm.Rank()
              << ", Metadata size = " << TSInfo.MetaEncodeBuffer->m_FixedSize << std::endl;
    std::cout << "key = " << key << std::endl;
    std::cout << "Printing the first 10 bytes of Metadata" << std::endl;
    char *data = reinterpret_cast<char *>(TSInfo.MetaEncodeBuffer->Data());
    for (int i = 0; i < 10; i++)
    {
        // std::cout << std::hex << std::setw(2) << std::setfill('0') <<
        // static_cast<int>(data[i]) << " ";
        std::cout << static_cast<int>(data[i]) << " ";
    }
    std::cout << std::endl;
    rc = daos_kv_put(oh, DAOS_TX_NONE, 0, key, TSInfo.MetaEncodeBuffer->m_FixedSize,
                     TSInfo.MetaEncodeBuffer->Data(), NULL);
    ASSERT(rc == 0, "daos_kv_put() failed with %d", rc);

    if (m_Parameters.AsyncWrite)
    {
        /* Start counting computation blocks between EndStep and next BeginStep
         * each time */
        {
            m_AsyncWriteLock.lock();
            m_ComputationBlockTimes.clear();
            m_ComputationBlocksLength = 0.0;
            m_ComputationBlockID = 0;
            m_AsyncWriteLock.unlock();
        }
    }

    m_Profiler.Stop("endstep");
    m_WriterStep++;
    m_EndStepEnd = Now();
    /* Seconds ts2 = Now() - m_EngineStart;
     std::cout << "END STEP ended at: " << ts2.count() << std::endl;*/
}

// PRIVATE
void DaosWriter::Init()
{
    m_BP5Serializer.m_Engine = this;
    m_RankMPI = m_Comm.Rank();
    InitParameters();
    InitAggregator();
    InitTransports();
    InitDAOS();
    InitBPBuffer();
}

void DaosWriter::InitParameters()
{
    ParseParams(m_IO, m_Parameters);
    m_WriteToBB = !(m_Parameters.BurstBufferPath.empty());
    m_DrainBB = m_WriteToBB && m_Parameters.BurstBufferDrain;

    unsigned int nproc = (unsigned int)m_Comm.Size();
    m_Parameters.NumAggregators = helper::SetWithinLimit(m_Parameters.NumAggregators, 0U, nproc);
    m_Parameters.NumSubFiles = helper::SetWithinLimit(m_Parameters.NumSubFiles, 0U, nproc);
    m_Parameters.AggregatorRatio = helper::SetWithinLimit(m_Parameters.AggregatorRatio, 0U, nproc);
    if (m_Parameters.NumAggregators == 0)
    {
        if (m_Parameters.AggregatorRatio > 0)
        {
            m_Parameters.NumAggregators =
                helper::SetWithinLimit(nproc / m_Parameters.AggregatorRatio, 0U, nproc);
        }
        else if (m_Parameters.NumSubFiles > 0)
        {
            m_Parameters.NumAggregators =
                helper::SetWithinLimit(m_Parameters.NumSubFiles, 0U, nproc);
        }
    }
    m_Parameters.NumSubFiles =
        helper::SetWithinLimit(m_Parameters.NumSubFiles, 0U, m_Parameters.NumAggregators);

    // Limiting to max 64MB page size
    m_Parameters.StripeSize = helper::SetWithinLimit(m_Parameters.StripeSize, 0U, 67108864U);
    if (m_Parameters.StripeSize == 0)
    {
        m_Parameters.StripeSize = 4096;
    }

    if (m_Parameters.DirectIO)
    {
        if (m_Parameters.DirectIOAlignBuffer == 0)
        {
            m_Parameters.DirectIOAlignBuffer = m_Parameters.DirectIOAlignOffset;
        }
        m_BP5Serializer.m_BufferBlockSize = m_Parameters.DirectIOAlignOffset;
        m_BP5Serializer.m_BufferAlign = m_Parameters.DirectIOAlignBuffer;
        if (m_Parameters.StripeSize % m_Parameters.DirectIOAlignOffset)
        {
            size_t k = m_Parameters.StripeSize / m_Parameters.DirectIOAlignOffset + 1;
            m_Parameters.StripeSize = k * m_Parameters.DirectIOAlignOffset;
        }
        if (m_Parameters.BufferChunkSize % m_Parameters.DirectIOAlignOffset)
        {
            size_t k = m_Parameters.BufferChunkSize / m_Parameters.DirectIOAlignOffset + 1;
            m_Parameters.BufferChunkSize = k * m_Parameters.DirectIOAlignOffset;
        }
    }

    m_BP5Serializer.m_StatsLevel = m_Parameters.StatsLevel;
}

uint64_t DaosWriter::CountStepsInMetadataIndex(format::BufferSTL &bufferSTL)
{
    const auto &buffer = bufferSTL.m_Buffer;
    size_t &position = bufferSTL.m_Position;

    if (buffer.size() < m_IndexHeaderSize)
    {
        m_AppendMetadataPos = 0;
        m_AppendMetaMetadataPos = 0;
        m_AppendMetadataIndexPos = 0;
        m_AppendDataPos.resize(m_Aggregator->m_NumAggregators,
                               0ULL); // safe bet
        return 0;
    }

    // Check endinanness
    position = m_EndianFlagPosition;
    const uint8_t endianness = helper::ReadValue<uint8_t>(buffer, position);
    bool IsLittleEndian = (endianness == 0) ? true : false;
    if (helper::IsLittleEndian() != IsLittleEndian)
    {
        std::string m = (IsLittleEndian ? "Little" : "Big");

        helper::Throw<std::runtime_error>("Engine", "DaosWriter", "CountStepsInMetadataIndex",
                                          "ADIOS2 BP5 Engine only supports appending with the same "
                                          "endianness. The existing file is " +
                                              m + "Endian");
    }

    // BP version
    position = m_BPVersionPosition;
    uint8_t Version = helper::ReadValue<uint8_t>(buffer, position, IsLittleEndian);
    if (Version != 5)
    {
        helper::Throw<std::runtime_error>("Engine", "DaosWriter", "CountStepsInMetadataIndex",
                                          "ADIOS2 BP5 Engine only supports bp format "
                                          "version 5, found " +
                                              std::to_string(Version) + " version");
    }

    // BP minor version
    position = m_BPMinorVersionPosition;
    uint8_t minorVersion = helper::ReadValue<uint8_t>(buffer, position, IsLittleEndian);
    if (minorVersion != m_BP5MinorVersion)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "DaosWriter", "CountStepsInMetadataIndex",
            "Current ADIOS2 BP5 Engine can only append to bp format 5." +
                std::to_string(m_BP5MinorVersion) + " but this file is 5." +
                std::to_string(minorVersion) + " version");
    }

    position = m_ColumnMajorFlagPosition;
    const uint8_t columnMajor = helper::ReadValue<uint8_t>(buffer, position, IsLittleEndian);
    const uint8_t NowColumnMajor = (m_IO.m_ArrayOrder == ArrayOrdering::ColumnMajor) ? 'y' : 'n';
    if (columnMajor != NowColumnMajor)
    {
        std::string m = (columnMajor == 'y' ? "column" : "row");
        helper::Throw<std::runtime_error>("Engine", "DaosWriter", "CountStepsInMetadataIndex",
                                          "ADIOS2 BP5 Engine only supports appending with the same "
                                          "column/row major settings as it was written."
                                          " Existing file is " +
                                              m + " major");
    }

    position = m_IndexHeaderSize; // after the header
    // Just count the steps first
    unsigned int availableSteps = 0;
    uint64_t nDataFiles = 0;
    while (position < buffer.size())
    {
        const unsigned char recordID =
            helper::ReadValue<unsigned char>(buffer, position, IsLittleEndian);
        position += sizeof(uint64_t); // recordLength

        switch (recordID)
        {
        case IndexRecord::WriterMapRecord: {
            m_AppendWriterCount = helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
            m_AppendAggregatorCount = helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
            m_AppendSubfileCount = helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
            if (m_AppendSubfileCount > nDataFiles)
            {
                nDataFiles = m_AppendSubfileCount;
            }
            // jump over writermap
            position += m_AppendWriterCount * sizeof(uint64_t);
            break;
        }
        case IndexRecord::StepRecord: {
            position += 2 * sizeof(uint64_t); // MetadataPos, MetadataSize
            const uint64_t FlushCount =
                helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
            // jump over the metadata positions
            position += sizeof(uint64_t) * m_AppendWriterCount * ((2 * FlushCount) + 1);
            availableSteps++;
            break;
        }
        }
    }

    unsigned int targetStep = 0;

    if (m_Parameters.AppendAfterSteps < 0)
    {
        // -1 means append after last step
        int s = (int)availableSteps + m_Parameters.AppendAfterSteps + 1;
        if (s < 0)
        {
            s = 0;
        }
        targetStep = static_cast<unsigned int>(s);
    }
    else
    {
        targetStep = static_cast<unsigned int>(m_Parameters.AppendAfterSteps);
    }
    if (targetStep > availableSteps)
    {
        targetStep = availableSteps;
    }

    m_AppendDataPos.resize(nDataFiles, 0ULL);

    if (!targetStep)
    {
        // append at 0 is like writing new file
        m_AppendMetadataPos = 0;
        m_AppendMetaMetadataPos = 0;
        m_AppendMetadataIndexPos = 0;
        return 0;
    }

    m_AppendMetadataPos = MaxSizeT; // size of header
    m_AppendMetaMetadataPos = MaxSizeT;
    m_AppendMetadataIndexPos = MaxSizeT;
    std::fill(m_AppendDataPos.begin(), m_AppendDataPos.end(), MaxSizeT);

    if (targetStep == availableSteps)
    {
        // append after existing steps
        return targetStep;
    }

    // append but not at 0 and not after existing steps
    // Read each record now completely to get offsets at step+1
    position = m_IndexHeaderSize;
    unsigned int currentStep = 0;
    std::vector<uint64_t> writerToFileMap;
    // reading one step beyond target to get correct offsets
    while (currentStep <= targetStep && position < buffer.size())
    {
        const unsigned char recordID =
            helper::ReadValue<unsigned char>(buffer, position, IsLittleEndian);
        position += sizeof(uint64_t); // recordLength

        switch (recordID)
        {
        case IndexRecord::WriterMapRecord: {
            m_AppendWriterCount = helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
            m_AppendAggregatorCount = helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
            m_AppendSubfileCount = helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);

            // Get the process -> subfile map
            writerToFileMap.clear();
            for (uint64_t i = 0; i < m_AppendWriterCount; i++)
            {
                const uint64_t subfileIdx =
                    helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
                writerToFileMap.push_back(subfileIdx);
            }
            break;
        }
        case IndexRecord::StepRecord: {
            m_AppendMetadataIndexPos =
                position - sizeof(unsigned char) - sizeof(uint64_t); // pos of RecordID
            const uint64_t MetadataPos =
                helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);
            position += sizeof(uint64_t); // MetadataSize
            const uint64_t FlushCount =
                helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian);

            m_AppendMetadataPos = static_cast<size_t>(MetadataPos);

            if (currentStep == targetStep)
            {
                // we need the very first (smallest) write position to each
                // subfile Offsets and sizes,  2*FlushCount + 1 per writer
                for (uint64_t i = 0; i < m_AppendWriterCount; i++)
                {
                    // first flush/write position will do
                    const size_t FirstDataPos = static_cast<size_t>(
                        helper::ReadValue<uint64_t>(buffer, position, IsLittleEndian));
                    position += sizeof(uint64_t) * 2 * FlushCount; // no need to read
                    /* std::cout << "Writer " << i << " subfile " <<
                       writerToFileMap[i]  << "  first data loc:" <<
                       FirstDataPos << std::endl; */
                    if (FirstDataPos < m_AppendDataPos[writerToFileMap[i]])
                    {
                        m_AppendDataPos[writerToFileMap[i]] = FirstDataPos;
                    }
                }
            }
            else
            {
                // jump over all data offsets in this step
                position += sizeof(uint64_t) * m_AppendWriterCount * (1 + 2 * FlushCount);
            }
            currentStep++;
            break;
        }
        }
    }
    return targetStep;
}

void DaosWriter::InitAggregator()
{
    // in BP5, aggregation is "always on", but processes may be alone, so
    // m_Aggregator.m_IsActive is always true
    // m_Aggregator.m_Comm.Rank() will always succeed (not abort)
    // m_Aggregator.m_SubFileIndex is always set

    if (m_Parameters.AggregationType == (int)AggregationType::EveryoneWrites ||
        m_Parameters.AggregationType == (int)AggregationType::EveryoneWritesSerial)
    {
        m_Parameters.NumSubFiles = m_Parameters.NumAggregators;
        m_AggregatorEveroneWrites.Init(m_Parameters.NumAggregators, m_Parameters.NumSubFiles,
                                       m_Comm);
        m_IAmDraining = m_AggregatorEveroneWrites.m_IsAggregator;
        m_IAmWritingData = true;
        DataWritingComm = &m_AggregatorEveroneWrites.m_Comm;
        m_Aggregator = static_cast<aggregator::MPIAggregator *>(&m_AggregatorEveroneWrites);
    }
    else
    {
        size_t numNodes = m_AggregatorTwoLevelShm.PreInit(m_Comm);
        (void)numNodes;
        m_AggregatorTwoLevelShm.Init(m_Parameters.NumAggregators, m_Parameters.NumSubFiles, m_Comm);

        /*std::cout << "Rank " << m_RankMPI << " aggr? "
                  << m_AggregatorTwoLevelShm.m_IsAggregator << " master? "
                  << m_AggregatorTwoLevelShm.m_IsMasterAggregator
                  << " aggr size = " << m_AggregatorTwoLevelShm.m_Size
                  << " rank = " << m_AggregatorTwoLevelShm.m_Rank
                  << " subfile = " << m_AggregatorTwoLevelShm.m_SubStreamIndex
                  << " type = " << m_Parameters.AggregationType << std::endl;*/

        m_IAmDraining = m_AggregatorTwoLevelShm.m_IsMasterAggregator;
        m_IAmWritingData = m_AggregatorTwoLevelShm.m_IsAggregator;
        DataWritingComm = &m_AggregatorTwoLevelShm.m_AggregatorChainComm;
        m_Aggregator = static_cast<aggregator::MPIAggregator *>(&m_AggregatorTwoLevelShm);
    }

    /* comm for Aggregators only.
     *  We are only interested in the chain of rank 0s
     */
    int color = m_Aggregator->m_Comm.Rank();
    m_CommAggregators = m_Comm.Split(color, 0, "creating level 2 chain of aggregators at Open");
}

void DaosWriter::InitTransports()
{
    if (m_IO.m_TransportsParameters.empty())
    {
        Params defaultTransportParameters;
        defaultTransportParameters["transport"] = "File";
        m_IO.m_TransportsParameters.push_back(defaultTransportParameters);
    }

    if (m_WriteToBB)
    {
        m_BBName = m_Parameters.BurstBufferPath + PathSeparator + m_Name;
    }
    else
    {
        m_BBName = m_Name;
    }
    /* From this point, engine writes to m_BBName, which points to either
        the BB file system if BB is turned on, or to the target file system.
        m_Name always points to the target file system, to which the drainer
        should write if BB is turned on
    */

    // Names passed to IO AddTransport option with key "Name"
    const std::vector<std::string> transportsNames =
        m_FileDataManager.GetFilesBaseNames(m_BBName, m_IO.m_TransportsParameters);

    // /path/name.bp.dir/name.bp.rank
    m_SubStreamNames = GetBPSubStreamNames(transportsNames, m_Aggregator->m_SubStreamIndex);

    if (m_IAmDraining)
    {
        // Only (master)aggregators will run draining processes
        if (m_DrainBB)
        {
            const std::vector<std::string> drainTransportNames =
                m_FileDataManager.GetFilesBaseNames(m_Name, m_IO.m_TransportsParameters);
            m_DrainSubStreamNames =
                GetBPSubStreamNames(drainTransportNames, m_Aggregator->m_SubStreamIndex);
            /* start up BB thread */
            //            m_FileDrainer.SetVerbose(
            //				     m_Parameters.BurstBufferVerbose,
            //				     m_Comm.Rank());
            m_FileDrainer.Start();
        }
    }

    /* Create the directories either on target or burst buffer if used */
    //    m_BP4Serializer.m_Profiler.Start("mkdir");

    if (m_Comm.Rank() == 0)
    {
        m_MetadataFileNames = GetBPMetadataFileNames(transportsNames);
        m_MetaMetadataFileNames = GetBPMetaMetadataFileNames(transportsNames);
        m_MetadataIndexFileNames = GetBPMetadataIndexFileNames(transportsNames);
    }
    m_FileMetadataManager.MkDirsBarrier(m_MetadataFileNames, m_IO.m_TransportsParameters,
                                        m_Parameters.NodeLocal || m_WriteToBB);
    /* Create the directories on burst buffer if used */
    if (m_DrainBB)
    {
        /* Create the directories on target anyway by main thread */
        m_FileDataManager.MkDirsBarrier(m_DrainSubStreamNames, m_IO.m_TransportsParameters,
                                        m_Parameters.NodeLocal);
    }

    /* Everyone opens its data file. Each aggregation chain opens
       one data file and does so in chain, not everyone at once */
    if (m_Parameters.AsyncOpen)
    {
        for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
        {
            m_IO.m_TransportsParameters[i]["asyncopen"] = "true";
        }
    }

    if (m_Parameters.DirectIO)
    {
        for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
        {
            m_IO.m_TransportsParameters[i]["DirectIO"] = "true";
        }
    }

    bool useProfiler = true;

    if (m_IAmWritingData)
    {
        m_FileDataManager.OpenFiles(m_SubStreamNames, m_OpenMode, m_IO.m_TransportsParameters,
                                    useProfiler, *DataWritingComm);
    }

    if (m_IAmDraining)
    {
        if (m_DrainBB)
        {
            for (const auto &name : m_DrainSubStreamNames)
            {
                m_FileDrainer.AddOperationOpen(name, m_OpenMode);
            }
        }
    }

    if (m_Comm.Rank() == 0)
    {
        // force turn off directio to metadata files
        for (size_t i = 0; i < m_IO.m_TransportsParameters.size(); ++i)
        {
            m_IO.m_TransportsParameters[i]["DirectIO"] = "false";
        }
        m_FileMetaMetadataManager.OpenFiles(m_MetaMetadataFileNames, m_OpenMode,
                                            m_IO.m_TransportsParameters, useProfiler);

        m_FileMetadataManager.OpenFiles(m_MetadataFileNames, m_OpenMode,
                                        m_IO.m_TransportsParameters, useProfiler);

        m_FileMetadataIndexManager.OpenFiles(m_MetadataIndexFileNames, m_OpenMode,
                                             m_IO.m_TransportsParameters, useProfiler);

        if (m_DrainBB)
        {
            const std::vector<std::string> drainTransportNames =
                m_FileDataManager.GetFilesBaseNames(m_Name, m_IO.m_TransportsParameters);
            m_DrainMetadataFileNames = GetBPMetadataFileNames(drainTransportNames);
            m_DrainMetadataIndexFileNames = GetBPMetadataIndexFileNames(drainTransportNames);

            for (const auto &name : m_DrainMetadataFileNames)
            {
                m_FileDrainer.AddOperationOpen(name, m_OpenMode);
            }
            for (const auto &name : m_DrainMetadataIndexFileNames)
            {
                m_FileDrainer.AddOperationOpen(name, m_OpenMode);
            }
        }
    }
}

void DaosWriter::InitDAOS()
{
    // Rank 0 - Connect to DAOS pool, and open container
    int rc;
    rc = gethostname(node, sizeof(node));
    ASSERT(rc == 0, "buffer for hostname too small");
    if (m_Comm.Rank() == 0)
    {
        /** connect to the just created DAOS pool */
        rc = daos_pool_connect(pool_label, DSS_PSETID,
                               // DAOS_PC_EX ,
                               DAOS_PC_RW /* read write access */, &poh /* returned pool handle */,
                               NULL /* returned pool info */, NULL /* event */);
        ASSERT(rc == 0, "pool connect failed with %d", rc);
    }

    /** share pool handle with peer tasks */
    daos_handle_share(&poh, DaosWriter::HANDLE_POOL);

    if (m_Comm.Rank() == 0)
    {
        /** open container */
        rc = daos_cont_open(poh, cont_label, DAOS_COO_RW, &coh, NULL, NULL);
        ASSERT(rc == 0, "container open failed with %d", rc);
    }

    /** share container handle with peer tasks */
    daos_handle_share(&coh, HANDLE_CO);

    if (m_Comm.Rank() == 0)
    {
        /** Open a DAOS KV object */
        rc = daos_obj_generate_oid(coh, &oid, DAOS_OT_KV_HASHED, OC_SX, 0, 0);
        ASSERT(rc == 0, "daos_obj_generate_oid failed with %d", rc);
    }

    // Rank 0 will broadcast the DAOS KV OID
    MPI_Bcast(&oid.hi, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&oid.lo, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    // Open KV object
    rc = daos_kv_open(coh, oid, 0, &oh, NULL);
    ASSERT(rc == 0, "daos_kv_open failed with %d", rc);
    FILE *fp = fopen("./share/oid.txt", "w");
    if (fp == NULL)
    {
        perror("fopen");
        exit(1);
    }
    fprintf(fp, "%" PRIu64 "\n%" PRIu64 "\n", oid.hi, oid.lo);
    fclose(fp);
}

/*generate the header for the metadata index file*/
void DaosWriter::MakeHeader(std::vector<char> &buffer, size_t &position, const std::string fileType,
                            const bool isActive)
{
    auto lf_CopyVersionChar = [](const std::string version, std::vector<char> &buffer,
                                 size_t &position) {
        helper::CopyToBuffer(buffer, position, version.c_str());
    };

    // auto &buffer = b.m_Buffer;
    // auto &position = b.m_Position;
    // auto &absolutePosition = b.m_AbsolutePosition;
    if (position > 0)
    {
        helper::Throw<std::invalid_argument>(
            "Engine", "DaosWriter", "MakeHeader",
            "BP4Serializer::MakeHeader can only be called for an empty "
            "buffer. This one for " +
                fileType + " already has content of " + std::to_string(position) + " bytes.");
    }

    if (buffer.size() < m_IndexHeaderSize)
    {
        buffer.resize(m_IndexHeaderSize);
    }

    const std::string majorVersion(std::to_string(ADIOS2_VERSION_MAJOR));
    const std::string minorVersion(std::to_string(ADIOS2_VERSION_MINOR));
    const std::string patchVersion(std::to_string(ADIOS2_VERSION_PATCH));

    // byte 0-31: Readable tag
    if (position != m_VersionTagPosition)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "DaosWriter", "MakeHeader",
            "ADIOS Coding ERROR in BP4Serializer::MakeHeader. Version Tag "
            "position mismatch");
    }
    std::string versionLongTag("ADIOS-BP v" + majorVersion + "." + minorVersion + "." +
                               patchVersion + " ");
    size_t maxTypeLen = m_VersionTagLength - versionLongTag.size();
    const std::string fileTypeStr = fileType.substr(0, maxTypeLen);
    versionLongTag += fileTypeStr;
    const size_t versionLongTagSize = versionLongTag.size();
    if (versionLongTagSize < m_VersionTagLength)
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), versionLongTagSize);
        position += m_VersionTagLength - versionLongTagSize;
    }
    else if (versionLongTagSize > m_VersionTagLength)
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), m_VersionTagLength);
    }
    else
    {
        helper::CopyToBuffer(buffer, position, versionLongTag.c_str(), m_VersionTagLength);
    }

    // byte 32-35: MAJOR MINOR PATCH Unused

    lf_CopyVersionChar(majorVersion, buffer, position);
    lf_CopyVersionChar(minorVersion, buffer, position);
    lf_CopyVersionChar(patchVersion, buffer, position);
    ++position;

    // Note: Reader does process and use bytes 36-38 in
    // BP4Deserialize.cpp::ParseMetadataIndex().
    // Order and position must match there.

    // byte 36: endianness
    if (position != m_EndianFlagPosition)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "DaosWriter", "MakeHeader",
            "ADIOS Coding ERROR in DaosWriter::MakeHeader. Endian Flag "
            "position mismatch");
    }
    const uint8_t endianness = helper::IsLittleEndian() ? 0 : 1;
    helper::CopyToBuffer(buffer, position, &endianness);

    // byte 37: BP Version 5
    if (position != m_BPVersionPosition)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "DaosWriter", "MakeHeader",
            "ADIOS Coding ERROR in DaosWriter::MakeHeader. BP Version "
            "position mismatch");
    }
    const uint8_t version = 5;
    helper::CopyToBuffer(buffer, position, &version);

    // byte 38: BP Minor version 1
    if (position != m_BPMinorVersionPosition)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "DaosWriter", "MakeHeader",
            "ADIOS Coding ERROR in DaosWriter::MakeHeader. BP Minor version "
            "position mismatch");
    }
    const uint8_t minorversion = m_BP5MinorVersion;
    helper::CopyToBuffer(buffer, position, &minorversion);

    // byte 39: Active flag (used in Index Table only)
    if (position != m_ActiveFlagPosition)
    {
        helper::Throw<std::runtime_error>(
            "Engine", "DaosWriter", "MakeHeader",
            "ADIOS Coding ERROR in DaosWriter::MakeHeader. Active Flag "
            "position mismatch");
    }
    const uint8_t activeFlag = (isActive ? 1 : 0);
    helper::CopyToBuffer(buffer, position, &activeFlag);

    // byte 40 columnMajor
    // write if data is column major in metadata and data
    const uint8_t columnMajor = (m_IO.m_ArrayOrder == ArrayOrdering::ColumnMajor) ? 'y' : 'n';
    helper::CopyToBuffer(buffer, position, &columnMajor);

    // byte 41-63: unused
    position += 23;
    // absolutePosition = position;
}

void DaosWriter::UpdateActiveFlag(const bool active)
{
    const char activeChar = (active ? '\1' : '\0');
    m_FileMetadataIndexManager.WriteFileAt(&activeChar, 1, m_ActiveFlagPosition);
    m_FileMetadataIndexManager.FlushFiles();
    m_FileMetadataIndexManager.SeekToFileEnd();
    if (m_DrainBB)
    {
        for (size_t i = 0; i < m_MetadataIndexFileNames.size(); ++i)
        {
            m_FileDrainer.AddOperationWriteAt(m_DrainMetadataIndexFileNames[i],
                                              m_ActiveFlagPosition, 1, &activeChar);
            m_FileDrainer.AddOperationSeekEnd(m_DrainMetadataIndexFileNames[i]);
        }
    }
}

void DaosWriter::InitBPBuffer()
{
    if (m_OpenMode == Mode::Append)
    {
        format::BufferSTL preMetadataIndex;
        size_t preMetadataIndexFileSize;

        if (m_Comm.Rank() == 0)
        {
            preMetadataIndexFileSize = m_FileMetadataIndexManager.GetFileSize(0);
            preMetadataIndex.m_Buffer.resize(preMetadataIndexFileSize);
            preMetadataIndex.m_Buffer.assign(preMetadataIndex.m_Buffer.size(), '\0');
            preMetadataIndex.m_Position = 0;
            m_FileMetadataIndexManager.ReadFile(preMetadataIndex.m_Buffer.data(),
                                                preMetadataIndexFileSize);
        }
        m_Comm.BroadcastVector(preMetadataIndex.m_Buffer);
        m_WriterStep = CountStepsInMetadataIndex(preMetadataIndex);

        // truncate and seek
        if (m_Aggregator->m_IsAggregator)
        {
            const size_t off = m_AppendDataPos[m_Aggregator->m_SubStreamIndex];
            if (off < MaxSizeT)
            {
                m_FileDataManager.Truncate(off);
                // Seek is needed since truncate does not seek.
                // SeekTo instead of SeetToFileEnd in case a transport
                // does not support actual truncate.
                m_FileDataManager.SeekTo(off);
                m_DataPos = off;
            }
            else
            {
                m_DataPos = m_FileDataManager.GetFileSize(0);
            }
        }

        if (m_Comm.Rank() == 0)
        {
            // Truncate existing metadata file
            if (m_AppendMetadataPos < MaxSizeT)
            {
                m_MetaDataPos = m_AppendMetadataPos;
                m_FileMetadataManager.Truncate(m_MetaDataPos);
                m_FileMetadataManager.SeekTo(m_MetaDataPos);
            }
            else
            {
                m_MetaDataPos = m_FileMetadataManager.GetFileSize(0);
                m_FileMetadataManager.SeekToFileEnd();
            }

            // Truncate existing meta-meta file
            if (m_AppendMetaMetadataPos < MaxSizeT)
            {
                m_FileMetaMetadataManager.Truncate(m_AppendMetaMetadataPos);
                m_FileMetaMetadataManager.SeekTo(m_AppendMetaMetadataPos);
            }
            else
            {
                m_FileMetadataIndexManager.SeekToFileEnd();
            }

            // Set the flag in the header of metadata index table to 1 again
            // to indicate a new run begins
            UpdateActiveFlag(true);

            // Truncate existing index file
            if (m_AppendMetadataIndexPos < MaxSizeT)
            {
                m_FileMetadataIndexManager.Truncate(m_AppendMetadataIndexPos);
                m_FileMetadataIndexManager.SeekTo(m_AppendMetadataIndexPos);
            }
            else
            {
                m_FileMetadataIndexManager.SeekToFileEnd();
            }
        }
        m_AppendDataPos.clear();
    }

    if (!m_WriterStep)
    {
        /* This is a new file or append at 0
         * Make headers in data buffer and metadata buffer (but do not write
         * them yet so that Open() can stay free of writing to disk)
         */
        if (m_Comm.Rank() == 0)
        {
            m_FileMetadataIndexManager.SeekToFileBegin();
            m_FileMetadataManager.SeekToFileBegin();
            m_FileMetaMetadataManager.SeekToFileBegin();
        }
        // last attempt to clean up datafile if called with append mode,
        // data existed but index was missing
        if (m_Aggregator->m_IsAggregator)
        {
            m_FileDataManager.SeekTo(0);
        }
    }

    if (m_Comm.Rank() == 0)
    {
        m_WriterDataPos.resize(m_Comm.Size());
    }

    if (!m_WriterStep || m_AppendWriterCount != static_cast<unsigned int>(m_Comm.Size()) ||
        m_AppendAggregatorCount != static_cast<unsigned int>(m_Aggregator->m_NumAggregators) ||
        m_AppendSubfileCount != static_cast<unsigned int>(m_Aggregator->m_SubStreams))
    {
        // new Writer Map is needed, generate now, write later
        const uint64_t a = static_cast<uint64_t>(m_Aggregator->m_SubStreamIndex);
        m_WriterSubfileMap = m_Comm.GatherValues(a, 0);
    }
}

void DaosWriter::EnterComputationBlock() noexcept
{
    if (m_Parameters.AsyncWrite && !m_BetweenStepPairs)
    {
        m_ComputationBlockStart = Now();
        {
            m_AsyncWriteLock.lock();
            m_InComputationBlock = true;
            m_AsyncWriteLock.unlock();
        }
    }
}

void DaosWriter::ExitComputationBlock() noexcept
{
    if (m_Parameters.AsyncWrite && m_InComputationBlock)
    {
        double t = Seconds(Now() - m_ComputationBlockStart).count();
        {
            m_AsyncWriteLock.lock();
            if (t > 0.1) // only register long enough intervals
            {
                m_ComputationBlockTimes.emplace_back(m_ComputationBlockID, t);
                m_ComputationBlocksLength += t;
            }
            m_InComputationBlock = false;
            ++m_ComputationBlockID;
            m_AsyncWriteLock.unlock();
        }
    }
}

void DaosWriter::FlushData(const bool isFinal)
{
    BufferV *DataBuf;
    if (m_Parameters.BufferVType == (int)BufferVType::MallocVType)
    {
        DataBuf = m_BP5Serializer.ReinitStepData(
            new MallocV("DaosWriter", false, m_BP5Serializer.m_BufferAlign,
                        m_BP5Serializer.m_BufferBlockSize, m_Parameters.InitialBufferSize,
                        m_Parameters.GrowthFactor),
            m_Parameters.AsyncWrite || m_Parameters.DirectIO);
    }
    else
    {
        DataBuf = m_BP5Serializer.ReinitStepData(
            new ChunkV("DaosWriter", false, m_BP5Serializer.m_BufferAlign,
                       m_BP5Serializer.m_BufferBlockSize, m_Parameters.BufferChunkSize),
            m_Parameters.AsyncWrite || m_Parameters.DirectIO);
    }

    auto databufsize = DataBuf->Size();
    WriteData(DataBuf);
    /* DataBuf is deleted in WriteData() */
    DataBuf = nullptr;

    m_ThisTimestepDataSize += databufsize;

    if (!isFinal)
    {
        size_t tmp[2];
        // aggregate start pos and data size to rank 0
        tmp[0] = m_StartDataPos;
        tmp[1] = databufsize;

        std::vector<size_t> RecvBuffer;
        if (m_Comm.Rank() == 0)
        {
            RecvBuffer.resize(m_Comm.Size() * 2);
        }
        m_Comm.GatherArrays(tmp, 2, RecvBuffer.data(), 0);
        if (m_Comm.Rank() == 0)
        {
            FlushPosSizeInfo.push_back(RecvBuffer);
        }
    }
}

void DaosWriter::Flush(const int transportIndex) {}

void DaosWriter::PerformDataWrite() { FlushData(false); }

void DaosWriter::DestructorClose(bool Verbose) noexcept
{
    if (Verbose)
    {
        std::cerr << "BP5 Writer \"" << m_Name << "\" Destroyed without a prior Close()."
                  << std::endl;
        std::cerr << "This may result in corrupt output." << std::endl;
    }
    // close metadata index file
    UpdateActiveFlag(false);
    m_IsOpen = false;
}

DaosWriter::~DaosWriter()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

void DaosWriter::DoClose(const int transportIndex)
{
    PERFSTUBS_SCOPED_TIMER("DaosWriter::Close");

    if ((m_WriterStep == 0) && !m_BetweenStepPairs)
    {
        /* never did begin step, do one now */
        BeginStep(StepMode::Update);
    }
    if (m_BetweenStepPairs)
    {
        EndStep();
    }

    TimePoint wait_start = Now();
    Seconds wait(0.0);
    if (m_WriteFuture.valid())
    {
        m_Profiler.Start("WaitOnAsync");
        m_AsyncWriteLock.lock();
        m_flagRush = true;
        m_AsyncWriteLock.unlock();
        m_WriteFuture.get();
        wait += Now() - wait_start;
        m_Profiler.Stop("WaitOnAsync");
    }

    m_FileDataManager.CloseFiles(transportIndex);
    // Delete files from temporary storage if draining was on

    if (m_Comm.Rank() == 0)
    {
        // close metadata file
        m_FileMetadataManager.CloseFiles();

        // close metametadata file
        m_FileMetaMetadataManager.CloseFiles();
    }

    if (m_Parameters.AsyncWrite)
    {
        // wait until all process' writing thread completes
        m_Profiler.Start("WaitOnAsync");
        wait_start = Now();
        m_Comm.Barrier();
        AsyncWriteDataCleanup();
        wait += Now() - wait_start;
        if (m_Comm.Rank() == 0 && m_Parameters.verbose > 0)
        {
            std::cout << "Close waited " << wait.count() << " seconds on async threads"
                      << std::endl;
        }
        m_Profiler.Stop("WaitOnAsync");
    }

    if (m_Comm.Rank() == 0)
    {
        if (m_Parameters.AsyncWrite)
        {
            WriteMetadataFileIndex(m_LatestMetaDataPos, m_LatestMetaDataSize);
        }
        // close metadata index file
        UpdateActiveFlag(false);
        m_FileMetadataIndexManager.CloseFiles();
    }

    FlushProfiler();
}

void DaosWriter::FlushProfiler()
{
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

    // m_Profiler.WriteOut(transportTypes, transportProfilers);

    const std::string lineJSON(m_Profiler.GetRankProfilingJSON(transportTypes, transportProfilers) +
                               ",\n");

    const std::vector<char> profilingJSON(m_Profiler.AggregateProfilingJSON(lineJSON));

    if (m_RankMPI == 0)
    {
        // std::cout << "write profiling file!" << std::endl;
        std::string profileFileName;
        if (m_DrainBB)
        {
            // auto bpTargetNames =
            // m_BP4Serializer.GetBPBaseNames({m_Name});
            std::vector<std::string> bpTargetNames = {m_Name};
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
            // auto bpBaseNames =
            // m_BP4Serializer.GetBPBaseNames({m_BBName});
            std::vector<std::string> bpBaseNames = {m_Name};
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

size_t DaosWriter::DebugGetDataBufferSize() const
{
    return m_BP5Serializer.DebugGetDataBufferSize();
}

void DaosWriter::PutCommon(VariableBase &variable, const void *values, bool sync)
{
    if (!m_BetweenStepPairs)
    {
        BeginStep(StepMode::Update);
    }

    // if the user buffer is allocated on the GPU always use sync mode
    auto memSpace = variable.GetMemorySpace(values);
    if (memSpace != MemorySpace::Host)
        sync = true;

    size_t *Shape = NULL;
    size_t *Start = NULL;
    size_t *Count = NULL;
    size_t DimCount = variable.m_Count.size();

    if (variable.m_ShapeID == ShapeID::GlobalArray)
    {
        Shape = variable.m_Shape.data();
        Count = variable.m_Count.data();
        Start = variable.m_Start.data();
    }
    else if (variable.m_ShapeID == ShapeID::LocalArray)
    {
        Count = variable.m_Count.data();
    }
    else if (variable.m_ShapeID == ShapeID::JoinedArray)
    {
        Count = variable.m_Count.data();
        Shape = variable.m_Shape.data();
    }

    size_t ObjSize;
    if (variable.m_Type == DataType::Struct)
    {
        ObjSize = variable.m_ElementSize;
    }
    else
    {
        ObjSize = helper::GetDataTypeSize(variable.m_Type);
    }

    if (!sync)
    {
        /* If arrays is small, force copying to internal buffer to aggregate
         * small writes */
        size_t n = helper::GetTotalSize(variable.m_Count) * ObjSize;
        if (n < m_Parameters.MinDeferredSize)
        {
            sync = true;
        }
    }

    if (!variable.m_MemoryCount.empty())
    {
        int DimCount = variable.m_Count.size();
        std::vector<size_t> ZeroDims(DimCount);
        // get a temporary span then fill with memselection now
        format::BufferV::BufferPos bp5span(0, 0, 0);
        m_BP5Serializer.Marshal((void *)&variable, variable.m_Name.c_str(), variable.m_Type,
                                variable.m_ElementSize, DimCount, Shape, Count, Start, nullptr,
                                false, &bp5span);
        void *ptr = m_BP5Serializer.GetPtr(bp5span.bufferIdx, bp5span.posInBuffer);

        const bool sourceRowMajor = helper::IsRowMajor(m_IO.m_HostLanguage);

        helper::NdCopy((const char *)values, helper::CoreDims(ZeroDims), variable.m_MemoryCount,
                       sourceRowMajor, false, (char *)ptr, variable.m_MemoryStart, variable.m_Count,
                       sourceRowMajor, false, ObjSize, helper::CoreDims(), helper::CoreDims(),
                       helper::CoreDims(), helper::CoreDims(), false /* safemode */, memSpace);
    }
    else
    {
        if (variable.m_Type == DataType::String)
        {
            std::string &source = *(std::string *)values;
            void *p = &(source[0]);
            m_BP5Serializer.Marshal((void *)&variable, variable.m_Name.c_str(), variable.m_Type,
                                    variable.m_ElementSize, DimCount, Shape, Count, Start, &p, sync,
                                    nullptr);
        }
        else
            m_BP5Serializer.Marshal((void *)&variable, variable.m_Name.c_str(), variable.m_Type,
                                    variable.m_ElementSize, DimCount, Shape, Count, Start, values,
                                    sync, nullptr);
    }
}

#define declare_type(T)                                                                            \
    void DaosWriter::DoPut(Variable<T> &variable, typename Variable<T>::Span &span,                \
                           const bool initialize, const T &value)                                  \
    {                                                                                              \
        PERFSTUBS_SCOPED_TIMER("DaosWriter::Put");                                                 \
        PutCommonSpan(variable, span, initialize, value);                                          \
    }

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

#define declare_type(T)                                                                            \
    void DaosWriter::DoPutSync(Variable<T> &variable, const T *data)                               \
    {                                                                                              \
        PutCommon(variable, data, true);                                                           \
    }                                                                                              \
    void DaosWriter::DoPutDeferred(Variable<T> &variable, const T *data)                           \
    {                                                                                              \
        PutCommon(variable, data, false);                                                          \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

#define declare_type(T, L)                                                                         \
    T *DaosWriter::DoBufferData_##L(const int bufferIdx, const size_t payloadPosition,             \
                                    const size_t bufferID) noexcept                                \
    {                                                                                              \
        return reinterpret_cast<T *>(m_BP5Serializer.GetPtr(bufferIdx, payloadPosition));          \
    }

ADIOS2_FOREACH_PRIMITVE_STDTYPE_2ARGS(declare_type)
#undef declare_type

void DaosWriter::DoPutStructSync(VariableStruct &variable, const void *data)
{
    PutCommon(variable, data, true);
}

void DaosWriter::DoPutStructDeferred(VariableStruct &variable, const void *data)
{
    PutCommon(variable, data, false);
}

void DaosWriter::daos_handle_share(daos_handle_t *hdl, int type)
{
    d_iov_t ghdl = {NULL, 0, 0};
    int rc;

    if (m_Comm.Rank() == 0)
    {
        /** fetch size of global handle */
        if (type == DaosWriter::HANDLE_POOL)
            rc = daos_pool_local2global(*hdl, &ghdl);
        else
            rc = daos_cont_local2global(*hdl, &ghdl);
        ASSERT(rc == 0, "local2global failed with %d", rc);
    }

    /** broadcast size of global handle to all peers */
    MPI_Bcast(&ghdl.iov_buf_len, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    /** allocate buffer for global pool handle */
    ghdl.iov_buf = malloc(ghdl.iov_buf_len);
    ghdl.iov_len = ghdl.iov_buf_len;

    if (m_Comm.Rank() == 0)
    {
        /** generate actual global handle to share with peer tasks */
        if (type == DaosWriter::HANDLE_POOL)
            rc = daos_pool_local2global(*hdl, &ghdl);
        else
            rc = daos_cont_local2global(*hdl, &ghdl);
        ASSERT(rc == 0, "local2global failed with %d", rc);
    }

    /** broadcast global handle to all peers */
    MPI_Bcast(ghdl.iov_buf, ghdl.iov_len, MPI_BYTE, 0, MPI_COMM_WORLD);

    if (m_Comm.Rank() != 0)
    {
        /** unpack global handle */
        if (type == DaosWriter::HANDLE_POOL)
        {
            /* NB: Only pool_global2local are different */
            rc = daos_pool_global2local(ghdl, hdl);
        }
        else
        {
            rc = daos_cont_global2local(poh, ghdl, hdl);
        }
        ASSERT(rc == 0, "global2local failed with %d", rc);
    }

    free(ghdl.iov_buf);

    MPI_Barrier(MPI_COMM_WORLD);
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
