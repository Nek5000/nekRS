/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP5Writer.cpp
 *
 */

#include "DaosWriter.h"
#include "DaosWriter.tcc"

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosFunctions.h" //CheckIndexRange
#include "adios2/toolkit/format/buffer/chunk/ChunkV.h"
#include "adios2/toolkit/format/buffer/malloc/MallocV.h"
#include "adios2/toolkit/transport/file/FileFStream.h"
#include <adios2-perfstubs-interface.h>

#include <algorithm> // max
#include <ctime>
#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

using namespace adios2::format;

DaosWriter::ComputationStatus DaosWriter::IsInComputationBlock(AsyncWriteInfo *info,
                                                               size_t &compBlockIdx)
{
    ComputationStatus compStatus = ComputationStatus::NotInComp_ExpectMore;
    size_t nExpectedBlocks = info->expectedComputationBlocks.size();

    if (compBlockIdx >= nExpectedBlocks)
    {
        compStatus = ComputationStatus::NoMoreComp;
    }
    else
    {
        bool inComp = false;
        size_t compBlockID = 0;
        // access variables modified by main thread to avoid data race
        info->lock->lock();
        compBlockID = *info->currentComputationBlockID;
        inComp = *info->inComputationBlock;
        info->lock->unlock();

        /* Track which computation block we are in */
        if (inComp)
        {
            while (compBlockIdx < nExpectedBlocks &&
                   info->expectedComputationBlocks[compBlockIdx].blockID < compBlockID)
            {
                ++compBlockIdx;
            }
            if (info->expectedComputationBlocks[compBlockIdx].blockID > compBlockID)
            {
                // the current computation block is a short one that was not
                // recorded
                compStatus = ComputationStatus::NotInComp_ExpectMore;
            }
            else
            {
                compStatus = ComputationStatus::InComp;
            }
        }
    }
    return compStatus;
}

void DaosWriter::AsyncWriteOwnData(AsyncWriteInfo *info, std::vector<core::iovec> &DataVec,
                                   const size_t totalsize, const bool seekOnFirstWrite)
{
    /* local variables to track variables modified by main thread */
    size_t compBlockIdx = 0; /* position in vector to get length */

    /* In a loop, write the data in smaller blocks */
    size_t nBlocks = DataVec.size();
    size_t wrote = 0;
    size_t block = 0;
    size_t temp_offset = 0;
    size_t max_size = std::max(1024 * 1024UL, totalsize / 100UL);

    bool firstWrite = seekOnFirstWrite;
    while (block < nBlocks)
    {
        bool doRush = false;
        bool doSleep = false;

        info->lock->lock();
        doRush = *info->flagRush;
        info->lock->unlock();

        if (!doRush)
        {
            ComputationStatus compStatus = IsInComputationBlock(info, compBlockIdx);

            /* Scheduling decisions:
               Cases:
               1. Not in a computation block AND we still expect more
               computation blocks down the line ==> Sleep
               2. In computation block ==> Write
               3. We are at the end of a computation block (how close??) AND we
               still expect more computation blocks down the line 3. ==> Sleep
               4. We are at the end of the LAST computation block ==> Write
               5. No more computation blocks expected ==> Write all at once
               6. Main thread set flagRush ==> Write all at once
               -- case 3 not handled yet properly
            */

            switch (compStatus)
            {
            case ComputationStatus::NotInComp_ExpectMore:
                // case 1
                doSleep = true;
                break;
            case ComputationStatus::NoMoreComp:
                // case 5
                doRush = true;
                break;
            default:
                // cases 2, 3, 4
                break;
            }
        }

        if (doRush)
        {
            auto vec = std::vector<adios2::core::iovec>(DataVec.begin() + block, DataVec.end());
            vec[0].iov_base = (const char *)DataVec[block].iov_base + temp_offset;
            vec[0].iov_len = DataVec[block].iov_len - temp_offset;
            size_t pos = MaxSizeT; // <==> no seek inside WriteFileAt
            if (firstWrite)
            {
                pos = info->startPos + wrote; // seek to pos
            }
            /*std::cout << "Async write on Rank " << info->rank_global
                      << " write the rest of  " << totalsize - wrote
                      << " bytes at pos " << pos << std::endl;*/

            info->tm->WriteFileAt(vec.data(), vec.size(), pos);

            break; /* Exit loop after this final write */
        }

        if (doSleep)
        {
            std::this_thread::sleep_for(core::Seconds(0.01));
            continue;
        }

        /* Write next batch of data */

        /* Get the next n bytes from the current block, current offset */
        size_t n = DataVec[block].iov_len - temp_offset;
        if (n > max_size)
        {
            n = max_size;
        }

        if (firstWrite)
        {
            info->tm->WriteFileAt((const char *)DataVec[block].iov_base + temp_offset, n,
                                  info->startPos);
            firstWrite = false;
        }
        else
        {
            info->tm->WriteFiles((const char *)DataVec[block].iov_base + temp_offset, n);
        }

        /* Have we processed the entire block or staying with it? */
        if (n + temp_offset < DataVec[block].iov_len)
        {
            temp_offset += n;
        }
        else
        {
            temp_offset = 0;
            ++block;
        }
        wrote += n;
    }
};

int DaosWriter::AsyncWriteThread_EveryoneWrites(AsyncWriteInfo *info)
{
    if (info->tokenChain)
    {
        if (info->rank_chain > 0)
        {
            info->tokenChain->RecvToken();
        }
    }

    std::vector<core::iovec> DataVec = info->Data->DataVec();
    const uint64_t mysize = info->Data->Size();
    AsyncWriteOwnData(info, DataVec, mysize, true);

    if (info->tokenChain)
    {
        uint64_t t = 1;
        info->tokenChain->SendToken(t);
        if (!info->rank_chain)
        {
            info->tokenChain->RecvToken();
        }
    }
    delete info->Data;
    return 1;
};

void DaosWriter::WriteData_EveryoneWrites_Async(format::BufferV *Data, bool SerializedWriters)
{

    const aggregator::MPIChain *a = dynamic_cast<aggregator::MPIChain *>(m_Aggregator);

    // new step writing starts at offset m_DataPos on aggregator
    // others will wait for the position to arrive from the rank below

    if (a->m_Comm.Rank() > 0)
    {
        a->m_Comm.Recv(&m_DataPos, 1, a->m_Comm.Rank() - 1, 0,
                       "Chain token in DaosWriter::WriteData_EveryoneWrites_Async");
    }

    // align to PAGE_SIZE
    m_DataPos += helper::PaddingToAlignOffset(m_DataPos, m_Parameters.StripeSize);
    m_StartDataPos = m_DataPos;

    if (a->m_Comm.Rank() < a->m_Comm.Size() - 1)
    {
        uint64_t nextWriterPos = m_DataPos + Data->Size();
        a->m_Comm.Isend(&nextWriterPos, 1, a->m_Comm.Rank() + 1, 0,
                        "Chain token in DaosWriter::WriteData_EveryoneWrites_Async");
    }

    m_DataPos += Data->Size();

    /* a->comm can span multiple nodes but we need comm inside a node
       when doing serialized aggregation */
    m_AsyncWriteInfo = new AsyncWriteInfo();
    m_AsyncWriteInfo->aggregator = nullptr;
    m_AsyncWriteInfo->rank_global = m_Comm.Rank();
    if (SerializedWriters)
    {
        m_AsyncWriteInfo->comm_chain = a->m_Comm.GroupByShm();
        m_AsyncWriteInfo->rank_chain = m_AsyncWriteInfo->comm_chain.Rank();
        m_AsyncWriteInfo->nproc_chain = m_AsyncWriteInfo->comm_chain.Size();
        m_AsyncWriteInfo->tokenChain = new shm::TokenChain<uint64_t>(&m_AsyncWriteInfo->comm_chain);
    }
    else
    {
        m_AsyncWriteInfo->comm_chain = helper::Comm(); // not needed
        m_AsyncWriteInfo->rank_chain = a->m_Comm.Rank();
        m_AsyncWriteInfo->nproc_chain = a->m_Comm.Size();
        m_AsyncWriteInfo->tokenChain = nullptr;
    }
    m_AsyncWriteInfo->tstart = m_EngineStart;
    m_AsyncWriteInfo->tm = &m_FileDataManager;
    m_AsyncWriteInfo->Data = Data;
    m_AsyncWriteInfo->startPos = m_StartDataPos;
    m_AsyncWriteInfo->totalSize = Data->Size();
    m_AsyncWriteInfo->deadline = m_ExpectedTimeBetweenSteps.count();
    m_AsyncWriteInfo->flagRush = &m_flagRush;
    m_AsyncWriteInfo->lock = &m_AsyncWriteLock;

    if (m_ComputationBlocksLength > 0.0 && m_Parameters.AsyncWrite == (int)AsyncWrite::Guided)
    {
        m_AsyncWriteInfo->inComputationBlock = &m_InComputationBlock;
        m_AsyncWriteInfo->computationBlocksLength = m_ComputationBlocksLength;
        if (m_AsyncWriteInfo->deadline < m_ComputationBlocksLength)
        {
            m_AsyncWriteInfo->deadline = m_ComputationBlocksLength;
        }
        m_AsyncWriteInfo->expectedComputationBlocks = m_ComputationBlockTimes; // copy!
        m_AsyncWriteInfo->currentComputationBlocks = &m_ComputationBlockTimes; // ptr!
        m_AsyncWriteInfo->currentComputationBlockID = &m_ComputationBlockID;

        /* Clear current block tracker now so that async thread does not get
        confused with the past info */
        m_ComputationBlockTimes.clear();
        m_ComputationBlocksLength = 0.0;
        m_ComputationBlockID = 0;
    }
    else
    {
        if (m_Parameters.AsyncWrite == (int)AsyncWrite::Naive)
        {
            m_AsyncWriteInfo->deadline = 0;
        }
        m_AsyncWriteInfo->inComputationBlock = nullptr;
        m_AsyncWriteInfo->computationBlocksLength = 0.0;
        m_AsyncWriteInfo->currentComputationBlocks = nullptr;
        m_AsyncWriteInfo->currentComputationBlockID = nullptr;
    }

    m_WriteFuture =
        std::async(std::launch::async, AsyncWriteThread_EveryoneWrites, m_AsyncWriteInfo);

    // At this point modifying Data in main thread is prohibited !!!

    if (a->m_Comm.Size() > 1)
    {
        // at the end, last rank sends back the final data pos to first rank
        // so it can update its data pos
        if (a->m_Comm.Rank() == a->m_Comm.Size() - 1)
        {
            a->m_Comm.Isend(&m_DataPos, 1, 0, 0,
                            "Final chain token in "
                            "DaosWriter::WriteData_EveryoneWrites_Async");
        }
        if (a->m_Comm.Rank() == 0)
        {
            a->m_Comm.Recv(&m_DataPos, 1, a->m_Comm.Size() - 1, 0,
                           "Chain token in DaosWriter::WriteData_EveryoneWrites_Async");
        }
    }
}

void DaosWriter::AsyncWriteDataCleanup_EveryoneWrites()
{
    if (m_AsyncWriteInfo->tokenChain)
    {
        delete m_AsyncWriteInfo->tokenChain;
    }
    delete m_AsyncWriteInfo;
    m_AsyncWriteInfo = nullptr;
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
