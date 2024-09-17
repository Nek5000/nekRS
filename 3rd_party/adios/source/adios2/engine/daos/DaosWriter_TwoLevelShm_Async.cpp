/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DaosWriter.cpp
 *
 */

#include "DaosWriter.h"

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosFunctions.h" //CheckIndexRange, PaddingToAlignOffset
#include "adios2/toolkit/format/buffer/chunk/ChunkV.h"
#include "adios2/toolkit/format/buffer/malloc/MallocV.h"
#include "adios2/toolkit/transport/file/FileFStream.h"
#include <adios2-perfstubs-interface.h>

#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

namespace adios2
{
namespace core
{
namespace engine
{

using namespace adios2::format;

/* Aggregator part of the async two level aggregation Guided version
   This process is the one writing to disk
*/
void DaosWriter::AsyncWriteThread_TwoLevelShm_Aggregator(AsyncWriteInfo *info)
{
    aggregator::MPIShmChain *a = dynamic_cast<aggregator::MPIShmChain *>(info->aggregator);
    uint64_t totalSize = info->totalSize;

    /* Write own data first */
    {
        std::vector<core::iovec> DataVec = info->Data->DataVec();
        const uint64_t mysize = info->Data->Size();
        info->tm->SeekTo(info->startPos);
        AsyncWriteOwnData(info, DataVec, mysize, false);
        totalSize -= mysize;
    }

    /* Write from shm until every non-aggr sent all data */
    std::vector<core::iovec> DataVec(1);
    size_t wrote = 0;
    while (wrote < totalSize)
    {
        /* Write the next shm block now */
        // potentially blocking call waiting on some non-aggr process
        aggregator::MPIShmChain::ShmDataBuffer *b = a->LockConsumerBuffer();
        // b->actual_size: how much we need to write
        DataVec[0].iov_base = b->buf;
        DataVec[0].iov_len = b->actual_size;
        AsyncWriteOwnData(info, DataVec, b->actual_size, false);
        wrote += b->actual_size;
        a->UnlockConsumerBuffer();
    }
}

/* Non-aggregator part of the async two level aggregation.
   This process passes data to Aggregator through SHM segment.
   tokenChain in caller ensures only one process (per aggregator chain)
   is running this function at a time
*/
void DaosWriter::AsyncWriteThread_TwoLevelShm_SendDataToAggregator(aggregator::MPIShmChain *a,
                                                                   format::BufferV *Data)
{
    /* In a loop, copy the local data into the shared memory, alternating
       between the two segments.
    */

    std::vector<core::iovec> DataVec = Data->DataVec();
    size_t nBlocks = DataVec.size();

    // size_t sent = 0;
    size_t block = 0;
    size_t temp_offset = 0;
    while (block < nBlocks)
    {
        // potentially blocking call waiting on Aggregator
        aggregator::MPIShmChain::ShmDataBuffer *b = a->LockProducerBuffer();
        // b->max_size: how much we can copy
        // b->actual_size: how much we actually copy
        b->actual_size = 0;
        while (true)
        {
            /* Copy n bytes from the current block, current offset to shm
                making sure to use up to shm_size bytes
            */
            size_t n = DataVec[block].iov_len - temp_offset;
            if (n > (b->max_size - b->actual_size))
            {
                n = b->max_size - b->actual_size;
            }
            std::memcpy(&b->buf[b->actual_size],
                        (const char *)DataVec[block].iov_base + temp_offset, n);
            b->actual_size += n;

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

            /* Have we reached the max allowed shm size ?*/
            if (b->actual_size >= b->max_size)
            {
                break;
            }
            if (block >= nBlocks)
            {
                break;
            }
        }
        // sent += b->actual_size;
        a->UnlockProducerBuffer();
    }
}

int DaosWriter::AsyncWriteThread_TwoLevelShm(AsyncWriteInfo *info)
{
    /* DO NOT use MPI in this separate thread, including destroying
       shm segments explicitely (a->DestroyShm) or implicitely (tokenChain) */
    Seconds ts = Now() - info->tstart;
    // std::cout << "ASYNC rank " << info->rank_global
    //          << " starts at: " << ts.count() << std::endl;
    aggregator::MPIShmChain *a = dynamic_cast<aggregator::MPIShmChain *>(info->aggregator);
    if (a->m_IsAggregator)
    {
        // Send token to first non-aggregator to start filling shm
        // Also informs next process its starting offset (for correct
        // metadata)
        uint64_t nextWriterPos = info->startPos + info->Data->Size();
        info->tokenChain->SendToken(nextWriterPos);
        AsyncWriteThread_TwoLevelShm_Aggregator(info);
        info->tokenChain->RecvToken();
    }
    else
    {
        // non-aggregators fill shared buffer in marching order
        // they also receive their starting offset this way
        uint64_t startPos = info->tokenChain->RecvToken();
        AsyncWriteThread_TwoLevelShm_SendDataToAggregator(a, info->Data);
        uint64_t nextWriterPos = startPos + info->Data->Size();
        info->tokenChain->SendToken(nextWriterPos);
    }
    delete info->Data;

    ts = Now() - info->tstart;
    /*std::cout << "ASYNC " << info->rank_global << " ended at: " << ts.count()
              << std::endl;*/
    return 1;
};

void DaosWriter::WriteData_TwoLevelShm_Async(format::BufferV *Data)
{
    aggregator::MPIShmChain *a = dynamic_cast<aggregator::MPIShmChain *>(m_Aggregator);

    // new step writing starts at offset m_DataPos on master aggregator
    // other aggregators to the same file will need to wait for the position
    // to arrive from the rank below

    // align to PAGE_SIZE (only valid on master aggregator at this point)
    m_DataPos += helper::PaddingToAlignOffset(m_DataPos, m_Parameters.StripeSize);

    // Each aggregator needs to know the total size they write
    // This calculation is valid on aggregators only
    std::vector<uint64_t> mySizes = a->m_Comm.GatherValues(Data->Size());
    uint64_t myTotalSize = 0;
    uint64_t maxSize = 0;
    for (auto s : mySizes)
    {
        myTotalSize += s;
        if (s > maxSize)
        {
            maxSize = s;
        }
    }

    if (a->m_Comm.Size() > 1)
    {
        size_t alignment_size = sizeof(max_align_t);
        if (m_Parameters.DirectIO)
        {
            alignment_size = m_Parameters.DirectIOAlignOffset;
        }
        a->CreateShm(static_cast<size_t>(maxSize), m_Parameters.MaxShmSize, alignment_size);
    }

    if (a->m_IsAggregator)
    {
        // In each aggregator chain, send from master down the line
        // these total sizes, so every aggregator knows where to start
        if (a->m_AggregatorChainComm.Rank() > 0)
        {
            a->m_AggregatorChainComm.Recv(
                &m_DataPos, 1, a->m_AggregatorChainComm.Rank() - 1, 0,
                "AggregatorChain token in DaosWriter::WriteData_TwoLevelShm");
            // align to PAGE_SIZE
            m_DataPos += helper::PaddingToAlignOffset(m_DataPos, m_Parameters.StripeSize);
        }
        m_StartDataPos = m_DataPos; // metadata needs this info
        if (a->m_AggregatorChainComm.Rank() < a->m_AggregatorChainComm.Size() - 1)
        {
            uint64_t nextWriterPos = m_DataPos + myTotalSize;
            a->m_AggregatorChainComm.Isend(&nextWriterPos, 1, a->m_AggregatorChainComm.Rank() + 1,
                                           0, "Chain token in DaosWriter::WriteData");
        }
        else if (a->m_AggregatorChainComm.Size() > 1)
        {
            // send back final position from last aggregator in file to master
            // aggregator
            uint64_t nextWriterPos = m_DataPos + myTotalSize;
            a->m_AggregatorChainComm.Isend(&nextWriterPos, 1, 0, 0,
                                           "Chain token in DaosWriter::WriteData");
        }

        // Master aggregator needs to know where the last writing ended by the
        // last aggregator in the chain, so that it can start from the correct
        // position at the next output step
        if (!a->m_AggregatorChainComm.Rank())
        {
            if (a->m_AggregatorChainComm.Size() > 1)
            {
                a->m_AggregatorChainComm.Recv(&m_DataPos, 1, a->m_AggregatorChainComm.Size() - 1, 0,
                                              "Chain token in DaosWriter::WriteData");
            }
            else
            {
                m_DataPos = m_StartDataPos + myTotalSize;
            }
        }
    }

    /*std::cout << "Rank " << m_Comm.Rank() << "  start data async "
              << " to subfile " << a->m_SubStreamIndex << " at pos "
              << m_StartDataPos << std::endl;*/

    m_AsyncWriteInfo = new AsyncWriteInfo();
    m_AsyncWriteInfo->aggregator = m_Aggregator;
    m_AsyncWriteInfo->rank_global = m_Comm.Rank();
    m_AsyncWriteInfo->rank_chain = a->m_Comm.Rank();
    m_AsyncWriteInfo->nproc_chain = a->m_Comm.Size();
    m_AsyncWriteInfo->comm_chain = helper::Comm(); // unused in this aggregation
    m_AsyncWriteInfo->tstart = m_EngineStart;
    m_AsyncWriteInfo->tokenChain = new shm::TokenChain<uint64_t>(&a->m_Comm);
    m_AsyncWriteInfo->tm = &m_FileDataManager;
    m_AsyncWriteInfo->Data = Data;
    m_AsyncWriteInfo->flagRush = &m_flagRush;
    m_AsyncWriteInfo->lock = &m_AsyncWriteLock;

    // Metadata collection needs m_StartDataPos correctly set on
    // every process before we call the async writing thread
    if (a->m_IsAggregator)
    {
        // Informs next process its starting offset (for correct metadata)
        uint64_t nextWriterPos = m_StartDataPos + Data->Size();
        m_AsyncWriteInfo->tokenChain->SendToken(nextWriterPos);
        m_AsyncWriteInfo->tokenChain->RecvToken();
    }
    else
    {
        // non-aggregators fill shared buffer in marching order
        // they also receive their starting offset this way
        m_StartDataPos = m_AsyncWriteInfo->tokenChain->RecvToken();
        uint64_t nextWriterPos = m_StartDataPos + Data->Size();
        m_AsyncWriteInfo->tokenChain->SendToken(nextWriterPos);
    }

    // Launch data writing thread, m_StartDataPos is valid
    // m_DataPos is already pointing to the end of the write, do not use here.
    m_AsyncWriteInfo->startPos = m_StartDataPos;
    m_AsyncWriteInfo->totalSize = myTotalSize;
    m_AsyncWriteInfo->deadline = m_ExpectedTimeBetweenSteps.count();

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

    m_WriteFuture = std::async(std::launch::async, AsyncWriteThread_TwoLevelShm, m_AsyncWriteInfo);

    /* At this point it is prohibited in the main thread
       - to modify Data, which will be deleted in the async thread any tiume
       - to use m_FileDataManager until next BeginStep, which is being used
         in the async thread to write data
    */
}

void DaosWriter::AsyncWriteDataCleanup_TwoLevelShm()
{
    aggregator::MPIShmChain *a =
        dynamic_cast<aggregator::MPIShmChain *>(m_AsyncWriteInfo->aggregator);
    if (a->m_Comm.Size() > 1)
    {
        a->DestroyShm();
    }
    delete m_AsyncWriteInfo->tokenChain;
    delete m_AsyncWriteInfo;
    m_AsyncWriteInfo = nullptr;
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
