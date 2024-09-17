/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * MPIShmChain.h
 *
 *  Created on: July 5, 2021
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */
#include "MPIShmChain.h"

#include "adios2/helper/adiosMemory.h" // PaddingToAlignOffset

#include <iostream>

namespace adios2
{
namespace aggregator
{

MPIShmChain::MPIShmChain() : MPIAggregator() {}

MPIShmChain::~MPIShmChain() { Close(); }

void MPIShmChain::Close()
{
    if (m_IsActive)
    {
        m_NodeComm.Free("free per-node comm in ~MPIShmChain()");
        m_OnePerNodeComm.Free("free chain of nodes in ~MPIShmChain()");
        m_AllAggregatorsComm.Free("free comm of all aggregators in ~MPIShmChain()");
        m_AggregatorChainComm.Free("free chains of aggregators in ~MPIShmChain()");
    }
    MPIAggregator::Close();
}

size_t MPIShmChain::PreInit(helper::Comm const &parentComm)
{
    /* Communicator connecting ranks on each Compute Node */
    m_NodeComm = parentComm.GroupByShm("creating per-node comm at Open");
    int NodeRank = m_NodeComm.Rank();

    /*
     *  Communicators connecting rank N of each node
     *  We are only interested in the chain of rank 0s
     */
    int color = (NodeRank ? 1 : 0);
    m_OnePerNodeComm = parentComm.Split(color, 0, "creating chain of nodes at Open");

    /* Number of nodes */
    if (!NodeRank)
    {
        m_NumNodes = static_cast<size_t>(m_OnePerNodeComm.Size());
    }
    m_NumNodes = m_NodeComm.BroadcastValue<size_t>(m_NumNodes, 0);
    PreInitCalled = true;
    return m_NumNodes;
}

void MPIShmChain::Init(const size_t numAggregators, const size_t subStreams,
                       helper::Comm const &parentComm)
{
    if (!PreInitCalled)
    {
        PreInit(parentComm);
    }

    // int AllRank = parentComm.Rank();
    // int AllSize = parentComm.Size();
    int NodeRank = m_NodeComm.Rank();
    size_t NodeSize = static_cast<size_t>(m_NodeComm.Size());

    /* Number of aggregators per node */
    size_t aggregatorPerNode = numAggregators / m_NumNodes;
    if (aggregatorPerNode == 0)
    {
        aggregatorPerNode = 1; /* default */
    }
    if (aggregatorPerNode > NodeSize)
    {
        aggregatorPerNode = NodeSize;
    }

    /* Create main communicator that splits the node comm into one or more
     * aggregator chains */
    float k = static_cast<float>(NodeSize) / static_cast<float>(aggregatorPerNode);
    float c = static_cast<float>(NodeRank) / k;
    int color = static_cast<int>(c);
    m_Comm = m_NodeComm.Split(color, 0, "creating aggregator groups at Open");
    m_Rank = m_Comm.Rank();
    m_Size = m_Comm.Size();
    if (m_Rank != 0)
    {
        m_IsAggregator = false;
        m_IsMasterAggregator = false;
    }

    /* Identify parent rank of aggregator process within each chain */
    if (!m_Rank)
    {
        m_AggregatorRank = parentComm.Rank();
    }
    m_AggregatorRank = m_Comm.BroadcastValue<int>(m_AggregatorRank, 0);

    /* Communicator for all Aggregators */
    color = (m_Rank ? 1 : 0);
    m_AllAggregatorsComm = parentComm.Split(color, 0, "creating comm of all aggregators at Open");

    /* Total number of aggregators */
    if (!NodeRank)
    {
        m_NumAggregators = static_cast<size_t>(m_AllAggregatorsComm.Size());
    }
    m_NumAggregators = m_NodeComm.BroadcastValue<size_t>(m_NumAggregators);

    /* Number of substreams */
    m_SubStreams = subStreams;
    if (m_SubStreams == 0)
    {
        m_SubStreams = m_NumAggregators; /* default */
    }
    if (m_SubStreams > m_NumAggregators)
    {
        m_SubStreams = m_NumAggregators;
    }

    if (!m_Rank)
    {
        k = static_cast<float>(m_NumAggregators) / static_cast<float>(m_SubStreams);
        /* 1.0 <= k <= m_NumAggregators */
        c = static_cast<float>(m_AllAggregatorsComm.Rank()) / k;
        m_SubStreamIndex = static_cast<int>(c);
    }
    m_SubStreamIndex = m_Comm.BroadcastValue<size_t>(m_SubStreamIndex);

    /* Create the communicator to connect aggregators writing to the same
     * substream */
    color = static_cast<int>(m_SubStreamIndex);
    m_AggregatorChainComm =
        m_AllAggregatorsComm.Split(color, 0, "creating chains of aggregators at Open");

    if (m_AggregatorChainComm.Rank() != 0)
    {
        m_IsMasterAggregator = false;
    }
    HandshakeStruct hsAC; // need persistence until complete
    if (m_AggregatorChainComm.Size() > 1)
    {
        HandshakeLinks_Start(m_AggregatorChainComm, hsAC);
    }

    m_IsActive = true;

    if (m_Comm.Size() > 1)
    {
        HandshakeStruct hs; // need persistence until complete
        HandshakeLinks_Start(m_Comm, hs);

        /* Create the shared memory segment */
        // CreateShm();

        HandshakeLinks_Complete(hs);
    }

    if (m_AggregatorChainComm.Size() > 1)
    {
        HandshakeLinks_Complete(hsAC);
    }
}

// PRIVATE
void MPIShmChain::HandshakeLinks_Start(helper::Comm &comm, HandshakeStruct &hs)
{
    int rank = comm.Rank();
    hs.sendToken = rank;
    if (rank < comm.Size() - 1) // send to next
    {
        hs.sendRequest = comm.Isend(&hs.sendToken, 1, rank + 1, 0,
                                    "Isend handshake with neighbor, MPIChain aggregator, at Open");
    }
    else // send to 0 to close the loop
    {
        hs.sendRequest = comm.Isend(&hs.sendToken, 1, 0, 0,
                                    "Isend handshake with rank 0, MPIChain aggregator, at Open");
    }

    if (comm.Rank() > 0) // receive from previous
    {
        hs.recvRequest = comm.Irecv(&hs.recvToken, 1, rank - 1, 0,
                                    "Irecv handshake with neighbor, MPIChain aggregator, at Open");
    }
    else // rank 0 receives from last
    {
        hs.recvRequest = comm.Irecv(&hs.recvToken, 1, comm.Size() - 1, 0,
                                    "Irecv handshake with neighbor, MPIChain aggregator, at Open");
    }
}

void MPIShmChain::HandshakeLinks_Complete(HandshakeStruct &hs)
{
    hs.recvRequest.Wait("Wait handshake with neighbor (recv), MPIChain "
                        "aggregator, at Open");
    hs.sendRequest.Wait("Wait handshake with neighbor (send), MPIChain "
                        "aggregator, at Open");
}

void MPIShmChain::CreateShm(size_t blocksize, const size_t maxsegmentsize,
                            const size_t alignment_size)
{
    if (!m_Comm.IsMPI())
    {
        helper::Throw<std::runtime_error>("Toolkit", "aggregator::mpi::MPIShmChain", "CreateShm",
                                          "called with a non-MPI communicator");
    }
    char *ptr;
    size_t structsize = sizeof(ShmSegment);
    structsize += helper::PaddingToAlignOffset(structsize, alignment_size);
    if (!m_Rank)
    {
        blocksize += helper::PaddingToAlignOffset(blocksize, alignment_size);
        size_t totalsize = structsize + 2 * blocksize;
        if (totalsize > maxsegmentsize)
        {
            // roll back and calculate sizes from maxsegmentsize
            totalsize = maxsegmentsize - alignment_size + 1;
            totalsize += helper::PaddingToAlignOffset(totalsize, alignment_size);
            blocksize = (totalsize - structsize) / 2 - alignment_size + 1;
            blocksize += helper::PaddingToAlignOffset(blocksize, alignment_size);
            totalsize = structsize + 2 * blocksize;
        }
        m_Win = m_Comm.Win_allocate_shared(totalsize, 1, &ptr);
    }
    else
    {
        m_Win = m_Comm.Win_allocate_shared(0, 1, &ptr);
        size_t shmsize;
        int disp_unit;
        m_Comm.Win_shared_query(m_Win, 0, &shmsize, &disp_unit, &ptr);
        blocksize = (shmsize - structsize) / 2;
    }
    m_Shm = reinterpret_cast<ShmSegment *>(ptr);
    m_ShmBufA = ptr + structsize;
    m_ShmBufB = m_ShmBufA + blocksize;

    if (!m_Rank)
    {
        m_Shm->producerBuffer = LastBufferUsed::None;
        m_Shm->consumerBuffer = LastBufferUsed::None;
        m_Shm->NumBuffersFull = 0;
        m_Shm->sdbA.buf = nullptr;
        m_Shm->sdbA.max_size = blocksize;
        m_Shm->sdbB.buf = nullptr;
        m_Shm->sdbB.max_size = blocksize;
    }
    /*std::cout << "Rank " << m_Rank << " shm = " << ptr
              << " bufA = " << static_cast<void *>(m_Shm->bufA)
              << " bufB = " << static_cast<void *>(m_Shm->bufB) << std::endl;*/
}

void MPIShmChain::DestroyShm() { m_Comm.Win_free(m_Win); }

/*
   The buffering strategy is the following.
   Assumptions: 1. Only one Producer (and one Consumer) is active at a time.

   The first Producer fills buffer A first then B and then is always
   alternating, blocking when Consumer is behind (NumBuffersFull == 2).
   The next Producer will continue with the alternating pattern where the
   previous Producer has finished.

   The Consumer is blocked until there is at least one buffer available.
   It takes buffer A at the first call, then it alternates between the two
   buffers.

   C++ atomic locks in Shm for IPC locking:
   - lockSegment is used to modify the m_Shm variables exclusively (short time),
   - lockA and lockB are used to give long term exclusive access to one buffer
   to be filled or consumed.

   The sleeping phases, to wait on the other party to catch up, are outside of
   the locking code areas.

   Note: the m_Shm->sdbX.buf pointers must be set on the local process every
   time, even tough it is stored on the shared memory segment, because the
   address of the segment is different on every process. Failing to set on the
   local process causes this pointer pointing to an invalid address (set on
   another process).

   Note: the sdbA and sdbB structs are stored on the shared memory segment
   because they contain 'actual_size' which is set on the Producer and used by
   the Consumer.

*/

MPIShmChain::ShmDataBuffer *MPIShmChain::LockProducerBuffer()
{
    MPIShmChain::ShmDataBuffer *sdb = nullptr;

    // Sleep until there is a buffer available at all
    while (m_Shm->NumBuffersFull == 2)
    {
        std::this_thread::sleep_for(std::chrono::duration<double>(0.00001));
    }

    m_Shm->lockSegment.lock();
    if (m_Shm->producerBuffer == LastBufferUsed::A)

    {
        m_Shm->producerBuffer = LastBufferUsed::B;
        sdb = &m_Shm->sdbB;
        // point to shm data buffer (in local process memory)
        sdb->buf = m_ShmBufB;
    }
    else // None or B
    {
        m_Shm->producerBuffer = LastBufferUsed::A;
        sdb = &m_Shm->sdbA;
        // point to shm data buffer (in local process memory)
        sdb->buf = m_ShmBufA;
    }
    m_Shm->lockSegment.unlock();

    // We determined we want a specific buffer
    // Now we need to get a lock on it in case consumer is using it
    if (m_Shm->producerBuffer == LastBufferUsed::A)
    {
        m_Shm->lockA.lock();
    }
    else
    {
        m_Shm->lockB.lock();
    }

    return sdb;
}

void MPIShmChain::UnlockProducerBuffer()
{
    m_Shm->lockSegment.lock();
    ++m_Shm->NumBuffersFull;
    m_Shm->lockSegment.unlock();

    if (m_Shm->producerBuffer == LastBufferUsed::A)
    {
        m_Shm->lockA.unlock();
    }
    else
    {
        m_Shm->lockB.unlock();
    }
}

MPIShmChain::ShmDataBuffer *MPIShmChain::LockConsumerBuffer()
{
    MPIShmChain::ShmDataBuffer *sdb = nullptr;

    // Sleep until there is at least one buffer filled
    while (m_Shm->NumBuffersFull < 1)
    {
        std::this_thread::sleep_for(std::chrono::duration<double>(0.00001));
    }
    // At this point we know buffer A has content or going to have content
    // when we successfully lock it

    m_Shm->lockSegment.lock();
    if (m_Shm->consumerBuffer == LastBufferUsed::A)

    {
        m_Shm->consumerBuffer = LastBufferUsed::B;
        sdb = &m_Shm->sdbB;
        // point to shm data buffer (in local process memory)
        sdb->buf = m_ShmBufB;
    }
    else // None or B
    {
        m_Shm->consumerBuffer = LastBufferUsed::A;
        sdb = &m_Shm->sdbA;
        // point to shm data buffer (in local process memory)
        sdb->buf = m_ShmBufA;
    }
    m_Shm->lockSegment.unlock();

    // We determined we want a specific buffer
    // Now we need to get a lock on it in case producer is using it
    if (m_Shm->consumerBuffer == LastBufferUsed::A)
    {
        m_Shm->lockA.lock();
    }
    else
    {
        m_Shm->lockB.lock();
    }

    return sdb;
}

void MPIShmChain::UnlockConsumerBuffer()
{
    m_Shm->lockSegment.lock();
    --m_Shm->NumBuffersFull;
    m_Shm->lockSegment.unlock();

    if (m_Shm->consumerBuffer == LastBufferUsed::A)
    {
        m_Shm->lockA.unlock();
    }
    else
    {
        m_Shm->lockB.unlock();
    }
}

} // end namespace aggregator
} // end namespace adios2
