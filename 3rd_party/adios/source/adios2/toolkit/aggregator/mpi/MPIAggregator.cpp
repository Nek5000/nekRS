/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * MPIAggregator.cpp
 *
 *  Created on: Feb 20, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "MPIAggregator.h"

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace aggregator
{

MPIAggregator::MPIAggregator() {}

MPIAggregator::~MPIAggregator() { Close(); }

void MPIAggregator::Init(const size_t numAggregators, const size_t subStreams,
                         helper::Comm const &parentComm)

{
    m_NumAggregators = numAggregators;
    m_SubStreams = subStreams;
}

void MPIAggregator::Close()
{
    if (m_IsActive)
    {
        m_Comm.Free("freeing aggregators comm at Close\n");
        m_IsActive = false;
    }
}

// PROTECTED
void MPIAggregator::InitComm(const size_t subStreams, helper::Comm const &parentComm)
{
    int parentRank = parentComm.Rank();
    int parentSize = parentComm.Size();

    const size_t process = static_cast<size_t>(parentRank);
    const size_t processes = static_cast<size_t>(parentSize);

    // Divide the processes into S=subStreams groups.
    const size_t q = processes / subStreams;
    const size_t r = processes % subStreams;

    // Groups [0,r) have size q+1.  Groups [r,S) have size q.
    const size_t firstInSmallGroups = r * (q + 1);

    // Within each group the first process becomes its consumer.
    if (process >= firstInSmallGroups)
    {
        m_SubStreamIndex = r + (process - firstInSmallGroups) / q;
        m_AggregatorRank = static_cast<int>(firstInSmallGroups + (m_SubStreamIndex - r) * q);
    }
    else
    {
        m_SubStreamIndex = process / (q + 1);
        m_AggregatorRank = static_cast<int>(m_SubStreamIndex * (q + 1));
    }

    m_Comm = parentComm.Split(m_AggregatorRank, parentRank,
                              "creating aggregators comm with split at Open");

    m_Rank = m_Comm.Rank();
    m_Size = m_Comm.Size();

    if (m_Rank != 0)
    {
        m_IsAggregator = false;
    }

    m_IsActive = true;
    m_SubStreams = subStreams;
}

void MPIAggregator::InitCommOnePerNode(helper::Comm const &parentComm)
{
    m_Comm = parentComm.GroupByShm("creating default aggregator setup at Open");
    m_Rank = m_Comm.Rank();
    m_Size = m_Comm.Size();

    if (m_Rank != 0)
    {
        m_IsAggregator = false;
    }

    m_IsActive = true;

    /* Determine number of aggregators (= nodes) */

    /*
     *  Communicators connecting rank N of each node
     *  We are only interested in the chain of rank 0s
     */
    int color = (m_Rank ? 1 : 0);
    helper::Comm onePerNodeComm =
        parentComm.Split(color, 0, "creating default aggregator setup at Open");

    if (!m_Rank)
    {
        m_SubStreamIndex = static_cast<size_t>(onePerNodeComm.Rank());
        m_SubStreams = static_cast<size_t>(onePerNodeComm.Size());
    }
    m_SubStreams = m_Comm.BroadcastValue<size_t>(m_SubStreams, 0);
    m_SubStreamIndex = m_Comm.BroadcastValue<size_t>(m_SubStreamIndex, 0);

    /* Identify parent rank of aggregator process within each group */
    if (!m_Rank)
    {
        m_AggregatorRank = parentComm.Rank();
    }
    m_AggregatorRank = m_Comm.BroadcastValue<int>(m_AggregatorRank, 0);
}

void MPIAggregator::HandshakeRank(const int rank)
{
    int message = -1;
    if (m_Rank == rank)
    {
        message = m_Rank;
    }

    m_Comm.Bcast(&message, 1, rank, "handshake with aggregator rank 0 at Open");
}

} // end namespace aggregator
} // end namespace adios2
