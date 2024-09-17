/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscWriterBase.cpp
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#include "SscWriterBase.h"
#include "adios2/helper/adiosMpiHandshake.h"
#include "adios2/helper/adiosString.h"

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

SscWriterBase::SscWriterBase(IO &io, const std::string &name, const Mode mode, MPI_Comm comm)
: m_Name(name), m_IO(io)
{

    helper::GetParameter(io.m_Parameters, "Verbose", m_Verbosity);
    helper::GetParameter(io.m_Parameters, "Threading", m_Threading);
    helper::GetParameter(io.m_Parameters, "OpenTimeoutSecs", m_OpenTimeoutSecs);

    int providedMpiMode;
    MPI_Query_thread(&providedMpiMode);
    if (m_Threading && providedMpiMode != MPI_THREAD_MULTIPLE)
    {
        m_Threading = false;
    }

    SyncMpiPattern(comm);
}

SscWriterBase::~SscWriterBase() = default;

void SscWriterBase::SyncMpiPattern(MPI_Comm comm)
{
    MPI_Group streamGroup;
    MPI_Group writerGroup;
    MPI_Comm readerComm;

    helper::HandshakeComm(m_Name, 'w', m_OpenTimeoutSecs, comm, streamGroup, writerGroup,
                          m_ReaderGroup, m_StreamComm, m_WriterComm, readerComm, m_Verbosity);

    MPI_Comm_rank(comm, &m_WriterRank);
    MPI_Comm_size(comm, &m_WriterSize);

    MPI_Comm_rank(m_StreamComm, &m_StreamRank);
    MPI_Comm_size(m_StreamComm, &m_StreamSize);

    int writerMasterStreamRank = -1;
    if (m_WriterRank == 0)
    {
        writerMasterStreamRank = m_StreamRank;
    }
    MPI_Allreduce(&writerMasterStreamRank, &m_WriterMasterStreamRank, 1, MPI_INT, MPI_MAX,
                  m_StreamComm);

    int readerMasterStreamRank = -1;
    MPI_Allreduce(&readerMasterStreamRank, &m_ReaderMasterStreamRank, 1, MPI_INT, MPI_MAX,
                  m_StreamComm);
}

}
}
}
}
