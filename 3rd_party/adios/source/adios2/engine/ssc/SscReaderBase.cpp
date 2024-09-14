/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderBase.cpp
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#include "SscReaderBase.h"
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

SscReaderBase::SscReaderBase(IO &io, const std::string &name, const Mode mode, MPI_Comm comm)
: m_Name(name), m_IO(io)
{
    helper::GetParameter(io.m_Parameters, "Verbose", m_Verbosity);
    helper::GetParameter(io.m_Parameters, "Threading", m_Threading);
    helper::GetParameter(io.m_Parameters, "OpenTimeoutSecs", m_OpenTimeoutSecs);

    SyncMpiPattern(comm);
}

SscReaderBase::~SscReaderBase() = default;

void SscReaderBase::SyncMpiPattern(MPI_Comm comm)
{

    MPI_Group streamGroup;
    MPI_Group readerGroup;
    MPI_Comm writerComm;

    helper::HandshakeComm(m_Name, 'r', m_OpenTimeoutSecs, comm, streamGroup, m_WriterGroup,
                          readerGroup, m_StreamComm, writerComm, m_ReaderComm);

    MPI_Comm_rank(comm, &m_ReaderRank);
    MPI_Comm_size(comm, &m_ReaderSize);
    MPI_Comm_rank(m_StreamComm, &m_StreamRank);
    MPI_Comm_size(m_StreamComm, &m_StreamSize);

    int writerMasterStreamRank = -1;
    MPI_Allreduce(&writerMasterStreamRank, &m_WriterMasterStreamRank, 1, MPI_INT, MPI_MAX,
                  m_StreamComm);

    int readerMasterStreamRank = -1;
    if (m_ReaderRank == 0)
    {
        readerMasterStreamRank = m_StreamRank;
    }
    MPI_Allreduce(&readerMasterStreamRank, &m_ReaderMasterStreamRank, 1, MPI_INT, MPI_MAX,
                  m_StreamComm);
}

}
}
}
}
