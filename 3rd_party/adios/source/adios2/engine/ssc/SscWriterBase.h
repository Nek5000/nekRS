/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscWriterBase.h
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCWRITERBASE_H_
#define ADIOS2_ENGINE_SSCWRITERBASE_H_

#include "SscHelper.h"
#include "adios2/core/IO.h"
#include <mpi.h>

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

class SscWriterBase
{

public:
    SscWriterBase(IO &io, const std::string &name, const Mode mode, MPI_Comm comm);
    virtual ~SscWriterBase();

    virtual StepStatus BeginStep(const StepMode mode, const float timeoutSeconds,
                                 const bool writerLocked) = 0;
    virtual size_t CurrentStep() = 0;
    virtual void PerformPuts() = 0;
    virtual void EndStep(const bool writerLocked) = 0;
    virtual void Close(const int transportIndex) = 0;

    virtual void PutDeferred(VariableBase &, const void *) = 0;

protected:
    void SyncMpiPattern(MPI_Comm comm);

    MPI_Group m_ReaderGroup;
    MPI_Comm m_StreamComm;
    MPI_Comm m_WriterComm;

    int m_StreamRank;
    int m_StreamSize;
    int m_WriterRank;
    int m_WriterSize;
    int m_WriterMasterStreamRank;
    int m_ReaderMasterStreamRank;

    ssc::Buffer m_Buffer;

    int64_t m_CurrentStep = -1;
    std::string m_Name;
    int m_Verbosity = 0;
    int m_OpenTimeoutSecs = 10;
    bool m_Threading = false;

    IO &m_IO;
};

}
}
}
}

#endif // ADIOS2_ENGINE_SSCWRITERBASE_H_
