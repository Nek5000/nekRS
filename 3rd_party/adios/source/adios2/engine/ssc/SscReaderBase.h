/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderBase.h
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCREADERBASE_H_
#define ADIOS2_ENGINE_SSCREADERBASE_H_

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

class SscReaderBase
{

public:
    SscReaderBase(IO &io, const std::string &name, const Mode mode, MPI_Comm comm);
    virtual ~SscReaderBase();

    virtual StepStatus BeginStep(const StepMode mode, const float timeoutSeconds,
                                 const bool readerLocked) = 0;
    virtual size_t CurrentStep() = 0;
    virtual void PerformGets() = 0;
    virtual void EndStep(const bool readerLocked) = 0;
    virtual void Close(const int transportIndex) = 0;

#define declare_type(T)                                                                            \
    virtual std::vector<typename Variable<T>::BPInfo> BlocksInfo(const Variable<T> &variable,      \
                                                                 const size_t step) const = 0;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    virtual std::vector<VariableStruct::BPInfo> BlocksInfo(const VariableStruct &variable,
                                                           const size_t step) const = 0;

    virtual void GetDeferred(VariableBase &, void *) = 0;

protected:
    void SyncMpiPattern(MPI_Comm comm);

    MPI_Group m_WriterGroup;
    MPI_Comm m_StreamComm;
    MPI_Comm m_ReaderComm;

    int m_StreamRank;
    int m_StreamSize;
    int m_ReaderRank;
    int m_ReaderSize;
    int m_WriterMasterStreamRank;
    int m_ReaderMasterStreamRank;

    std::unordered_map<std::string, StructDefinition> m_StructDefinitions;

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

#endif // ADIOS2_ENGINE_SSCREADERBASE_H_
