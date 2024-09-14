/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderNaive.h
 *
 *  Created on: Mar 7, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCREADERNAIVE_H_
#define ADIOS2_ENGINE_SSCREADERNAIVE_H_

#include "SscReaderBase.h"
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

class SscReaderNaive : public SscReaderBase
{

public:
    SscReaderNaive(IO &io, const std::string &name, const Mode mode, MPI_Comm comm);
    ~SscReaderNaive() = default;

    StepStatus BeginStep(const StepMode mode, const float timeoutSeconds,
                         const bool readerLocked) final;
    size_t CurrentStep() final;
    void PerformGets() final;
    void EndStep(const bool readerLocked) final;
    void Close(const int transportIndex) final;

#define declare_type(T)                                                                            \
    std::vector<typename Variable<T>::BPInfo> BlocksInfo(const Variable<T> &variable,              \
                                                         const size_t step) const final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    std::vector<VariableStruct::BPInfo> BlocksInfo(const VariableStruct &variable,
                                                   const size_t step) const final;

    void GetDeferred(VariableBase &, void *) final;

private:
    template <typename T>
    std::vector<typename Variable<T>::BPInfo> BlocksInfoCommon(const Variable<T> &variable,
                                                               const size_t step) const;

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);

    std::unordered_map<std::string, ssc::BlockVec> m_BlockMap;
};

}
}
}
}

#endif // ADIOS2_ENGINE_SSCREADERNAIVE_H_
