/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscWriterNaive.h
 *
 *  Created on: Mar 7, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCWRITERNAIVE_H_
#define ADIOS2_ENGINE_SSCWRITERNAIVE_H_

#include "SscWriterBase.h"
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

class SscWriterNaive : public SscWriterBase
{

public:
    SscWriterNaive(IO &io, const std::string &name, const Mode mode, MPI_Comm comm);
    ~SscWriterNaive() = default;

    StepStatus BeginStep(const StepMode mode, const float timeoutSeconds,
                         const bool writerLocked) final;
    size_t CurrentStep() final;
    void PerformPuts() final;
    void EndStep(const bool writerLocked) final;
    void Close(const int transportIndex) final;

    void PutDeferred(VariableBase &, const void *) final;

private:
    template <class T>
    void PutDeferredCommon(Variable<T> &variable, const T *values);

    ssc::BlockVec m_Metadata;
};

}
}
}
}

#endif // ADIOS2_ENGINE_SSCWRITENAIVE_H_
