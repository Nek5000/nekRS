/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscWriterGeneric.h
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCWRITERGENERIC_H_
#define ADIOS2_ENGINE_SSCWRITERGENERIC_H_

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

class SscWriterGeneric : public SscWriterBase
{

public:
    SscWriterGeneric(IO &io, const std::string &name, const Mode mode, MPI_Comm comm);
    ~SscWriterGeneric() = default;

    StepStatus BeginStep(const StepMode mode, const float timeoutSeconds,
                         const bool writerLocked) final;
    size_t CurrentStep() final;
    void PerformPuts() final;
    void EndStep(const bool writerLocked) final;
    void Close(const int transportIndex) final;

    void PutDeferred(VariableBase &, const void *) final;

    std::unordered_map<std::string, StructDefinition> m_StructDefinitions;

private:
    MPI_Win m_MpiWin;
    std::thread m_EndStepThread;
    ssc::BlockVecVec m_GlobalWritePattern;
    ssc::BlockVecVec m_GlobalReadPattern;
    std::vector<MPI_Request> m_MpiRequests;
    ssc::RankPosMap m_AllSendingReaderRanks;

    bool m_WriterDefinitionsLocked = false;
    bool m_ReaderSelectionsLocked = false;

    template <class T>
    void PutDeferredCommon(Variable<T> &variable, const T *values);
    void SyncWritePattern(bool finalStep = false);
    void SyncReadPattern();
    void EndStepFirst();
    void EndStepConsequentFixed();
    void EndStepConsequentFlexible();
    void CalculatePosition(ssc::BlockVecVec &writerMapVec, ssc::BlockVecVec &readerMapVec,
                           const int writerRank, ssc::RankPosMap &allOverlapRanks);
};

}
}
}
}

#endif // ADIOS2_ENGINE_SSCWRITERGENERIC_H_
