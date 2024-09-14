/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderGeneric.h
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCREADERGENERIC_H_
#define ADIOS2_ENGINE_SSCREADERGENERIC_H_

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

class SscReaderGeneric : public SscReaderBase
{

public:
    SscReaderGeneric(IO &io, const std::string &name, const Mode mode, MPI_Comm comm);
    ~SscReaderGeneric() = default;

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
    bool m_StepBegun = false;
    bool m_WriterDefinitionsLocked = false;
    bool m_ReaderSelectionsLocked = false;
    std::thread m_EndStepThread;
    StepStatus m_StepStatus;
    std::vector<MPI_Request> m_MpiRequests;
    ssc::RankPosMap m_AllReceivingWriterRanks;
    ssc::BlockVecVec m_GlobalWritePattern;
    ssc::BlockVec m_LocalReadPattern;
    ssc::Buffer m_GlobalWritePatternBuffer;
    MPI_Win m_MpiWin;

    bool SyncWritePattern();
    void SyncReadPattern();
    void BeginStepConsequentFixed();
    void BeginStepFlexible(StepStatus &status);
    void EndStepFixed();
    void EndStepFirstFlexible();
    void EndStepConsequentFlexible();
    void CalculatePosition(ssc::BlockVecVec &mapVec, ssc::RankPosMap &allOverlapRanks);

    template <typename T>
    std::vector<typename Variable<T>::BPInfo> BlocksInfoCommon(const Variable<T> &variable,
                                                               const size_t step) const;

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);

    void GetDeferredDeltaCommon(VariableBase &variable, void *data);
};

}
}
}
}

#endif // ADIOS2_ENGINE_SSCREADERGENERIC_H_
