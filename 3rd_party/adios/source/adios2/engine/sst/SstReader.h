/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SstReader.h
 *
 *  Created on: Aug 17, 2017
 *      Author: Greg Eisenhauer
 */

#ifndef ADIOS2_ENGINE_SST_SSTREADER_H_
#define ADIOS2_ENGINE_SST_SSTREADER_H_

#include "adios2/toolkit/sst/sst.h"

#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/format/bp/bp3/BP3Deserializer.h"
#include "adios2/toolkit/format/bp5/BP5Deserializer.h"

namespace adios2
{
namespace core
{
namespace engine
{

class SstReader : public Engine
{

public:
    /**
     * Constructor for sst engine Reader
     * @param adios
     * @param name
     * @param accessMode
     * @param comm
     * @param method
     * @param nthreads
     */
    SstReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    virtual ~SstReader();

    StepStatus BeginStep();
    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0);
    size_t CurrentStep() const final;
    void EndStep();
    void PerformGets();
    void Flush(const int transportIndex = -1) final;
    MinVarInfo *MinBlocksInfo(const VariableBase &, const size_t Step) const;
    bool VarShape(const VariableBase &Var, const size_t Step, Dims &Shape) const;
    bool VariableMinMax(const VariableBase &, const size_t Step, MinMaxStruct &MinMax);

private:
    template <class T>
    void ReadVariableBlocksRequests(Variable<T> &variable, std::vector<void *> &sstReadHandlers,
                                    std::vector<std::vector<char>> &buffers);

    template <class T>
    void ReadVariableBlocksFill(Variable<T> &variable, std::vector<std::vector<char>> &buffers,
                                size_t &iter);

    template <class T>
    void SstBPPerformGets();
    void BP5PerformGets();
    void Init();
    SstStream m_Input;
    SstMarshalMethod m_WriterMarshalMethod;
    int m_WriterIsRowMajor;
    bool m_DefinitionsNotified = false;

    /* --- Used only with BP marshaling --- */
    SstFullMetadata m_CurrentStepMetaData = NULL;
    format::BP3Deserializer *m_BP3Deserializer;
    format::BP5Deserializer *m_BP5Deserializer = nullptr;
    /* --- Used only with BP marshaling --- */

    struct _SstParams Params;

    std::unordered_map<const VariableBase *, std::unique_ptr<MinVarInfo>> m_InfoMap;
#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;                                                  \
                                                                                                   \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> DoAllStepsBlocksInfo(              \
        const Variable<T> &variable) const final;                                                  \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> DoBlocksInfo(const Variable<T> &variable,            \
                                                           const size_t step) const final;

    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoGetStructSync(VariableStruct &, void *);
    void DoGetStructDeferred(VariableStruct &, void *);

    /**
     * Called if destructor is called on an open engine.  Should warn or take
     * any non-complex measure that might help recover.
     */
    void DestructorClose(bool Verbose) noexcept final;

    void DoClose(const int transportIndex = -1) final;

    template <class T>
    void GetSyncCommon(Variable<T> &variable, T *data);

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_SST_SSTREADER_H_ */
