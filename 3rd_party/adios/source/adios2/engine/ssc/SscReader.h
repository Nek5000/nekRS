/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReader.h
 *
 *  Created on: Nov 1, 2018
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCREADER_H_
#define ADIOS2_ENGINE_SSCREADER_H_

#include "SscReaderBase.h"
#include "adios2/core/Engine.h"

namespace adios2
{
namespace core
{
namespace engine
{

class SscReader : public Engine
{
public:
    SscReader(IO &adios, const std::string &name, const Mode mode, helper::Comm comm);
    ~SscReader();
    StepStatus BeginStep(StepMode stepMode = StepMode::Read,
                         const float timeoutSeconds = std::numeric_limits<float>::max()) final;
    void PerformGets() final;
    size_t CurrentStep() const final;
    void EndStep() final;

private:
#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;                                                  \
    std::vector<typename Variable<T>::BPInfo> DoBlocksInfo(const Variable<T> &variable,            \
                                                           const size_t step) const final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoGetStructSync(VariableStruct &, void *) final;
    void DoGetStructDeferred(VariableStruct &, void *) final;
    std::vector<VariableStruct::BPInfo> DoBlocksInfoStruct(const VariableStruct &variable,
                                                           const size_t step) const final;

    void DoClose(const int transportIndex = -1) final;

    /**
     * Called if destructor is called on an open engine.  Should warn or take
     * any non-complex measure that might help recover.
     */
    void DestructorClose(bool Verbose) noexcept final{};

    int m_Verbosity = 0;
    std::string m_EngineMode = "generic";
    std::shared_ptr<ssc::SscReaderBase> m_EngineInstance;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif // ADIOS2_ENGINE_SSCREADER_H_
