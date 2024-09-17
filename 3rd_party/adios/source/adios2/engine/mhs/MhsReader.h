/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * MhsReader.h
 *
 *  Created on: Aug 04, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#ifndef ADIOS2_ENGINE_MHSREADER_H_
#define ADIOS2_ENGINE_MHSREADER_H_

#include "adios2/core/Engine.h"
#include "adios2/operator/compress/CompressSirius.h"

namespace adios2
{
namespace core
{
namespace engine
{

class MhsReader : public Engine
{
public:
    MhsReader(IO &adios, const std::string &name, const Mode mode, helper::Comm comm);
    virtual ~MhsReader();

    StepStatus BeginStep(StepMode mode = StepMode::Read, const float timeoutSeconds = -1.0) final;
    size_t CurrentStep() const final;
    void PerformGets() final;
    void EndStep() final;

private:
    std::vector<IO *> m_SubIOs;
    std::vector<Engine *> m_SubEngines;
    std::shared_ptr<compress::CompressSirius> m_SiriusCompressor;
    int m_Tiers;

#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoClose(const int transportIndex = -1);

    /**
     * Called if destructor is called on an open engine.  Should warn or take
     * any non-complex measure that might help recover.
     */
    void DestructorClose(bool Verbose) noexcept final{};

    template <class T>
    void GetSyncCommon(Variable<T> &variable, T *data);

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_MHSREADER_H_ */
