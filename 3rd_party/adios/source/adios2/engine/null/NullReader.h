/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * NullReader.h
 *
 *  Created on: 23 March 2022
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#ifndef ADIOS2_ENGINE_NULL2_NULLREADER_H_
#define ADIOS2_ENGINE_NULL2_NULLREADER_H_

#include <memory>

#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosComm.h"

namespace adios2
{
namespace core
{
namespace engine
{

class NullReader : public core::Engine
{

public:
    NullReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    virtual ~NullReader();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) override;
    size_t CurrentStep() const override;
    void EndStep() override;
    void PerformGets() override;

protected:
#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) override;                                                   \
    void DoGetDeferred(Variable<T> &, T *) override;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoClose(const int transportIndex) override;

    template <class T>
    void GetSyncCommon(Variable<T> &variable, T *data);

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);

protected:
    void DestructorClose(bool Verbose) noexcept final;

private:
    struct NullReaderImpl;
    std::unique_ptr<NullReaderImpl> Impl;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_NULL2_NULLREADER_H_ */
