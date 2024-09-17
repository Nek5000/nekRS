/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * NullWriter.h
 *
 *  Created on: 16 Apr 19
 *      Author: Chuck Atkins chuck.atkins@kitware.com
 */

#ifndef ADIOS2_ENGINE_NULL2_NULLWRITER_H_
#define ADIOS2_ENGINE_NULL2_NULLWRITER_H_

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

class NullWriter : public core::Engine
{

public:
    NullWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    virtual ~NullWriter();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) override;
    size_t CurrentStep() const override;
    void EndStep() override;
    void PerformPuts() override;
    void Flush(const int transportIndex = -1) override;

protected:
    // Put
#define declare_type(T)                                                                            \
    void DoPut(Variable<T> &variable, typename Variable<T>::Span &span, const bool initialize,     \
               const T &value) override;

    ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

#define declare_type(T)                                                                            \
    void DoPutSync(Variable<T> &, const T *) override;                                             \
    void DoPutDeferred(Variable<T> &, const T *) override;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoClose(const int transportIndex) override;

protected:
    void DestructorClose(bool Verbose) noexcept final{};

private:
    struct NullWriterImpl;
    std::unique_ptr<NullWriterImpl> Impl;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_NULL2_NULLWRITER_H_ */
