/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * PluginEngine.h Support for an engine implemented outside libadios2
 *
 *  Created on: July 5, 2021
 *      Author: Chuck Atkins <chuck.atkins@kitware.com>
 *              Caitlin Ross <caitlin.ross@kitware.com>
 */

#ifndef ADIOS2_ENGINE_PLUGIN_PLUGINENGINE_H_
#define ADIOS2_ENGINE_PLUGIN_PLUGINENGINE_H_

#include "PluginEngineInterface.h"

#include <memory> // for unique_ptr
#include <string> // for string

#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/core/Variable.h"
#include "adios2/helper/adiosComm.h"

namespace adios2
{
namespace plugin
{

/** A front-end wrapper for an engine implemented outside of libadios2 */
class PluginEngine : public core::Engine
{
public:
    PluginEngine(core::IO &io, const std::string &name, const Mode mode, helper::Comm comm);
    virtual ~PluginEngine();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = 0.f) override;
    void PerformPuts() override;
    void PerformGets() override;
    void EndStep() override;

protected:
#define declare(T)                                                                                 \
    void DoPutSync(core::Variable<T> &, const T *) override;                                       \
    void DoPutDeferred(core::Variable<T> &, const T *) override;                                   \
    void DoGetSync(core::Variable<T> &, T *) override;                                             \
    void DoGetDeferred(core::Variable<T> &, T *) override;

    ADIOS2_FOREACH_STDTYPE_1ARG(declare)
#undef declare

    void DoClose(const int transportIndex = -1) override;

    // probably should go to plugin if this is to be used for non-trivial
    // engines
    void DestructorClose(bool Verbose) noexcept final{};

private:
    struct Impl;
    std::unique_ptr<Impl> m_Impl;
};

} // end namespace plugin
} // end namespace adios2

#endif /* ADIOS2_ENGINE_PLUGIN_PLUGINENGINE_H_ */
