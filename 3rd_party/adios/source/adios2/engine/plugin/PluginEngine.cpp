/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * PluginEngine.cpp
 *
 *  Created on: July 5, 2021
 *      Author: Chuck Atkins <chuck.atkins@kitware.com>
 *              Caitlin Ross <caitlin.ross@kitware.com>
 */

#include "PluginEngine.h"
#include "PluginEngineInterface.h"

#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>

#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosPluginManager.h"

namespace adios2
{
namespace plugin
{

/******************************************************************************/

struct PluginEngine::Impl
{
    PluginManager::EngineCreateFun m_HandleCreate;
    PluginManager::EngineDestroyFun m_HandleDestroy;
    PluginEngineInterface *m_Plugin = nullptr;
};

/******************************************************************************/

PluginEngine::PluginEngine(core::IO &io, const std::string &name, const Mode mode,
                           helper::Comm comm)
: Engine("Plugin", io, name, mode, comm.Duplicate()), m_Impl(new Impl)
{
    auto pluginNameIt = m_IO.m_Parameters.find("PluginName");
    if (pluginNameIt == m_IO.m_Parameters.end())
    {
        helper::Throw<std::runtime_error>("Plugins", "PluginEngine", "PluginEngine",
                                          "PluginName must be specified in the engine parameters");
    }

    auto pluginLibIt = m_IO.m_Parameters.find("PluginLibrary");
    if (pluginLibIt == m_IO.m_Parameters.end())
    {
        helper::Throw<std::runtime_error>(
            "Plugins", "PluginEngine", "PluginEngine",
            "PluginLibrary must be specified in the engine parameters");
    }

    auto &pluginManager = PluginManager::GetInstance();
    pluginManager.SetParameters(m_IO.m_Parameters);
    pluginManager.LoadPlugin(pluginNameIt->second, pluginLibIt->second);
    m_Impl->m_HandleCreate = pluginManager.GetEngineCreateFun(pluginNameIt->second);
    m_Impl->m_HandleDestroy = pluginManager.GetEngineDestroyFun(pluginNameIt->second);
    m_Impl->m_Plugin = m_Impl->m_HandleCreate(io, pluginNameIt->second, mode, comm.Duplicate());
    m_IsOpen = true;
}

PluginEngine::~PluginEngine() { m_Impl->m_HandleDestroy(m_Impl->m_Plugin); }

StepStatus PluginEngine::BeginStep(StepMode mode, const float timeoutSeconds)
{
    return m_Impl->m_Plugin->BeginStep(mode, timeoutSeconds);
}

void PluginEngine::PerformPuts() { m_Impl->m_Plugin->PerformPuts(); }

void PluginEngine::PerformGets() { m_Impl->m_Plugin->PerformGets(); }

void PluginEngine::EndStep() { m_Impl->m_Plugin->EndStep(); }

#define declare(T)                                                                                 \
    void PluginEngine::DoPutSync(core::Variable<T> &variable, const T *values)                     \
    {                                                                                              \
        m_Impl->m_Plugin->DoPutSync(variable, values);                                             \
    }                                                                                              \
    void PluginEngine::DoPutDeferred(core::Variable<T> &variable, const T *values)                 \
    {                                                                                              \
        m_Impl->m_Plugin->DoPutDeferred(variable, values);                                         \
    }                                                                                              \
    void PluginEngine::DoGetSync(core::Variable<T> &variable, T *values)                           \
    {                                                                                              \
        m_Impl->m_Plugin->DoGetSync(variable, values);                                             \
    }                                                                                              \
    void PluginEngine::DoGetDeferred(core::Variable<T> &variable, T *values)                       \
    {                                                                                              \
        m_Impl->m_Plugin->DoGetDeferred(variable, values);                                         \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare)
#undef declare

void PluginEngine::DoClose(const int transportIndex) { m_Impl->m_Plugin->Close(transportIndex); }

} // end namespace plugin
} // end namespace adios2
