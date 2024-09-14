/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosPluginManager.h
 *
 * Created on: Dec 14, 2021
 *     Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#ifndef ADIOS2_HELPER_PLUGINMANAGER_H
#define ADIOS2_HELPER_PLUGINMANAGER_H

#include "adios2/engine/plugin/PluginEngineInterface.h"
#include "adios2/operator/plugin/PluginOperatorInterface.h"

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

namespace adios2
{
namespace plugin
{

class PluginManager
{
public:
    using EngineCreatePtr = std::add_pointer<PluginEngineInterface *(
        core::IO &, const std::string &, const Mode, helper::Comm)>::type;
    using EngineDestroyPtr = std::add_pointer<void(PluginEngineInterface *)>::type;
    using EngineCreateFun = std::function<typename std::remove_pointer<EngineCreatePtr>::type>;
    using EngineDestroyFun = std::function<typename std::remove_pointer<EngineDestroyPtr>::type>;

    using OperatorCreatePtr = std::add_pointer<PluginOperatorInterface *(const Params &)>::type;
    using OperatorDestroyPtr = std::add_pointer<void(PluginOperatorInterface *)>::type;
    using OperatorCreateFun = std::function<std::remove_pointer<OperatorCreatePtr>::type>;
    using OperatorDestroyFun = std::function<std::remove_pointer<OperatorDestroyPtr>::type>;

    static PluginManager &GetInstance();

    void SetParameters(const Params &params);

    /**
     * Attempts to load a single plugin specified by pluginName and
     * pluginLibrary.
     */
    bool LoadPlugin(const std::string &pluginName, const std::string &pluginLibrary);

    EngineCreateFun GetEngineCreateFun(const std::string &name);
    EngineDestroyFun GetEngineDestroyFun(const std::string &name);

    OperatorCreateFun GetOperatorCreateFun(const std::string &name);
    OperatorDestroyFun GetOperatorDestroyFun(const std::string &name);

private:
    PluginManager();
    PluginManager(const PluginManager &) = delete;
    PluginManager &operator=(const PluginManager &) = delete;
    virtual ~PluginManager();

    static void CreateInstance();

    bool OpenPlugin(const std::string &pluginName, const std::string &pluginLibrary,
                    const std::string &pluginPath);

    static PluginManager *m_Instance;
    static bool m_Destroyed;

    struct Impl;
    std::unique_ptr<Impl> m_Impl;
};

} // end namespace plugin
} // end namespace adios2

#endif /* ADIOS2_HELPER_PLUGINMANAGER_H */
