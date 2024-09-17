/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosPluginManager.cpp
 *
 * Created on: Dec 14, 2021
 *     Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#include "adios2/helper/adiosPluginManager.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/IO.h"
#include "adios2/engine/plugin/PluginEngineInterface.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosDynamicBinder.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"

#include <adios2sys/SystemTools.hxx>

#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace adios2
{
namespace plugin
{

namespace
{
const std::string pluginEnvVarName = "ADIOS2_PLUGIN_PATH";

struct EnginePluginInfo
{
    std::string m_LibraryName;
    std::unique_ptr<helper::DynamicBinder> m_Binder;
    PluginManager::EngineCreateFun m_HandleCreate;
    PluginManager::EngineDestroyFun m_HandleDestroy;
};

struct OperatorPluginInfo
{
    std::string m_LibraryName;
    std::unique_ptr<helper::DynamicBinder> m_Binder;
    PluginManager::OperatorCreateFun m_HandleCreate;
    PluginManager::OperatorDestroyFun m_HandleDestroy;
};

} // end anon namespace

struct PluginManager::Impl
{
    std::unordered_map<std::string, EnginePluginInfo> m_EngineRegistry;
    std::unordered_map<std::string, OperatorPluginInfo> m_OperatorRegistry;
    adios2::Params m_Parameters;
    int m_Verbosity = 0;
};

PluginManager *PluginManager::m_Instance = nullptr;
bool PluginManager::m_Destroyed = false;

PluginManager::PluginManager() : m_Impl(new Impl) {}

PluginManager::~PluginManager()
{
    m_Instance = nullptr;
    m_Destroyed = true;
}

PluginManager &PluginManager::GetInstance()
{
    if (!m_Instance)
    {
        if (m_Destroyed)
        {
            throw std::runtime_error("Dead reference to PluginManager singleton");
        }
        else
        {
            CreateInstance();
        }
    }
    return *m_Instance;
}

void PluginManager::CreateInstance()
{
    static PluginManager theInstance;
    m_Instance = &theInstance;
}

void PluginManager::SetParameters(const Params &params)
{
    helper::GetParameter(params, "verbose", m_Impl->m_Verbosity);
}

bool PluginManager::LoadPlugin(const std::string &pluginName, const std::string &pluginLibrary)
{
    if (m_Impl->m_EngineRegistry.find(pluginName) != m_Impl->m_EngineRegistry.end() ||
        m_Impl->m_OperatorRegistry.find(pluginName) != m_Impl->m_OperatorRegistry.end())
    {
        return true;
    }

    std::string allPluginPaths;
    adios2sys::SystemTools::GetEnv(pluginEnvVarName, allPluginPaths);
    if (allPluginPaths.empty())
    {
        return OpenPlugin(pluginName, pluginLibrary, "");
    }

#ifdef _WIN32
    char platform_separator = ';';
#else
    char platform_separator = ':';
#endif

    auto pathsSplit =
        adios2sys::SystemTools::SplitString(allPluginPaths, platform_separator, false);

    bool loaded = false;
    auto pathIt = pathsSplit.begin();
    while (pathIt != pathsSplit.end() && !loaded)
    {
        try
        {
            loaded = OpenPlugin(pluginName, pluginLibrary, *pathIt);
        }
        catch (std::exception &e)
        {
            // this is not necessarily an error, because you could have
            // multiple paths in ADIOS2_PLUGIN_PATH variable
            helper::Log("Plugins", "PluginManager", "LoadPlugin",
                        std::string("OpenPlugin failed: ") + e.what(), 5, m_Impl->m_Verbosity,
                        helper::LogMode::INFO);
            loaded = false;
        }
        ++pathIt;
    }

    if (!loaded)
    {
        helper::Log("Plugins", "PluginManager", "LoadPlugin",
                    "The plugin " + pluginLibrary +
                        " could not be loaded."
                        " Double check ADIOS2_PLUGIN_PATH or plugin library "
                        "name is correct."
                        "\nADIOS2_PLUGIN_PATH: " +
                        allPluginPaths,
                    helper::LogMode::WARNING);
    }

    return loaded;
}

bool PluginManager::OpenPlugin(const std::string &pluginName, const std::string &pluginLibrary,
                               const std::string &pluginPath)
{
    helper::Log("Plugins", "PluginManager", "OpenPlugin",
                "Attempting to open plugin " + pluginLibrary + " at path " + pluginPath, 5,
                m_Impl->m_Verbosity, helper::LogMode::INFO);
    std::unique_ptr<helper::DynamicBinder> binder(
        new helper::DynamicBinder(pluginLibrary, pluginPath));
    if (auto createHandle = binder->GetSymbol("EngineCreate"))
    {
        // we have an engine plugin
        EnginePluginInfo plugin;
        plugin.m_LibraryName = pluginLibrary;
        plugin.m_HandleCreate = reinterpret_cast<EngineCreatePtr>(createHandle);
        if (!plugin.m_HandleCreate)
        {
            helper::Throw<std::runtime_error>("Plugins", "PluginManager", "OpenPlugin",
                                              "Unable to locate EngineCreate"
                                              " symbol in library " +
                                                  pluginLibrary);
        }

        plugin.m_HandleDestroy =
            reinterpret_cast<EngineDestroyPtr>(binder->GetSymbol("EngineDestroy"));
        if (!plugin.m_HandleDestroy)
        {
            helper::Throw<std::runtime_error>("Plugins", "PluginManager", "OpenPlugin",
                                              "Unable to locate EngineDestroy"
                                              " symbol in library " +
                                                  pluginLibrary);
        }
        plugin.m_Binder = std::move(binder);
        m_Impl->m_EngineRegistry[pluginName] = std::move(plugin);
        helper::Log("Plugins", "PluginManager", "OpenPlugin",
                    "Engine Plugin " + pluginName + " successfully opened", 5, m_Impl->m_Verbosity,
                    helper::LogMode::INFO);
        return true;
    }
    else if (auto createHandle = binder->GetSymbol("OperatorCreate"))
    {
        // should be an operator plugin
        OperatorPluginInfo plugin;
        plugin.m_LibraryName = pluginLibrary;
        plugin.m_HandleCreate = reinterpret_cast<OperatorCreatePtr>(createHandle);
        if (!plugin.m_HandleCreate)
        {
            helper::Throw<std::runtime_error>("Plugins", "PluginManager", "OpenPlugin",
                                              "Unable to locate OperatorCreate"
                                              " symbol in library " +
                                                  pluginLibrary);
        }

        plugin.m_HandleDestroy =
            reinterpret_cast<OperatorDestroyPtr>(binder->GetSymbol("OperatorDestroy"));
        if (!plugin.m_HandleDestroy)
        {
            helper::Throw<std::runtime_error>("Plugins", "PluginManager", "OpenPlugin",
                                              "Unable to locate OperatorDestroy"
                                              " symbol in library " +
                                                  pluginLibrary);
        }
        plugin.m_Binder = std::move(binder);
        m_Impl->m_OperatorRegistry[pluginName] = std::move(plugin);
        helper::Log("Plugins", "PluginManager", "OpenPlugin",
                    "Operator Plugin " + pluginName + " successfully opened", 5,
                    m_Impl->m_Verbosity, helper::LogMode::INFO);
        return true;
    }
    helper::Throw<std::runtime_error>("Plugins", "PluginManager", "OpenPlugin",
                                      "Unable to locate Create/Destroy symbols in library " +
                                          pluginLibrary);
    return false;
}

PluginManager::EngineCreateFun PluginManager::GetEngineCreateFun(const std::string &name)
{
    auto pluginIt = m_Impl->m_EngineRegistry.find(name);
    if (pluginIt == m_Impl->m_EngineRegistry.end())
    {
        helper::Throw<std::runtime_error>("Plugins", "PluginManager", "GetEngineCreateFun",
                                          "Couldn't find engine plugin named " + name);
    }

    return pluginIt->second.m_HandleCreate;
}

PluginManager::EngineDestroyFun PluginManager::GetEngineDestroyFun(const std::string &name)
{
    auto pluginIt = m_Impl->m_EngineRegistry.find(name);
    if (pluginIt == m_Impl->m_EngineRegistry.end())
    {
        helper::Throw<std::runtime_error>("Plugins", "PluginManager", "GetEngineDestroyFun",
                                          "Couldn't find engine plugin named " + name);
    }

    return pluginIt->second.m_HandleDestroy;
}

PluginManager::OperatorCreateFun PluginManager::GetOperatorCreateFun(const std::string &name)
{
    auto pluginIt = m_Impl->m_OperatorRegistry.find(name);
    if (pluginIt == m_Impl->m_OperatorRegistry.end())
    {
        helper::Throw<std::runtime_error>("Plugins", "PluginManager", "GetOperatorCreateFun",
                                          "Couldn't find operator plugin named " + name);
    }

    return pluginIt->second.m_HandleCreate;
}

PluginManager::OperatorDestroyFun PluginManager::GetOperatorDestroyFun(const std::string &name)
{
    auto pluginIt = m_Impl->m_OperatorRegistry.find(name);
    if (pluginIt == m_Impl->m_OperatorRegistry.end())
    {
        helper::Throw<std::runtime_error>("Plugins", "PluginManager", "GetOperatorDestroyFun",
                                          "Couldn't find operator plugin named " + name);
    }

    return pluginIt->second.m_HandleDestroy;
}

} // end namespace plugin
} // end namespace adios2
