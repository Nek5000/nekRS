/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * PluginOperator.cpp
 *
 *  Created on: Dec 7, 2021
 *      Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#include "PluginOperator.h"
#include "PluginOperatorInterface.h"

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>

#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosPluginManager.h"
#include "adios2/helper/adiosString.h"

#include <adios2sys/SystemTools.hxx>

namespace adios2
{
namespace plugin
{

/******************************************************************************/

struct PluginOperator::Impl
{
    Params m_PluginParams;
    PluginManager::OperatorCreateFun m_HandleCreate;
    PluginManager::OperatorDestroyFun m_HandleDestroy;
    PluginOperatorInterface *m_Plugin = nullptr;
    int m_Verbosity = 0;
};

/******************************************************************************/

PluginOperator::PluginOperator(const Params &parameters)
: Operator("plugin", PLUGIN_INTERFACE, "plugin", parameters), m_Impl(new Impl)
{
    helper::GetParameter(m_Parameters, "verbose", m_Impl->m_Verbosity);

    // It is possible for 'pluginname' and 'pluginlibrary' to not be in
    // m_Parameters (e.g., when using bpls, which will call InverseOperate).
    // PluginName and PluginLibrary which will be saved in header info in the
    // Operate() call, so that in this situation, the plugin can be loaded
    // to perform InverseOperate.
    auto pluginNameIt = m_Parameters.find("pluginname");
    auto pluginLibIt = m_Parameters.find("pluginlibrary");
    if (pluginNameIt != m_Parameters.end() && pluginLibIt != m_Parameters.end())
    {
        m_Impl->m_PluginParams["PluginName"] = pluginNameIt->second;
        m_Impl->m_PluginParams["PluginLibrary"] = pluginLibIt->second;
        PluginInit(pluginNameIt->second, pluginLibIt->second);
    }
}

PluginOperator::~PluginOperator() { m_Impl->m_HandleDestroy(m_Impl->m_Plugin); }

void PluginOperator::PluginInit(const std::string &pluginName, const std::string &pluginLibrary)
{
    if (m_Impl->m_Plugin)
    {
        return;
    }

    auto &pluginManager = PluginManager::GetInstance();
    pluginManager.SetParameters(m_Parameters);
    pluginManager.LoadPlugin(pluginName, pluginLibrary);

    m_Impl->m_HandleCreate = pluginManager.GetOperatorCreateFun(pluginName);
    m_Impl->m_HandleDestroy = pluginManager.GetOperatorDestroyFun(pluginName);
    m_Impl->m_Plugin = m_Impl->m_HandleCreate(m_Parameters);
}

size_t PluginOperator::GetEstimatedSize(const size_t ElemCount, const size_t ElemSize,
                                        const size_t ndims, const size_t *dims) const
{
    // Need to calculate the size of the header written by Operate and then add our plugin's size.
    constexpr size_t commonHeaderSize =
        sizeof(m_TypeEnum) + sizeof(std::uint8_t) + sizeof(std::uint16_t);

    auto &pp = m_Impl->m_PluginParams;
    // Want to use std::transform_reduce but C++11
    size_t paramsSize = 1; // for the number of parameters
    for (auto &&p : pp)
    {
        // Need length and string for key and values.
        paramsSize += p.first.size() + p.second.size() + 2;
    }

    // Plugin's estimate of size so it doesn't need to know about headers.
    auto implSize = m_Impl->m_Plugin->GetEstimatedSize(ElemCount, ElemSize, ndims, dims);

    return commonHeaderSize + paramsSize + implSize;
}

size_t PluginOperator::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                               const DataType type, char *bufferOut)
{
    // handle common header first
    size_t offset = 0;
    const uint8_t bufferVersion = 1;
    MakeCommonHeader(bufferOut, offset, bufferVersion);

    // plugin specific header, this enables InverseOperate to work in situations
    // where user parameters aren't passed such as bpls
    PutParameters(bufferOut, offset, m_Impl->m_PluginParams);

    // add offset to the bufferOut pointer, so that way the plugin doesn't
    // need to know anything about the plugin header.
    size_t pluginSize =
        m_Impl->m_Plugin->Operate(dataIn, blockStart, blockCount, type, bufferOut + offset);
    return offset + pluginSize;
}

size_t PluginOperator::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t offset = 4; // skip 4 bytes for the common header

    // now handle plugin specific header
    m_Impl->m_PluginParams = GetParameters(bufferIn, offset);

    auto paramPluginNameIt = m_Impl->m_PluginParams.find("PluginName");
    if (paramPluginNameIt == m_Impl->m_PluginParams.end())
    {
        helper::Throw<std::runtime_error>("Plugins", "PluginOperator", "InverseOperate",
                                          "PluginName could not be found in the plugin header");
    }
    std::string pluginName = paramPluginNameIt->second;

    auto paramPluginLibraryIt = m_Impl->m_PluginParams.find("PluginLibrary");
    if (paramPluginLibraryIt == m_Impl->m_PluginParams.end())
    {
        helper::Throw<std::runtime_error>("Plugins", "PluginOperator", "InverseOperate",
                                          "PluginLibrary could not be found in the plugin header");
    }
    const std::string &pluginLibrary = paramPluginLibraryIt->second;

    // now set up the plugin if it hasn't already
    PluginInit(pluginName, pluginLibrary);

    // add offset to bufferIn, so plugin doesn't have to worry about plugin
    // header or common header
    size_t pluginSize =
        m_Impl->m_Plugin->InverseOperate(bufferIn + offset, sizeIn - offset, dataOut);
    return pluginSize;
}

bool PluginOperator::IsDataTypeValid(const DataType type) const
{
    return m_Impl->m_Plugin->IsDataTypeValid(type);
}

} // end namespace plugin
} // end namespace adios2
