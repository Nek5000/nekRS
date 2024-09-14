/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * PluginEngineInterface.h Engines using the plugin interface should derive from
 * this class.
 *
 *  Created on: July 5, 2021
 *      Author: Chuck Atkins <chuck.atkins@kitware.com>
 *              Caitlin Ross <caitlin.ross@kitware.com>
 */

#ifndef ADIOS2_ENGINE_PLUGIN_PLUGINENGINEINTERFACE_H_
#define ADIOS2_ENGINE_PLUGIN_PLUGINENGINEINTERFACE_H_

/// \cond EXCLUDE_FROM_DOXYGEN
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"

namespace adios2
{
namespace plugin
{

/** An engine interface to be used by the plugin infrastructure */
class PluginEngineInterface : public core::Engine
{
    // Give the plugin engine access to everything
    friend class PluginEngine;

public:
    PluginEngineInterface(core::IO &io, const std::string &name, const Mode mode,
                          helper::Comm comm);
    virtual ~PluginEngineInterface() = default;
};

} // end namespace plugin
} // end namespace adios2

#endif /* ADIOS2_ENGINE_PLUGIN_PLUGINENGINEINTERFACE_H_ */
