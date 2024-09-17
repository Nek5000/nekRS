/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * PluginEngineInterface.cxx
 *
 *  Created on: July 5, 2021
 *      Author: Chuck Atkins <chuck.atkins@kitware.com>
 *              Caitlin Ross <caitlin.ross@kitware.com>
 */

#include "PluginEngineInterface.h"

namespace adios2
{
namespace plugin
{

PluginEngineInterface::PluginEngineInterface(core::IO &io, const std::string &name, const Mode mode,
                                             helper::Comm comm)
: Engine("PluginInterface", io, name, mode, std::move(comm))
{
}

} // end namespace plugin
} // end namespace adios2
