/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * PluginOperatorInterface.h Operators using the plugin interface should derive
 * from this class.
 *
 *  Created on: Dec 7, 2021
 *      Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#ifndef ADIOS2_OPERATOR_PLUGIN_PLUGINOPERATORINTERFACE_H_
#define ADIOS2_OPERATOR_PLUGIN_PLUGINOPERATORINTERFACE_H_

/// \cond EXCLUDE_FROM_DOXYGEN
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/IO.h"
#include "adios2/core/Operator.h"
#include "adios2/helper/adiosComm.h"

namespace adios2
{
namespace plugin
{

/** An operator interface to be used by the plugin infrastructure */
class PluginOperatorInterface : public core::Operator
{
    // Give the plugin operator access to everything
    friend class PluginOperator;

public:
    PluginOperatorInterface(const Params &parameters);
    virtual ~PluginOperatorInterface() = default;
};

} // end namespace plugin
} // end namespace adios2

#endif /* ADIOS2_OPERATOR_PLUGIN_PLUGINOPERATORINTERFACE_H_ */
