/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosXML.h basic XML parsing functionality for ADIOS config file schema
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_HELPER_ADIOSXML_H_
#define ADIOS2_HELPER_ADIOSXML_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <map>
#include <memory> //std::shared_ptr
#include <string>
#include <utility> //std::pair
/// \endcond

#include "adios2/core/ADIOS.h"
#include "adios2/core/IO.h"

namespace adios2
{
namespace helper
{

void ParseConfigXMLIO(core::ADIOS &adios, const std::string &configFileXML,
                      const std::string &configFileContents, core::IO &io,
                      std::unordered_map<std::string, std::pair<std::string, Params>> &operators);

std::string
ParseConfigXML(core::ADIOS &adios, const std::string &configFile,
               std::map<std::string, core::IO> &ios,
               std::unordered_map<std::string, std::pair<std::string, Params>> &operators);

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSXML_H_ */
