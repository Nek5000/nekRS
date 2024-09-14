/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosYAML.h basic YAML parsing functionality for ADIOS config file schema
 *
 *  Created on: Oct 24, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_HELPER_ADIOSYAML_H_
#define ADIOS2_HELPER_ADIOSYAML_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <map>
#include <memory> //std::shared_ptr
#include <string>
#include <utility> //std::pair
/// \endcond

#include "adios2/common/ADIOSTypes.h" // UserOptions
#include "adios2/core/ADIOS.h"
#include "adios2/core/IO.h"

namespace adios2
{
namespace helper
{
void ParseConfigYAMLIO(core::ADIOS &adios, const std::string &configFileYAML,
                       const std::string &configFileContents, core::IO &io);

std::string ParseConfigYAML(core::ADIOS &adios, const std::string &configFileYAML,
                            std::map<std::string, core::IO> &ios);

void ParseUserOptionsFile(Comm &comm, const std::string &configFileYAML, UserOptions &options,
                          std::string &homePath);

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSXML_H_ */
