/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosNetwork.h
 * Common Network functions needed by staging and streaming engines
 *
 *  Created on: March 22, 2019
 *      Author: Jason Wang
 */

#ifndef ADIOS2_HELPER_ADIOSNETWORK_H_
#define ADIOS2_HELPER_ADIOSNETWORK_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <string>
#include <vector>
/// \endcond

namespace adios2
{
namespace helper
{

class Comm;

/**
 * returns the (a) fully qualified domain name of the current machine.
 * the result is "Unknown_Host_Name" if fqdn is not found
 */
std::string GetFQDN() noexcept;

/**
 * returns a hostname from an FQDN but not the login* or batch* name,
 * instead the second part from such names
 * e.g. (login01.summit.ornl.gov -> summit)
 */
std::string GetClusterName() noexcept;

#ifndef _WIN32
// The function implementations use <nlohmann/json.hpp> which does not
// work with all of the compilers we support.  Provide these functions
// only when one of the engines that needs them is enabled.
#if defined(ADIOS2_HAVE_DATAMAN) || defined(ADIOS2_HAVE_TABLE)

/**
 * returns a vector of strings with all available IP addresses on the node
 * @return vector of strings
 */
std::vector<std::string> AvailableIpAddresses() noexcept;

void HandshakeWriter(Comm const &comm, size_t &appID, std::vector<std::string> &fullAddresses,
                     const std::string &name, const std::string &engineName, const int basePort,
                     const int channelsPerRank, const int maxRanksPerNode = 100,
                     const int maxAppsPerNode = 10);

void HandshakeReader(Comm const &comm, size_t &appID, std::vector<std::string> &fullAddresses,
                     const std::string &name, const std::string &engineName);

#endif // ADIOS2_HAVE_DATAMAN || ADIOS2_HAVE_TABLE
#endif // _WIN32

} // end namespace helper
} // end namespace adios2

#endif // ADIOS2_HELPER_ADIOSNETWORK_H_
