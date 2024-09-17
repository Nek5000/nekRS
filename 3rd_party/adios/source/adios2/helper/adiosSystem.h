/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosSystem.h  system related functions using std or POSIX,
 * we might wait for C++17 filesystem
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_HELPER_ADIOSSYSTEM_H_
#define ADIOS2_HELPER_ADIOSSYSTEM_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <string>
#include <vector>
/// \endcond

#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"

namespace adios2
{
namespace helper
{
/**
 * Creates a chain of directories using POSIX systems calls (stat, mkdir),
 * Verifies if directory exists before creating a new one. Permissions are 777
 * for now
 * @param fullPath /full/path/for/directory
 * @return true: directory exists, false: failed to create or access directory
 */
bool CreateDirectory(const std::string &fullPath) noexcept;
/**
 * Check if system is little endian
 * @return true: little endian, false: big endian
 */
bool IsLittleEndian() noexcept;

/**
 * returns a string with current local time and date information from std::ctime
 * @return string from char* std::ctime
 */
std::string LocalTimeDate() noexcept;

/**
 * Support for language bindings, identify if data is row-major (C, C++) or not
 * (Fortran, R)
 * @param hostLanguage input host language
 * @return true: row-major, false: column-major
 */
bool IsRowMajor(const std::string hostLanguage) noexcept;

/**
 * Support for language bindings, identify if data is zero-indexed (C, C++) or
 * not (Fortran, R)
 * @param hostLanguage input host language
 * @return true: zero-indexed, false: one-indexed
 */
bool IsZeroIndexed(const std::string hostLanguage) noexcept;

/**
 * Function to be called inside a catch(...) block to convert C++ exceptions to
 * error integers. Used by C, Fortran bindings.
 * @param function function name to be printed in std::cerr
 * @return error integers mapped to exceptions:
 * <pre>
 * 	none                   : 0
 *  std::invalid_argument  : 1
 *  std::system_error      : 2
 *  std::runtime_exception : 3
 *  std::exception         : 4
 * </pre>
 */
int ExceptionToError(const std::string &function);

bool IsHDF5File(const std::string &name, core::IO &io, helper::Comm &comm,
                const std::vector<Params> &transportsParameters) noexcept;
char BPVersion(const std::string &name, helper::Comm &comm,
               const std::vector<Params> &transportsParameters) noexcept;

/** Return the number of available hardware threads on the node.
 * It might return 0 if the detection does not work
 */
unsigned int NumHardwareThreadsPerNode();

/** Attempt to raise the limit of number of files opened at once
 *  Return: the limit on number of open files
 */
size_t RaiseLimitNoFile();

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSSYSTEM_H_ */
