/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Utils.h
 *
 *  Created on: Oct 24, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef SOURCE_UTILS_UTILS_H_
#define SOURCE_UTILS_UTILS_H_

#include <string>
#include <vector>

#include "adios2/common/ADIOSTypes.h"

namespace adios2
{

class Utils
{
public:
    const std::string m_Name;

    Utils(const std::string name, int argc, char *argv[]);

    virtual ~Utils() = default;

    virtual void Run() = 0;

protected:
    const std::vector<std::string> m_Arguments;
    Params m_Parameters;

    virtual void ParseArguments() = 0;
    virtual void ProcessParameters() = 0;
    virtual void PrintUsage() const noexcept = 0;
    virtual void PrintExamples() const noexcept = 0;
    virtual void SetParameters(const std::string argument, const bool isLong) = 0;
};

} /// end namespace adios2

#endif /* SOURCE_UTILS_UTILS_H_ */
