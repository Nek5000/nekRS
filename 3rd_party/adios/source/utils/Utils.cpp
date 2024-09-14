/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Utils.cpp
 *
 *  Created on: Oct 24, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "Utils.h"

namespace adios2
{

Utils::Utils(const std::string name, int argc, char *argv[])
: m_Name(name), m_Arguments(argv, argv + argc)
{
}

} // end namespace adios2
