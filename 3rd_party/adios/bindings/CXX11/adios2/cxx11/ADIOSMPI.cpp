/* Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOSMPI.cpp : MPI constructors of public ADIOS class for C++11 bindings
 */

#include "ADIOS.h"

#include "adios2/core/ADIOS.h"

#include "adios2/helper/adiosCommMPI.h"

namespace adios2
{
ADIOS::ADIOS(const std::string &configFile, MPI_Comm comm)
: m_ADIOS(std::make_shared<core::ADIOS>(configFile, helper::CommDupMPI(comm), "C++"))
{
}

ADIOS::ADIOS(MPI_Comm comm) : ADIOS("", comm) {}

ADIOS::ADIOS(const std::string &configFile, MPI_Comm comm, const std::string &hostLanguage)
: m_ADIOS(std::make_shared<core::ADIOS>(configFile, helper::CommDupMPI(comm), hostLanguage))
{
}

} // end namespace adios2
