/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2fstreamMPI.cpp : MPI interfaces to adios2::fstream C++11 bindings
 */

#include "ADIOS2fstream.h"
#include "ADIOS2fstream.tcc"

#include "adios2/helper/adiosCommMPI.h"

namespace adios2
{

fstream::fstream(const std::string &name, const openmode mode, MPI_Comm comm,
                 const std::string engineType)
: m_Stream(std::make_shared<core::Stream>(name, ToMode(mode), helper::CommDupMPI(comm), engineType,
                                          "C++"))
{
}

fstream::fstream(const std::string &name, const openmode mode, MPI_Comm comm,
                 const std::string &configFile, const std::string ioInConfigFile)
: m_Stream(std::make_shared<core::Stream>(name, ToMode(mode), helper::CommDupMPI(comm), configFile,
                                          ioInConfigFile, "C++"))
{
}

void fstream::open(const std::string &name, const openmode mode, MPI_Comm comm,
                   const std::string engineType)
{
    CheckOpen(name);
    m_Stream = std::make_shared<core::Stream>(name, ToMode(mode), helper::CommDupMPI(comm),
                                              engineType, "C++");
}

void fstream::open(const std::string &name, const openmode mode, MPI_Comm comm,
                   const std::string configFile, const std::string ioInConfigFile)
{
    CheckOpen(name);
    m_Stream = std::make_shared<core::Stream>(name, ToMode(mode), helper::CommDupMPI(comm),
                                              configFile, ioInConfigFile, "C++");
}

} // end namespace adios2
