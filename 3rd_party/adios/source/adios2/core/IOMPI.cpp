/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOMPI.cpp : MPI-specific IO engine factory registration
 */

#include "IO.h"
#include "adios2/engine/ssc/SscReader.h"
#include "adios2/engine/ssc/SscWriter.h"
#include "adios2/helper/adiosCommMPI.h"

#ifdef ADIOS2_HAVE_DATASPACES // external dependencies
#include "adios2/engine/dataspaces/DataSpacesReader.h"
#include "adios2/engine/dataspaces/DataSpacesWriter.h"
#endif

namespace adios2
{

#if defined(ADIOS2_HAVE_HDF5_PARALLEL)
namespace interop
{
void RegisterHDF5Common_MPI_API();
} // namespace interop
#endif

namespace core
{

namespace
{
template <typename T>
std::shared_ptr<Engine> MakeEngineMPI(IO &io, const std::string &name, const Mode mode,
                                      helper::Comm comm)
{
    if (!comm.IsMPI())
    {
        helper::Throw<std::invalid_argument>("Core", "IOMPI", "MakeEngineMPI",
                                             "A MPI-only engine cannot be used with a "
                                             "communicator that is not MPI-based.");
    }
    return IO::MakeEngine<T>(io, name, mode, std::move(comm));
}
}

void RegisterMPIEngines()
{
    IO::RegisterEngine("ssc", IO::EngineFactoryEntry{MakeEngineMPI<engine::SscReader>,
                                                     MakeEngineMPI<engine::SscWriter>});
#ifdef ADIOS2_HAVE_DATASPACES
    IO::RegisterEngine("dataspaces",
                       IO::EngineFactoryEntry{MakeEngineMPI<engine::DataSpacesReader>,
                                              MakeEngineMPI<engine::DataSpacesWriter>});
#endif
#if defined(ADIOS2_HAVE_HDF5_PARALLEL)
    interop::RegisterHDF5Common_MPI_API();
#endif
}

} // end namespace core
} // end namespace adios2
