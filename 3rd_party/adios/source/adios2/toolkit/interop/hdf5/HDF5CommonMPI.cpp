/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5CommonMPI.cpp : MPI-specific HDF5 engine functionality
 */

#include "adios2/helper/adiosCommMPI.h"

#include "adios2/toolkit/interop/hdf5/HDF5Common.h"

#include <mutex>

namespace adios2
{

namespace interop
{

namespace
{

bool init(helper::Comm const &comm, hid_t id, int *rank, int *size)
{
    MPI_Comm mpiComm = helper::CommAsMPI(comm);
    if (mpiComm == MPI_COMM_NULL)
    {
        return false;
    }
    MPI_Comm_rank(mpiComm, rank);
    MPI_Comm_size(mpiComm, size);
    if (*size != 1)
    {
        H5Pset_fapl_mpio(id, mpiComm, MPI_INFO_NULL);
    }
    return true;
}

herr_t set_dxpl_mpio(hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode)
{
    return H5Pset_dxpl_mpio(dxpl_id, xfer_mode);
}

HDF5Common::MPI_API HDF5Common_MPI_API_Impl = {&init, &set_dxpl_mpio};

} // end anonymous namespace

extern std::mutex HDF5Common_MPI_API_Mutex;
extern HDF5Common::MPI_API const *HDF5Common_MPI_API;

void RegisterHDF5Common_MPI_API()
{
    std::lock_guard<std::mutex> guard(HDF5Common_MPI_API_Mutex);
    HDF5Common_MPI_API = &HDF5Common_MPI_API_Impl;
}

} // namespace interop
} // end namespace adios2
