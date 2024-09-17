/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_hdf5_a.cpp
 *
 * Write output with sequential HDF5, one file per process, one separate set per
 * timestep
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <hdf5.h>

IO::IO(const Settings &s, MPI_Comm comm) : m_outputfilename{s.outputfile} {}

IO::~IO() {}

void IO::write(int step, const HeatTransfer &ht, const Settings &s, MPI_Comm comm)
{
    m_outputfilename = MakeFilename(s.outputfile, ".h5", s.rank, step);

    // for time measurements, let's synchronize the processes
    MPI_Barrier(comm);
    double time_start = MPI_Wtime();

    hsize_t dims[2] = {static_cast<hsize_t>(s.ndx), static_cast<hsize_t>(s.ndy)};

    hid_t space = H5Screate_simple(2, dims, NULL);
    hid_t file = H5Fcreate(m_outputfilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset =
        H5Dcreate(file, "T", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ht.data_noghost().data());

    H5Dclose(dset);
    H5Sclose(space);
    H5Fclose(file);

    MPI_Barrier(comm);
    double total_time = MPI_Wtime() - time_start;
    uint64_t adios_totalsize = 2 * sizeof(int) + 2 * s.ndx * s.ndy * sizeof(double);
    uint64_t sizeMB = adios_totalsize * s.nproc / 1024 / 1024 / 1024; // size in MB
    double mbs = static_cast<double>(sizeMB) / total_time;
    if (s.rank == 0)
        std::cout << "Step " << step << ": " << m_outputfilename << " " << sizeMB << " "
                  << total_time << "" << mbs << std::endl;
}
