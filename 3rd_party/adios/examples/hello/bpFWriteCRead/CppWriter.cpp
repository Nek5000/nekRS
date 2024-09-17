/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CppWriter.cpp
 *
 *  Created on: Jan 29, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <adios2.h>
#include <mpi.h>

#include <iostream>

int main(int argc, char *argv[])
{
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const size_t nx = 4;
    const size_t ny = 3;
    std::vector<float> data(nx * ny);

    for (size_t i = 0; i < nx; ++i)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            data[i * ny + j] = static_cast<float>(rank * nx * ny + i * ny + j);
        }
    }

    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io = adios.DeclareIO("CppWriter");

    adios2::Variable<float> bpFloats = io.DefineVariable<float>(
        "data2D", {size * nx, ny}, {rank * nx, 0}, {nx, ny}, adios2::ConstantDims);

    adios2::Engine engine = io.Open("CppWriter.bp", adios2::Mode::Write);
    engine.BeginStep();
    engine.Put(bpFloats, data.data());
    engine.EndStep();
    engine.Close();

    MPI_Finalize();
}
