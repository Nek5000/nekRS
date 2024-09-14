/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "restart.h"

#include <iostream>

static bool firstCkpt = true;

using DataType = std::complex<double>;

void WriteCkpt(MPI_Comm comm, const int step, const Settings &settings, const GrayScott &sim,
               adios2::IO io)
{
    int rank, nproc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);
    std::cout << "checkpoint at step " << step << " create file " << settings.checkpoint_output
              << std::endl;
    adios2::Engine writer = io.Open(settings.checkpoint_output, adios2::Mode::Write);
    if (writer)
    {
        adios2::Variable<DataType> var_uv;
        adios2::Variable<int> var_step;

        if (firstCkpt)
        {
            size_t X = sim.size_x + 2;
            size_t Y = sim.size_y + 2;
            size_t Z = sim.size_z + 2;
            size_t R = static_cast<size_t>(rank);
            size_t N = static_cast<size_t>(nproc);

            var_uv = io.DefineVariable<DataType>("UV", {N, X, Y, Z}, {R, 0, 0, 0}, {1, X, Y, Z});

            var_step = io.DefineVariable<int>("step");
            firstCkpt = false;
        }
        else
        {
            var_uv = io.InquireVariable<DataType>("UV");
            var_step = io.InquireVariable<int>("step");
        }

        writer.Put<int>(var_step, &step);
        const DataType *ptr = reinterpret_cast<const DataType *>(sim.d_ghost().data());
        writer.Put<DataType>(var_uv, ptr);

        writer.Close();
    }
}

int ReadRestart(MPI_Comm comm, const Settings &settings, GrayScott &sim, adios2::IO io)
{
    int step = 0;
    int rank, nproc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);
    if (!rank)
    {
        std::cout << "restart from file " << settings.restart_input << std::endl;
    }
    adios2::Engine reader = io.Open(settings.restart_input, adios2::Mode::ReadRandomAccess);
    if (reader)
    {
        adios2::Variable<int> var_step = io.InquireVariable<int>("step");
        adios2::Variable<DataType> var_uv = io.InquireVariable<DataType>("UV");
        size_t X = sim.size_x + 2;
        size_t Y = sim.size_y + 2;
        size_t Z = sim.size_z + 2;
        size_t R = static_cast<size_t>(rank);
        std::vector<GrayScott::MemLayout> uv(X * Y * Z);

        var_uv.SetSelection({{R, 0, 0, 0}, {1, X, Y, Z}});
        reader.Get<int>(var_step, step);
        reader.Get<DataType>(var_uv, reinterpret_cast<DataType *>(uv.data()));
        reader.Close();

        if (!rank)
        {
            std::cout << "restart from step " << step << std::endl;
        }
        sim.restart(uv);
    }
    else
    {
        std::cout << "    failed to open file " << std::endl;
    }
    return step;
}
