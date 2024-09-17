/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestUtilsChangingShape.cpp :
 *
 *  Created on: Dec 17, 2018
 *      Author: Norbert Podhorszki, Keichi Takahashi
 */

#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

std::string engineName; // comes from command line

int main(int argc, char **argv)
{
    // Each process would write a 4x2 array and all processes would
    // form a 2D 4 * (NumberOfProcess * Nx) matrix where Nx is 2 here

    const std::string fname("TestUtilsChangingShape.bp");
    const int nsteps = 10;
    // Number of rows
    const size_t Nx = 8;

    int rank = 0, nproc = 1;

#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    // Writer
    adios2::IO outIO = adios.DeclareIO("Output");

    if (!engineName.empty())
    {
        outIO.SetEngine(engineName);
    }

    adios2::Engine writer = outIO.Open(fname, adios2::Mode::Write);

    adios2::Dims shape{static_cast<size_t>(nproc), Nx};
    adios2::Dims start{static_cast<size_t>(rank * Nx), 0};
    adios2::Dims count{1, Nx};
    // variable that changes shape every step
    auto var_ch = outIO.DefineVariable<double>("ChangingShapeVar", shape, start, count);
    // variable that is written every other steps but not changing shape
    auto var_alt = outIO.DefineVariable<double>("AlternatingStepsVar", shape, start, count);
    // variable that is written every other steps AND is changing shape
    auto var_altch =
        outIO.DefineVariable<double>("AlternatingStepsAndChangingShapeVar", shape, start, count);
    // Other ("normal") variables
    auto var_fixed = outIO.DefineVariable<double>("FixedShapeVar", shape, start, count);
    auto var_single = outIO.DefineVariable<double>("SingleStepVar", shape, start, count);

    std::vector<double> buf(Nx + nsteps / 2 + 1, 0.0);

    if (!rank)
    {
        std::cout << "Writing:" << std::endl;
    }
    for (size_t i = 0; i < nsteps; i++)
    {
        for (size_t j = 0; j < buf.size(); j++)
        {
            buf[j] = i + rank / 10.0 + static_cast<double>(j) / 100.0;
        }

        writer.BeginStep();

        if (i > 0)
        {
            shape[1] = Nx - nsteps / 2 + i;
            count[1] = shape[1];
            var_ch.SetShape(shape);
            var_ch.SetSelection({start, count});
            var_altch.SetShape(shape);
            var_altch.SetSelection({start, count});
        }

        if (!rank)
        {
            std::cout << "Step " << i << " shape (" << var_ch.Shape()[0] << ", "
                      << var_ch.Shape()[1] << ")" << std::endl;
        }

        if (i == 0)
        {
            writer.Put(var_single, buf.data());
        }

        writer.Put(var_ch, buf.data());
        writer.Put(var_fixed, buf.data());

        if (i % 2 == 0)
        {
            writer.Put(var_alt, buf.data());
            writer.Put(var_altch, buf.data());
        }

        writer.EndStep();
    }

    writer.Close();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
