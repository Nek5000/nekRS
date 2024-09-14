/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * variablesShapes_hl.cpp : adios2 high-level API example to write and read
 *                   supported Variables shapes using stepping (streaming) mode
 *
 *  Created on: Nov 14, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <cstddef>   //std::size_t
#include <iostream>  // std::cout
#include <limits>    // std::numeric_limits
#include <numeric>   //std::iota
#include <stdexcept> //std::exception

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

void writer(const std::size_t nx, const std::size_t nsteps, const int rank, const int size)
{
    auto lf_compute = [](const std::size_t step, const std::size_t nx,
                         const int rank) -> std::vector<float> {
        const float value = static_cast<float>(step + rank * nx);
        std::vector<float> array(nx);
        std::iota(array.begin(), array.end(), value);
        return array;
    };

    // You can define variables according to:
    // type <T>: string, uint8_t, int8_t8, ... , float, double
    // shape: global : value or array
    //        local  : value (return a global array), array

    // global shape -> this is the physical dimension across MPI processes
    const adios2::Dims shape = {static_cast<std::size_t>(size * nx)};

    // local start for rank offset -> this is the local origin for the rank
    // domain
    const adios2::Dims start = {static_cast<std::size_t>(rank * nx)};

    // local count -> this is the local size from the local start for the rank
    // domain
    const adios2::Dims count = {nx};

    // adios2::Dims is an alias to std::vector<std::size_t>
    // helps remember the inputs to adios2 functions DefineVariable (write) and
    // SetSelection (read) make sure you always pass std::size_t types

#if ADIOS2_USE_MPI
    adios2::fstream out("variables-shapes_hl.bp", adios2::fstream::out, MPI_COMM_WORLD);
#else
    adios2::fstream out("variables-shapes_hl.bp", adios2::fstream::out);
#endif

    for (size_t step = 0; step < nsteps; ++step)
    {
        // this part mimics the compute portion in an application
        const std::vector<float> array = lf_compute(step, nx, rank);

        // ADIOS2 I/O portion

        // minimize global and local values footprint, by only one rank writing
        // the variables
        if (rank == 0)
        {
            // Global value changing over steps
            out.write("Step", static_cast<uint64_t>(step));

            if (step == 0)
            {
                // Constant Global value
                out.write("GlobalValueString", std::string("ADIOS2 Basics Variable Example"));

                // Constant Local value
                out.write("LocalValueInt32", static_cast<int32_t>(rank), adios2::LocalValue);
            }
        }

        // for this example all ranks write a global and a local array
        out.write("GlobalArray", array.data(), shape, start, count);
        out.write("LocalArray", array.data(), {}, {}, count);

        out.end_step();
    }
    out.close();
}

void reader(const int rank, const int size)
{
    auto lf_ArrayToString = [](const std::vector<float> &array) -> std::string {
        std::string contents = "{ ";
        for (const float value : array)
        {
            contents += std::to_string(static_cast<int>(value)) + " ";
        }
        contents += "}";
        return contents;
    };

// all ranks opening the bp file have access to the entire metadata
#if ADIOS2_USE_MPI
    adios2::fstream in("variables-shapes_hl.bp", adios2::fstream::in, MPI_COMM_WORLD);
#else
    adios2::fstream in("variables-shapes_hl.bp", adios2::fstream::in);
#endif

    // reading in streaming mode, supported by all engines
    // similar to std::getline in std::fstream
    adios2::fstep inStep;
    while (adios2::getstep(in, inStep))
    {
        const std::size_t currentStep = inStep.current_step();

        const std::vector<uint64_t> steps = inStep.read<uint64_t>("Step");
        if (!steps.empty() && rank == 0)
        {
            std::cout << "Found Step " << steps.front() << " in currentStep " << currentStep
                      << "\n";
        }

        const std::vector<std::string> globalValueString =
            inStep.read<std::string>("GlobalValueString");
        if (!globalValueString.empty() && rank == 0)
        {
            std::cout << "Found GlobalValueString " << globalValueString.front()
                      << " in currentStep " << currentStep << "\n";
        }

        const std::vector<int32_t> ranks = inStep.read<int32_t>("Ranks");
        if (!ranks.empty() && rank == 0)
        {
            std::cout << "Found rank " << ranks.front() << " in currentStep " << currentStep
                      << "\n";
        }

        const std::vector<float> globalArray = inStep.read<float>("GlobalArray");
        if (!globalArray.empty() && rank == 0)
        {
            std::cout << "Found globalArray " << lf_ArrayToString(globalArray) + " in currentStep "
                      << currentStep << "\n";
        }

        // default reads block 0
        const std::vector<float> localArray = inStep.read<float>("LocalArray");
        if (!localArray.empty() && rank == 0)
        {
            std::cout << "Found localArray " << lf_ArrayToString(localArray) + " in currentStep "
                      << currentStep << "\n";
        }
        // indicate end of adios2 operations for this step
        in.end_step();
    }
    in.close();
}

int main(int argc, char *argv[])
{
#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
#endif
    int rank = 0;
    int size = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    try
    {
        constexpr std::size_t nx = 10;
        constexpr std::size_t nsteps = 3;

        writer(nx, nsteps, rank, size);
        reader(rank, size);
    }
    catch (std::exception &e)
    {
        std::cout << "ERROR: ADIOS2 exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
