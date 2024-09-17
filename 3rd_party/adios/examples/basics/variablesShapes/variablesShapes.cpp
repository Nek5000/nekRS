/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * * variablesShapes.cpp : adios2 low-level API example to write and read
 *                          supported Variables shapes using stepping
 * (streaming) mode
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

void writer(adios2::ADIOS &adios, const std::size_t nx, const std::size_t nsteps, const int rank,
            const int size)
{
    auto lf_compute = [](const std::size_t step, const std::size_t nx,
                         const int rank) -> std::vector<float> {
        const float value = static_cast<float>(step + rank * nx);
        std::vector<float> array(nx);
        std::iota(array.begin(), array.end(), value);
        return array;
    };

    adios2::IO io = adios.DeclareIO("variables-shapes_writer");

    // You can define variables according to:
    // type <T>: string, uint8_t, int8_t8, ... , float, double
    // shape: global : value or array
    //        local  : value (return a global array), array

    /********** GLOBAL VALUE *************/
    // string variables are always of value type, can't pass dimensions
    adios2::Variable<std::string> varGlobalValueString =
        io.DefineVariable<std::string>("GlobalValueString");

    // global value can change on each step. Example: Step
    adios2::Variable<uint64_t> varStep = io.DefineVariable<uint64_t>("Step");

    /********** GLOBAL ARRAYS *************/
    // For a regular 1D decomposition:

    // 0*nx      1*nx     2*nx      3*nx       shape (4*nx)
    //--------//-------//--------//---------//
    //    nx      nx        nx        nx

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

    // global array with dimensions (shape, start, count)
    // last argument indicates to adios2 that dimensions won't change
    adios2::Variable<float> varGlobalArray =
        io.DefineVariable<float>("GlobalArray", shape, start, count, adios2::ConstantDims);

    /********** LOCAL VALUE **************/
    // Independent values per rank at write, but presented as a global array at
    // read Example: store current rank, but presented as an array of ranks
    adios2::Variable<int32_t> varLocalValueInt32 =
        io.DefineVariable<int32_t>("Rank", {adios2::LocalValueDim});

    /********** LOCAL ARRAY **************/
    // Independent values and dimensions per rank, there is no notion of
    // "continuity", each start from 0-origin to count. Example: mesh
    // nodes-per-rank
    adios2::Variable<float> varLocalArray =
        io.DefineVariable<float>("LocalArray", {}, {}, count, adios2::ConstantDims);

    adios2::Engine writer = io.Open("variables-shapes.bp", adios2::Mode::Write);

    for (size_t step = 0; step < nsteps; ++step)
    {
        // this part mimics the compute portion in an application
        const std::vector<float> array = lf_compute(step, nx, rank);

        // ADIOS2 I/O portion

        // BeginStep/EndStep is the streaming API -> supported by all engines
        writer.BeginStep();

        // minimize global and local values footprint, by only one rank putting
        // the variables
        if (rank == 0)
        {
            // Global value changing over steps
            writer.Put(varStep, static_cast<uint64_t>(step));

            if (step == 0)
            {
                // Global absolute value
                writer.Put(varGlobalValueString, std::string("ADIOS2 Basics Variable Example"));
                // Local absolute value
                writer.Put(varLocalValueInt32, static_cast<int32_t>(rank));
            }
        }

        // for this example all ranks put a global and a local array
        writer.Put(varGlobalArray, array.data());
        writer.Put(varLocalArray, array.data());
        writer.EndStep();
    }
    writer.Close();
}

void reader(adios2::ADIOS &adios, const int rank, const int size)
{
    adios2::IO io = adios.DeclareIO("variables-shapes_reader");
    // all ranks opening the bp file have access to the entire metadata
    adios2::Engine reader = io.Open("variables-shapes.bp", adios2::Mode::Read);

    // reading in streaming mode
    while (reader.BeginStep() != adios2::StepStatus::EndOfStream)
    {
        // scope between BeginStep and EndStep is only for the current step
        const size_t currentStep = reader.CurrentStep();

        // Typical flow: InquireVariable
        adios2::Variable<uint64_t> varStep = io.InquireVariable<uint64_t>("Step");
        uint64_t step = std::numeric_limits<uint64_t>::max();
        // check Variable existence
        if (varStep)
        {
            if (rank == 0)
            {
                // variable objects are "printable" reporting Name and Type
                std::cout << "Found Global Value " << varStep << " in step " << currentStep << "\n";
                // output: Found Global Value Variable<uint64_t>(Name: "Step")
                // in step 0
            }
            reader.Get(varStep, step);
        }

        // GlobalValueString
        adios2::Variable<std::string> varGlobalValueString =
            io.InquireVariable<std::string>("GlobalValueString");
        std::string globalValueString;
        // check Variable existence and Get
        if (varGlobalValueString)
        {
            if (rank == 0)
            {
                std::cout << "Found Global Value " << varGlobalValueString << " in step "
                          << currentStep << "\n";
            }
            reader.Get(varGlobalValueString, globalValueString);
        }

        // Global Arrays at read from local values at write
        adios2::Variable<int32_t> varRanks = io.InquireVariable<int32_t>("Ranks");
        std::vector<int32_t> ranks;
        if (varRanks)
        {
            if (rank == 0)
            {
                std::cout << "Found Global Array " << varRanks << " in step " << currentStep
                          << "\n";
            }
            // passing a vector convenience: adios2 would resize it
            // automatically
            reader.Get(varRanks, ranks);
        }

        // Global Array
        adios2::Variable<float> varGlobalArray = io.InquireVariable<float>("GlobalArray");
        std::vector<float> globalArray;
        if (varGlobalArray)
        {
            if (rank == 0)
            {
                std::cout << "Found GlobalArray " << varGlobalArray << " in step " << currentStep
                          << "\n";
            }
            reader.Get(varGlobalArray, globalArray);
        }

        // Local Array
        adios2::Variable<float> varLocalArray = io.InquireVariable<float>("LocalArray");
        std::vector<float> localArray;
        if (varLocalArray)
        {
            // local arrays require an extra step to select the block of
            // interest (0 is default) we only select block 0 in this example
            varLocalArray.SetBlockSelection(0);
            if (rank == 0)
            {
                std::cout << "Found LocalArray " << varLocalArray << " in step " << currentStep
                          << "\n";
            }
            reader.Get(varLocalArray, localArray);
        }

        // since all Get calls are "deferred" all the data would be populated at
        // EndStep
        reader.EndStep();

        // data is available

        if (rank == 0)
        {
            std::cout << "\n";
        }
    }

    reader.Close();
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
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        constexpr std::size_t nx = 10;
        constexpr std::size_t nsteps = 3;

        writer(adios, nx, nsteps, rank, size);
        reader(adios, rank, size);
    }
    catch (const std::exception &e)
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
