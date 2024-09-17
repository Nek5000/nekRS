/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <chrono>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <thread>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

class TestOnDemandMPI : public ::testing::Test
{
public:
    TestOnDemandMPI() = default;
};

void DoWriter(adios2::IO io)
{

    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    unsigned int myStart = (int)0;
    unsigned int myCount = (int)Nx;
    adios2::Dims shape{static_cast<unsigned int>(Nx)};
    adios2::Dims start{static_cast<unsigned int>(myStart)};
    adios2::Dims count{static_cast<unsigned int>(myCount)};

    auto var = io.DefineVariable<double>("r64", shape, start, count);
    auto stepvar = io.DefineVariable<size_t>("Step");

    // Create the Engine
    io.SetEngine(engine);
    engineParams["StepDistributionMode"] = "OnDemand";
    io.SetParameters(engineParams);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    for (size_t step = 0; step < (size_t)NSteps; step++)
    {

        // Generate test data for each process uniquely
        std::vector<double> data_forward;

        generateSimpleForwardData(data_forward, (int)step, myStart, myCount, (int)Nx);

        engine.BeginStep();

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces
        adios2::Box<adios2::Dims> sel({myStart}, {myCount});

        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count
        const adios2::Mode sync = GlobalWriteMode;

        var.SetSelection(sel);
        engine.Put(var, data_forward.data(), sync);
        engine.Put(stepvar, step);
        engine.EndStep();
    }
    // Close the file
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    engine.Close();
    int steps = 0;
    int received_steps = 0;
    MPI_Reduce(&steps, &received_steps, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    EXPECT_EQ(NSteps, received_steps);
}

void DoReader(adios2::IO io, int Rank)
{
    io.SetEngine(engine);
    adios2::Engine sstReader = io.Open(fname, adios2::Mode::Read);

    int get_count;
    const std::size_t my_start = 0;
    const adios2::Dims pos_start{my_start};
    const adios2::Dims count{Nx};
    const adios2::Box<adios2::Dims> sel(pos_start, count);

    int steps = 0;
    std::vector<double> myFloats;
    while (sstReader.BeginStep() == adios2::StepStatus::OK)
    {
        adios2::Variable<size_t> stepVar = io.InquireVariable<size_t>("Step");
        adios2::Variable<double> floatVar = io.InquireVariable<double>("r64");
        size_t writerStep;
        sstReader.Get(floatVar, myFloats);
        sstReader.Get(stepVar, writerStep);
        //	    std::cout << "Reader " << Rank << " got writerStep " <<
        // writerStep << std::endl;
        sstReader.EndStep();
        steps += 1;
    }
    sstReader.Close();
    std::cout << "Reader " << Rank << " got " << steps << " steps " << std::endl;
    MPI_Reduce(&steps, &get_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}

// ADIOS2 COMMON write
TEST_F(TestOnDemandMPI, ADIOS2OnDemandMPI)
{
    // form a mpiSize * Nx 1D array
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Init without MPI.  MPI only used to coordinate non-MPI actors
    adios2::ADIOS adios;

    adios2::IO io = adios.DeclareIO("TestIO");

    if (mpiRank == 0)
    {
        DoWriter(io);
    }
    else
    {

        DoReader(io, mpiRank);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(nullptr, nullptr);

    int result;
    ::testing::InitGoogleTest(&argc, argv);

    NSteps = 100;

    ParseArgs(argc, argv);

    result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
