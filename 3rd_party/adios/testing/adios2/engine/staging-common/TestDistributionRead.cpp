/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <chrono>
#include <cstdint>
#include <cstring>
#include <stdint.h> /* SIZE_MAX */
#include <thread>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::this_thread::sleep_for;

class CommonReadTest : public ::testing::Test
{
public:
    CommonReadTest() = default;
};

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

// ADIOS2 Common read
TEST_F(CommonReadTest, ADIOS2CommonRead1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
    MPI_Comm_size(testComm, &mpiSize);
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    /* we won't share IOs on the reader side */
    adios2::IO io1 = adios.DeclareIO("TestIO");
    adios2::IO io2 = adios.DeclareIO("TestIO2");

    // Create the Engines
    io1.SetEngine(engine);
    io1.SetParameters(engineParams);
    io2.SetEngine(engine);
    io2.SetParameters(engineParams);

    adios2::Engine engine1 = io1.Open(fname, adios2::Mode::Read);

    std::string varname1 = "r64";

    size_t first_step = SIZE_MAX;
    size_t total_steps = 0;
    while (engine1.BeginStep() == adios2::StepStatus::OK)
    {
        size_t writerSize;
        size_t step;
        auto var1 = io1.InquireVariable<double>(varname1);

        EXPECT_TRUE(var1);
        ASSERT_EQ(var1.ShapeID(), adios2::ShapeID::GlobalArray);

        auto step_var = io1.InquireVariable<size_t>("Step");

        EXPECT_TRUE(step_var);
        ASSERT_EQ(step_var.ShapeID(), adios2::ShapeID::GlobalValue);

        /* take the first size as something that gives us writer size */
        writerSize = var1.Shape()[0] / 10;

        long unsigned int myStart = (long unsigned int)(writerSize * Nx / mpiSize) * mpiRank;
        long unsigned int myLength = (long unsigned int)((writerSize * Nx + mpiSize - 1) / mpiSize);

        if (myStart + myLength > writerSize * Nx)
        {
            myLength = (long unsigned int)(writerSize + 1) * (int)Nx - myStart;
        }
        const adios2::Dims start{myStart};
        const adios2::Dims count{myLength};
        std::vector<double> in_R64_1;

        const adios2::Box<adios2::Dims> sel(start, count);

        var1.SetSelection(sel);

        in_R64_1.resize(myLength);
        engine1.Get(var1, in_R64_1.data());
        engine1.Get(step_var, &step);
        engine1.EndStep();

        int result =
            validateSimpleForwardData(in_R64_1, (int)step, myStart, myLength, writerSize * Nx);
        if (first_step == SIZE_MAX)
        {
            first_step = step;
        }
        if (result != 0)
        {
            std::cout << "Read Data Validation failed on node " << mpiRank << " timestep " << step
                      << std::endl;
        }
        EXPECT_EQ(result, 0);
        total_steps++;
    }
    if (RoundRobin)
    {
        if (first_step == 0)
        {
            EXPECT_EQ(total_steps, 4);
        }
        else
        {
            EXPECT_EQ(total_steps, 3);
        }
    }
    else if (OnDemand)
    {
        std::cout << " Total steps received is " << total_steps << std::endl;
    }
    else
    {
        EXPECT_EQ(first_step, 0);
        EXPECT_EQ(total_steps, 10);
    }
    // Close the file
    engine1.Close();
}

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    MPI_Init(nullptr, nullptr);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 2;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &testComm);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);

    ParseArgs(argc, argv);

    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
#ifdef CRAY_MPICH_VERSION
    MPI_Barrier(MPI_COMM_WORLD);
#else
    MPI_Finalize();
#endif
#endif

    return result;
}
