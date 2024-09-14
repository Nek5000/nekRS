/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

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
    adios2::IO io = adios.DeclareIO("TestIO");

    // Create the Engines
    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Read);

    unsigned int t = 0;

    std::vector<std::time_t> write_times;

    std::string varname = "r64";

    while (engine.BeginStep() == adios2::StepStatus::OK)
    {
        size_t writerSize;

        auto var = io.InquireVariable<double>(varname);

        EXPECT_TRUE(var);
        ASSERT_EQ(var.ShapeID(), adios2::ShapeID::GlobalArray);

        /* take the first size as something that gives us writer size */
        writerSize = var.Shape()[0] / 10;

        long unsigned int myStart = (long unsigned int)(writerSize * Nx / mpiSize) * mpiRank;
        long unsigned int myLength = (long unsigned int)(((writerSize)*Nx + mpiSize) / mpiSize);

        if (myStart + myLength > writerSize * Nx)
        {
            myLength = (long unsigned int)writerSize * (int)Nx - myStart;
        }
        const adios2::Dims start{myStart};
        const adios2::Dims count{myLength};
        std::vector<double> in_R64;

        const adios2::Box<adios2::Dims> sel(start, count);

        var.SetSelection(sel);

        in_R64.resize(myLength);
        double VarMin = var.Min();
        double VarMax = var.Max();
        engine.Get(var, in_R64.data());
        engine.EndStep();

        int result = validateSimpleForwardData(in_R64, t, myStart, myLength, writerSize * Nx);

        EXPECT_EQ(VarMin, t * 100);
        EXPECT_EQ(VarMax, t * 100 + writerSize * 10 - 1);

        if (result != 0)
        {
            std::cout << "Read Data Validation failed on node " << mpiRank << " timestep " << t
                      << std::endl;
        }
        EXPECT_EQ(result, 0);

        ++t;
    }

    // Close the file
    engine.Close();
}

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);

    ParseArgs(argc, argv);

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engine == "SST" || engine == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;
    MPI_Init_thread(nullptr, nullptr, thread_support_level, &provided);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 2;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &testComm);
#endif

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
