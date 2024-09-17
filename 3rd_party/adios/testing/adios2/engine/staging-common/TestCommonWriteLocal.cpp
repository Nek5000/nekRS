/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>
#include <ctime>

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <thread>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

class CommonWriteTest : public ::testing::Test
{
public:
    CommonWriteTest() = default;
};

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

// ADIOS2 COMMON write
TEST_F(CommonWriteTest, ADIOS2CommonWrite)
{
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
    adios2::IO io = adios.DeclareIO("TestIO");

    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().

    /* we'll be dropping the size for some vars */
    if (VaryingDataSize)
    {
        assert(Nx > static_cast<size_t>(NSteps + mpiSize));
    }

    {
        adios2::Dims count{static_cast<unsigned int>(Nx)};
        adios2::Dims time_shape{static_cast<unsigned int>(mpiSize)};
        adios2::Dims time_start{static_cast<unsigned int>(mpiRank)};
        adios2::Dims time_count{1};

        (void)io.DefineVariable<float>("r32", {}, {}, count);
        (void)io.DefineVariable<int64_t>("time", time_shape, time_start, time_count);
    }

    // Create the Engine
    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    data_R32.resize(Nx);
    for (int step = 0; step < NSteps; ++step)
    {
        engine.BeginStep();
        auto var_r32 = io.InquireVariable<float>("r32");
        auto var_time = io.InquireVariable<int64_t>("time");

        const adios2::Mode sync = adios2::Mode::Deferred;

        for (int index = 0; index < LocalCount; index++)
        {
            for (size_t i = 0; i < Nx; i++)
            {
                int64_t j = (mpiRank * LocalCount + index) * Nx * 10 + step;
                data_R32[i] = (float)j + 10 * i;
            }

            engine.Put(var_r32, data_R32.data(), sync);
        }
        // Advance to the next time step
        std::time_t localtime = std::time(NULL);
        engine.Put(var_time, (int64_t *)&localtime);
        engine.EndStep();
        std::this_thread::sleep_for(
            std::chrono::milliseconds(DelayMS)); /* sleep for DelayMS milliseconds */
    }

    // Close the file
    engine.Close();
}

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);

    DelayMS = 0; // zero default for common writer

    ParseArgs(argc, argv);

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engine == "SST" || engine == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, thread_support_level, &provided);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 1;
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
