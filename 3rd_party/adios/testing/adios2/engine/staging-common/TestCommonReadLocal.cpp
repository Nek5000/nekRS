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

    int TimeGapDetected = 0;
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

    // Create the Engine
    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Read);

    unsigned int t = 0;

    std::vector<std::time_t> write_times;

    while (engine.BeginStep() == adios2::StepStatus::OK)
    {
        const size_t currentStep = engine.CurrentStep();
        EXPECT_EQ(currentStep, static_cast<size_t>(t));

        size_t writerSize;

        auto var_time = io.InquireVariable<int64_t>("time");
        EXPECT_TRUE(var_time);
        ASSERT_EQ(var_time.ShapeID(), adios2::ShapeID::GlobalArray);

        writerSize = var_time.Shape()[0];

        int rankToRead = mpiRank;
        if (writerSize < static_cast<size_t>(mpiSize))
        {
            rankToRead = mpiRank % writerSize;
        }

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);

        long unsigned int hisLength = (long unsigned int)Nx;

        for (int index = 0; index < LocalCount; index++)
        {
            int indexToRead = rankToRead * LocalCount;
            ASSERT_EQ(engine.BlocksInfo(var_r32, currentStep).at(indexToRead).Count[0], hisLength);
        }

        const adios2::Dims start_time{0};
        const adios2::Dims count_time{1};
        const adios2::Box<adios2::Dims> sel_time(start_time, count_time);
        var_time.SetSelection(sel_time);
        auto BI = engine.BlocksInfo(var_r32, currentStep);
        in_R32_blocks.resize(BI.size());
        {
            for (size_t index = 0; index < BI.size(); index++)
            {
                in_R32_blocks[index].resize(hisLength);
                var_r32.SetBlockSelection(index);
                engine.Get(var_r32, in_R32_blocks[index].data());
            }
        }
        std::time_t write_time;
        engine.Get(var_time, (int64_t *)&write_time);
        engine.EndStep();

        int result = 0;
        for (size_t index = 0; index < BI.size(); index++)
        {
            for (size_t i = 0; i < Nx; i++)
            {
                int64_t j = index * Nx * 10 + t;
                float expected = (float)j + 10 * i;
                if (in_R32_blocks[index][i] != expected)
                {
                    std::cout << "Expected " << expected << ", got " << in_R32_blocks[index][i]
                              << " for in_R32[" << i << "][" << index << "[" << i << "], timestep "
                              << t << std::endl;
                    result++;
                }
            }
        }

        if (result != 0)
        {
            std::cout << "Read Data Validation failed on node " << mpiRank << " timestep " << t
                      << std::endl;
        }
        EXPECT_EQ(result, 0);
        write_times.push_back(write_time);
        ++t;
    }

    if ((write_times.size() > 1) && ((write_times.back() - write_times.front()) > 1))
    {
        TimeGapDetected++;
    }

    if (!IgnoreTimeGap)
    {
        if (TimeGapExpected)
        {
            EXPECT_TRUE(TimeGapDetected);
        }
        else
        {
            EXPECT_FALSE(TimeGapDetected);
        }
    }

    EXPECT_EQ(t, NSteps);

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

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
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
