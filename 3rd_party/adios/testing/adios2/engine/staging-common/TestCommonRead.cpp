/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <thread>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

class CommonReadTest : public ::testing::Test
{
public:
    CommonReadTest() = default;
};

typedef std::chrono::duration<double> Seconds;

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

    adios2::Engine engine;

    auto ts = std::chrono::steady_clock::now();
    if (ExpectOpenTimeout)
    {
        ASSERT_ANY_THROW(engine = io.Open(fname, adios2::Mode::Read));
        return;
    }
    else
    {
        engine = io.Open(fname, adios2::Mode::Read);
    }
    Seconds timeOpen = std::chrono::steady_clock::now() - ts;
    EXPECT_TRUE(engine);

    size_t t = 0;

    std::vector<std::time_t> write_times;
    std::vector<Seconds> begin_times;
    std::vector<adios2::StepStatus> begin_statuses;

    while (true)
    {
        if (EarlyExit && (t == static_cast<size_t>(NSteps)))
        {
            break;
        }
        ts = std::chrono::steady_clock::now();
        adios2::StepStatus status = engine.BeginStep();
        auto av = io.AvailableVariables();

        Seconds timeBeginStep = std::chrono::steady_clock::now() - ts;
        begin_statuses.push_back(status);
        begin_times.push_back(timeBeginStep);

        if (status != adios2::StepStatus::OK)
        {
            break;
        }
        const size_t currentStep = engine.CurrentStep();
        EXPECT_EQ(currentStep, t);

        size_t writerSize = 0;

        auto scalar_r64 = io.InquireVariable<double>("scalar_r64");
        if (!NoData)
            EXPECT_TRUE(scalar_r64);
        else
            EXPECT_FALSE(scalar_r64);

        auto var_time = io.InquireVariable<int64_t>("time");
        if (NoData)
        {
            EXPECT_FALSE(var_time);
        }
        else
        {
            EXPECT_TRUE(var_time);
            ASSERT_EQ(var_time.ShapeID(), adios2::ShapeID::GlobalArray);
            writerSize = var_time.Shape()[0];
        }

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        if (NoData)
        {
            EXPECT_FALSE(var_time);
        }
        else
        {
            EXPECT_TRUE(var_i8);
            ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            /* take the first size as something that gives us writer size */
            Nx = var_i8.Shape()[0] / writerSize;
        }

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        if (NoData)
        {
            EXPECT_FALSE(var_time);
        }
        else
        {
            EXPECT_TRUE(var_i16);
            ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i16.Shape()[0], writerSize * Nx);
        }

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        if (NoData)
        {
            EXPECT_FALSE(var_time);
        }
        else
        {
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], writerSize * Nx);
        }

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        if (NoData)
        {
            EXPECT_FALSE(var_time);
        }
        else
        {
            EXPECT_TRUE(var_i64);
            ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i64.Shape()[0], writerSize * Nx);
        }

        auto var_r32 = io.InquireVariable<float>("r32");
        if (NoData)
        {
            EXPECT_FALSE(var_time);
        }
        else
        {
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], writerSize * Nx);
        }

        auto var_r64 = io.InquireVariable<double>("r64");
        if (NoData)
        {
            EXPECT_FALSE(var_time);
        }
        else
        {
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Shape()[0], writerSize * Nx);
        }

        auto var_c32 = io.InquireVariable<std::complex<float>>("c32");
        auto var_c64 = io.InquireVariable<std::complex<double>>("c64");
        auto var_r64_2d = io.InquireVariable<double>("r64_2d");
        auto var_r64_2d_rev = io.InquireVariable<double>("r64_2d_rev");
        if (var_c32 && !NoData)
        {
            EXPECT_TRUE(var_c32);
            ASSERT_EQ(var_c32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_c32.Shape()[0], writerSize * Nx);

            EXPECT_TRUE(var_c64);
            ASSERT_EQ(var_c64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_c64.Shape()[0], writerSize * Nx);

            EXPECT_TRUE(var_r64_2d);
            ASSERT_EQ(var_r64_2d.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64_2d.Shape()[0], writerSize * Nx);
            ASSERT_EQ(var_r64_2d.Shape()[1], 2);

            EXPECT_TRUE(var_r64_2d_rev);
            ASSERT_EQ(var_r64_2d_rev.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64_2d_rev.Shape()[0], 2);
            ASSERT_EQ(var_r64_2d_rev.Shape()[1], writerSize * Nx);
        }
        else
        {
            EXPECT_FALSE(var_c32);
            EXPECT_FALSE(var_c64);
            EXPECT_FALSE(var_r64_2d);
            EXPECT_FALSE(var_r64_2d_rev);
        }

        if (!NoData)
        {
            const std::vector<adios2::Variable<int8_t>::Info> i8Info =
                engine.BlocksInfo(var_i8, engine.CurrentStep());
            const std::vector<adios2::Variable<int16_t>::Info> i16Info =
                engine.BlocksInfo(var_i16, engine.CurrentStep());
            const std::vector<adios2::Variable<int32_t>::Info> i32Info =
                engine.BlocksInfo(var_i32, engine.CurrentStep());
            const std::vector<adios2::Variable<int64_t>::Info> i64Info =
                engine.BlocksInfo(var_i64, engine.CurrentStep());
            const std::vector<adios2::Variable<float>::Info> r32Info =
                engine.BlocksInfo(var_r32, engine.CurrentStep());
            const std::vector<adios2::Variable<double>::Info> r64Info =
                engine.BlocksInfo(var_r64, engine.CurrentStep());
            EXPECT_EQ(i8Info.size(), writerSize);
            EXPECT_EQ(i16Info.size(), writerSize);
            EXPECT_TRUE((i32Info.size() == writerSize) || (i32Info.size() == writerSize * 3));
            EXPECT_EQ(i64Info.size(), writerSize);
            EXPECT_EQ(r32Info.size(), writerSize);
            EXPECT_EQ(r64Info.size(), writerSize);
            for (size_t i = 0; i < writerSize; ++i)
            {
                EXPECT_FALSE(i8Info[0].IsValue);
                EXPECT_FALSE(i16Info[0].IsValue);
                EXPECT_FALSE(i32Info[0].IsValue);
                EXPECT_FALSE(i64Info[0].IsValue);
                EXPECT_FALSE(r32Info[0].IsValue);
                EXPECT_FALSE(r64Info[0].IsValue);
            }
        }

        if (var_c32)
        {
            const std::vector<adios2::Variable<std::complex<float>>::Info> c32Info =
                engine.BlocksInfo(var_c32, engine.CurrentStep());
            const std::vector<adios2::Variable<std::complex<double>>::Info> c64Info =
                engine.BlocksInfo(var_c64, engine.CurrentStep());
            EXPECT_EQ(c32Info.size(), writerSize);
            EXPECT_EQ(c64Info.size(), writerSize);
            for (size_t i = 0; i < writerSize; ++i)
            {
                EXPECT_FALSE(c32Info[0].IsValue);
                EXPECT_FALSE(c64Info[0].IsValue);
            }
        }

        long unsigned int myStart = (long unsigned int)(writerSize * Nx / mpiSize) * mpiRank;
        long unsigned int myLength = (long unsigned int)((writerSize * Nx + mpiSize - 1) / mpiSize);

        if (myStart + myLength > writerSize * Nx)
        {
            myLength = (long unsigned int)writerSize * (int)Nx - myStart;
        }
        const adios2::Dims start{myStart};
        const adios2::Dims count{myLength};
        const adios2::Dims start2{myStart, 0};
        const adios2::Dims count2{myLength, 2};
        const adios2::Dims start3{0, myStart};
        const adios2::Dims count3{2, myLength};
        const adios2::Dims start_time{0};
        const adios2::Dims count_time{1};

        const adios2::Box<adios2::Dims> sel(start, count);
        const adios2::Box<adios2::Dims> sel2(start2, count2);
        const adios2::Box<adios2::Dims> sel3(start3, count3);
        const adios2::Box<adios2::Dims> sel_time(start_time, count_time);
        std::time_t write_time;

        if (!NoData)
        {
            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);

            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);
            if (var_c32)
                var_c32.SetSelection(sel);
            if (var_c64)
                var_c64.SetSelection(sel);
            if (var_r64_2d)
                var_r64_2d.SetSelection(sel2);
            if (var_r64_2d_rev)
                var_r64_2d_rev.SetSelection(sel3);

            var_time.SetSelection(sel_time);

            in_I8.resize(myLength);
            in_I16.resize(myLength);
            in_I32.resize(myLength);
            in_I64.resize(myLength);
            in_R32.resize(myLength);
            in_R64.resize(myLength);
            in_C32.resize(myLength);
            in_C64.resize(myLength);
            in_R64_2d.resize(myLength * 2);
            in_R64_2d_rev.resize(myLength * 2);
            if (!NoData)
            {
                engine.Get(var_i8, in_I8.data(), GlobalReadMode);
                engine.Get(var_i16, in_I16.data(), GlobalReadMode);
                engine.Get(var_i32, in_I32.data(), GlobalReadMode);
                engine.Get(var_i64, in_I64.data(), GlobalReadMode);

                engine.Get(scalar_r64, in_scalar_R64, GlobalReadMode);

                engine.Get(var_r32, in_R32.data(), GlobalReadMode);
                engine.Get(var_r64, in_R64.data(), GlobalReadMode);
                if (!mpiRank)
                    engine.Get(var_time, (int64_t *)&write_time, GlobalReadMode);
            }
            else
            {
                if (NoDataNode != -1)
                {
                    // someone wrote everything, get our part
                    engine.Get(var_r64, in_R64.data(), GlobalReadMode);
                }
            }
            if (var_c32)
                engine.Get(var_c32, in_C32.data(), GlobalReadMode);
            if (var_c64)
                engine.Get(var_c64, in_C64.data(), GlobalReadMode);
            if (var_r64_2d)
                engine.Get(var_r64_2d, in_R64_2d.data(), GlobalReadMode);
            if (var_r64_2d_rev)
                engine.Get(var_r64_2d_rev, in_R64_2d_rev.data(), GlobalReadMode);
            if (LockGeometry)
            {
                // we'll never change our data decomposition
                engine.LockReaderSelections();
            }
        }
        if (!NoData && (GlobalReadMode == adios2::Mode::Sync))
        {
            // go ahead and test data now, it should be valid
            int result = validateCommonTestData(myStart, myLength, t, !var_c32);
            if (result != 0)
            {
                std::cout << "Read Data Validation failed on node " << mpiRank << " timestep " << t
                          << std::endl;
            }
            EXPECT_EQ(result, 0);
        }
        engine.EndStep();

        if (!NoData)
        {
            int result = validateCommonTestData(myStart, myLength, t, !var_c32);
            if (result != 0)
            {
                std::cout << "Read Data Validation failed on node " << mpiRank << " timestep " << t
                          << std::endl;
            }
            EXPECT_EQ(result, 0);
            if (AdvancingAttrs)
            {
                /* we only succeed if every attribute from every prior step is
                 * there, but not the next few */
                for (size_t step = 0; step <= currentStep + 2; step++)
                {
                    const std::string r64_Single =
                        std::string("r64_PerStep_") + std::to_string(step);
                    auto attr_r64 = io.InquireAttribute<double>(r64_Single);
                    std::cout << "Testing for attribute " << r64_Single << std::endl;
                    if (step <= currentStep)
                    {
                        EXPECT_TRUE(attr_r64);
                        ASSERT_EQ(attr_r64.Data().size() == 1, true);
                        ASSERT_EQ(attr_r64.Type(), adios2::GetType<double>());
                        ASSERT_EQ(attr_r64.Data().front(), (double)(step * 10.0));
                    }
                    else
                    {
                        // The file engines let attributes appear early, so only
                        // enforce non-appearance if that changes.
                        // EXPECT_FALSE(attr_r64);
                    }
                }
            }
            if (!mpiRank)
                write_times.push_back(write_time);
        }
        else
        {
            if (NoDataNode != -1)
            {
                EXPECT_EQ(validateCommonTestDataR64(myStart, myLength, t, !var_c32), 0);
            }
        }
        ++t;
    }

    if (!mpiRank)
    {
        std::cout << "Reader Open took " << std::fixed << std::setprecision(9) << timeOpen.count()
                  << " seconds" << std::endl;
        for (size_t i = 0; i < begin_times.size(); ++i)
        {
            std::cout << "Reader BeginStep t = " << i << " had status = " << begin_statuses[i]
                      << " after " << std::fixed << std::setprecision(9) << begin_times[i].count()
                      << " seconds" << std::endl;
        }
    }

    EXPECT_EQ(t, static_cast<size_t>(NSteps));
    if (!NoData && !mpiRank)
    {
        if ((write_times.back() - write_times.front()) > 1)
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
    }

    // Close the file
    if (!DontClose)
    {
        engine.Close();
    }
    if (TestVarDestruction)
    {
        // Engine close should delete was was created by the reader, not other
        // vars
        EXPECT_FALSE(io.InquireVariable<std::complex<float>>("c32"));
    }
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
