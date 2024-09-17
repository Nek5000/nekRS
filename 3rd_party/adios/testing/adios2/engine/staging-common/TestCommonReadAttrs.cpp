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

    auto ts = std::chrono::steady_clock::now();
    adios2::Engine engine = io.Open(fname, adios2::Mode::Read);
    Seconds timeOpen = std::chrono::steady_clock::now() - ts;
    EXPECT_TRUE(engine);

    const std::string zero = std::to_string(0);
    const std::string s1_Single = std::string("s1_Single_") + zero;
    const std::string s1_Array = std::string("s1_Array_") + zero;
    const std::string i8_Single = std::string("i8_Single_") + zero;
    const std::string i16_Single = std::string("i16_Single_") + zero;
    const std::string i32_Single = std::string("i32_Single_") + zero;
    const std::string i64_Single = std::string("i64_Single_") + zero;
    const std::string u8_Single = std::string("u8_Single_") + zero;
    const std::string u16_Single = std::string("u16_Single_") + zero;
    const std::string u32_Single = std::string("u32_Single_") + zero;
    const std::string u64_Single = std::string("u64_Single_") + zero;
    const std::string r32_Single = std::string("r32_Single_") + zero;
    const std::string r64_Single = std::string("r64_Single_") + zero;

    unsigned int t = 0;

    std::vector<std::time_t> write_times;
    std::vector<Seconds> begin_times;
    std::vector<adios2::StepStatus> begin_statuses;

    while (true)
    {
        ts = std::chrono::steady_clock::now();
        adios2::StepStatus status = engine.BeginStep();
        Seconds timeBeginStep = std::chrono::steady_clock::now() - ts;
        begin_statuses.push_back(status);
        begin_times.push_back(timeBeginStep);

        if (status != adios2::StepStatus::OK)
        {
            break;
        }
        const size_t currentStep = engine.CurrentStep();
        EXPECT_EQ(currentStep, static_cast<size_t>(t));

        size_t writerSize;

        generateCommonTestData((int)0, mpiRank, mpiSize, (int)Nx, (int)Nx);
        auto attr_s1 = io.InquireAttribute<std::string>(s1_Single);
        //        auto attr_s1a = io.InquireAttribute<std::string>(s1_Array);
        auto attr_i8 = io.InquireAttribute<int8_t>(i8_Single);
        auto attr_i16 = io.InquireAttribute<int16_t>(i16_Single);
        auto attr_i32 = io.InquireAttribute<int32_t>(i32_Single);
        auto attr_i64 = io.InquireAttribute<int64_t>(i64_Single);

        auto attr_r32 = io.InquireAttribute<float>(r32_Single);
        auto attr_r64 = io.InquireAttribute<double>(r64_Single);

        EXPECT_TRUE(attr_s1);
        ASSERT_EQ(attr_s1.Name(), s1_Single);
        ASSERT_EQ(attr_s1.Data().size() == 1, true);
        ASSERT_EQ(attr_s1.Type(), adios2::GetType<std::string>());
        ASSERT_EQ(attr_s1.Data().front(), data_S1);

        // EXPECT_TRUE(attr_s1a);
        // ASSERT_EQ(attr_s1a.Name(), s1_Array);
        // ASSERT_EQ(attr_s1a.Data().size() == 1, true);
        // ASSERT_EQ(attr_s1a.Type(), adios2::GetType<std::string>());
        // ASSERT_EQ(attr_s1a.Data()[0], currentTestData.S1array[0]);

        EXPECT_TRUE(attr_i8);
        ASSERT_EQ(attr_i8.Name(), i8_Single);
        ASSERT_EQ(attr_i8.Data().size() == 1, true);
        ASSERT_EQ(attr_i8.Type(), adios2::GetType<int8_t>());
        ASSERT_EQ(attr_i8.Data().front(), data_I8.front());

        EXPECT_TRUE(attr_i16);
        ASSERT_EQ(attr_i16.Name(), i16_Single);
        ASSERT_EQ(attr_i16.Data().size() == 1, true);
        ASSERT_EQ(attr_i16.Type(), adios2::GetType<int16_t>());
        ASSERT_EQ(attr_i16.Data().front(), data_I16.front());

        EXPECT_TRUE(attr_i32);
        ASSERT_EQ(attr_i32.Name(), i32_Single);
        ASSERT_EQ(attr_i32.Data().size() == 1, true);
        ASSERT_EQ(attr_i32.Type(), adios2::GetType<int32_t>());
        ASSERT_EQ(attr_i32.Data().front(), data_I32.front());

        EXPECT_TRUE(attr_i64);
        ASSERT_EQ(attr_i64.Name(), i64_Single);
        ASSERT_EQ(attr_i64.Data().size() == 1, true);
        ASSERT_EQ(attr_i64.Type(), adios2::GetType<int64_t>());
        ASSERT_EQ(attr_i64.Data().front(), data_I64.front());

        EXPECT_TRUE(attr_r32);
        ASSERT_EQ(attr_r32.Name(), r32_Single);
        ASSERT_EQ(attr_r32.Data().size() == 1, true);
        ASSERT_EQ(attr_r32.Type(), adios2::GetType<float>());
        ASSERT_EQ(attr_r32.Data().front(), data_R32.front());

        EXPECT_TRUE(attr_r64);
        ASSERT_EQ(attr_r64.Name(), r64_Single);
        ASSERT_EQ(attr_r64.Data().size() == 1, true);
        ASSERT_EQ(attr_r64.Type(), adios2::GetType<double>());
        if (ModifiableAttributes)
        {
            ASSERT_EQ(attr_r64.Data().front(), (double)3.14159 + (double)t);
        }
        else
        {
            ASSERT_EQ(attr_r64.Data().front(), (double)3.14159);
        }

        auto scalar_r64 = io.InquireVariable<double>("scalar_r64");
        EXPECT_TRUE(scalar_r64);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        /* must be a multiple of Nx */
        ASSERT_EQ(var_i8.Shape()[0] % Nx, 0);

        /* take the first size as something that gives us writer size */
        writerSize = var_i8.Shape()[0] / 10;

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Shape()[0], writerSize * Nx);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Shape()[0], writerSize * Nx);

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Shape()[0], writerSize * Nx);

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Shape()[0], writerSize * Nx);

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Shape()[0], writerSize * Nx);

        auto var_c32 = io.InquireVariable<std::complex<float>>("c32");
        auto var_c64 = io.InquireVariable<std::complex<double>>("c64");
        auto var_r64_2d = io.InquireVariable<double>("r64_2d");
        auto var_r64_2d_rev = io.InquireVariable<double>("r64_2d_rev");
        if (var_c32)
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
        auto var_time = io.InquireVariable<int64_t>("time");
        EXPECT_TRUE(var_time);
        ASSERT_EQ(var_time.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_time.Shape()[0], writerSize);

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
        const adios2::Dims start_time{myStart};
        const adios2::Dims count_time{1};

        const adios2::Box<adios2::Dims> sel(start, count);
        const adios2::Box<adios2::Dims> sel2(start2, count2);
        const adios2::Box<adios2::Dims> sel3(start3, count3);
        const adios2::Box<adios2::Dims> sel_time(start_time, count_time);

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
        engine.Get(var_i8, in_I8.data());
        engine.Get(var_i16, in_I16.data());
        engine.Get(var_i32, in_I32.data());
        engine.Get(var_i64, in_I64.data());

        engine.Get(scalar_r64, in_scalar_R64);

        engine.Get(var_r32, in_R32.data());
        engine.Get(var_r64, in_R64.data());
        if (var_c32)
            engine.Get(var_c32, in_C32.data());
        if (var_c64)
            engine.Get(var_c64, in_C64.data());
        if (var_r64_2d)
            engine.Get(var_r64_2d, in_R64_2d.data());
        if (var_r64_2d_rev)
            engine.Get(var_r64_2d_rev, in_R64_2d_rev.data());
        std::time_t write_time;
        engine.Get(var_time, (int64_t *)&write_time);
        engine.EndStep();

        int result = validateCommonTestData(myStart, myLength, t, !var_c32);
        if (result != 0)
        {
            std::cout << "Read Data Validation failed on node " << mpiRank << " timestep " << t
                      << std::endl;
        }
        EXPECT_EQ(result, 0);
        write_times.push_back(write_time);
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
