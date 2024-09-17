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

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPWriteReadAsStreamTestADIOS2_Threads : public ::testing::Test
{
public:
    BPWriteReadAsStreamTestADIOS2_Threads() = default;

    SmallTestData m_TestData;
};

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

TEST_F(BPWriteReadAsStreamTestADIOS2_Threads, ADIOS2BPWriteRead1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadAsStream_Threads1D8_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadAsStream_Threads1D8.bp");
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
            const adios2::Dims count{Nx};

            io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.SetParameter("Threads", "2");
        io.AddTransport("file");

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();
            bpWriter.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter.Put<float>("r32", currentTestData.R32.data());
            bpWriter.Put<double>("r64", currentTestData.R64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        io.SetParameter("Threads", "2");
        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        std::string IString;
        std::array<int8_t, Nx> I8;
        std::array<int16_t, Nx> I16;
        std::array<int32_t, Nx> I32;
        std::array<int64_t, Nx> I64;
        std::array<uint8_t, Nx> U8;
        std::array<uint16_t, Nx> U16;
        std::array<uint32_t, Nx> U32;
        std::array<uint64_t, Nx> U64;
        std::array<float, Nx> R32;
        std::array<double, Nx> R64;

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        unsigned int t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx);

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx);

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx);

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx);

            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            EXPECT_TRUE(var_u8);
            ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx);

            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            EXPECT_TRUE(var_u16);
            ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx);

            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            EXPECT_TRUE(var_u32);
            ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx);

            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            EXPECT_TRUE(var_u64);
            ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx);

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);

            var_u8.SetSelection(sel);
            var_u16.SetSelection(sel);
            var_u32.SetSelection(sel);
            var_u64.SetSelection(sel);

            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            bpReader.Get(var_i8, I8.data());
            bpReader.Get(var_i16, I16.data());
            bpReader.Get(var_i32, I32.data());
            bpReader.Get(var_i64, I64.data());

            bpReader.Get(var_u8, U8.data());
            bpReader.Get(var_u16, U16.data());
            bpReader.Get(var_u32, U32.data());
            bpReader.Get(var_u64, U64.data());

            bpReader.Get(var_r32, R32.data());
            bpReader.Get(var_r64, R64.data());

            bpReader.EndStep();

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

TEST_F(BPWriteReadAsStreamTestADIOS2_Threads, ADIOS2BPWriteRead2D2x4)
{
    // Each process would write a 2x4 array and all processes would
    // form a 2D 2 * (numberOfProcess*Nx) matrix where Nx is 4 here

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 4;

    // Number of rows
    const std::size_t Ny = 2;

    // Number of steps
    const std::size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadAsStream_Threads2D2x4Test_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadAsStream_Threads2D2x4Test.bp");
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetParameter("Threads", "2");
        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
            const adios2::Dims count{Ny, Nx};

            io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.SetParameter("Threads", "2");
        io.AddTransport("file");

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();
            bpWriter.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter.Put<float>("r32", currentTestData.R32.data());
            bpWriter.Put<double>("r64", currentTestData.R64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        io.SetParameter("Threads", "2");
        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        unsigned int t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i8.Shape()[0], Ny);
            ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i16.Shape()[0], Ny);
            ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], Ny);
            ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i64.Shape()[0], Ny);
            ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            EXPECT_TRUE(var_u8);
            ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u8.Shape()[0], Ny);
            ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            EXPECT_TRUE(var_u16);
            ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u16.Shape()[0], Ny);
            ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            EXPECT_TRUE(var_u32);
            ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u32.Shape()[0], Ny);
            ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            EXPECT_TRUE(var_u64);
            ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u64.Shape()[0], Ny);
            ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], Ny);
            ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Shape()[0], Ny);
            ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);

            var_u8.SetSelection(sel);
            var_u16.SetSelection(sel);
            var_u32.SetSelection(sel);
            var_u64.SetSelection(sel);

            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            bpReader.Get(var_i8, I8.data());
            bpReader.Get(var_i16, I16.data());
            bpReader.Get(var_i32, I32.data());
            bpReader.Get(var_i64, I64.data());

            bpReader.Get(var_u8, U8.data());
            bpReader.Get(var_u16, U16.data());
            bpReader.Get(var_u32, U32.data());
            bpReader.Get(var_u64, U64.data());

            bpReader.Get(var_r32, R32.data());
            bpReader.Get(var_r64, R64.data());

            bpReader.PerformGets();
            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            ++t;
        }
        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

//******************************************************************************
// 2D 4x2 test data
//******************************************************************************

TEST_F(BPWriteReadAsStreamTestADIOS2_Threads, ADIOS2BPWriteRead2D4x2)
{
    // Each process would write a 4x2 array and all processes would
    // form a 2D 4 * (NumberOfProcess * Nx) matrix where Nx is 2 here

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 2;
    // Number of cols
    const std::size_t Ny = 4;

    // Number of steps
    const std::size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadAsStream_Threads2D4x2Test_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadAsStream_Threads2D4x2Test.bp");
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetParameter("Threads", "2");
        // Declare 2D variables (4 * (NumberOfProcess * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            adios2::Dims shape{static_cast<unsigned int>(Ny),
                               static_cast<unsigned int>(mpiSize * Nx)};
            adios2::Dims start{static_cast<unsigned int>(0),
                               static_cast<unsigned int>(mpiRank * Nx)};
            adios2::Dims count{static_cast<unsigned int>(Ny), static_cast<unsigned int>(Nx)};

            io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.SetParameter("Threads", "2");
        io.AddTransport("file");

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();
            bpWriter.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter.Put<float>("r32", currentTestData.R32.data());
            bpWriter.Put<double>("r64", currentTestData.R64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        io.SetParameter("Threads", "2");
        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        // If the size of the array is smaller than the data
        // the result is weird... double and uint64_t would get
        // completely garbage data
        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        unsigned int t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i8.Shape()[0], Ny);
            ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i16.Shape()[0], Ny);
            ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], Ny);
            ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i64.Shape()[0], Ny);
            ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            EXPECT_TRUE(var_u8);
            ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u8.Shape()[0], Ny);
            ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            EXPECT_TRUE(var_u16);
            ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u16.Shape()[0], Ny);
            ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            EXPECT_TRUE(var_u32);
            ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u32.Shape()[0], Ny);
            ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            EXPECT_TRUE(var_u64);
            ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_u64.Shape()[0], Ny);
            ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], Ny);
            ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Shape()[0], Ny);
            ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);

            var_u8.SetSelection(sel);
            var_u16.SetSelection(sel);
            var_u32.SetSelection(sel);
            var_u64.SetSelection(sel);

            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            bpReader.Get(var_i8, I8.data());
            bpReader.Get(var_i16, I16.data());
            bpReader.Get(var_i32, I32.data());
            bpReader.Get(var_i64, I64.data());

            bpReader.Get(var_u8, U8.data());
            bpReader.Get(var_u16, U16.data());
            bpReader.Get(var_u32, U32.data());
            bpReader.Get(var_u64, U64.data());

            bpReader.Get(var_r32, R32.data());
            bpReader.Get(var_r64, R64.data());

            bpReader.PerformGets();

            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            ++t;
        }
        EXPECT_EQ(t, NSteps);
        bpReader.Close();
    }
}

/* This test fails with BP5
TEST_F(BPWriteReadAsStreamTestADIOS2_Threads, ADIOS2BPWriteRead1D8MissingPerformGets)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("ADIOS2BPWriteReadAsStream1D8MissingPerformGetsTest.bp");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
            const adios2::Dims count{Nx};

            io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.SetParameter("Threads", "2");
        io.AddTransport("file");

        // QUESTION: It seems that BPFilterWriter cannot overwrite existing
        // files
        // Ex. if you tune Nx and NSteps, the test would fail. But if you clear
        // the cache in
        // ${adios2Build}/testing/adios2/engine/bp/ADIOS2BPWriteRead1D8.bp.dir,
        // then it works
        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();
            bpWriter.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter.Put<float>("r32", currentTestData.R32.data());
            bpWriter.Put<double>("r64", currentTestData.R64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        io.SetParameter("Threads", "2");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx);

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx);

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx);

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx);

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx);

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx);

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx);

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

        std::string IString;
        std::array<int8_t, Nx> I8;
        std::array<int16_t, Nx> I16;
        std::array<int32_t, Nx> I32;
        std::array<int64_t, Nx> I64;
        std::array<uint8_t, Nx> U8;
        std::array<uint16_t, Nx> U16;
        std::array<uint32_t, Nx> U32;
        std::array<uint64_t, Nx> U64;
        std::array<float, Nx> R32;
        std::array<double, Nx> R64;

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        unsigned int t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            bpReader.Get(var_i8, I8.data());
            bpReader.Get(var_i16, I16.data());
            bpReader.Get(var_i32, I32.data());
            bpReader.Get(var_i64, I64.data());

            bpReader.Get(var_u8, U8.data());
            bpReader.Get(var_u16, U16.data());
            bpReader.Get(var_u32, U32.data());
            bpReader.Get(var_u64, U64.data());

            bpReader.Get(var_r32, R32.data());
            bpReader.Get(var_r64, R64.data());

            bpReader.EndStep();

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}
*/

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

/* This test fails with BP5
TEST_F(BPWriteReadAsStreamTestADIOS2_Threads, ADIOS2BPWriteRead2D2x4MissingPerformGets)
{
    // Each process would write a 2x4 array and all processes would
    // form a 2D 2 * (numberOfProcess*Nx) matrix where Nx is 4 here
    const std::string fname("ADIOS2BPWriteReadAsStream2D2x4MissingPerformGetsTest.bp");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 4;

    // Number of rows
    const std::size_t Ny = 2;

    // Number of steps
    const std::size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
            const adios2::Dims count{Ny, Nx};

            io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.SetParameter("Threads", "2");
        io.AddTransport("file");

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();
            bpWriter.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter.Put<float>("r32", currentTestData.R32.data());
            bpWriter.Put<double>("r64", currentTestData.R64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        io.SetParameter("Threads", "2");
        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], Ny);
        ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], Ny);
        ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], Ny);
        ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], Ny);
        ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], Ny);
        ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], Ny);
        ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], Ny);
        ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], Ny);
        ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], Ny);
        ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        unsigned int t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            bpReader.Get(var_i8, I8.data());
            bpReader.Get(var_i16, I16.data());
            bpReader.Get(var_i32, I32.data());
            bpReader.Get(var_i64, I64.data());

            bpReader.Get(var_u8, U8.data());
            bpReader.Get(var_u16, U16.data());
            bpReader.Get(var_u32, U32.data());
            bpReader.Get(var_u64, U64.data());

            bpReader.Get(var_r32, R32.data());
            bpReader.Get(var_r64, R64.data());

            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            ++t;
        }
        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}
*/

//******************************************************************************
// 2D 4x2 test data
//******************************************************************************

/* This test fails with BP5
TEST_F(BPWriteReadAsStreamTestADIOS2_Threads, ADIOS2BPWriteRead2D4x2MissingPerformGets)
{
    // Each process would write a 4x2 array and all processes would
    // form a 2D 4 * (NumberOfProcess * Nx) matrix where Nx is 2 here
    const std::string fname("ADIOS2BPWriteReadAsStream2D4x2MissingPerformGetsTest.bp");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 2;
    // Number of cols
    const std::size_t Ny = 4;

    // Number of steps
    const std::size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        // Declare 2D variables (4 * (NumberOfProcess * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            adios2::Dims shape{static_cast<unsigned int>(Ny),
                               static_cast<unsigned int>(mpiSize * Nx)};
            adios2::Dims start{static_cast<unsigned int>(0),
                               static_cast<unsigned int>(mpiRank * Nx)};
            adios2::Dims count{static_cast<unsigned int>(Ny), static_cast<unsigned int>(Nx)};

            io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.SetParameter("Threads", "2");
        io.AddTransport("file");

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();
            bpWriter.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter.Put<float>("r32", currentTestData.R32.data());
            bpWriter.Put<double>("r64", currentTestData.R64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        io.SetParameter("Threads", "2");
        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], Ny);
        ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], Ny);
        ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], Ny);
        ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], Ny);
        ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], Ny);
        ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], Ny);
        ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], Ny);
        ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], Ny);
        ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], Ny);
        ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        // If the size of the array is smaller than the data
        // the result is weird... double and uint64_t would get
        // completely garbage data
        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        unsigned int t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            bpReader.Get(var_i8, I8.data());
            bpReader.Get(var_i16, I16.data());
            bpReader.Get(var_i32, I32.data());
            bpReader.Get(var_i64, I64.data());

            bpReader.Get(var_u8, U8.data());
            bpReader.Get(var_u16, U16.data());
            bpReader.Get(var_u32, U32.data());
            bpReader.Get(var_u64, U64.data());

            bpReader.Get(var_r32, R32.data());
            bpReader.Get(var_r64, R64.data());

            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            ++t;
        }
        EXPECT_EQ(t, NSteps);
        bpReader.Close();
    }
}
*/

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);

    if (argc > 1)
    {
        engineName = std::string(argv[1]);
    }
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
