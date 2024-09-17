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

class BPWriteFlushRead : public ::testing::Test
{
public:
    BPWriteFlushRead() = default;

    SmallTestData m_TestData;
    SmallTestData m_OriginalData1D;
    SmallTestData m_OriginalData2D;
};

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

TEST_F(BPWriteFlushRead, ADIOS2BPWrite1D2D)
{
    int mpiRank = 0, mpiSize = 1;

    const size_t Nx1D = 10;

    // Number of rows
    const std::size_t Nx2D = 4;

    // Number of rows
    const std::size_t Ny2D = 2;

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
        adios2::IO io1D = adios.DeclareIO("Flush1D");
        adios2::IO io2D = adios.DeclareIO("Flush2D");
        if (!engineName.empty())
        {
            io1D.SetEngine(engineName);
            io2D.SetEngine(engineName);
        }

        io1D.SetParameter("FlushStepsCount", std::to_string(NSteps + 1));
        io2D.SetParameter("FlushStepsCount", std::to_string(NSteps + 1));

        // io1D Variables
        {
            const adios2::Dims shape{static_cast<size_t>(Nx1D * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx1D * mpiRank)};
            const adios2::Dims count{Nx1D};

            io1D.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        // io2D variables
        {
            const adios2::Dims shape{Ny2D, static_cast<size_t>(Nx2D * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx2D)};
            const adios2::Dims count{Ny2D, Nx2D};

            io2D.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

#if ADIOS2_USE_MPI
        adios2::Engine bpWriter1D = io1D.Open("Flush1D_MPI.bp", adios2::Mode::Write);
        adios2::Engine bpWriter2D = io2D.Open("Flush2D_MPI.bp", adios2::Mode::Write);
#else
        adios2::Engine bpWriter1D = io1D.Open("Flush1D.bp", adios2::Mode::Write);
        adios2::Engine bpWriter2D = io2D.Open("Flush2D.bp", adios2::Mode::Write);
#endif
        for (size_t step = 0; step < NSteps / 2; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter1D.CurrentStep(), step);
            EXPECT_EQ(bpWriter2D.CurrentStep(), step);

            bpWriter1D.BeginStep();
            bpWriter1D.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter1D.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter1D.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter1D.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter1D.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter1D.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter1D.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter1D.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter1D.Put<float>("r32", currentTestData.R32.data());
            bpWriter1D.Put<double>("r64", currentTestData.R64.data());
            bpWriter1D.EndStep();

            bpWriter2D.BeginStep();
            bpWriter2D.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter2D.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter2D.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter2D.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter2D.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter2D.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter2D.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter2D.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter2D.Put<float>("r32", currentTestData.R32.data());
            bpWriter2D.Put<double>("r64", currentTestData.R64.data());
            bpWriter2D.EndStep();
        }

        adios.FlushAll(); // checkpoint for 1D and 2D writers

        // Partial 1D read
        {
            adios2::IO io = adios.DeclareIO("ReadIO1");

            if (!engineName.empty())
            {
                io.SetEngine(engineName);
            }

#if ADIOS2_USE_MPI
            adios2::Engine bpReader = io.Open("Flush1D_MPI.bp", adios2::Mode::Read);
#else
            adios2::Engine bpReader = io.Open("Flush1D.bp", adios2::Mode::Read);
#endif
            unsigned int t = 0;

            while (bpReader.BeginStep(adios2::StepMode::Read, 0.0) == adios2::StepStatus::OK)
            {
                auto var_i8 = io.InquireVariable<int8_t>("i8");
                EXPECT_TRUE(var_i8);
                ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i8.Steps(), NSteps / 2);
                ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx1D);

                auto var_i16 = io.InquireVariable<int16_t>("i16");
                EXPECT_TRUE(var_i16);
                ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i16.Steps(), NSteps / 2);
                ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx1D);

                auto var_i32 = io.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i32.Steps(), NSteps / 2);
                ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx1D);

                auto var_i64 = io.InquireVariable<int64_t>("i64");
                EXPECT_TRUE(var_i64);
                ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i64.Steps(), NSteps / 2);
                ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx1D);

                auto var_u8 = io.InquireVariable<uint8_t>("u8");
                EXPECT_TRUE(var_u8);
                ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u8.Steps(), NSteps / 2);
                ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx1D);

                auto var_u16 = io.InquireVariable<uint16_t>("u16");
                EXPECT_TRUE(var_u16);
                ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u16.Steps(), NSteps / 2);
                ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx1D);

                auto var_u32 = io.InquireVariable<uint32_t>("u32");
                EXPECT_TRUE(var_u32);
                ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u32.Steps(), NSteps / 2);
                ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx1D);

                auto var_u64 = io.InquireVariable<uint64_t>("u64");
                EXPECT_TRUE(var_u64);
                ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u64.Steps(), NSteps / 2);
                ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx1D);

                auto var_r32 = io.InquireVariable<float>("r32");
                EXPECT_TRUE(var_r32);
                ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r32.Steps(), NSteps / 2);
                ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx1D);

                auto var_r64 = io.InquireVariable<double>("r64");
                EXPECT_TRUE(var_r64);
                ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r64.Steps(), NSteps / 2);
                ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx1D);

                std::string IString;
                std::array<int8_t, Nx1D> I8;
                std::array<int16_t, Nx1D> I16;
                std::array<int32_t, Nx1D> I32;
                std::array<int64_t, Nx1D> I64;
                std::array<uint8_t, Nx1D> U8;
                std::array<uint16_t, Nx1D> U16;
                std::array<uint32_t, Nx1D> U32;
                std::array<uint64_t, Nx1D> U64;
                std::array<float, Nx1D> R32;
                std::array<double, Nx1D> R64;

                const adios2::Dims start{mpiRank * Nx1D};
                const adios2::Dims count{Nx1D};

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

                const size_t currentStep = bpReader.CurrentStep();
                EXPECT_EQ(currentStep, static_cast<size_t>(t));

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

                UpdateSmallTestData(m_OriginalData1D, static_cast<int>(t), mpiRank, mpiSize);

                for (size_t i = 0; i < Nx1D; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    EXPECT_EQ(I8[i], m_OriginalData1D.I8[i]) << msg;
                    EXPECT_EQ(I16[i], m_OriginalData1D.I16[i]) << msg;
                    EXPECT_EQ(I32[i], m_OriginalData1D.I32[i]) << msg;
                    EXPECT_EQ(I64[i], m_OriginalData1D.I64[i]) << msg;
                    EXPECT_EQ(U8[i], m_OriginalData1D.U8[i]) << msg;
                    EXPECT_EQ(U16[i], m_OriginalData1D.U16[i]) << msg;
                    EXPECT_EQ(U32[i], m_OriginalData1D.U32[i]) << msg;
                    EXPECT_EQ(U64[i], m_OriginalData1D.U64[i]) << msg;
                    EXPECT_EQ(R32[i], m_OriginalData1D.R32[i]) << msg;
                    EXPECT_EQ(R64[i], m_OriginalData1D.R64[i]) << msg;
                }
                ++t;
            }

            EXPECT_EQ(t, NSteps / 2);

            bpReader.Close();
        }

        // Partial 2D Read
        {
            adios2::IO io = adios.DeclareIO("ReadIO2");

            if (!engineName.empty())
            {
                io.SetEngine(engineName);
            }

#if ADIOS2_USE_MPI
            adios2::Engine bpReader = io.Open("Flush2D_MPI.bp", adios2::Mode::Read);
#else
            adios2::Engine bpReader = io.Open("Flush2D.bp", adios2::Mode::Read);
#endif
            unsigned int t = 0;

            while (bpReader.BeginStep(adios2::StepMode::Read, 0.0) == adios2::StepStatus::OK)
            {
                auto var_i8 = io.InquireVariable<int8_t>("i8");
                EXPECT_TRUE(var_i8);
                ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i8.Steps(), NSteps / 2);
                ASSERT_EQ(var_i8.Shape()[0], Ny2D);
                ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i16 = io.InquireVariable<int16_t>("i16");
                EXPECT_TRUE(var_i16);
                ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i16.Steps(), NSteps / 2);
                ASSERT_EQ(var_i16.Shape()[0], Ny2D);
                ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i32 = io.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i32.Steps(), NSteps / 2);
                ASSERT_EQ(var_i32.Shape()[0], Ny2D);
                ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i64 = io.InquireVariable<int64_t>("i64");
                EXPECT_TRUE(var_i64);
                ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i64.Steps(), NSteps / 2);
                ASSERT_EQ(var_i64.Shape()[0], Ny2D);
                ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u8 = io.InquireVariable<uint8_t>("u8");
                EXPECT_TRUE(var_u8);
                ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u8.Steps(), NSteps / 2);
                ASSERT_EQ(var_u8.Shape()[0], Ny2D);
                ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u16 = io.InquireVariable<uint16_t>("u16");
                EXPECT_TRUE(var_u16);
                ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u16.Steps(), NSteps / 2);
                ASSERT_EQ(var_u16.Shape()[0], Ny2D);
                ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u32 = io.InquireVariable<uint32_t>("u32");
                EXPECT_TRUE(var_u32);
                ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u32.Steps(), NSteps / 2);
                ASSERT_EQ(var_u32.Shape()[0], Ny2D);
                ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u64 = io.InquireVariable<uint64_t>("u64");
                EXPECT_TRUE(var_u64);
                ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u64.Steps(), NSteps / 2);
                ASSERT_EQ(var_u64.Shape()[0], Ny2D);
                ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_r32 = io.InquireVariable<float>("r32");
                EXPECT_TRUE(var_r32);
                ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r32.Steps(), NSteps / 2);
                ASSERT_EQ(var_r32.Shape()[0], Ny2D);
                ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_r64 = io.InquireVariable<double>("r64");
                EXPECT_TRUE(var_r64);
                ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r64.Steps(), NSteps / 2);
                ASSERT_EQ(var_r64.Shape()[0], Ny2D);
                ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                std::array<int8_t, Nx2D * Ny2D> I8;
                std::array<int16_t, Nx2D * Ny2D> I16;
                std::array<int32_t, Nx2D * Ny2D> I32;
                std::array<int64_t, Nx2D * Ny2D> I64;
                std::array<uint8_t, Nx2D * Ny2D> U8;
                std::array<uint16_t, Nx2D * Ny2D> U16;
                std::array<uint32_t, Nx2D * Ny2D> U32;
                std::array<uint64_t, Nx2D * Ny2D> U64;
                std::array<float, Nx2D * Ny2D> R32;
                std::array<double, Nx2D * Ny2D> R64;

                const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx2D)};
                const adios2::Dims count{Ny2D, Nx2D};

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

                const size_t currentStep = bpReader.CurrentStep();
                EXPECT_EQ(currentStep, static_cast<size_t>(t));

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

                // Generate test data for each rank uniquely
                UpdateSmallTestData(m_OriginalData2D, static_cast<int>(t), mpiRank, mpiSize);

                for (size_t i = 0; i < Nx2D * Ny2D; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    EXPECT_EQ(I8[i], m_OriginalData2D.I8[i]) << msg;
                    EXPECT_EQ(I16[i], m_OriginalData2D.I16[i]) << msg;
                    EXPECT_EQ(I32[i], m_OriginalData2D.I32[i]) << msg;
                    EXPECT_EQ(I64[i], m_OriginalData2D.I64[i]) << msg;
                    EXPECT_EQ(U8[i], m_OriginalData2D.U8[i]) << msg;
                    EXPECT_EQ(U16[i], m_OriginalData2D.U16[i]) << msg;
                    EXPECT_EQ(U32[i], m_OriginalData2D.U32[i]) << msg;
                    EXPECT_EQ(U64[i], m_OriginalData2D.U64[i]) << msg;
                    EXPECT_EQ(R32[i], m_OriginalData2D.R32[i]) << msg;
                    EXPECT_EQ(R64[i], m_OriginalData2D.R64[i]) << msg;
                }
                ++t;
            }
            EXPECT_EQ(t, NSteps / 2);

            bpReader.Close();
        }

        bpWriter1D.Close();
        bpWriter2D.Close();
    }
}

TEST_F(BPWriteFlushRead, ADIOS2BPWrite1D2Dstdio)
{
    int mpiRank = 0, mpiSize = 1;

    const size_t Nx1D = 10;

    // Number of rows
    const std::size_t Nx2D = 4;

    // Number of rows
    const std::size_t Ny2D = 2;

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
        adios2::IO io1D = adios.DeclareIO("Flush1D");
        adios2::IO io2D = adios.DeclareIO("Flush2D");

        if (!engineName.empty())
        {
            io1D.SetEngine(engineName);
            io2D.SetEngine(engineName);
        }

        io1D.SetParameter("FlushStepsCount", std::to_string(NSteps + 1));
        io2D.SetParameter("FlushStepsCount", std::to_string(NSteps + 1));

        io1D.AddTransport("File", {{"Library", "stdio"}});
        io2D.AddTransport("File", {{"Library", "stdio"}});

        // io1D Variables
        {
            const adios2::Dims shape{static_cast<size_t>(Nx1D * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx1D * mpiRank)};
            const adios2::Dims count{Nx1D};

            io1D.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        // io2D variables
        {
            const adios2::Dims shape{Ny2D, static_cast<size_t>(Nx2D * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx2D)};
            const adios2::Dims count{Ny2D, Nx2D};

            io2D.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

#if ADIOS2_USE_MPI
        adios2::Engine bpWriter1D = io1D.Open("Flush1Dstdio_MPI.bp", adios2::Mode::Write);
        adios2::Engine bpWriter2D = io2D.Open("Flush2Dstdio_MPI.bp", adios2::Mode::Write);
#else
        adios2::Engine bpWriter1D = io1D.Open("Flush1Dstdio.bp", adios2::Mode::Write);
        adios2::Engine bpWriter2D = io2D.Open("Flush2Dstdio.bp", adios2::Mode::Write);
#endif

        for (size_t step = 0; step < NSteps / 2; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter1D.CurrentStep(), step);
            EXPECT_EQ(bpWriter2D.CurrentStep(), step);

            bpWriter1D.BeginStep();
            bpWriter1D.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter1D.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter1D.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter1D.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter1D.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter1D.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter1D.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter1D.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter1D.Put<float>("r32", currentTestData.R32.data());
            bpWriter1D.Put<double>("r64", currentTestData.R64.data());
            bpWriter1D.EndStep();

            bpWriter2D.BeginStep();
            bpWriter2D.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter2D.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter2D.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter2D.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter2D.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter2D.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter2D.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter2D.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter2D.Put<float>("r32", currentTestData.R32.data());
            bpWriter2D.Put<double>("r64", currentTestData.R64.data());
            bpWriter2D.EndStep();
        }

        adios.FlushAll(); // checkpoint for 1D and 2D writers

        // Partial 1D read
        {
            adios2::IO io = adios.DeclareIO("ReadIO1");

            if (!engineName.empty())
            {
                io.SetEngine(engineName);
            }

#if ADIOS2_USE_MPI
            adios2::Engine bpReader = io.Open("Flush1Dstdio_MPI.bp", adios2::Mode::Read);
#else
            adios2::Engine bpReader = io.Open("Flush1Dstdio.bp", adios2::Mode::Read);
#endif
            unsigned int t = 0;

            while (bpReader.BeginStep(adios2::StepMode::Read, 0.0) == adios2::StepStatus::OK)
            {
                auto var_i8 = io.InquireVariable<int8_t>("i8");
                EXPECT_TRUE(var_i8);
                ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i8.Steps(), NSteps / 2);
                ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx1D);

                auto var_i16 = io.InquireVariable<int16_t>("i16");
                EXPECT_TRUE(var_i16);
                ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i16.Steps(), NSteps / 2);
                ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx1D);

                auto var_i32 = io.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i32.Steps(), NSteps / 2);
                ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx1D);

                auto var_i64 = io.InquireVariable<int64_t>("i64");
                EXPECT_TRUE(var_i64);
                ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i64.Steps(), NSteps / 2);
                ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx1D);

                auto var_u8 = io.InquireVariable<uint8_t>("u8");
                EXPECT_TRUE(var_u8);
                ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u8.Steps(), NSteps / 2);
                ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx1D);

                auto var_u16 = io.InquireVariable<uint16_t>("u16");
                EXPECT_TRUE(var_u16);
                ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u16.Steps(), NSteps / 2);
                ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx1D);

                auto var_u32 = io.InquireVariable<uint32_t>("u32");
                EXPECT_TRUE(var_u32);
                ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u32.Steps(), NSteps / 2);
                ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx1D);

                auto var_u64 = io.InquireVariable<uint64_t>("u64");
                EXPECT_TRUE(var_u64);
                ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u64.Steps(), NSteps / 2);
                ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx1D);

                auto var_r32 = io.InquireVariable<float>("r32");
                EXPECT_TRUE(var_r32);
                ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r32.Steps(), NSteps / 2);
                ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx1D);

                auto var_r64 = io.InquireVariable<double>("r64");
                EXPECT_TRUE(var_r64);
                ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r64.Steps(), NSteps / 2);
                ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx1D);

                std::string IString;
                std::array<int8_t, Nx1D> I8;
                std::array<int16_t, Nx1D> I16;
                std::array<int32_t, Nx1D> I32;
                std::array<int64_t, Nx1D> I64;
                std::array<uint8_t, Nx1D> U8;
                std::array<uint16_t, Nx1D> U16;
                std::array<uint32_t, Nx1D> U32;
                std::array<uint64_t, Nx1D> U64;
                std::array<float, Nx1D> R32;
                std::array<double, Nx1D> R64;

                const adios2::Dims start{mpiRank * Nx1D};
                const adios2::Dims count{Nx1D};

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

                const size_t currentStep = bpReader.CurrentStep();
                EXPECT_EQ(currentStep, static_cast<size_t>(t));

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

                UpdateSmallTestData(m_OriginalData1D, static_cast<int>(t), mpiRank, mpiSize);

                for (size_t i = 0; i < Nx1D; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    EXPECT_EQ(I8[i], m_OriginalData1D.I8[i]) << msg;
                    EXPECT_EQ(I16[i], m_OriginalData1D.I16[i]) << msg;
                    EXPECT_EQ(I32[i], m_OriginalData1D.I32[i]) << msg;
                    EXPECT_EQ(I64[i], m_OriginalData1D.I64[i]) << msg;
                    EXPECT_EQ(U8[i], m_OriginalData1D.U8[i]) << msg;
                    EXPECT_EQ(U16[i], m_OriginalData1D.U16[i]) << msg;
                    EXPECT_EQ(U32[i], m_OriginalData1D.U32[i]) << msg;
                    EXPECT_EQ(U64[i], m_OriginalData1D.U64[i]) << msg;
                    EXPECT_EQ(R32[i], m_OriginalData1D.R32[i]) << msg;
                    EXPECT_EQ(R64[i], m_OriginalData1D.R64[i]) << msg;
                }
                ++t;
            }

            EXPECT_EQ(t, NSteps / 2);

            bpReader.Close();
        }

        // Partial 2D Read
        {
            adios2::IO io = adios.DeclareIO("ReadIO2");

            if (!engineName.empty())
            {
                io.SetEngine(engineName);
            }

#if ADIOS2_USE_MPI
            adios2::Engine bpReader = io.Open("Flush2Dstdio_MPI.bp", adios2::Mode::Read);
#else
            adios2::Engine bpReader = io.Open("Flush2Dstdio.bp", adios2::Mode::Read);
#endif

            unsigned int t = 0;

            while (bpReader.BeginStep(adios2::StepMode::Read, 0.0) == adios2::StepStatus::OK)
            {
                auto var_i8 = io.InquireVariable<int8_t>("i8");
                EXPECT_TRUE(var_i8);
                ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i8.Steps(), NSteps / 2);
                ASSERT_EQ(var_i8.Shape()[0], Ny2D);
                ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i16 = io.InquireVariable<int16_t>("i16");
                EXPECT_TRUE(var_i16);
                ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i16.Steps(), NSteps / 2);
                ASSERT_EQ(var_i16.Shape()[0], Ny2D);
                ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i32 = io.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i32.Steps(), NSteps / 2);
                ASSERT_EQ(var_i32.Shape()[0], Ny2D);
                ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i64 = io.InquireVariable<int64_t>("i64");
                EXPECT_TRUE(var_i64);
                ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i64.Steps(), NSteps / 2);
                ASSERT_EQ(var_i64.Shape()[0], Ny2D);
                ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u8 = io.InquireVariable<uint8_t>("u8");
                EXPECT_TRUE(var_u8);
                ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u8.Steps(), NSteps / 2);
                ASSERT_EQ(var_u8.Shape()[0], Ny2D);
                ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u16 = io.InquireVariable<uint16_t>("u16");
                EXPECT_TRUE(var_u16);
                ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u16.Steps(), NSteps / 2);
                ASSERT_EQ(var_u16.Shape()[0], Ny2D);
                ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u32 = io.InquireVariable<uint32_t>("u32");
                EXPECT_TRUE(var_u32);
                ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u32.Steps(), NSteps / 2);
                ASSERT_EQ(var_u32.Shape()[0], Ny2D);
                ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u64 = io.InquireVariable<uint64_t>("u64");
                EXPECT_TRUE(var_u64);
                ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u64.Steps(), NSteps / 2);
                ASSERT_EQ(var_u64.Shape()[0], Ny2D);
                ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_r32 = io.InquireVariable<float>("r32");
                EXPECT_TRUE(var_r32);
                ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r32.Steps(), NSteps / 2);
                ASSERT_EQ(var_r32.Shape()[0], Ny2D);
                ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_r64 = io.InquireVariable<double>("r64");
                EXPECT_TRUE(var_r64);
                ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r64.Steps(), NSteps / 2);
                ASSERT_EQ(var_r64.Shape()[0], Ny2D);
                ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                std::array<int8_t, Nx2D * Ny2D> I8;
                std::array<int16_t, Nx2D * Ny2D> I16;
                std::array<int32_t, Nx2D * Ny2D> I32;
                std::array<int64_t, Nx2D * Ny2D> I64;
                std::array<uint8_t, Nx2D * Ny2D> U8;
                std::array<uint16_t, Nx2D * Ny2D> U16;
                std::array<uint32_t, Nx2D * Ny2D> U32;
                std::array<uint64_t, Nx2D * Ny2D> U64;
                std::array<float, Nx2D * Ny2D> R32;
                std::array<double, Nx2D * Ny2D> R64;

                const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx2D)};
                const adios2::Dims count{Ny2D, Nx2D};

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

                const size_t currentStep = bpReader.CurrentStep();
                EXPECT_EQ(currentStep, static_cast<size_t>(t));

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

                // Generate test data for each rank uniquely
                UpdateSmallTestData(m_OriginalData2D, static_cast<int>(t), mpiRank, mpiSize);

                for (size_t i = 0; i < Nx2D * Ny2D; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    EXPECT_EQ(I8[i], m_OriginalData2D.I8[i]) << msg;
                    EXPECT_EQ(I16[i], m_OriginalData2D.I16[i]) << msg;
                    EXPECT_EQ(I32[i], m_OriginalData2D.I32[i]) << msg;
                    EXPECT_EQ(I64[i], m_OriginalData2D.I64[i]) << msg;
                    EXPECT_EQ(U8[i], m_OriginalData2D.U8[i]) << msg;
                    EXPECT_EQ(U16[i], m_OriginalData2D.U16[i]) << msg;
                    EXPECT_EQ(U32[i], m_OriginalData2D.U32[i]) << msg;
                    EXPECT_EQ(U64[i], m_OriginalData2D.U64[i]) << msg;
                    EXPECT_EQ(R32[i], m_OriginalData2D.R32[i]) << msg;
                    EXPECT_EQ(R64[i], m_OriginalData2D.R64[i]) << msg;
                }
                ++t;
            }
            EXPECT_EQ(t, NSteps / 2);

            bpReader.Close();
        }

        bpWriter1D.Close();
        bpWriter2D.Close();
    }
}

TEST_F(BPWriteFlushRead, ADIOS2BPWrite1D2Dfstream)
{
    int mpiRank = 0, mpiSize = 1;

    const size_t Nx1D = 10;

    // Number of rows
    const std::size_t Nx2D = 4;

    // Number of rows
    const std::size_t Ny2D = 2;

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
        adios2::IO io1D = adios.DeclareIO("Flush1D");
        adios2::IO io2D = adios.DeclareIO("Flush2D");

        if (!engineName.empty())
        {
            io1D.SetEngine(engineName);
            io2D.SetEngine(engineName);
        }

        io1D.SetParameter("FlushStepsCount", std::to_string(NSteps + 1));
        io2D.SetParameter("FlushStepsCount", std::to_string(NSteps + 1));

        io1D.AddTransport("File", {{"Library", "fstream"}});
        io2D.AddTransport("File", {{"Library", "fstream"}});

        // io1D Variables
        {
            const adios2::Dims shape{static_cast<size_t>(Nx1D * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx1D * mpiRank)};
            const adios2::Dims count{Nx1D};

            io1D.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io1D.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io1D.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

        // io2D variables
        {
            const adios2::Dims shape{Ny2D, static_cast<size_t>(Nx2D * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx2D)};
            const adios2::Dims count{Ny2D, Nx2D};

            io2D.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);

            io2D.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
            io2D.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        }

#if ADIOS2_USE_MPI
        adios2::Engine bpWriter1D = io1D.Open("Flush1Dfstream_MPI.bp", adios2::Mode::Write);
        adios2::Engine bpWriter2D = io2D.Open("Flush2Dfstream_MPI.bp", adios2::Mode::Write);
#else
        adios2::Engine bpWriter1D = io1D.Open("Flush1Dfstream.bp", adios2::Mode::Write);
        adios2::Engine bpWriter2D = io2D.Open("Flush2Dfstream.bp", adios2::Mode::Write);
#endif

        for (size_t step = 0; step < NSteps / 2; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter1D.CurrentStep(), step);
            EXPECT_EQ(bpWriter2D.CurrentStep(), step);

            bpWriter1D.BeginStep();
            bpWriter1D.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter1D.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter1D.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter1D.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter1D.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter1D.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter1D.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter1D.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter1D.Put<float>("r32", currentTestData.R32.data());
            bpWriter1D.Put<double>("r64", currentTestData.R64.data());
            bpWriter1D.EndStep();

            bpWriter2D.BeginStep();
            bpWriter2D.Put<int8_t>("i8", currentTestData.I8.data());
            bpWriter2D.Put<int16_t>("i16", currentTestData.I16.data());
            bpWriter2D.Put<int32_t>("i32", currentTestData.I32.data());
            bpWriter2D.Put<int64_t>("i64", currentTestData.I64.data());
            bpWriter2D.Put<uint8_t>("u8", currentTestData.U8.data());
            bpWriter2D.Put<uint16_t>("u16", currentTestData.U16.data());
            bpWriter2D.Put<uint32_t>("u32", currentTestData.U32.data());
            bpWriter2D.Put<uint64_t>("u64", currentTestData.U64.data());
            bpWriter2D.Put<float>("r32", currentTestData.R32.data());
            bpWriter2D.Put<double>("r64", currentTestData.R64.data());
            bpWriter2D.EndStep();
        }

        adios.FlushAll(); // checkpoint for 1D and 2D writers

        // Partial 1D read
        {
            adios2::IO io = adios.DeclareIO("ReadIO1");

            if (!engineName.empty())
            {
                io.SetEngine(engineName);
            }

#if ADIOS2_USE_MPI
            adios2::Engine bpReader = io.Open("Flush1Dfstream_MPI.bp", adios2::Mode::Read);
#else
            adios2::Engine bpReader = io.Open("Flush1Dfstream.bp", adios2::Mode::Read);
#endif

            unsigned int t = 0;

            while (bpReader.BeginStep(adios2::StepMode::Read, 0.0) == adios2::StepStatus::OK)
            {
                auto var_i8 = io.InquireVariable<int8_t>("i8");
                EXPECT_TRUE(var_i8);
                ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i8.Steps(), NSteps / 2);
                ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx1D);

                auto var_i16 = io.InquireVariable<int16_t>("i16");
                EXPECT_TRUE(var_i16);
                ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i16.Steps(), NSteps / 2);
                ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx1D);

                auto var_i32 = io.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i32.Steps(), NSteps / 2);
                ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx1D);

                auto var_i64 = io.InquireVariable<int64_t>("i64");
                EXPECT_TRUE(var_i64);
                ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i64.Steps(), NSteps / 2);
                ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx1D);

                auto var_u8 = io.InquireVariable<uint8_t>("u8");
                EXPECT_TRUE(var_u8);
                ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u8.Steps(), NSteps / 2);
                ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx1D);

                auto var_u16 = io.InquireVariable<uint16_t>("u16");
                EXPECT_TRUE(var_u16);
                ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u16.Steps(), NSteps / 2);
                ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx1D);

                auto var_u32 = io.InquireVariable<uint32_t>("u32");
                EXPECT_TRUE(var_u32);
                ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u32.Steps(), NSteps / 2);
                ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx1D);

                auto var_u64 = io.InquireVariable<uint64_t>("u64");
                EXPECT_TRUE(var_u64);
                ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u64.Steps(), NSteps / 2);
                ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx1D);

                auto var_r32 = io.InquireVariable<float>("r32");
                EXPECT_TRUE(var_r32);
                ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r32.Steps(), NSteps / 2);
                ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx1D);

                auto var_r64 = io.InquireVariable<double>("r64");
                EXPECT_TRUE(var_r64);
                ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r64.Steps(), NSteps / 2);
                ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx1D);

                std::string IString;
                std::array<int8_t, Nx1D> I8;
                std::array<int16_t, Nx1D> I16;
                std::array<int32_t, Nx1D> I32;
                std::array<int64_t, Nx1D> I64;
                std::array<uint8_t, Nx1D> U8;
                std::array<uint16_t, Nx1D> U16;
                std::array<uint32_t, Nx1D> U32;
                std::array<uint64_t, Nx1D> U64;
                std::array<float, Nx1D> R32;
                std::array<double, Nx1D> R64;

                const adios2::Dims start{mpiRank * Nx1D};
                const adios2::Dims count{Nx1D};

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

                const size_t currentStep = bpReader.CurrentStep();
                EXPECT_EQ(currentStep, static_cast<size_t>(t));

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

                UpdateSmallTestData(m_OriginalData1D, static_cast<int>(t), mpiRank, mpiSize);

                for (size_t i = 0; i < Nx1D; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    EXPECT_EQ(I8[i], m_OriginalData1D.I8[i]) << msg;
                    EXPECT_EQ(I16[i], m_OriginalData1D.I16[i]) << msg;
                    EXPECT_EQ(I32[i], m_OriginalData1D.I32[i]) << msg;
                    EXPECT_EQ(I64[i], m_OriginalData1D.I64[i]) << msg;
                    EXPECT_EQ(U8[i], m_OriginalData1D.U8[i]) << msg;
                    EXPECT_EQ(U16[i], m_OriginalData1D.U16[i]) << msg;
                    EXPECT_EQ(U32[i], m_OriginalData1D.U32[i]) << msg;
                    EXPECT_EQ(U64[i], m_OriginalData1D.U64[i]) << msg;
                    EXPECT_EQ(R32[i], m_OriginalData1D.R32[i]) << msg;
                    EXPECT_EQ(R64[i], m_OriginalData1D.R64[i]) << msg;
                }
                ++t;
            }

            EXPECT_EQ(t, NSteps / 2);

            bpReader.Close();
        }

        // Partial 2D Read
        {
            adios2::IO io = adios.DeclareIO("ReadIO2");

            if (!engineName.empty())
            {
                io.SetEngine(engineName);
            }

#if ADIOS2_USE_MPI
            adios2::Engine bpReader = io.Open("Flush2Dfstream_MPI.bp", adios2::Mode::Read);
#else
            adios2::Engine bpReader = io.Open("Flush2Dfstream.bp", adios2::Mode::Read);
#endif

            unsigned int t = 0;

            while (bpReader.BeginStep(adios2::StepMode::Read, 0.0) == adios2::StepStatus::OK)
            {
                auto var_i8 = io.InquireVariable<int8_t>("i8");
                EXPECT_TRUE(var_i8);
                ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i8.Steps(), NSteps / 2);
                ASSERT_EQ(var_i8.Shape()[0], Ny2D);
                ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i16 = io.InquireVariable<int16_t>("i16");
                EXPECT_TRUE(var_i16);
                ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i16.Steps(), NSteps / 2);
                ASSERT_EQ(var_i16.Shape()[0], Ny2D);
                ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i32 = io.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i32.Steps(), NSteps / 2);
                ASSERT_EQ(var_i32.Shape()[0], Ny2D);
                ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_i64 = io.InquireVariable<int64_t>("i64");
                EXPECT_TRUE(var_i64);
                ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i64.Steps(), NSteps / 2);
                ASSERT_EQ(var_i64.Shape()[0], Ny2D);
                ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u8 = io.InquireVariable<uint8_t>("u8");
                EXPECT_TRUE(var_u8);
                ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u8.Steps(), NSteps / 2);
                ASSERT_EQ(var_u8.Shape()[0], Ny2D);
                ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u16 = io.InquireVariable<uint16_t>("u16");
                EXPECT_TRUE(var_u16);
                ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u16.Steps(), NSteps / 2);
                ASSERT_EQ(var_u16.Shape()[0], Ny2D);
                ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u32 = io.InquireVariable<uint32_t>("u32");
                EXPECT_TRUE(var_u32);
                ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u32.Steps(), NSteps / 2);
                ASSERT_EQ(var_u32.Shape()[0], Ny2D);
                ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_u64 = io.InquireVariable<uint64_t>("u64");
                EXPECT_TRUE(var_u64);
                ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u64.Steps(), NSteps / 2);
                ASSERT_EQ(var_u64.Shape()[0], Ny2D);
                ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_r32 = io.InquireVariable<float>("r32");
                EXPECT_TRUE(var_r32);
                ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r32.Steps(), NSteps / 2);
                ASSERT_EQ(var_r32.Shape()[0], Ny2D);
                ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                auto var_r64 = io.InquireVariable<double>("r64");
                EXPECT_TRUE(var_r64);
                ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r64.Steps(), NSteps / 2);
                ASSERT_EQ(var_r64.Shape()[0], Ny2D);
                ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx2D));

                std::array<int8_t, Nx2D * Ny2D> I8;
                std::array<int16_t, Nx2D * Ny2D> I16;
                std::array<int32_t, Nx2D * Ny2D> I32;
                std::array<int64_t, Nx2D * Ny2D> I64;
                std::array<uint8_t, Nx2D * Ny2D> U8;
                std::array<uint16_t, Nx2D * Ny2D> U16;
                std::array<uint32_t, Nx2D * Ny2D> U32;
                std::array<uint64_t, Nx2D * Ny2D> U64;
                std::array<float, Nx2D * Ny2D> R32;
                std::array<double, Nx2D * Ny2D> R64;

                const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx2D)};
                const adios2::Dims count{Ny2D, Nx2D};

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

                const size_t currentStep = bpReader.CurrentStep();
                EXPECT_EQ(currentStep, static_cast<size_t>(t));

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

                // Generate test data for each rank uniquely
                UpdateSmallTestData(m_OriginalData2D, static_cast<int>(t), mpiRank, mpiSize);

                for (size_t i = 0; i < Nx2D * Ny2D; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    EXPECT_EQ(I8[i], m_OriginalData2D.I8[i]) << msg;
                    EXPECT_EQ(I16[i], m_OriginalData2D.I16[i]) << msg;
                    EXPECT_EQ(I32[i], m_OriginalData2D.I32[i]) << msg;
                    EXPECT_EQ(I64[i], m_OriginalData2D.I64[i]) << msg;
                    EXPECT_EQ(U8[i], m_OriginalData2D.U8[i]) << msg;
                    EXPECT_EQ(U16[i], m_OriginalData2D.U16[i]) << msg;
                    EXPECT_EQ(U32[i], m_OriginalData2D.U32[i]) << msg;
                    EXPECT_EQ(U64[i], m_OriginalData2D.U64[i]) << msg;
                    EXPECT_EQ(R32[i], m_OriginalData2D.R32[i]) << msg;
                    EXPECT_EQ(R64[i], m_OriginalData2D.R64[i]) << msg;
                }
                ++t;
            }
            EXPECT_EQ(t, NSteps / 2);

            bpReader.Close();
        }

        bpWriter1D.Close();
        bpWriter2D.Close();
    }
}

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
