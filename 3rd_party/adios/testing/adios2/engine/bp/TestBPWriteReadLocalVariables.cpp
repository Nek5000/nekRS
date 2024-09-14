/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <numeric>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPWriteReadLocalVariables : public ::testing::Test
{
public:
    BPWriteReadLocalVariables() = default;

    SmallTestData m_TestData;
};

TEST_F(BPWriteReadLocalVariables, ADIOS2BPWriteReadLocal1D)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadLocal1D_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadLocal1D.bp");
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{};
            const adios2::Dims start{};
            const adios2::Dims count{Nx};

            io.DefineVariable<int32_t>("stepsGlobalValue");
            io.DefineVariable<std::string>("stepsGlobalValueString");

            io.DefineVariable<int32_t>("ranksLocalValue", {adios2::LocalValueDim});
            io.DefineVariable<std::string>("ranksLocalValueString", {adios2::LocalValueDim});

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

            io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                   adios2::ConstantDims);
            io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                    adios2::ConstantDims);
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();

            const int32_t step32 = static_cast<int32_t>(step);
            bpWriter.Put<int32_t>("stepsGlobalValue", step32);
            bpWriter.Put<std::string>("stepsGlobalValueString", std::to_string(step));

            bpWriter.Put<int32_t>("ranksLocalValue", mpiRank);
            bpWriter.Put<std::string>("ranksLocalValueString", std::to_string(mpiRank));

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
            bpWriter.Put<std::complex<float>>("cr32", currentTestData.CR32.data());
            bpWriter.Put<std::complex<double>>("cr64", currentTestData.CR64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    // if (mpiRank == 0)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

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
        std::array<std::complex<float>, Nx> CR32;
        std::array<std::complex<double>, Nx> CR64;

        size_t t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            auto var_StepsGlobalValue = io.InquireVariable<int32_t>("stepsGlobalValue");
            auto var_StepsGlobalValueString =
                io.InquireVariable<std::string>("stepsGlobalValueString");
            auto var_RanksLocalValue = io.InquireVariable<int32_t>("ranksLocalValue");
            auto var_RanksLocalValueString =
                io.InquireVariable<std::string>("ranksLocalValueString");

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            auto var_i16 = io.InquireVariable<int16_t>("i16");
            auto var_i32 = io.InquireVariable<int32_t>("i32");
            auto var_i64 = io.InquireVariable<int64_t>("i64");
            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            auto var_r32 = io.InquireVariable<float>("r32");
            auto var_r64 = io.InquireVariable<double>("r64");
            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");

            // Global value
            EXPECT_TRUE(var_StepsGlobalValue);
            EXPECT_EQ(var_StepsGlobalValue.ShapeID(), adios2::ShapeID::GlobalValue);
            EXPECT_EQ(var_StepsGlobalValue.Shape().size(), 0);
            EXPECT_EQ(var_StepsGlobalValue.Min(), static_cast<int32_t>(currentStep));
            EXPECT_EQ(var_StepsGlobalValue.Max(), static_cast<int32_t>(currentStep));
            int32_t stepsGlobalValueData;
            bpReader.Get(var_StepsGlobalValue, stepsGlobalValueData, adios2::Mode::Sync);
            EXPECT_EQ(stepsGlobalValueData, static_cast<int32_t>(currentStep));

            EXPECT_TRUE(var_StepsGlobalValueString);
            EXPECT_EQ(var_StepsGlobalValueString.ShapeID(), adios2::ShapeID::GlobalValue);
            EXPECT_EQ(var_StepsGlobalValueString.Shape().size(), 0);
            std::string stepsGlobalValueStringDataString;
            bpReader.Get(var_StepsGlobalValueString, stepsGlobalValueStringDataString,
                         adios2::Mode::Sync);
            EXPECT_EQ(stepsGlobalValueStringDataString, std::to_string(currentStep));

            // Local values will be read as 1D GlobalArrays
            EXPECT_TRUE(var_RanksLocalValue);
            EXPECT_EQ(var_RanksLocalValue.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_RanksLocalValue.Shape().size(), 1);
            EXPECT_EQ(var_RanksLocalValue.Shape()[0], mpiSize);
            EXPECT_EQ(var_RanksLocalValue.Min(), 0);
            EXPECT_EQ(var_RanksLocalValue.Max(), mpiSize - 1);
            std::vector<int32_t> rankLocalValueData;
            bpReader.Get(var_RanksLocalValue, rankLocalValueData);
            EXPECT_EQ(rankLocalValueData.size(), mpiSize);
            for (size_t r = 0; r < rankLocalValueData.size(); ++r)
            {
                EXPECT_EQ(rankLocalValueData[r], static_cast<int32_t>(r));
            }

            if ((engineName != "BP3") && (engineName != "BP4") && (mpiSize > 1))
            {
                // This test exposes a difference in behaviour between
                // BP3/4 and BP5.  BP5 resets the variables to their
                // default state at each BeginStep where prior engines
                // do not.  Therefore setting the block selection in
                // this portion causes failures in future timesteps in
                // prior engines, but not in BP5.  So we'll just test
                // BP5 for the moment.
                for (size_t r = 0; r < rankLocalValueData.size(); ++r)
                {
                    int32_t val;
                    var_RanksLocalValue.SetBlockSelection(r);
                    bpReader.Get(var_RanksLocalValue, &val);
                    EXPECT_EQ(val, static_cast<int32_t>(r));
                }
            }

            bpReader.Get(var_RanksLocalValue, rankLocalValueData);

            EXPECT_TRUE(var_RanksLocalValueString);
            EXPECT_EQ(var_RanksLocalValueString.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_RanksLocalValueString.Shape().size(), 1);
            EXPECT_EQ(var_RanksLocalValueString.Shape()[0], mpiSize);
            std::vector<std::string> rankLocalValueDataString;
            bpReader.Get(var_RanksLocalValueString, rankLocalValueDataString, adios2::Mode::Sync);
            EXPECT_EQ(rankLocalValueDataString.size(), mpiSize);
            for (size_t r = 0; r < rankLocalValueDataString.size(); ++r)
            {
                EXPECT_EQ(rankLocalValueDataString[r], std::to_string(r));
            }

            EXPECT_TRUE(var_i8);
            EXPECT_TRUE(var_i16);
            EXPECT_TRUE(var_i32);
            EXPECT_TRUE(var_i64);
            EXPECT_TRUE(var_u8);
            EXPECT_TRUE(var_u16);
            EXPECT_TRUE(var_u32);
            EXPECT_TRUE(var_u64);
            EXPECT_TRUE(var_cr32);
            EXPECT_TRUE(var_cr64);

            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i8.Shape().size(), 0);
            EXPECT_EQ(var_i8.Start().size(), 0);
            EXPECT_EQ(var_i8.Count()[0], Nx);
            //
            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i16.Shape().size(), 0);
            EXPECT_EQ(var_i16.Start().size(), 0);
            EXPECT_EQ(var_i16.Count()[0], Nx);
            //
            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i32.Shape().size(), 0);
            EXPECT_EQ(var_i32.Start().size(), 0);
            EXPECT_EQ(var_i32.Count()[0], Nx);

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i64.Shape().size(), 0);
            EXPECT_EQ(var_i64.Start().size(), 0);
            EXPECT_EQ(var_i64.Count()[0], Nx);

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r32.Shape().size(), 0);
            EXPECT_EQ(var_r32.Start().size(), 0);
            EXPECT_EQ(var_r32.Count()[0], Nx);

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r64.Shape().size(), 0);
            EXPECT_EQ(var_r64.Start().size(), 0);
            EXPECT_EQ(var_r64.Count()[0], Nx);

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr32.Shape().size(), 0);
            EXPECT_EQ(var_cr32.Start().size(), 0);
            EXPECT_EQ(var_cr32.Count()[0], Nx);

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr64.Shape().size(), 0);
            EXPECT_EQ(var_cr64.Start().size(), 0);
            EXPECT_EQ(var_cr64.Count()[0], Nx);

            for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
            {
                var_i8.SetBlockSelection(b);
                var_i16.SetBlockSelection(b);
                var_i32.SetBlockSelection(b);
                var_i64.SetBlockSelection(b);

                var_u8.SetBlockSelection(b);
                var_u16.SetBlockSelection(b);
                var_u32.SetBlockSelection(b);
                var_u64.SetBlockSelection(b);

                var_r32.SetBlockSelection(b);
                var_r64.SetBlockSelection(b);

                var_cr32.SetBlockSelection(b);
                var_cr64.SetBlockSelection(b);

                EXPECT_EQ(var_i8.BlockID(), b);
                EXPECT_EQ(var_i16.BlockID(), b);
                EXPECT_EQ(var_i32.BlockID(), b);
                EXPECT_EQ(var_i64.BlockID(), b);
                EXPECT_EQ(var_u8.BlockID(), b);
                EXPECT_EQ(var_u16.BlockID(), b);
                EXPECT_EQ(var_u32.BlockID(), b);
                EXPECT_EQ(var_u64.BlockID(), b);
                EXPECT_EQ(var_r32.BlockID(), b);
                EXPECT_EQ(var_r64.BlockID(), b);
                EXPECT_EQ(var_cr32.BlockID(), b);
                EXPECT_EQ(var_cr64.BlockID(), b);

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

                bpReader.Get(var_cr32, CR32.data());
                bpReader.Get(var_cr64, CR64.data());

                bpReader.PerformGets();

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                for (size_t i = 0; i < Nx; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    if (var_i8)
                    {
                        EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                    }
                    if (var_i16)
                    {
                        EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                    }
                    if (var_i32)
                    {
                        EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                    }
                    if (var_i64)
                    {
                        EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                    }

                    if (var_u8)
                    {
                        EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                    }
                    if (var_u16)
                    {
                        EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                    }
                    if (var_u32)
                    {
                        EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                    }
                    if (var_u64)
                    {
                        EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                    }
                    if (var_r32)
                    {
                        EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                    }
                    if (var_r64)
                    {
                        EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
                    }

                    if (var_cr32)
                    {
                        EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                    }
                    if (var_cr64)
                    {
                        EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
                    }
                }
            } // block loop

            bpReader.EndStep();
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariables, ADIOS2BPWriteReadLocal2D2x4)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    const size_t Nx = 4;
    const size_t Ny = 2;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadLocal2D2x4_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadLocal2D2x4.bp");
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{};
            const adios2::Dims start{};
            const adios2::Dims count{Ny, Nx};

            io.DefineVariable<int32_t>("stepsGlobalValue");
            io.DefineVariable<std::string>("stepsGlobalValueString");

            io.DefineVariable<int32_t>("ranksLocalValue", {adios2::LocalValueDim});
            io.DefineVariable<std::string>("ranksLocalValueString", {adios2::LocalValueDim});

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

            io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                   adios2::ConstantDims);
            io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                    adios2::ConstantDims);
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();

            const int32_t step32 = static_cast<int32_t>(step);
            bpWriter.Put<int32_t>("stepsGlobalValue", step32);
            bpWriter.Put<std::string>("stepsGlobalValueString", std::to_string(step));

            bpWriter.Put<int32_t>("ranksLocalValue", mpiRank);
            bpWriter.Put<std::string>("ranksLocalValueString", std::to_string(mpiRank));

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
            bpWriter.Put<std::complex<float>>("cr32", currentTestData.CR32.data());
            bpWriter.Put<std::complex<double>>("cr64", currentTestData.CR64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    // if (mpiRank == 0)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        std::string IString;
        std::array<int8_t, Ny * Nx> I8;
        std::array<int16_t, Ny * Nx> I16;
        std::array<int32_t, Ny * Nx> I32;
        std::array<int64_t, Ny * Nx> I64;
        std::array<uint8_t, Ny * Nx> U8;
        std::array<uint16_t, Ny * Nx> U16;
        std::array<uint32_t, Ny * Nx> U32;
        std::array<uint64_t, Ny * Nx> U64;
        std::array<float, Ny * Nx> R32;
        std::array<double, Ny * Nx> R64;
        std::array<std::complex<float>, Ny * Nx> CR32;
        std::array<std::complex<double>, Ny * Nx> CR64;

        size_t t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            auto var_StepsGlobalValue = io.InquireVariable<int32_t>("stepsGlobalValue");
            auto var_StepsGlobalValueString =
                io.InquireVariable<std::string>("stepsGlobalValueString");
            auto var_RanksLocalValue = io.InquireVariable<int32_t>("ranksLocalValue");
            auto var_RanksLocalValueString =
                io.InquireVariable<std::string>("ranksLocalValueString");

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            auto var_i16 = io.InquireVariable<int16_t>("i16");
            auto var_i32 = io.InquireVariable<int32_t>("i32");
            auto var_i64 = io.InquireVariable<int64_t>("i64");
            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            auto var_r32 = io.InquireVariable<float>("r32");
            auto var_r64 = io.InquireVariable<double>("r64");
            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");

            // Global value
            EXPECT_TRUE(var_StepsGlobalValue);
            EXPECT_EQ(var_StepsGlobalValue.ShapeID(), adios2::ShapeID::GlobalValue);
            EXPECT_EQ(var_StepsGlobalValue.Shape().size(), 0);
            EXPECT_EQ(var_StepsGlobalValue.Min(), static_cast<int32_t>(currentStep));
            EXPECT_EQ(var_StepsGlobalValue.Max(), static_cast<int32_t>(currentStep));
            int32_t stepsGlobalValueData;
            bpReader.Get(var_StepsGlobalValue, stepsGlobalValueData, adios2::Mode::Sync);
            EXPECT_EQ(stepsGlobalValueData, static_cast<int32_t>(currentStep));

            EXPECT_TRUE(var_StepsGlobalValueString);
            EXPECT_EQ(var_StepsGlobalValueString.ShapeID(), adios2::ShapeID::GlobalValue);
            EXPECT_EQ(var_StepsGlobalValueString.Shape().size(), 0);
            std::string stepsGlobalValueStringDataString;
            bpReader.Get(var_StepsGlobalValueString, stepsGlobalValueStringDataString,
                         adios2::Mode::Sync);
            EXPECT_EQ(stepsGlobalValueStringDataString, std::to_string(currentStep));

            // Local values are read as 1D Global Arrays
            EXPECT_TRUE(var_RanksLocalValue);
            EXPECT_EQ(var_RanksLocalValue.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_RanksLocalValue.Shape()[0], mpiSize);
            EXPECT_EQ(var_RanksLocalValue.Min(), 0);
            EXPECT_EQ(var_RanksLocalValue.Max(), mpiSize - 1);
            std::vector<int32_t> rankLocalValueData;
            bpReader.Get(var_RanksLocalValue, rankLocalValueData);
            EXPECT_EQ(rankLocalValueData.size(), mpiSize);
            for (size_t r = 0; r < rankLocalValueData.size(); ++r)
            {
                EXPECT_EQ(rankLocalValueData[r], static_cast<int32_t>(r));
            }

            EXPECT_TRUE(var_RanksLocalValueString);
            EXPECT_EQ(var_RanksLocalValueString.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_RanksLocalValue.Shape()[0], mpiSize);
            std::vector<std::string> rankLocalValueDataString;
            bpReader.Get(var_RanksLocalValueString, rankLocalValueDataString, adios2::Mode::Sync);
            EXPECT_EQ(rankLocalValueData.size(), mpiSize);
            for (size_t r = 0; r < rankLocalValueData.size(); ++r)
            {
                EXPECT_EQ(rankLocalValueDataString[r], std::to_string(r));
            }

            EXPECT_TRUE(var_i8);
            EXPECT_TRUE(var_i16);
            EXPECT_TRUE(var_i32);
            EXPECT_TRUE(var_i64);
            EXPECT_TRUE(var_u8);
            EXPECT_TRUE(var_u16);
            EXPECT_TRUE(var_u32);
            EXPECT_TRUE(var_u64);
            EXPECT_TRUE(var_cr32);
            EXPECT_TRUE(var_cr64);

            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i8.Shape().size(), 0);
            EXPECT_EQ(var_i8.Start().size(), 0);
            EXPECT_EQ(var_i8.Count()[0], Ny);
            EXPECT_EQ(var_i8.Count()[1], Nx);
            //
            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i16.Shape().size(), 0);
            EXPECT_EQ(var_i16.Start().size(), 0);
            EXPECT_EQ(var_i16.Count()[0], Ny);
            EXPECT_EQ(var_i16.Count()[1], Nx);
            //
            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i32.Shape().size(), 0);
            EXPECT_EQ(var_i32.Start().size(), 0);
            EXPECT_EQ(var_i32.Count()[0], Ny);
            EXPECT_EQ(var_i32.Count()[1], Nx);

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i64.Shape().size(), 0);
            EXPECT_EQ(var_i64.Start().size(), 0);
            EXPECT_EQ(var_i64.Count()[0], Ny);
            EXPECT_EQ(var_i64.Count()[1], Nx);

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r32.Shape().size(), 0);
            EXPECT_EQ(var_r32.Start().size(), 0);
            EXPECT_EQ(var_r32.Count()[0], Ny);
            EXPECT_EQ(var_r32.Count()[1], Nx);

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r64.Shape().size(), 0);
            EXPECT_EQ(var_r64.Start().size(), 0);
            EXPECT_EQ(var_r64.Count()[0], Ny);
            EXPECT_EQ(var_r64.Count()[1], Nx);

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr32.Shape().size(), 0);
            EXPECT_EQ(var_cr32.Start().size(), 0);
            EXPECT_EQ(var_cr32.Count()[0], Ny);
            EXPECT_EQ(var_cr32.Count()[1], Nx);

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr64.Shape().size(), 0);
            EXPECT_EQ(var_cr64.Start().size(), 0);
            EXPECT_EQ(var_cr64.Count()[0], Ny);
            EXPECT_EQ(var_cr64.Count()[1], Nx);

            for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
            {
                var_i8.SetBlockSelection(b);
                var_i16.SetBlockSelection(b);
                var_i32.SetBlockSelection(b);
                var_i64.SetBlockSelection(b);

                var_u8.SetBlockSelection(b);
                var_u16.SetBlockSelection(b);
                var_u32.SetBlockSelection(b);
                var_u64.SetBlockSelection(b);

                var_r32.SetBlockSelection(b);
                var_r64.SetBlockSelection(b);

                var_cr32.SetBlockSelection(b);
                var_cr64.SetBlockSelection(b);

                EXPECT_EQ(var_i8.BlockID(), b);
                EXPECT_EQ(var_i16.BlockID(), b);
                EXPECT_EQ(var_i32.BlockID(), b);
                EXPECT_EQ(var_i64.BlockID(), b);
                EXPECT_EQ(var_u8.BlockID(), b);
                EXPECT_EQ(var_u16.BlockID(), b);
                EXPECT_EQ(var_u32.BlockID(), b);
                EXPECT_EQ(var_u64.BlockID(), b);
                EXPECT_EQ(var_r32.BlockID(), b);
                EXPECT_EQ(var_r64.BlockID(), b);
                EXPECT_EQ(var_cr32.BlockID(), b);
                EXPECT_EQ(var_cr64.BlockID(), b);

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

                bpReader.Get(var_cr32, CR32.data());
                bpReader.Get(var_cr64, CR64.data());

                bpReader.PerformGets();

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                for (size_t i = 0; i < Nx * Ny; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    if (var_i8)
                    {
                        EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                    }
                    if (var_i16)
                    {
                        EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                    }
                    if (var_i32)
                    {
                        EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                    }
                    if (var_i64)
                    {
                        EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                    }

                    if (var_u8)
                    {
                        EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                    }
                    if (var_u16)
                    {
                        EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                    }
                    if (var_u32)
                    {
                        EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                    }
                    if (var_u64)
                    {
                        EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                    }
                    if (var_r32)
                    {
                        EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                    }
                    if (var_r64)
                    {
                        EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
                    }

                    if (var_cr32)
                    {
                        EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                    }
                    if (var_cr64)
                    {
                        EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
                    }
                }
            } // block loop

            bpReader.EndStep();
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariables, ADIOS2BPWriteReadLocal2D4x2)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    const size_t Nx = 2;
    const size_t Ny = 4;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadLocal2D4x2_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadLocal2D4x2.bp");
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{};
            const adios2::Dims start{};
            const adios2::Dims count{Ny, Nx};

            io.DefineVariable<int32_t>("stepsGlobalValue");
            io.DefineVariable<std::string>("stepsGlobalValueString");

            io.DefineVariable<int32_t>("ranksLocalValue", {adios2::LocalValueDim});
            io.DefineVariable<std::string>("ranksLocalValueString", {adios2::LocalValueDim});

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

            io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                   adios2::ConstantDims);
            io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                    adios2::ConstantDims);
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();
            const int32_t step32 = static_cast<int32_t>(step);
            bpWriter.Put<int32_t>("stepsGlobalValue", step32);
            bpWriter.Put<std::string>("stepsGlobalValueString", std::to_string(step));

            bpWriter.Put<int32_t>("ranksLocalValue", mpiRank);
            bpWriter.Put<std::string>("ranksLocalValueString", std::to_string(mpiRank));

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
            bpWriter.Put<std::complex<float>>("cr32", currentTestData.CR32.data());
            bpWriter.Put<std::complex<double>>("cr64", currentTestData.CR64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    // if (mpiRank == 0)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        std::string IString;
        std::array<int8_t, Ny * Nx> I8;
        std::array<int16_t, Ny * Nx> I16;
        std::array<int32_t, Ny * Nx> I32;
        std::array<int64_t, Ny * Nx> I64;
        std::array<uint8_t, Ny * Nx> U8;
        std::array<uint16_t, Ny * Nx> U16;
        std::array<uint32_t, Ny * Nx> U32;
        std::array<uint64_t, Ny * Nx> U64;
        std::array<float, Ny * Nx> R32;
        std::array<double, Ny * Nx> R64;
        std::array<std::complex<float>, Ny * Nx> CR32;
        std::array<std::complex<double>, Ny * Nx> CR64;

        size_t t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            auto var_StepsGlobalValue = io.InquireVariable<int32_t>("stepsGlobalValue");
            auto var_StepsGlobalValueString =
                io.InquireVariable<std::string>("stepsGlobalValueString");
            auto var_RanksLocalValue = io.InquireVariable<int32_t>("ranksLocalValue");
            auto var_RanksLocalValueString =
                io.InquireVariable<std::string>("ranksLocalValueString");

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            auto var_i16 = io.InquireVariable<int16_t>("i16");
            auto var_i32 = io.InquireVariable<int32_t>("i32");
            auto var_i64 = io.InquireVariable<int64_t>("i64");
            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            auto var_r32 = io.InquireVariable<float>("r32");
            auto var_r64 = io.InquireVariable<double>("r64");
            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");

            // Global value
            EXPECT_TRUE(var_StepsGlobalValue);
            EXPECT_EQ(var_StepsGlobalValue.ShapeID(), adios2::ShapeID::GlobalValue);
            EXPECT_EQ(var_StepsGlobalValue.Shape().size(), 0);
            EXPECT_EQ(var_StepsGlobalValue.Min(), static_cast<int32_t>(currentStep));
            EXPECT_EQ(var_StepsGlobalValue.Max(), static_cast<int32_t>(currentStep));
            int32_t stepsGlobalValueData;
            bpReader.Get(var_StepsGlobalValue, stepsGlobalValueData, adios2::Mode::Sync);
            EXPECT_EQ(stepsGlobalValueData, static_cast<int32_t>(currentStep));

            EXPECT_TRUE(var_StepsGlobalValueString);
            EXPECT_EQ(var_StepsGlobalValueString.ShapeID(), adios2::ShapeID::GlobalValue);
            EXPECT_EQ(var_StepsGlobalValueString.Shape().size(), 0);
            std::string stepsGlobalValueStringDataString;
            bpReader.Get(var_StepsGlobalValueString, stepsGlobalValueStringDataString,
                         adios2::Mode::Sync);
            EXPECT_EQ(stepsGlobalValueStringDataString, std::to_string(currentStep));

            // Local values are read as 1D Global Arrays
            EXPECT_TRUE(var_RanksLocalValue);
            EXPECT_EQ(var_RanksLocalValue.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_RanksLocalValue.Shape().size(), 1);
            EXPECT_EQ(var_RanksLocalValue.Shape()[0], mpiSize);
            EXPECT_EQ(var_RanksLocalValue.Min(), 0);
            EXPECT_EQ(var_RanksLocalValue.Max(), mpiSize - 1);
            std::vector<int32_t> rankLocalValueData;
            bpReader.Get(var_RanksLocalValue, rankLocalValueData);
            EXPECT_EQ(rankLocalValueData.size(), mpiSize);
            for (size_t r = 0; r < rankLocalValueData.size(); ++r)
            {
                EXPECT_EQ(rankLocalValueData[r], static_cast<int32_t>(r));
            }

            EXPECT_TRUE(var_RanksLocalValueString);
            EXPECT_EQ(var_RanksLocalValueString.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_RanksLocalValue.Shape().size(), 1);
            EXPECT_EQ(var_RanksLocalValue.Shape()[0], mpiSize);
            std::vector<std::string> rankLocalValueDataString;
            bpReader.Get(var_RanksLocalValueString, rankLocalValueDataString, adios2::Mode::Sync);
            EXPECT_EQ(rankLocalValueData.size(), mpiSize);
            for (size_t r = 0; r < rankLocalValueData.size(); ++r)
            {
                EXPECT_EQ(rankLocalValueDataString[r], std::to_string(r));
            }

            EXPECT_TRUE(var_i8);
            EXPECT_TRUE(var_i16);
            EXPECT_TRUE(var_i32);
            EXPECT_TRUE(var_i64);
            EXPECT_TRUE(var_u8);
            EXPECT_TRUE(var_u16);
            EXPECT_TRUE(var_u32);
            EXPECT_TRUE(var_u64);
            EXPECT_TRUE(var_cr32);
            EXPECT_TRUE(var_cr64);

            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i8.Shape().size(), 0);
            EXPECT_EQ(var_i8.Start().size(), 0);
            EXPECT_EQ(var_i8.Count()[0], Ny);
            EXPECT_EQ(var_i8.Count()[1], Nx);
            //
            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i16.Shape().size(), 0);
            EXPECT_EQ(var_i16.Start().size(), 0);
            EXPECT_EQ(var_i16.Count()[0], Ny);
            EXPECT_EQ(var_i16.Count()[1], Nx);
            //
            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i32.Shape().size(), 0);
            EXPECT_EQ(var_i32.Start().size(), 0);
            EXPECT_EQ(var_i32.Count()[0], Ny);
            EXPECT_EQ(var_i32.Count()[1], Nx);

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i64.Shape().size(), 0);
            EXPECT_EQ(var_i64.Start().size(), 0);
            EXPECT_EQ(var_i64.Count()[0], Ny);
            EXPECT_EQ(var_i64.Count()[1], Nx);

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r32.Shape().size(), 0);
            EXPECT_EQ(var_r32.Start().size(), 0);
            EXPECT_EQ(var_r32.Count()[0], Ny);
            EXPECT_EQ(var_r32.Count()[1], Nx);

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r64.Shape().size(), 0);
            EXPECT_EQ(var_r64.Start().size(), 0);
            EXPECT_EQ(var_r64.Count()[0], Ny);
            EXPECT_EQ(var_r64.Count()[1], Nx);

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr32.Shape().size(), 0);
            EXPECT_EQ(var_cr32.Start().size(), 0);
            EXPECT_EQ(var_cr32.Count()[0], Ny);
            EXPECT_EQ(var_cr32.Count()[1], Nx);

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr64.Shape().size(), 0);
            EXPECT_EQ(var_cr64.Start().size(), 0);
            EXPECT_EQ(var_cr64.Count()[0], Ny);
            EXPECT_EQ(var_cr64.Count()[1], Nx);

            for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
            {
                var_i8.SetBlockSelection(b);
                var_i16.SetBlockSelection(b);
                var_i32.SetBlockSelection(b);
                var_i64.SetBlockSelection(b);

                var_u8.SetBlockSelection(b);
                var_u16.SetBlockSelection(b);
                var_u32.SetBlockSelection(b);
                var_u64.SetBlockSelection(b);

                var_r32.SetBlockSelection(b);
                var_r64.SetBlockSelection(b);

                var_cr32.SetBlockSelection(b);
                var_cr64.SetBlockSelection(b);

                EXPECT_EQ(var_i8.BlockID(), b);
                EXPECT_EQ(var_i16.BlockID(), b);
                EXPECT_EQ(var_i32.BlockID(), b);
                EXPECT_EQ(var_i64.BlockID(), b);
                EXPECT_EQ(var_u8.BlockID(), b);
                EXPECT_EQ(var_u16.BlockID(), b);
                EXPECT_EQ(var_u32.BlockID(), b);
                EXPECT_EQ(var_u64.BlockID(), b);
                EXPECT_EQ(var_r32.BlockID(), b);
                EXPECT_EQ(var_r64.BlockID(), b);
                EXPECT_EQ(var_cr32.BlockID(), b);
                EXPECT_EQ(var_cr64.BlockID(), b);

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

                bpReader.Get(var_cr32, CR32.data());
                bpReader.Get(var_cr64, CR64.data());

                bpReader.PerformGets();

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                for (size_t i = 0; i < Nx * Ny; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    if (var_i8)
                    {
                        EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                    }
                    if (var_i16)
                    {
                        EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                    }
                    if (var_i32)
                    {
                        EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                    }
                    if (var_i64)
                    {
                        EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                    }

                    if (var_u8)
                    {
                        EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                    }
                    if (var_u16)
                    {
                        EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                    }
                    if (var_u32)
                    {
                        EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                    }
                    if (var_u64)
                    {
                        EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                    }
                    if (var_r32)
                    {
                        EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                    }
                    if (var_r64)
                    {
                        EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
                    }

                    if (var_cr32)
                    {
                        EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                    }
                    if (var_cr64)
                    {
                        EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
                    }
                }
            } // block loop

            bpReader.EndStep();
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariables, ADIOS2BPWriteReadLocal1DAllSteps)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadLocal1DAllSteps_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadLocal1DAllSteps.bp");
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{};
            const adios2::Dims start{};
            const adios2::Dims count{Nx};

            io.DefineVariable<int32_t>("stepsGlobalValue");
            io.DefineVariable<std::string>("stepsGlobalValueString");

            io.DefineVariable<int32_t>("ranksLocalValue", {adios2::LocalValueDim});
            io.DefineVariable<std::string>("ranksLocalValueString", {adios2::LocalValueDim});

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

            io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                   adios2::ConstantDims);
            io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                    adios2::ConstantDims);
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();

            const int32_t step32 = static_cast<int32_t>(step);
            bpWriter.Put<int32_t>("stepsGlobalValue", step32);
            bpWriter.Put<std::string>("stepsGlobalValueString", std::to_string(step));

            bpWriter.Put<int32_t>("ranksLocalValue", mpiRank);
            bpWriter.Put<std::string>("ranksLocalValueString", std::to_string(mpiRank));

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
            bpWriter.Put<std::complex<float>>("cr32", currentTestData.CR32.data());
            bpWriter.Put<std::complex<double>>("cr64", currentTestData.CR64.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    // if (mpiRank == 0)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        std::vector<int8_t> I8;
        std::vector<int16_t> I16;
        std::vector<int32_t> I32;
        std::vector<int64_t> I64;
        std::vector<uint8_t> U8;
        std::vector<uint16_t> U16;
        std::vector<uint32_t> U32;
        std::vector<uint64_t> U64;
        std::vector<float> R32;
        std::vector<double> R64;
        std::vector<std::complex<float>> CR32;
        std::vector<std::complex<double>> CR64;

        size_t t = 0;

        auto var_StepsGlobalValue = io.InquireVariable<int32_t>("stepsGlobalValue");
        auto var_StepsGlobalValueString = io.InquireVariable<std::string>("stepsGlobalValueString");
        auto var_RanksLocalValue = io.InquireVariable<int32_t>("ranksLocalValue");
        auto var_RanksLocalValueString = io.InquireVariable<std::string>("ranksLocalValueString");

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        auto var_i16 = io.InquireVariable<int16_t>("i16");
        auto var_i32 = io.InquireVariable<int32_t>("i32");
        auto var_i64 = io.InquireVariable<int64_t>("i64");
        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        auto var_r32 = io.InquireVariable<float>("r32");
        auto var_r64 = io.InquireVariable<double>("r64");
        auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
        auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");

        EXPECT_EQ(var_StepsGlobalValue.Min(), 0);
        EXPECT_EQ(var_StepsGlobalValue.Max(), static_cast<int32_t>(NSteps - 1));

        // Read all steps at once
        var_StepsGlobalValue.SetStepSelection({0, var_StepsGlobalValue.Steps()});
        var_StepsGlobalValueString.SetStepSelection({0, var_StepsGlobalValueString.Steps()});

        var_RanksLocalValue.SetStepSelection({0, var_RanksLocalValue.Steps()});
        var_RanksLocalValueString.SetStepSelection({0, var_RanksLocalValueString.Steps()});

        var_i8.SetStepSelection({0, var_i8.Steps()});
        var_i16.SetStepSelection({0, var_i16.Steps()});
        var_i32.SetStepSelection({0, var_i32.Steps()});
        var_i64.SetStepSelection({0, var_i64.Steps()});

        var_u8.SetStepSelection({0, var_u8.Steps()});
        var_u16.SetStepSelection({0, var_u16.Steps()});
        var_u32.SetStepSelection({0, var_u32.Steps()});
        var_u64.SetStepSelection({0, var_u64.Steps()});

        var_r32.SetStepSelection({0, var_r32.Steps()});
        var_r64.SetStepSelection({0, var_r64.Steps()});
        var_cr32.SetStepSelection({0, var_cr32.Steps()});
        var_cr64.SetStepSelection({0, var_cr64.Steps()});

        for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
        {
            var_i8.SetBlockSelection(b);
            var_i16.SetBlockSelection(b);
            var_i32.SetBlockSelection(b);
            var_i64.SetBlockSelection(b);

            var_u8.SetBlockSelection(b);
            var_u16.SetBlockSelection(b);
            var_u32.SetBlockSelection(b);
            var_u64.SetBlockSelection(b);

            var_r32.SetBlockSelection(b);
            var_r64.SetBlockSelection(b);

            var_cr32.SetBlockSelection(b);
            var_cr64.SetBlockSelection(b);

            bpReader.Get(var_i8, I8);
            bpReader.Get(var_i16, I16);
            bpReader.Get(var_i32, I32);
            bpReader.Get(var_i64, I64);

            bpReader.Get(var_u8, U8);
            bpReader.Get(var_u16, U16);
            bpReader.Get(var_u32, U32);
            bpReader.Get(var_u64, U64);

            bpReader.Get(var_r32, R32);
            bpReader.Get(var_r64, R64);

            bpReader.Get(var_cr32, CR32);
            bpReader.Get(var_cr64, CR64);
            bpReader.PerformGets();

            for (size_t s = 0; s < NSteps; ++s)
            {
                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(s), static_cast<int>(b), mpiSize);

                for (size_t i = 0; i < Nx; ++i)
                {
                    std::stringstream ss;
                    ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    if (var_i8)
                    {
                        ASSERT_EQ(I8[s * Nx + i], currentTestData.I8[i]) << msg;
                    }
                    if (var_i16)
                    {
                        ASSERT_EQ(I16[s * Nx + i], currentTestData.I16[i]) << msg;
                    }
                    if (var_i32)
                    {
                        EXPECT_EQ(I32[s * Nx + i], currentTestData.I32[i]) << msg;
                    }
                    if (var_i64)
                    {
                        EXPECT_EQ(I64[s * Nx + i], currentTestData.I64[i]) << msg;
                    }

                    if (var_u8)
                    {
                        EXPECT_EQ(U8[s * Nx + i], currentTestData.U8[i]) << msg;
                    }
                    if (var_u16)
                    {
                        EXPECT_EQ(U16[s * Nx + i], currentTestData.U16[i]) << msg;
                    }
                    if (var_u32)
                    {
                        EXPECT_EQ(U32[s * Nx + i], currentTestData.U32[i]) << msg;
                    }
                    if (var_u64)
                    {
                        EXPECT_EQ(U64[s * Nx + i], currentTestData.U64[i]) << msg;
                    }
                    if (var_r32)
                    {
                        EXPECT_EQ(R32[s * Nx + i], currentTestData.R32[i]) << msg;
                    }
                    if (var_r64)
                    {
                        EXPECT_EQ(R64[s * Nx + i], currentTestData.R64[i]) << msg;
                    }

                    if (var_cr32)
                    {
                        EXPECT_EQ(CR32[s * Nx + i], currentTestData.CR32[i]) << msg;
                    }
                    if (var_cr64)
                    {
                        EXPECT_EQ(CR64[s * Nx + i], currentTestData.CR64[i]) << msg;
                    }
                }
            }
        }

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariables, ADIOS2BPWriteReadLocal1DBlockInfo)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteReadLocal1DBlockInfo_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadLocal1DBlockInfo.bp");
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{};
            const adios2::Dims start{};
            const adios2::Dims count{Nx};

            io.DefineVariable<int32_t>("stepsGlobalValue");
            io.DefineVariable<std::string>("stepsGlobalValueString");

            io.DefineVariable<int32_t>("ranksLocalValue", {adios2::LocalValueDim});
            io.DefineVariable<std::string>("ranksLocalValueString", {adios2::LocalValueDim});

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

            io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                   adios2::ConstantDims);
            io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                    adios2::ConstantDims);
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.BeginStep();

            const int32_t step32 = static_cast<int32_t>(step);
            bpWriter.Put<int32_t>("stepsGlobalValue", step32);
            bpWriter.Put<std::string>("stepsGlobalValueString", std::to_string(step));

            bpWriter.Put<int32_t>("ranksLocalValue", mpiRank);
            bpWriter.Put<std::string>("ranksLocalValueString", std::to_string(mpiRank));

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
            bpWriter.Put<std::complex<float>>("cr32", currentTestData.CR32.data());
            bpWriter.Put<std::complex<double>>("cr64", currentTestData.CR64.data());
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

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_StepsGlobalValue = io.InquireVariable<int32_t>("stepsGlobalValue");
        auto var_StepsGlobalValueString = io.InquireVariable<std::string>("stepsGlobalValueString");
        auto var_RanksLocalValue = io.InquireVariable<int32_t>("ranksLocalValue");
        auto var_RanksLocalValueString = io.InquireVariable<std::string>("ranksLocalValueString");

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        auto var_i16 = io.InquireVariable<int16_t>("i16");
        auto var_i32 = io.InquireVariable<int32_t>("i32");
        auto var_i64 = io.InquireVariable<int64_t>("i64");
        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        auto var_r32 = io.InquireVariable<float>("r32");
        auto var_r64 = io.InquireVariable<double>("r64");
        auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
        auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");

        (void)var_StepsGlobalValue;
        (void)var_StepsGlobalValueString;
        (void)var_RanksLocalValue;
        (void)var_RanksLocalValueString;
        (void)var_i8;
        (void)var_i16;
        (void)var_i32;
        (void)var_i64;
        (void)var_u8;
        (void)var_u16;
        (void)var_u32;
        (void)var_u64;
        (void)var_r32;
        (void)var_r64;
        (void)var_cr32;
        (void)var_cr64;

        for (size_t s = 0; s < NSteps; ++s)
        {
            const std::vector<typename adios2::Variable<int32_t>::Info> i16_blocks =
                bpReader.BlocksInfo(var_RanksLocalValue, s);

            EXPECT_EQ(i16_blocks.size(), mpiSize);

            size_t i = 0;
            for (const auto &i16_block : i16_blocks)
            {
                std::cout << "Block " << i << "\n";

                size_t j = 0;
                for (const size_t c : i16_block.Count)
                {
                    std::cout << j << ": " << c << "\n";
                    ++j;
                }
                ++i;
                std::cout << "\n";
            }
        }

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariables, ADIOS2BPWriteReadLocal1DSubFile)
{
    const std::string fname("BPWriteReadLocal1DSubFile.bp");

    int mpiRank = 0;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif

    const size_t Nx0 = 3 + static_cast<size_t>(mpiRank);
    const size_t Nx1 = Nx0 + 1;
    // data for block 0 and block 1 per rank
    std::vector<std::vector<int32_t>> data(2);
    data[0] = std::vector<int32_t>(Nx0);
    data[1] = std::vector<int32_t>(Nx1);

    const int32_t startBlock0 = 3 * mpiRank + 1;
    const int32_t startBlock1 = 3 * mpiRank + 15;

    std::iota(data[0].begin(), data[0].end(), startBlock0);
    std::iota(data[1].begin(), data[1].end(), startBlock1);

    /* This is a test for only BP3 */
    if (engineName != "BP3")
    {
        return;
    }

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetEngine("BP3");
        io.SetParameter("AggregatorRatio", "1");
        const adios2::Dims shape{};
        const adios2::Dims start{};
        const adios2::Dims count{Nx0};

        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count, false);

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);
        bpWriter.Put(var_i32, data[0].data());

        var_i32.SetSelection({{}, {Nx1}});
        bpWriter.Put(var_i32, data[1].data());
        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        // each rank opens a local subfile
        const std::string subFileName = fname + ".dir/" + fname + "." + std::to_string(mpiRank);

// requires using MPI_COMM_SELF so metadata won't be shared when reading
// each subfile independently
#if ADIOS2_USE_MPI
        adios2::Engine bpReader =
            io.Open(subFileName, adios2::Mode::ReadRandomAccess, MPI_COMM_SELF);
#else
        adios2::Engine bpReader = io.Open(subFileName, adios2::Mode::ReadRandomAccess);
#endif

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
        EXPECT_EQ(var_i32.Steps(), 1);
        EXPECT_EQ(var_i32.Shape().size(), 0);
        EXPECT_EQ(var_i32.Start().size(), 0);

        std::vector<adios2::Variable<int32_t>::Info> info = bpReader.BlocksInfo(var_i32, 0);

        EXPECT_EQ(info[0].Count[0], Nx0) << " rank= " << mpiRank;
        ASSERT_EQ(info[1].Count[0], Nx1) << " rank= " << mpiRank;

        std::vector<int32_t> inI32;

        for (size_t b = 0; b < 2; ++b)
        {
            const size_t blockNx = (b == 0) ? Nx0 : Nx1;

            var_i32.SetBlockSelection(b);
            EXPECT_EQ(var_i32.BlockID(), b);

            std::stringstream ss;
            ss << "rank= " << mpiRank << " block=" << b << " ";

            EXPECT_EQ(var_i32.Count()[0], blockNx) << ss.str();
            EXPECT_EQ(var_i32.SelectionSize(), blockNx) << ss.str();

            bpReader.Get(var_i32, inI32, adios2::Mode::Sync);

            EXPECT_EQ(inI32.size(), blockNx);

            for (size_t i = 0; i < blockNx; ++i)
            {
                ss << i << ": " << inI32[i] << " " << data[b][i] << " ";
                EXPECT_EQ(inI32[i], data[b][i]) << ss.str();
            }
        } // block loop

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariables, ADIOS2BPWriteReadLocal2DChangeCount)
{

    int mpiRank = 0;
    int mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRLocal2DChangeCount_" + engineName + "_MPI.bp");
#else
    const std::string fname("BPWRLocal2DChangeCount_" + engineName + ".bp");
#endif

    const size_t Nx0 = static_cast<size_t>(std::pow(2 - mpiRank, 2)) + 1;
    const size_t Nx1 = 5;

    // data for block 0 and block 1 per rank
    std::vector<float> data(Nx0 * Nx1);
    const float startBlock = static_cast<float>(3 * mpiRank + 1);

    std::iota(data.begin(), data.end(), startBlock);

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        const adios2::Dims shape{};
        const adios2::Dims start{};
        const adios2::Dims count{Nx0, Nx1};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);
        bpWriter.Put(var_r32, data.data());
        bpWriter.Close();
    }
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    {
        adios2::IO io = adios.DeclareIO("ReaderIO");
#if ADIOS2_USE_MPI
        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess, MPI_COMM_SELF);
#else
        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);
#endif

        auto var_r32 = io.InquireVariable<float>("r32");
        std::map<size_t, std::vector<float>> dataIn;

        for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
        {
            var_r32.SetBlockSelection(b);
            bpReader.Get(var_r32, dataIn[b]);
        }
        bpReader.Close();

        for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
        {
            for (size_t i = 0; i < dataIn[b].size(); ++i)
            {
                EXPECT_EQ(dataIn[b][i], static_cast<float>(3 * b + 1 + i));
            }
        }
    }
}

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
