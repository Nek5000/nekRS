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

class BPWriteReadLocalVariablesSel : public ::testing::Test
{
public:
    BPWriteReadLocalVariablesSel() = default;

    SmallTestData m_TestData;
};

TEST_F(BPWriteReadLocalVariablesSel, BPWriteReadLocal1DSel)
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
    const std::string fname("BPWriteReadLocal1DSel_MPI.bp");
#else
    const std::string fname("BPWriteReadLocal1DSel.bp");
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

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        std::string IString;

        std::vector<int8_t> I8(Nx);
        std::vector<int16_t> I16(Nx);
        std::vector<int32_t> I32(Nx);
        std::vector<int64_t> I64(Nx);
        std::vector<uint8_t> U8(Nx);
        std::vector<uint16_t> U16(Nx);
        std::vector<uint32_t> U32(Nx);
        std::vector<uint64_t> U64(Nx);
        std::vector<float> R32(Nx);
        std::vector<double> R64(Nx);
        std::vector<std::complex<float>> CR32(Nx);
        std::vector<std::complex<double>> CR64(Nx);

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
            EXPECT_EQ(stepsGlobalValueData, currentStep);

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
            EXPECT_EQ(var_RanksLocalValue.Shape()[0], static_cast<size_t>(mpiSize));
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

            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i16.Shape().size(), 0);

            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i32.Shape().size(), 0);

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i64.Shape().size(), 0);

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r32.Shape().size(), 0);

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r64.Shape().size(), 0);

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr32.Shape().size(), 0);

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr64.Shape().size(), 0);

            // loop blocks
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

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                // loop selections
                for (size_t s = 0; s < Nx; ++s)
                {
                    var_i8.SetSelection({{s}, {Nx - s}});
                    var_i16.SetSelection({{s}, {Nx - s}});
                    var_i32.SetSelection({{s}, {Nx - s}});
                    var_i64.SetSelection({{s}, {Nx - s}});

                    var_u8.SetSelection({{s}, {Nx - s}});
                    var_u16.SetSelection({{s}, {Nx - s}});
                    var_u32.SetSelection({{s}, {Nx - s}});
                    var_u64.SetSelection({{s}, {Nx - s}});

                    var_r32.SetSelection({{s}, {Nx - s}});
                    var_r64.SetSelection({{s}, {Nx - s}});

                    var_cr32.SetSelection({{s}, {Nx - s}});
                    var_cr64.SetSelection({{s}, {Nx - s}});

                    EXPECT_EQ(var_i8.Start()[0], s);
                    EXPECT_EQ(var_i8.Count()[0], Nx - s);

                    EXPECT_EQ(var_i16.Start()[0], s);
                    EXPECT_EQ(var_i16.Count()[0], Nx - s);

                    EXPECT_EQ(var_i32.Start()[0], s);
                    EXPECT_EQ(var_i32.Count()[0], Nx - s);

                    EXPECT_EQ(var_i64.Start()[0], s);
                    EXPECT_EQ(var_i64.Count()[0], Nx - s);

                    EXPECT_EQ(var_u8.Start()[0], s);
                    EXPECT_EQ(var_u8.Count()[0], Nx - s);

                    EXPECT_EQ(var_u16.Start()[0], s);
                    EXPECT_EQ(var_u16.Count()[0], Nx - s);

                    EXPECT_EQ(var_u32.Start()[0], s);
                    EXPECT_EQ(var_u32.Count()[0], Nx - s);

                    EXPECT_EQ(var_u64.Start()[0], s);
                    EXPECT_EQ(var_u64.Count()[0], Nx - s);

                    EXPECT_EQ(var_r32.Start()[0], s);
                    EXPECT_EQ(var_r32.Count()[0], Nx - s);

                    EXPECT_EQ(var_r64.Start()[0], s);
                    EXPECT_EQ(var_r64.Count()[0], Nx - s);

                    EXPECT_EQ(var_cr32.Start()[0], s);
                    EXPECT_EQ(var_cr32.Count()[0], Nx - s);

                    EXPECT_EQ(var_cr64.Start()[0], s);
                    EXPECT_EQ(var_cr64.Count()[0], Nx - s);

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

                    for (size_t i = s; i < Nx; ++i)
                    {
                        std::stringstream ss;
                        ss << "t=" << t << " s=" << s << " rank=" << mpiRank;
                        std::string msg = ss.str();

                        if (var_i8)
                        {
                            EXPECT_EQ(I8[i - s], currentTestData.I8[i]) << msg;
                        }

                        if (var_i16)
                        {
                            EXPECT_EQ(I16[i - s], currentTestData.I16[i]) << msg;
                        }
                        if (var_i32)
                        {
                            EXPECT_EQ(I32[i - s], currentTestData.I32[i]) << msg;
                        }
                        if (var_i64)
                        {
                            EXPECT_EQ(I64[i - s], currentTestData.I64[i]) << msg;
                        }

                        if (var_u8)
                        {
                            EXPECT_EQ(U8[i - s], currentTestData.U8[i]) << msg;
                        }
                        if (var_u16)
                        {
                            EXPECT_EQ(U16[i - s], currentTestData.U16[i]) << msg;
                        }
                        if (var_u32)
                        {
                            EXPECT_EQ(U32[i - s], currentTestData.U32[i]) << msg;
                        }
                        if (var_u64)
                        {
                            EXPECT_EQ(U64[i - s], currentTestData.U64[i]) << msg;
                        }
                        if (var_r32)
                        {
                            EXPECT_EQ(R32[i - s], currentTestData.R32[i]) << msg;
                        }
                        if (var_r64)
                        {
                            EXPECT_EQ(R64[i - s], currentTestData.R64[i]) << msg;
                        }

                        if (var_cr32)
                        {
                            EXPECT_EQ(CR32[i - s], currentTestData.CR32[i]) << msg;
                        }
                        if (var_cr64)
                        {
                            EXPECT_EQ(CR64[i - s], currentTestData.CR64[i]) << msg;
                        }
                    }

                } // selection loop
            }     // block loop

            bpReader.EndStep();
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariablesSel, BPWriteReadLocal2D2x4Sel)
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
    const std::string fname("BPWriteReadLocal2D2x4Sel_MPI.bp");
#else
    const std::string fname("BPWriteReadLocal2D2x4Sel.bp");
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
        std::vector<int8_t> I8(Ny * Nx);
        std::vector<int16_t> I16(Ny * Nx);
        std::vector<int32_t> I32(Ny * Nx);
        std::vector<int64_t> I64(Ny * Nx);
        std::vector<uint8_t> U8(Ny * Nx);
        std::vector<uint16_t> U16(Ny * Nx);
        std::vector<uint32_t> U32(Ny * Nx);
        std::vector<uint64_t> U64(Ny * Nx);
        std::vector<float> R32(Ny * Nx);
        std::vector<double> R64(Ny * Nx);
        std::vector<std::complex<float>> CR32(Ny * Nx);
        std::vector<std::complex<double>> CR64(Ny * Nx);

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
            EXPECT_EQ(stepsGlobalValueData, currentStep);

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

            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i16.Shape().size(), 0);

            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i32.Shape().size(), 0);

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i64.Shape().size(), 0);

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r32.Shape().size(), 0);

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r64.Shape().size(), 0);

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr32.Shape().size(), 0);

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr64.Shape().size(), 0);

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

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                // selection loop
                for (size_t j = 0; j < Ny; ++j)
                {
                    for (size_t i = 0; i < Nx; ++i)
                    {
                        var_i8.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_i16.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_i32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_i64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        var_u8.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_u16.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_u32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_u64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        var_r32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_r64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        var_cr32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_cr64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        EXPECT_EQ(var_i8.Count()[0], Ny - j);
                        EXPECT_EQ(var_i8.Count()[1], Nx - i);
                        EXPECT_EQ(var_i8.Start()[0], j);
                        EXPECT_EQ(var_i8.Start()[1], i);

                        EXPECT_EQ(var_i16.Count()[0], Ny - j);
                        EXPECT_EQ(var_i16.Count()[1], Nx - i);
                        EXPECT_EQ(var_i16.Start()[0], j);
                        EXPECT_EQ(var_i16.Start()[1], i);

                        EXPECT_EQ(var_i32.Count()[0], Ny - j);
                        EXPECT_EQ(var_i32.Count()[1], Nx - i);
                        EXPECT_EQ(var_i32.Start()[0], j);
                        EXPECT_EQ(var_i32.Start()[1], i);

                        EXPECT_EQ(var_i64.Count()[0], Ny - j);
                        EXPECT_EQ(var_i64.Count()[1], Nx - i);
                        EXPECT_EQ(var_i64.Start()[0], j);
                        EXPECT_EQ(var_i64.Start()[1], i);

                        EXPECT_EQ(var_u8.Count()[0], Ny - j);
                        EXPECT_EQ(var_u8.Count()[1], Nx - i);
                        EXPECT_EQ(var_u8.Start()[0], j);
                        EXPECT_EQ(var_u8.Start()[1], i);

                        EXPECT_EQ(var_u16.Count()[0], Ny - j);
                        EXPECT_EQ(var_u16.Count()[1], Nx - i);
                        EXPECT_EQ(var_u16.Start()[0], j);
                        EXPECT_EQ(var_u16.Start()[1], i);

                        EXPECT_EQ(var_u32.Count()[0], Ny - j);
                        EXPECT_EQ(var_u32.Count()[1], Nx - i);
                        EXPECT_EQ(var_u32.Start()[0], j);
                        EXPECT_EQ(var_u32.Start()[1], i);

                        EXPECT_EQ(var_u64.Count()[0], Ny - j);
                        EXPECT_EQ(var_u64.Count()[1], Nx - i);
                        EXPECT_EQ(var_u64.Start()[0], j);
                        EXPECT_EQ(var_u64.Start()[1], i);

                        EXPECT_EQ(var_r32.Count()[0], Ny - j);
                        EXPECT_EQ(var_r32.Count()[1], Nx - i);
                        EXPECT_EQ(var_r32.Start()[0], j);
                        EXPECT_EQ(var_r32.Start()[1], i);

                        EXPECT_EQ(var_r64.Count()[0], Ny - j);
                        EXPECT_EQ(var_r64.Count()[1], Nx - i);
                        EXPECT_EQ(var_r64.Start()[0], j);
                        EXPECT_EQ(var_r64.Start()[1], i);

                        EXPECT_EQ(var_cr32.Count()[0], Ny - j);
                        EXPECT_EQ(var_cr32.Count()[1], Nx - i);
                        EXPECT_EQ(var_cr32.Start()[0], j);
                        EXPECT_EQ(var_cr32.Start()[1], i);

                        EXPECT_EQ(var_cr64.Count()[0], Ny - j);
                        EXPECT_EQ(var_cr64.Count()[1], Nx - i);
                        EXPECT_EQ(var_cr64.Start()[0], j);
                        EXPECT_EQ(var_cr64.Start()[1], i);

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

                        for (size_t q = j; q < Ny; ++q)
                        {
                            for (size_t p = i; p < Nx; ++p)
                            {
                                std::stringstream ss;
                                ss << "t=" << t << " q=" << q << " p=" << p << " rank=" << mpiRank;
                                std::string msg = ss.str();

                                const size_t indexSel = (q - j) * (Nx - i) + (p - i);
                                const size_t indexBlock = q * Nx + p;

                                if (var_i8)
                                {
                                    EXPECT_EQ(I8[indexSel], currentTestData.I8[indexBlock]) << msg;
                                }

                                if (var_i16)
                                {
                                    EXPECT_EQ(I16[indexSel], currentTestData.I16[indexBlock])
                                        << msg;
                                }

                                if (var_i32)
                                {
                                    EXPECT_EQ(I32[indexSel], currentTestData.I32[indexBlock])
                                        << msg;
                                }

                                if (var_i64)
                                {
                                    EXPECT_EQ(I64[indexSel], currentTestData.I64[indexBlock])
                                        << msg;
                                }

                                if (var_u8)
                                {
                                    EXPECT_EQ(U8[indexSel], currentTestData.U8[indexBlock]) << msg;
                                }

                                if (var_u16)
                                {
                                    EXPECT_EQ(U16[indexSel], currentTestData.U16[indexBlock])
                                        << msg;
                                }

                                if (var_u32)
                                {
                                    EXPECT_EQ(U32[indexSel], currentTestData.U32[indexBlock])
                                        << msg;
                                }

                                if (var_u64)
                                {
                                    EXPECT_EQ(U64[indexSel], currentTestData.U64[indexBlock])
                                        << msg;
                                }

                                if (var_r32)
                                {
                                    EXPECT_EQ(R32[indexSel], currentTestData.R32[indexBlock])
                                        << msg;
                                }

                                if (var_r64)
                                {
                                    EXPECT_EQ(R64[indexSel], currentTestData.R64[indexBlock])
                                        << msg;
                                }

                                if (var_cr32)
                                {
                                    EXPECT_EQ(CR32[indexSel], currentTestData.CR32[indexBlock])
                                        << msg;
                                }

                                if (var_cr64)
                                {
                                    EXPECT_EQ(CR64[indexSel], currentTestData.CR64[indexBlock])
                                        << msg;
                                }
                            }
                        }
                    } // i selection loop
                }     // j selection loop
            }         // block loop

            bpReader.EndStep();
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariablesSel, BPWriteReadLocal2D4x2Sel)
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
    const std::string fname("BPWriteReadLocal2D4x2Sel_MPI.bp");
#else
    const std::string fname("BPWriteReadLocal2D4x2Sel.bp");
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
        std::vector<int8_t> I8(Ny * Nx);
        std::vector<int16_t> I16(Ny * Nx);
        std::vector<int32_t> I32(Ny * Nx);
        std::vector<int64_t> I64(Ny * Nx);
        std::vector<uint8_t> U8(Ny * Nx);
        std::vector<uint16_t> U16(Ny * Nx);
        std::vector<uint32_t> U32(Ny * Nx);
        std::vector<uint64_t> U64(Ny * Nx);
        std::vector<float> R32(Ny * Nx);
        std::vector<double> R64(Ny * Nx);
        std::vector<std::complex<float>> CR32(Ny * Nx);
        std::vector<std::complex<double>> CR64(Ny * Nx);

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
            EXPECT_EQ(stepsGlobalValueData, currentStep);

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

            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i16.Shape().size(), 0);

            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i32.Shape().size(), 0);

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i64.Shape().size(), 0);

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r32.Shape().size(), 0);

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r64.Shape().size(), 0);

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr32.Shape().size(), 0);

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr64.Shape().size(), 0);

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

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                // selection loop
                for (size_t j = 0; j < Ny; ++j)
                {
                    for (size_t i = 0; i < Nx; ++i)
                    {
                        var_i8.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_i16.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_i32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_i64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        var_u8.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_u16.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_u32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_u64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        var_r32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_r64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        var_cr32.SetSelection({{j, i}, {Ny - j, Nx - i}});
                        var_cr64.SetSelection({{j, i}, {Ny - j, Nx - i}});

                        EXPECT_EQ(var_i8.Count()[0], Ny - j);
                        EXPECT_EQ(var_i8.Count()[1], Nx - i);
                        EXPECT_EQ(var_i8.Start()[0], j);
                        EXPECT_EQ(var_i8.Start()[1], i);

                        EXPECT_EQ(var_i16.Count()[0], Ny - j);
                        EXPECT_EQ(var_i16.Count()[1], Nx - i);
                        EXPECT_EQ(var_i16.Start()[0], j);
                        EXPECT_EQ(var_i16.Start()[1], i);

                        EXPECT_EQ(var_i32.Count()[0], Ny - j);
                        EXPECT_EQ(var_i32.Count()[1], Nx - i);
                        EXPECT_EQ(var_i32.Start()[0], j);
                        EXPECT_EQ(var_i32.Start()[1], i);

                        EXPECT_EQ(var_i64.Count()[0], Ny - j);
                        EXPECT_EQ(var_i64.Count()[1], Nx - i);
                        EXPECT_EQ(var_i64.Start()[0], j);
                        EXPECT_EQ(var_i64.Start()[1], i);

                        EXPECT_EQ(var_u8.Count()[0], Ny - j);
                        EXPECT_EQ(var_u8.Count()[1], Nx - i);
                        EXPECT_EQ(var_u8.Start()[0], j);
                        EXPECT_EQ(var_u8.Start()[1], i);

                        EXPECT_EQ(var_u16.Count()[0], Ny - j);
                        EXPECT_EQ(var_u16.Count()[1], Nx - i);
                        EXPECT_EQ(var_u16.Start()[0], j);
                        EXPECT_EQ(var_u16.Start()[1], i);

                        EXPECT_EQ(var_u32.Count()[0], Ny - j);
                        EXPECT_EQ(var_u32.Count()[1], Nx - i);
                        EXPECT_EQ(var_u32.Start()[0], j);
                        EXPECT_EQ(var_u32.Start()[1], i);

                        EXPECT_EQ(var_u64.Count()[0], Ny - j);
                        EXPECT_EQ(var_u64.Count()[1], Nx - i);
                        EXPECT_EQ(var_u64.Start()[0], j);
                        EXPECT_EQ(var_u64.Start()[1], i);

                        EXPECT_EQ(var_r32.Count()[0], Ny - j);
                        EXPECT_EQ(var_r32.Count()[1], Nx - i);
                        EXPECT_EQ(var_r32.Start()[0], j);
                        EXPECT_EQ(var_r32.Start()[1], i);

                        EXPECT_EQ(var_r64.Count()[0], Ny - j);
                        EXPECT_EQ(var_r64.Count()[1], Nx - i);
                        EXPECT_EQ(var_r64.Start()[0], j);
                        EXPECT_EQ(var_r64.Start()[1], i);

                        EXPECT_EQ(var_cr32.Count()[0], Ny - j);
                        EXPECT_EQ(var_cr32.Count()[1], Nx - i);
                        EXPECT_EQ(var_cr32.Start()[0], j);
                        EXPECT_EQ(var_cr32.Start()[1], i);

                        EXPECT_EQ(var_cr64.Count()[0], Ny - j);
                        EXPECT_EQ(var_cr64.Count()[1], Nx - i);
                        EXPECT_EQ(var_cr64.Start()[0], j);
                        EXPECT_EQ(var_cr64.Start()[1], i);

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

                        for (size_t q = j; q < Ny; ++q)
                        {
                            for (size_t p = i; p < Nx; ++p)
                            {
                                std::stringstream ss;
                                ss << "t=" << t << " q=" << q << " p=" << p << " rank=" << mpiRank;
                                std::string msg = ss.str();

                                const size_t indexSel = (q - j) * (Nx - i) + (p - i);
                                const size_t indexBlock = q * Nx + p;

                                if (var_i8)
                                {
                                    EXPECT_EQ(I8[indexSel], currentTestData.I8[indexBlock]) << msg;
                                }

                                if (var_i16)
                                {
                                    EXPECT_EQ(I16[indexSel], currentTestData.I16[indexBlock])
                                        << msg;
                                }

                                if (var_i32)
                                {
                                    EXPECT_EQ(I32[indexSel], currentTestData.I32[indexBlock])
                                        << msg;
                                }

                                if (var_i64)
                                {
                                    EXPECT_EQ(I64[indexSel], currentTestData.I64[indexBlock])
                                        << msg;
                                }

                                if (var_u8)
                                {
                                    EXPECT_EQ(U8[indexSel], currentTestData.U8[indexBlock]) << msg;
                                }

                                if (var_u16)
                                {
                                    EXPECT_EQ(U16[indexSel], currentTestData.U16[indexBlock])
                                        << msg;
                                }

                                if (var_u32)
                                {
                                    EXPECT_EQ(U32[indexSel], currentTestData.U32[indexBlock])
                                        << msg;
                                }

                                if (var_u64)
                                {
                                    EXPECT_EQ(U64[indexSel], currentTestData.U64[indexBlock])
                                        << msg;
                                }

                                if (var_r32)
                                {
                                    EXPECT_EQ(R32[indexSel], currentTestData.R32[indexBlock])
                                        << msg;
                                }

                                if (var_r64)
                                {
                                    EXPECT_EQ(R64[indexSel], currentTestData.R64[indexBlock])
                                        << msg;
                                }

                                if (var_cr32)
                                {
                                    EXPECT_EQ(CR32[indexSel], currentTestData.CR32[indexBlock])
                                        << msg;
                                }

                                if (var_cr64)
                                {
                                    EXPECT_EQ(CR64[indexSel], currentTestData.CR64[indexBlock])
                                        << msg;
                                }
                            }
                        }
                    } // i selection loop
                }     // j selection loop
            }         // block loop

            bpReader.EndStep();
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadLocalVariablesSel, BPWriteReadLocal1DAllStepsSel)
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
    const std::string fname("BPWriteReadLocal1DAllStepsSel_MPI.bp");
#else
    const std::string fname("BPWriteReadLocal1DAllStepsSel.bp");
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

            for (size_t i = 0; i < Nx; ++i)
            {
                var_i8.SetSelection({{i}, {Nx - i}});
                var_i16.SetSelection({{i}, {Nx - i}});
                var_i32.SetSelection({{i}, {Nx - i}});
                var_i64.SetSelection({{i}, {Nx - i}});

                var_u8.SetSelection({{i}, {Nx - i}});
                var_u16.SetSelection({{i}, {Nx - i}});
                var_u32.SetSelection({{i}, {Nx - i}});
                var_u64.SetSelection({{i}, {Nx - i}});

                var_r32.SetSelection({{i}, {Nx - i}});
                var_r64.SetSelection({{i}, {Nx - i}});

                var_cr32.SetSelection({{i}, {Nx - i}});
                var_cr64.SetSelection({{i}, {Nx - i}});

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

                    for (size_t j = i; j < Nx; ++j)
                    {
                        std::stringstream ss;
                        ss << "t=" << t << " j=" << j << " rank=" << mpiRank;
                        std::string msg = ss.str();

                        if (var_i8)
                        {
                            EXPECT_EQ(I8[s * (Nx - i) + j - i], currentTestData.I8[j]) << msg;
                        }

                        if (var_i16)
                        {
                            EXPECT_EQ(I16[s * (Nx - i) + j - i], currentTestData.I16[j]) << msg;
                        }

                        if (var_i32)
                        {
                            EXPECT_EQ(I32[s * (Nx - i) + j - i], currentTestData.I32[j]) << msg;
                        }

                        if (var_i64)
                        {
                            EXPECT_EQ(I64[s * (Nx - i) + j - i], currentTestData.I64[j]) << msg;
                        }

                        if (var_u8)
                        {
                            EXPECT_EQ(U8[s * (Nx - i) + j - i], currentTestData.U8[j]) << msg;
                        }

                        if (var_u16)
                        {
                            EXPECT_EQ(U16[s * (Nx - i) + j - i], currentTestData.U16[j]) << msg;
                        }

                        if (var_u32)
                        {
                            EXPECT_EQ(U32[s * (Nx - i) + j - i], currentTestData.U32[j]) << msg;
                        }

                        if (var_u64)
                        {
                            EXPECT_EQ(U64[s * (Nx - i) + j - i], currentTestData.U64[j]) << msg;
                        }

                        if (var_r32)
                        {
                            EXPECT_EQ(R32[s * (Nx - i) + j - i], currentTestData.R32[j]) << msg;
                        }

                        if (var_r64)
                        {
                            EXPECT_EQ(R64[s * (Nx - i) + j - i], currentTestData.R64[j]) << msg;
                        }

                        if (var_cr32)
                        {
                            EXPECT_EQ(CR32[s * (Nx - i) + j - i], currentTestData.CR32[j]) << msg;
                        }

                        if (var_cr64)
                        {
                            EXPECT_EQ(CR64[s * (Nx - i) + j - i], currentTestData.CR64[j]) << msg;
                        }
                    }
                } // steps loop
            }     // selection loop
        }

        bpReader.Close();
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
