/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <algorithm> //std::min_element, std::max_element
#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPWriteReadMultiblockTest : public ::testing::Test
{
public:
    BPWriteReadMultiblockTest() = default;

    SmallTestData m_TestData;
};

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

TEST_F(BPWriteReadMultiblockTest, ADIOS2BPWriteReadMultiblock1D8)
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
    const std::string fname("ADIOS2BPWriteReadMultiblock1D8_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadMultiblock1D8.bp");
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    // Make a 1D selection to describe the local dimensions of the
    // variable we write and its offsets in the global spaces
    const adios2::Box<adios2::Dims> sel1({mpiRank * Nx}, {Nx / 2});
    const adios2::Box<adios2::Dims> sel2({mpiRank * Nx + Nx / 2}, {Nx - Nx / 2});

    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
            const adios2::Dims count{Nx};

            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            EXPECT_TRUE(var_i8);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            EXPECT_TRUE(var_i16);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            EXPECT_TRUE(var_i64);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            EXPECT_TRUE(var_u8);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            EXPECT_TRUE(var_u16);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            EXPECT_TRUE(var_u32);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            EXPECT_TRUE(var_u64);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            EXPECT_TRUE(var_r32);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
            EXPECT_TRUE(var_r64);
            auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count);
            EXPECT_TRUE(var_cr32);
            auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count);
            EXPECT_TRUE(var_cr64);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("BPFile");
        }

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
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto var_iString = io.InquireVariable<std::string>("iString");
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

            bpWriter.BeginStep();
            bpWriter.Put(var_iString, currentTestData.S1);

            var_i8.SetSelection(sel1);
            bpWriter.Put(var_i8, currentTestData.I8.data());
            var_i8.SetSelection(sel2);
            bpWriter.Put(var_i8, currentTestData.I8.data() + Nx / 2);

            var_i16.SetSelection(sel1);
            bpWriter.Put(var_i16, currentTestData.I16.data());
            var_i16.SetSelection(sel2);
            bpWriter.Put(var_i16, currentTestData.I16.data() + Nx / 2);

            var_i32.SetSelection(sel1);
            bpWriter.Put(var_i32, currentTestData.I32.data());
            var_i32.SetSelection(sel2);
            bpWriter.Put(var_i32, currentTestData.I32.data() + Nx / 2);

            var_i64.SetSelection(sel1);
            bpWriter.Put(var_i64, currentTestData.I64.data());
            var_i64.SetSelection(sel2);
            bpWriter.Put(var_i64, currentTestData.I64.data() + Nx / 2);

            var_u8.SetSelection(sel1);
            bpWriter.Put(var_u8, currentTestData.U8.data());
            var_u8.SetSelection(sel2);
            bpWriter.Put(var_u8, currentTestData.U8.data() + Nx / 2);

            var_u16.SetSelection(sel1);
            bpWriter.Put(var_u16, currentTestData.U16.data());
            var_u16.SetSelection(sel2);
            bpWriter.Put(var_u16, currentTestData.U16.data() + Nx / 2);

            var_u32.SetSelection(sel1);
            bpWriter.Put(var_u32, currentTestData.U32.data());
            var_u32.SetSelection(sel2);
            bpWriter.Put(var_u32, currentTestData.U32.data() + Nx / 2);

            var_u64.SetSelection(sel1);
            bpWriter.Put(var_u64, currentTestData.U64.data());
            var_u64.SetSelection(sel2);
            bpWriter.Put(var_u64, currentTestData.U64.data() + Nx / 2);

            var_r32.SetSelection(sel1);
            bpWriter.Put(var_r32, currentTestData.R32.data());
            var_r32.SetSelection(sel2);
            bpWriter.Put(var_r32, currentTestData.R32.data() + Nx / 2);

            var_r64.SetSelection(sel1);
            bpWriter.Put(var_r64, currentTestData.R64.data());
            var_r64.SetSelection(sel2);
            bpWriter.Put(var_r64, currentTestData.R64.data() + Nx / 2);

            var_cr32.SetSelection(sel1);
            bpWriter.Put(var_cr32, currentTestData.CR32.data());
            var_cr32.SetSelection(sel2);
            bpWriter.Put(var_cr32, currentTestData.CR32.data() + Nx / 2);

            var_cr64.SetSelection(sel1);
            bpWriter.Put(var_cr64, currentTestData.CR64.data());
            var_cr64.SetSelection(sel2);
            bpWriter.Put(var_cr64, currentTestData.CR64.data() + Nx / 2);

            bpWriter.EndStep();
        }

        // Close the file
        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        EXPECT_EQ(var_iString.Shape().size(), 0);
        EXPECT_EQ(var_iString.Steps(), NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i8.Steps(), NSteps);
        EXPECT_EQ(var_i8.Shape()[0], mpiSize * Nx);

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i16.Steps(), NSteps);
        EXPECT_EQ(var_i16.Shape()[0], mpiSize * Nx);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i32.Steps(), NSteps);
        EXPECT_EQ(var_i32.Shape()[0], mpiSize * Nx);

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i64.Steps(), NSteps);
        EXPECT_EQ(var_i64.Shape()[0], mpiSize * Nx);

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u8.Steps(), NSteps);
        EXPECT_EQ(var_u8.Shape()[0], mpiSize * Nx);

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u16.Steps(), NSteps);
        EXPECT_EQ(var_u16.Shape()[0], mpiSize * Nx);

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u32.Steps(), NSteps);
        EXPECT_EQ(var_u32.Shape()[0], mpiSize * Nx);

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u64.Steps(), NSteps);
        EXPECT_EQ(var_u64.Shape()[0], mpiSize * Nx);

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_r32.Steps(), NSteps);
        EXPECT_EQ(var_r32.Shape()[0], mpiSize * Nx);

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_r64.Steps(), NSteps);
        EXPECT_EQ(var_r64.Shape()[0], mpiSize * Nx);

        auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
        EXPECT_TRUE(var_cr32);
        EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_cr32.Steps(), NSteps);
        EXPECT_EQ(var_cr32.Shape()[0], mpiSize * Nx);

        auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
        EXPECT_TRUE(var_cr64);
        EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_cr64.Steps(), NSteps);
        EXPECT_EQ(var_cr64.Shape()[0], mpiSize * Nx);

        // TODO: other types

        SmallTestData testData;

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

        const auto iStringAllInfo = bpReader.AllStepsBlocksInfo(var_iString);

        const auto i8AllInfo = bpReader.AllStepsBlocksInfo(var_i8);
        const auto i16AllInfo = bpReader.AllStepsBlocksInfo(var_i16);
        const auto i32AllInfo = bpReader.AllStepsBlocksInfo(var_i32);
        const auto i64AllInfo = bpReader.AllStepsBlocksInfo(var_i64);

        const auto u8AllInfo = bpReader.AllStepsBlocksInfo(var_u8);
        const auto u16AllInfo = bpReader.AllStepsBlocksInfo(var_u16);
        const auto u32AllInfo = bpReader.AllStepsBlocksInfo(var_u32);
        const auto u64AllInfo = bpReader.AllStepsBlocksInfo(var_u64);

        const auto r32AllInfo = bpReader.AllStepsBlocksInfo(var_r32);
        const auto r64AllInfo = bpReader.AllStepsBlocksInfo(var_r64);

        const auto cr32AllInfo = bpReader.AllStepsBlocksInfo(var_cr32);
        const auto cr64AllInfo = bpReader.AllStepsBlocksInfo(var_cr64);

        EXPECT_EQ(iStringAllInfo.size(), NSteps);
        EXPECT_EQ(i8AllInfo.size(), NSteps);
        EXPECT_EQ(i16AllInfo.size(), NSteps);
        EXPECT_EQ(i32AllInfo.size(), NSteps);
        EXPECT_EQ(i64AllInfo.size(), NSteps);
        EXPECT_EQ(u8AllInfo.size(), NSteps);
        EXPECT_EQ(u16AllInfo.size(), NSteps);
        EXPECT_EQ(u32AllInfo.size(), NSteps);
        EXPECT_EQ(u64AllInfo.size(), NSteps);
        EXPECT_EQ(r32AllInfo.size(), NSteps);
        EXPECT_EQ(r64AllInfo.size(), NSteps);
        EXPECT_EQ(cr32AllInfo.size(), NSteps);
        EXPECT_EQ(cr64AllInfo.size(), NSteps);

        for (size_t CurrentStep = 0; CurrentStep < NSteps; CurrentStep++)
        {
            const std::vector<adios2::Variable<std::string>::Info> &iStringInfo =
                iStringAllInfo.at(t);

            const std::vector<adios2::Variable<int8_t>::Info> &i8Info = i8AllInfo.at(t);
            const std::vector<adios2::Variable<int16_t>::Info> &i16Info = i16AllInfo.at(t);
            const std::vector<adios2::Variable<int32_t>::Info> &i32Info = i32AllInfo.at(t);
            const std::vector<adios2::Variable<int64_t>::Info> &i64Info = i64AllInfo.at(t);

            const std::vector<adios2::Variable<uint8_t>::Info> &u8Info = u8AllInfo.at(t);
            const std::vector<adios2::Variable<uint16_t>::Info> &u16Info = u16AllInfo.at(t);
            const std::vector<adios2::Variable<uint32_t>::Info> &u32Info = u32AllInfo.at(t);
            const std::vector<adios2::Variable<uint64_t>::Info> &u64Info = u64AllInfo.at(t);

            const std::vector<adios2::Variable<float>::Info> &r32Info = r32AllInfo.at(t);
            const std::vector<adios2::Variable<double>::Info> &r64Info = r64AllInfo.at(t);

            const std::vector<adios2::Variable<std::complex<float>>::Info> &cr32Info =
                cr32AllInfo.at(t);
            const std::vector<adios2::Variable<std::complex<double>>::Info> &cr64Info =
                cr64AllInfo.at(t);

            EXPECT_EQ(iStringInfo.size(), mpiSize);
            EXPECT_EQ(i8Info.size(), 2 * mpiSize);
            EXPECT_EQ(i16Info.size(), 2 * mpiSize);
            EXPECT_EQ(i32Info.size(), 2 * mpiSize);
            EXPECT_EQ(i64Info.size(), 2 * mpiSize);
            EXPECT_EQ(u8Info.size(), 2 * mpiSize);
            EXPECT_EQ(u16Info.size(), 2 * mpiSize);
            EXPECT_EQ(u32Info.size(), 2 * mpiSize);
            EXPECT_EQ(u64Info.size(), 2 * mpiSize);
            EXPECT_EQ(r32Info.size(), 2 * mpiSize);
            EXPECT_EQ(r64Info.size(), 2 * mpiSize);
            EXPECT_EQ(cr32Info.size(), 2 * mpiSize);
            EXPECT_EQ(cr64Info.size(), 2 * mpiSize);

            // String
            for (size_t i = 0; i < static_cast<size_t>(mpiSize); ++i)
            {
                EXPECT_TRUE(iStringInfo[i].IsValue);
                EXPECT_EQ(iStringInfo[i].Value, "Testing ADIOS2 String type");
            }

            for (size_t i = 0; i < 2 * static_cast<size_t>(mpiSize); ++i)
            {
                EXPECT_FALSE(i8Info[0].IsValue);
                EXPECT_EQ(i8Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(i16Info[0].IsValue);
                EXPECT_EQ(i16Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(i32Info[0].IsValue);
                EXPECT_EQ(i32Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(i64Info[0].IsValue);
                EXPECT_EQ(i64Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(u8Info[0].IsValue);
                EXPECT_EQ(u8Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(u16Info[0].IsValue);
                EXPECT_EQ(u16Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(u32Info[0].IsValue);
                EXPECT_EQ(u32Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(u64Info[0].IsValue);
                EXPECT_EQ(u64Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(r32Info[0].IsValue);
                EXPECT_EQ(r32Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(r64Info[0].IsValue);
                EXPECT_EQ(r64Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(cr32Info[0].IsValue);
                EXPECT_EQ(cr32Info[i].Count[0], Nx / 2);

                EXPECT_FALSE(cr64Info[0].IsValue);
                EXPECT_EQ(cr64Info[i].Count[0], Nx / 2);

                const size_t inRank = i / 2;
                int8_t i8Min, i8Max;
                int16_t i16Min, i16Max;
                int32_t i32Min, i32Max;
                int64_t i64Min, i64Max;
                uint8_t u8Min, u8Max;
                uint16_t u16Min, u16Max;
                uint32_t u32Min, u32Max;
                uint64_t u64Min, u64Max;
                float r32Min, r32Max;
                double r64Min, r64Max;
                std::complex<float> cr32Min, cr32Max;
                std::complex<double> cr64Min, cr64Max;

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(t), static_cast<int>(inRank), mpiSize);

                if (i % 2 == 0)
                {
                    EXPECT_EQ(i8Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(i16Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(i32Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(i64Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(u8Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(u16Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(u32Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(u64Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(r32Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(r64Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(cr32Info[i].Start[0], inRank * Nx);
                    EXPECT_EQ(cr64Info[i].Start[0], inRank * Nx);

                    i8Min = *std::min_element(currentTestData.I8.begin(),
                                              currentTestData.I8.begin() + Nx / 2);
                    i8Max = *std::max_element(currentTestData.I8.begin(),
                                              currentTestData.I8.begin() + Nx / 2);
                    i16Min = *std::min_element(currentTestData.I16.begin(),
                                               currentTestData.I16.begin() + Nx / 2);
                    i16Max = *std::max_element(currentTestData.I16.begin(),
                                               currentTestData.I16.begin() + Nx / 2);
                    i32Min = *std::min_element(currentTestData.I32.begin(),
                                               currentTestData.I32.begin() + Nx / 2);
                    i32Max = *std::max_element(currentTestData.I32.begin(),
                                               currentTestData.I32.begin() + Nx / 2);
                    i64Min = *std::min_element(currentTestData.I64.begin(),
                                               currentTestData.I64.begin() + Nx / 2);
                    i64Max = *std::max_element(currentTestData.I64.begin(),
                                               currentTestData.I64.begin() + Nx / 2);
                    u8Min = *std::min_element(currentTestData.U8.begin(),
                                              currentTestData.U8.begin() + Nx / 2);
                    u8Max = *std::max_element(currentTestData.U8.begin(),
                                              currentTestData.U8.begin() + Nx / 2);
                    u16Min = *std::min_element(currentTestData.U16.begin(),
                                               currentTestData.U16.begin() + Nx / 2);
                    u16Max = *std::max_element(currentTestData.U16.begin(),
                                               currentTestData.U16.begin() + Nx / 2);
                    u32Min = *std::min_element(currentTestData.U32.begin(),
                                               currentTestData.U32.begin() + Nx / 2);
                    u32Max = *std::max_element(currentTestData.U32.begin(),
                                               currentTestData.U32.begin() + Nx / 2);
                    u64Min = *std::min_element(currentTestData.U64.begin(),
                                               currentTestData.U64.begin() + Nx / 2);
                    u64Max = *std::max_element(currentTestData.U64.begin(),
                                               currentTestData.U64.begin() + Nx / 2);
                    r32Min = *std::min_element(currentTestData.R32.begin(),
                                               currentTestData.R32.begin() + Nx / 2);
                    r32Max = *std::max_element(currentTestData.R32.begin(),
                                               currentTestData.R32.begin() + Nx / 2);
                    r64Min = *std::min_element(currentTestData.R64.begin(),
                                               currentTestData.R64.begin() + Nx / 2);
                    r64Max = *std::max_element(currentTestData.R64.begin(),
                                               currentTestData.R64.begin() + Nx / 2);

                    cr32Min = currentTestData.CR32.front();
                    cr32Max = currentTestData.CR32.front();
                    for (auto it = currentTestData.CR32.begin();
                         it != currentTestData.CR32.begin() + Nx / 2; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr32Min))
                        {
                            cr32Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr32Max))
                        {
                            cr32Max = *it;
                        }
                    }

                    cr64Min = currentTestData.CR64.front();
                    cr64Max = currentTestData.CR64.front();
                    for (auto it = currentTestData.CR64.begin();
                         it != currentTestData.CR64.begin() + Nx / 2; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr64Min))
                        {
                            cr64Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr64Max))
                        {
                            cr64Max = *it;
                        }
                    }
                }
                else
                {
                    EXPECT_EQ(i8Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(i16Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(i32Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(i64Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(u8Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(u16Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(u32Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(u64Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(r32Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(r64Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(cr32Info[i].Start[0], inRank * Nx + Nx / 2);
                    EXPECT_EQ(cr64Info[i].Start[0], inRank * Nx + Nx / 2);

                    i8Min = *std::min_element(currentTestData.I8.begin() + Nx / 2,
                                              currentTestData.I8.begin() + Nx);
                    i8Max = *std::max_element(currentTestData.I8.begin() + Nx / 2,
                                              currentTestData.I8.begin() + Nx);

                    i16Min = *std::min_element(currentTestData.I16.begin() + Nx / 2,
                                               currentTestData.I16.begin() + Nx);
                    i16Max = *std::max_element(currentTestData.I16.begin() + Nx / 2,
                                               currentTestData.I16.begin() + Nx);

                    i32Min = *std::min_element(currentTestData.I32.begin() + Nx / 2,
                                               currentTestData.I32.begin() + Nx);
                    i32Max = *std::max_element(currentTestData.I32.begin() + Nx / 2,
                                               currentTestData.I32.begin() + Nx);

                    i64Min = *std::min_element(currentTestData.I64.begin() + Nx / 2,
                                               currentTestData.I64.begin() + Nx);
                    i64Max = *std::max_element(currentTestData.I64.begin() + Nx / 2,
                                               currentTestData.I64.begin() + Nx);

                    u8Min = *std::min_element(currentTestData.U8.begin() + Nx / 2,
                                              currentTestData.U8.begin() + Nx);
                    u8Max = *std::max_element(currentTestData.U8.begin() + Nx / 2,
                                              currentTestData.U8.begin() + Nx);

                    u16Min = *std::min_element(currentTestData.U16.begin() + Nx / 2,
                                               currentTestData.U16.begin() + Nx);
                    u16Max = *std::max_element(currentTestData.U16.begin() + Nx / 2,
                                               currentTestData.U16.begin() + Nx);

                    u32Min = *std::min_element(currentTestData.U32.begin() + Nx / 2,
                                               currentTestData.U32.begin() + Nx);
                    u32Max = *std::max_element(currentTestData.U32.begin() + Nx / 2,
                                               currentTestData.U32.begin() + Nx);

                    u64Min = *std::min_element(currentTestData.U64.begin() + Nx / 2,
                                               currentTestData.U64.begin() + Nx);
                    u64Max = *std::max_element(currentTestData.U64.begin() + Nx / 2,
                                               currentTestData.U64.begin() + Nx);

                    r32Min = *std::min_element(currentTestData.R32.begin() + Nx / 2,
                                               currentTestData.R32.begin() + Nx);
                    r32Max = *std::max_element(currentTestData.R32.begin() + Nx / 2,
                                               currentTestData.R32.begin() + Nx);

                    r64Min = *std::min_element(currentTestData.R64.begin() + Nx / 2,
                                               currentTestData.R64.begin() + Nx);
                    r64Max = *std::max_element(currentTestData.R64.begin() + Nx / 2,
                                               currentTestData.R64.begin() + Nx);

                    cr32Min = currentTestData.CR32[Nx / 2];
                    cr32Max = currentTestData.CR32[Nx / 2];
                    for (auto it = currentTestData.CR32.begin() + Nx / 2;
                         it != currentTestData.CR32.begin() + Nx; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr32Min))
                        {
                            cr32Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr32Max))
                        {
                            cr32Max = *it;
                        }
                    }

                    cr64Min = currentTestData.CR64[Nx / 2];
                    cr64Max = currentTestData.CR64[Nx / 2];
                    for (auto it = currentTestData.CR64.begin() + Nx / 2;
                         it != currentTestData.CR64.begin() + Nx; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr64Min))
                        {
                            cr64Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr64Max))
                        {
                            cr64Max = *it;
                        }
                    }
                }

                EXPECT_EQ(i8Info[i].Min, i8Min);
                EXPECT_EQ(i8Info[i].Max, i8Max);
                EXPECT_EQ(i16Info[i].Min, i16Min);
                EXPECT_EQ(i16Info[i].Max, i16Max);
                EXPECT_EQ(i32Info[i].Min, i32Min);
                EXPECT_EQ(i32Info[i].Max, i32Max);
                EXPECT_EQ(i64Info[i].Min, i64Min);
                EXPECT_EQ(i64Info[i].Max, i64Max);

                EXPECT_EQ(u8Info[i].Min, u8Min);
                EXPECT_EQ(u8Info[i].Max, u8Max);
                EXPECT_EQ(u16Info[i].Min, u16Min);
                EXPECT_EQ(u16Info[i].Max, u16Max);
                EXPECT_EQ(u32Info[i].Min, u32Min);
                EXPECT_EQ(u32Info[i].Max, u32Max);
                EXPECT_EQ(u64Info[i].Min, u64Min);
                EXPECT_EQ(u64Info[i].Max, u64Max);

                EXPECT_EQ(r32Info[i].Min, r32Min);
                EXPECT_EQ(r32Info[i].Max, r32Max);
                EXPECT_EQ(r64Info[i].Min, r64Min);
                EXPECT_EQ(r64Info[i].Max, r64Max);
            }

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            var_iString.SetStepSelection({CurrentStep, 1});
            bpReader.Get(var_iString, IString);

            var_i8.SetStepSelection({CurrentStep, 1});
            var_i8.SetSelection(sel1);
            bpReader.Get(var_i8, I8.data());
            var_i8.SetSelection(sel2);
            bpReader.Get(var_i8, I8.data() + Nx / 2);

            var_i16.SetStepSelection({CurrentStep, 1});
            var_i16.SetSelection(sel1);
            bpReader.Get(var_i16, I16.data());
            var_i16.SetSelection(sel2);
            bpReader.Get(var_i16, I16.data() + Nx / 2);

            var_i32.SetStepSelection({CurrentStep, 1});
            var_i32.SetSelection(sel1);
            bpReader.Get(var_i32, I32.data());
            var_i32.SetSelection(sel2);
            bpReader.Get(var_i32, I32.data() + Nx / 2);

            var_i64.SetStepSelection({CurrentStep, 1});
            var_i64.SetSelection(sel1);
            bpReader.Get(var_i64, I64.data());
            var_i64.SetSelection(sel2);
            bpReader.Get(var_i64, I64.data() + Nx / 2);

            var_u8.SetStepSelection({CurrentStep, 1});
            var_u8.SetSelection(sel1);
            bpReader.Get(var_u8, U8.data());
            var_u8.SetSelection(sel2);
            bpReader.Get(var_u8, U8.data() + Nx / 2);

            var_u16.SetStepSelection({CurrentStep, 1});
            var_u16.SetSelection(sel1);
            bpReader.Get(var_u16, U16.data());
            var_u16.SetSelection(sel2);
            bpReader.Get(var_u16, U16.data() + Nx / 2);

            var_u32.SetStepSelection({CurrentStep, 1});
            var_u32.SetSelection(sel1);
            bpReader.Get(var_u32, U32.data());
            var_u32.SetSelection(sel2);
            bpReader.Get(var_u32, U32.data() + Nx / 2);

            var_u64.SetStepSelection({CurrentStep, 1});
            var_u64.SetSelection(sel1);
            bpReader.Get(var_u64, U64.data());
            var_u64.SetSelection(sel2);
            bpReader.Get(var_u64, U64.data() + Nx / 2);

            var_r32.SetStepSelection({CurrentStep, 1});
            var_r32.SetSelection(sel1);
            bpReader.Get(var_r32, R32.data());
            var_r32.SetSelection(sel2);
            bpReader.Get(var_r32, R32.data() + Nx / 2);

            var_r64.SetStepSelection({CurrentStep, 1});
            var_r64.SetSelection(sel1);
            bpReader.Get(var_r64, R64.data());
            var_r64.SetSelection(sel2);
            bpReader.Get(var_r64, R64.data() + Nx / 2);

            var_cr32.SetStepSelection({CurrentStep, 1});
            var_cr32.SetSelection(sel1);
            bpReader.Get(var_cr32, CR32.data());
            var_cr32.SetSelection(sel2);
            bpReader.Get(var_cr32, CR32.data() + Nx / 2);

            var_cr64.SetStepSelection({CurrentStep, 1});
            var_cr64.SetSelection(sel1);
            bpReader.Get(var_cr64, CR64.data());
            var_cr64.SetSelection(sel2);
            bpReader.Get(var_cr64, CR64.data() + Nx / 2);

            EXPECT_EQ(IString, currentTestData.S1);

            bpReader.PerformGets();
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
                EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
            }

            ++t;
        }
        bpReader.Close();
    }
}

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

TEST_F(BPWriteReadMultiblockTest, ADIOS2BPWriteReadMultiblock2D2x4)
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
    const std::string fname("ADIOS2BPWriteReadMultiblock2D2x4Test_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadMultiblock2D2x4Test.bp");
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    // Make a 2D selection to describe the local dimensions of the
    // variable we write and its offsets in the global spaces
    const adios2::Box<adios2::Dims> sel1({0, static_cast<size_t>(mpiRank * Nx)}, {Ny / 2, Nx});

    const adios2::Box<adios2::Dims> sel2({Ny / 2, static_cast<size_t>(mpiRank * Nx)},
                                         {Ny - Ny / 2, Nx});

    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
            const adios2::Dims count{Ny, Nx};

            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            EXPECT_TRUE(var_i8);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            EXPECT_TRUE(var_i16);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            EXPECT_TRUE(var_i64);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            EXPECT_TRUE(var_u8);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            EXPECT_TRUE(var_u16);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            EXPECT_TRUE(var_u32);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            EXPECT_TRUE(var_u64);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            EXPECT_TRUE(var_r32);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
            EXPECT_TRUE(var_r64);
            auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count);
            EXPECT_TRUE(var_cr32);
            auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count);
            EXPECT_TRUE(var_cr64);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("BPFile");
        }
        io.AddTransport("file");

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto var_iString = io.InquireVariable<std::string>("iString");
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

            bpWriter.BeginStep();
            bpWriter.Put(var_iString, currentTestData.S1);

            var_i8.SetSelection(sel1);
            bpWriter.Put(var_i8, currentTestData.I8.data());
            var_i8.SetSelection(sel2);
            bpWriter.Put(var_i8, currentTestData.I8.data() + Ny * Nx / 2);

            var_i16.SetSelection(sel1);
            bpWriter.Put(var_i16, currentTestData.I16.data());
            var_i16.SetSelection(sel2);
            bpWriter.Put(var_i16, currentTestData.I16.data() + Ny * Nx / 2);

            var_i32.SetSelection(sel1);
            bpWriter.Put(var_i32, currentTestData.I32.data());
            var_i32.SetSelection(sel2);
            bpWriter.Put(var_i32, currentTestData.I32.data() + Ny * Nx / 2);

            var_i64.SetSelection(sel1);
            bpWriter.Put(var_i64, currentTestData.I64.data());
            var_i64.SetSelection(sel2);
            bpWriter.Put(var_i64, currentTestData.I64.data() + Ny * Nx / 2);

            var_u8.SetSelection(sel1);
            bpWriter.Put(var_u8, currentTestData.U8.data());
            var_u8.SetSelection(sel2);
            bpWriter.Put(var_u8, currentTestData.U8.data() + Ny * Nx / 2);

            var_u16.SetSelection(sel1);
            bpWriter.Put(var_u16, currentTestData.U16.data());
            var_u16.SetSelection(sel2);
            bpWriter.Put(var_u16, currentTestData.U16.data() + Ny * Nx / 2);

            var_u32.SetSelection(sel1);
            bpWriter.Put(var_u32, currentTestData.U32.data());
            var_u32.SetSelection(sel2);
            bpWriter.Put(var_u32, currentTestData.U32.data() + Ny * Nx / 2);

            var_u64.SetSelection(sel1);
            bpWriter.Put(var_u64, currentTestData.U64.data());
            var_u64.SetSelection(sel2);
            bpWriter.Put(var_u64, currentTestData.U64.data() + Ny * Nx / 2);

            var_r32.SetSelection(sel1);
            bpWriter.Put(var_r32, currentTestData.R32.data());
            var_r32.SetSelection(sel2);
            bpWriter.Put(var_r32, currentTestData.R32.data() + Ny * Nx / 2);

            var_r64.SetSelection(sel1);
            bpWriter.Put(var_r64, currentTestData.R64.data());
            var_r64.SetSelection(sel2);
            bpWriter.Put(var_r64, currentTestData.R64.data() + Ny * Nx / 2);

            var_cr32.SetSelection(sel1);
            bpWriter.Put(var_cr32, currentTestData.CR32.data());
            var_cr32.SetSelection(sel2);
            bpWriter.Put(var_cr32, currentTestData.CR32.data() + Ny * Nx / 2);

            var_cr64.SetSelection(sel1);
            bpWriter.Put(var_cr64, currentTestData.CR64.data());
            var_cr64.SetSelection(sel2);
            bpWriter.Put(var_cr64, currentTestData.CR64.data() + Ny * Nx / 2);

            bpWriter.EndStep();
        }

        // Close the file
        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        EXPECT_EQ(var_iString.Shape().size(), 0);
        EXPECT_EQ(var_iString.Steps(), NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i8.Steps(), NSteps);
        EXPECT_EQ(var_i8.Shape()[0], Ny);
        EXPECT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i16.Steps(), NSteps);
        EXPECT_EQ(var_i16.Shape()[0], Ny);
        EXPECT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i32.Steps(), NSteps);
        EXPECT_EQ(var_i32.Shape()[0], Ny);
        EXPECT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i64.Steps(), NSteps);
        EXPECT_EQ(var_i64.Shape()[0], Ny);
        EXPECT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u8.Steps(), NSteps);
        EXPECT_EQ(var_u8.Shape()[0], Ny);
        EXPECT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u16.Steps(), NSteps);
        EXPECT_EQ(var_u16.Shape()[0], Ny);
        EXPECT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u32.Steps(), NSteps);
        EXPECT_EQ(var_u32.Shape()[0], Ny);
        EXPECT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u64.Steps(), NSteps);
        EXPECT_EQ(var_u64.Shape()[0], Ny);
        EXPECT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_r32.Steps(), NSteps);
        EXPECT_EQ(var_r32.Shape()[0], Ny);
        EXPECT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_r64.Steps(), NSteps);
        EXPECT_EQ(var_r64.Shape()[0], Ny);
        EXPECT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
        EXPECT_TRUE(var_cr32);
        EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_cr32.Steps(), NSteps);
        EXPECT_EQ(var_cr32.Shape()[0], Ny);
        EXPECT_EQ(var_cr32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
        EXPECT_TRUE(var_cr64);
        EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_cr64.Steps(), NSteps);
        EXPECT_EQ(var_cr64.Shape()[0], Ny);
        EXPECT_EQ(var_cr64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::string IString;
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
        std::array<std::complex<float>, Nx * Ny> CR32;
        std::array<std::complex<double>, Nx * Ny> CR64;

        size_t t = 0;

        for (size_t CurrentStep = 0; CurrentStep < NSteps; CurrentStep++)
        {
            var_iString.SetStepSelection({CurrentStep, 1});
            bpReader.Get(var_iString, IString);

            var_i8.SetStepSelection({CurrentStep, 1});
            var_i8.SetSelection(sel1);
            bpReader.Get(var_i8, I8.data());
            var_i8.SetSelection(sel2);
            bpReader.Get(var_i8, I8.data() + Ny * Nx / 2);

            var_i16.SetStepSelection({CurrentStep, 1});
            var_i16.SetSelection(sel1);
            bpReader.Get(var_i16, I16.data());
            var_i16.SetSelection(sel2);
            bpReader.Get(var_i16, I16.data() + Ny * Nx / 2);

            var_i32.SetStepSelection({CurrentStep, 1});
            var_i32.SetSelection(sel1);
            bpReader.Get(var_i32, I32.data());
            var_i32.SetSelection(sel2);
            bpReader.Get(var_i32, I32.data() + Ny * Nx / 2);

            var_i64.SetStepSelection({CurrentStep, 1});
            var_i64.SetSelection(sel1);
            bpReader.Get(var_i64, I64.data());
            var_i64.SetSelection(sel2);
            bpReader.Get(var_i64, I64.data() + Ny * Nx / 2);

            var_u8.SetStepSelection({CurrentStep, 1});
            var_u8.SetSelection(sel1);
            bpReader.Get(var_u8, U8.data());
            var_u8.SetSelection(sel2);
            bpReader.Get(var_u8, U8.data() + Ny * Nx / 2);

            var_u16.SetStepSelection({CurrentStep, 1});
            var_u16.SetSelection(sel1);
            bpReader.Get(var_u16, U16.data());
            var_u16.SetSelection(sel2);
            bpReader.Get(var_u16, U16.data() + Ny * Nx / 2);

            var_u32.SetStepSelection({CurrentStep, 1});
            var_u32.SetSelection(sel1);
            bpReader.Get(var_u32, U32.data());
            var_u32.SetSelection(sel2);
            bpReader.Get(var_u32, U32.data() + Ny * Nx / 2);

            var_u64.SetStepSelection({CurrentStep, 1});
            var_u64.SetSelection(sel1);
            bpReader.Get(var_u64, U64.data());
            var_u64.SetSelection(sel2);
            bpReader.Get(var_u64, U64.data() + Ny * Nx / 2);

            var_r32.SetStepSelection({CurrentStep, 1});
            var_r32.SetSelection(sel1);
            bpReader.Get(var_r32, R32.data());
            var_r32.SetSelection(sel2);
            bpReader.Get(var_r32, R32.data() + Ny * Nx / 2);

            var_r64.SetStepSelection({CurrentStep, 1});
            var_r64.SetSelection(sel1);
            bpReader.Get(var_r64, R64.data());
            var_r64.SetSelection(sel2);
            bpReader.Get(var_r64, R64.data() + Ny * Nx / 2);

            var_cr32.SetStepSelection({CurrentStep, 1});
            var_cr32.SetSelection(sel1);
            bpReader.Get(var_cr32, CR32.data());
            var_cr32.SetSelection(sel2);
            bpReader.Get(var_cr32, CR32.data() + Ny * Nx / 2);

            var_cr64.SetStepSelection({CurrentStep, 1});
            var_cr64.SetSelection(sel1);
            bpReader.Get(var_cr64, CR64.data());
            var_cr64.SetSelection(sel2);
            bpReader.Get(var_cr64, CR64.data() + Ny * Nx / 2);

            bpReader.PerformGets();

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString, currentTestData.S1);

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
                EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
            }

            ++t;
        }
        bpReader.Close();
    }
}

//******************************************************************************
// 2D 4x2 test data
//******************************************************************************

TEST_F(BPWriteReadMultiblockTest, ADIOS2BPWriteReadMultiblock2D4x2)
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
    const std::string fname("ADIOS2BPWriteReadMultiblock2D4x2Test_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteReadMultiblock2D4x2Test.bp");
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    // Make a 2D selection to describe the local dimensions of the
    // variable we write and its offsets in the global spaces
    const adios2::Box<adios2::Dims> sel1({0, static_cast<size_t>(mpiRank * Nx)}, {Ny / 2, Nx});

    const adios2::Box<adios2::Dims> sel2({Ny / 2, static_cast<size_t>(mpiRank * Nx)},
                                         {Ny - Ny / 2, Nx});

    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 2D variables (4 * (NumberOfProcess * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            adios2::Dims shape{static_cast<size_t>(Ny), static_cast<size_t>(mpiSize * Nx)};
            adios2::Dims start{static_cast<size_t>(0), static_cast<size_t>(mpiRank * Nx)};
            adios2::Dims count{static_cast<size_t>(Ny), static_cast<size_t>(Nx)};
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            EXPECT_TRUE(var_i8);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            EXPECT_TRUE(var_i16);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            EXPECT_TRUE(var_i64);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            EXPECT_TRUE(var_u8);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            EXPECT_TRUE(var_u16);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            EXPECT_TRUE(var_u32);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            EXPECT_TRUE(var_u64);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            EXPECT_TRUE(var_r32);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
            EXPECT_TRUE(var_r64);
            auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count);
            EXPECT_TRUE(var_cr32);
            auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count);
            EXPECT_TRUE(var_cr64);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("BPFile");
        }

        io.AddTransport("file");
        io.SetParameter("AggregatorRatio", "1");

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
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

            bpWriter.BeginStep();

            var_i8.SetSelection(sel1);
            bpWriter.Put(var_i8, currentTestData.I8.data());
            var_i8.SetSelection(sel2);
            bpWriter.Put(var_i8, currentTestData.I8.data() + Ny * Nx / 2);

            var_i16.SetSelection(sel1);
            bpWriter.Put(var_i16, currentTestData.I16.data());
            var_i16.SetSelection(sel2);
            bpWriter.Put(var_i16, currentTestData.I16.data() + Ny * Nx / 2);

            var_i32.SetSelection(sel1);
            bpWriter.Put(var_i32, currentTestData.I32.data());
            var_i32.SetSelection(sel2);
            bpWriter.Put(var_i32, currentTestData.I32.data() + Ny * Nx / 2);

            var_i64.SetSelection(sel1);
            bpWriter.Put(var_i64, currentTestData.I64.data());
            var_i64.SetSelection(sel2);
            bpWriter.Put(var_i64, currentTestData.I64.data() + Ny * Nx / 2);

            var_u8.SetSelection(sel1);
            bpWriter.Put(var_u8, currentTestData.U8.data());
            var_u8.SetSelection(sel2);
            bpWriter.Put(var_u8, currentTestData.U8.data() + Ny * Nx / 2);

            var_u16.SetSelection(sel1);
            bpWriter.Put(var_u16, currentTestData.U16.data());
            var_u16.SetSelection(sel2);
            bpWriter.Put(var_u16, currentTestData.U16.data() + Ny * Nx / 2);

            var_u32.SetSelection(sel1);
            bpWriter.Put(var_u32, currentTestData.U32.data());
            var_u32.SetSelection(sel2);
            bpWriter.Put(var_u32, currentTestData.U32.data() + Ny * Nx / 2);

            var_u64.SetSelection(sel1);
            bpWriter.Put(var_u64, currentTestData.U64.data());
            var_u64.SetSelection(sel2);
            bpWriter.Put(var_u64, currentTestData.U64.data() + Ny * Nx / 2);

            var_r32.SetSelection(sel1);
            bpWriter.Put(var_r32, currentTestData.R32.data());
            var_r32.SetSelection(sel2);
            bpWriter.Put(var_r32, currentTestData.R32.data() + Ny * Nx / 2);

            var_r64.SetSelection(sel1);
            bpWriter.Put(var_r64, currentTestData.R64.data());
            var_r64.SetSelection(sel2);
            bpWriter.Put(var_r64, currentTestData.R64.data() + Ny * Nx / 2);

            var_cr32.SetSelection(sel1);
            bpWriter.Put(var_cr32, currentTestData.CR32.data());
            var_cr32.SetSelection(sel2);
            bpWriter.Put(var_cr32, currentTestData.CR32.data() + Ny * Nx / 2);

            var_cr64.SetSelection(sel1);
            bpWriter.Put(var_cr64, currentTestData.CR64.data());
            var_cr64.SetSelection(sel2);
            bpWriter.Put(var_cr64, currentTestData.CR64.data() + Ny * Nx / 2);

            bpWriter.EndStep();
        }

        // Close the file
        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i8.Steps(), NSteps);
        EXPECT_EQ(var_i8.Shape()[0], Ny);
        EXPECT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i16.Steps(), NSteps);
        EXPECT_EQ(var_i16.Shape()[0], Ny);
        EXPECT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i32.Steps(), NSteps);
        EXPECT_EQ(var_i32.Shape()[0], Ny);
        EXPECT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i64.Steps(), NSteps);
        EXPECT_EQ(var_i64.Shape()[0], Ny);
        EXPECT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u8.Steps(), NSteps);
        EXPECT_EQ(var_u8.Shape()[0], Ny);
        EXPECT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u16.Steps(), NSteps);
        EXPECT_EQ(var_u16.Shape()[0], Ny);
        EXPECT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u32.Steps(), NSteps);
        EXPECT_EQ(var_u32.Shape()[0], Ny);
        EXPECT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_u64.Steps(), NSteps);
        EXPECT_EQ(var_u64.Shape()[0], Ny);
        EXPECT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_r32.Steps(), NSteps);
        EXPECT_EQ(var_r32.Shape()[0], Ny);
        EXPECT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_r64.Steps(), NSteps);
        EXPECT_EQ(var_r64.Shape()[0], Ny);
        EXPECT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
        EXPECT_TRUE(var_cr32);
        EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_cr32.Steps(), NSteps);
        EXPECT_EQ(var_cr32.Shape()[0], Ny);
        EXPECT_EQ(var_cr32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
        EXPECT_TRUE(var_cr64);
        EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_cr64.Steps(), NSteps);
        EXPECT_EQ(var_cr64.Shape()[0], Ny);
        EXPECT_EQ(var_cr64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

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
        std::array<std::complex<float>, Nx * Ny> CR32;
        std::array<std::complex<double>, Nx * Ny> CR64;

        size_t t = 0;

        for (size_t CurrentStep = 0; CurrentStep < NSteps; CurrentStep++)
        {
            const std::vector<adios2::Variable<int8_t>::Info> i8Info =
                bpReader.BlocksInfo(var_i8, CurrentStep);
            const std::vector<adios2::Variable<int16_t>::Info> i16Info =
                bpReader.BlocksInfo(var_i16, CurrentStep);
            const std::vector<adios2::Variable<int32_t>::Info> i32Info =
                bpReader.BlocksInfo(var_i32, CurrentStep);
            const std::vector<adios2::Variable<int64_t>::Info> i64Info =
                bpReader.BlocksInfo(var_i64, CurrentStep);
            const std::vector<adios2::Variable<uint8_t>::Info> u8Info =
                bpReader.BlocksInfo(var_u8, CurrentStep);
            const std::vector<adios2::Variable<uint16_t>::Info> u16Info =
                bpReader.BlocksInfo(var_u16, CurrentStep);
            const std::vector<adios2::Variable<uint32_t>::Info> u32Info =
                bpReader.BlocksInfo(var_u32, CurrentStep);
            const std::vector<adios2::Variable<uint64_t>::Info> u64Info =
                bpReader.BlocksInfo(var_u64, CurrentStep);
            const std::vector<adios2::Variable<float>::Info> r32Info =
                bpReader.BlocksInfo(var_r32, CurrentStep);
            const std::vector<adios2::Variable<double>::Info> r64Info =
                bpReader.BlocksInfo(var_r64, CurrentStep);

            const std::vector<adios2::Variable<std::complex<float>>::Info> cr32Info =
                bpReader.BlocksInfo(var_cr32, CurrentStep);
            const std::vector<adios2::Variable<std::complex<double>>::Info> cr64Info =
                bpReader.BlocksInfo(var_cr64, CurrentStep);

            EXPECT_EQ(i8Info.size(), 2 * mpiSize);
            EXPECT_EQ(i16Info.size(), 2 * mpiSize);
            EXPECT_EQ(i32Info.size(), 2 * mpiSize);
            EXPECT_EQ(i64Info.size(), 2 * mpiSize);
            EXPECT_EQ(u8Info.size(), 2 * mpiSize);
            EXPECT_EQ(u16Info.size(), 2 * mpiSize);
            EXPECT_EQ(u32Info.size(), 2 * mpiSize);
            EXPECT_EQ(u64Info.size(), 2 * mpiSize);
            EXPECT_EQ(r32Info.size(), 2 * mpiSize);
            EXPECT_EQ(r64Info.size(), 2 * mpiSize);
            EXPECT_EQ(cr32Info.size(), 2 * mpiSize);
            EXPECT_EQ(cr64Info.size(), 2 * mpiSize);

            for (size_t i = 0; i < 2 * static_cast<size_t>(mpiSize); ++i)
            {
                EXPECT_FALSE(i8Info[0].IsValue);
                EXPECT_FALSE(i16Info[0].IsValue);
                EXPECT_FALSE(i32Info[0].IsValue);
                EXPECT_FALSE(i64Info[0].IsValue);
                EXPECT_FALSE(u8Info[0].IsValue);
                EXPECT_FALSE(u16Info[0].IsValue);
                EXPECT_FALSE(u32Info[0].IsValue);
                EXPECT_FALSE(u64Info[0].IsValue);
                EXPECT_FALSE(r32Info[0].IsValue);
                EXPECT_FALSE(r64Info[0].IsValue);
                EXPECT_FALSE(cr32Info[0].IsValue);
                EXPECT_FALSE(cr64Info[0].IsValue);

                const size_t inRank = i / 2;
                int8_t i8Min, i8Max;
                int16_t i16Min, i16Max;
                int32_t i32Min, i32Max;
                int64_t i64Min, i64Max;
                uint8_t u8Min, u8Max;
                uint16_t u16Min, u16Max;
                uint32_t u32Min, u32Max;
                uint64_t u64Min, u64Max;
                float r32Min, r32Max;
                double r64Min, r64Max;
                std::complex<float> cr32Min, cr32Max;
                std::complex<double> cr64Min, cr64Max;

                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(t), static_cast<int>(inRank), mpiSize);

                if (i % 2 == 0)
                {
                    EXPECT_EQ(i8Info[i].Start[0], 0);
                    EXPECT_EQ(i8Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i8Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(i8Info[i].Count[1], Nx);
                    EXPECT_EQ(i8Info[i].WriterID, inRank);

                    EXPECT_EQ(i16Info[i].Start[0], 0);
                    EXPECT_EQ(i16Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i16Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(i16Info[i].Count[1], Nx);
                    EXPECT_EQ(i16Info[i].WriterID, inRank);

                    EXPECT_EQ(i32Info[i].Start[0], 0);
                    EXPECT_EQ(i32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i32Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(i32Info[i].Count[1], Nx);
                    EXPECT_EQ(i32Info[i].WriterID, inRank);

                    EXPECT_EQ(i64Info[i].Start[0], 0);
                    EXPECT_EQ(i64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i64Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(i64Info[i].Count[1], Nx);
                    EXPECT_EQ(i64Info[i].WriterID, inRank);

                    EXPECT_EQ(u8Info[i].Start[0], 0);
                    EXPECT_EQ(u8Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u8Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(u8Info[i].Count[1], Nx);
                    EXPECT_EQ(u8Info[i].WriterID, inRank);

                    EXPECT_EQ(u16Info[i].Start[0], 0);
                    EXPECT_EQ(u16Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u16Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(u16Info[i].Count[1], Nx);
                    EXPECT_EQ(u16Info[i].WriterID, inRank);

                    EXPECT_EQ(u32Info[i].Start[0], 0);
                    EXPECT_EQ(u32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u32Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(u32Info[i].Count[1], Nx);
                    EXPECT_EQ(u32Info[i].WriterID, inRank);

                    EXPECT_EQ(u64Info[i].Start[0], 0);
                    EXPECT_EQ(u64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u64Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(u64Info[i].Count[1], Nx);
                    EXPECT_EQ(u64Info[i].WriterID, inRank);

                    EXPECT_EQ(r32Info[i].Start[0], 0);
                    EXPECT_EQ(r32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(r32Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(r32Info[i].Count[1], Nx);
                    EXPECT_EQ(r32Info[i].WriterID, inRank);

                    EXPECT_EQ(r64Info[i].Start[0], 0);
                    EXPECT_EQ(r64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(r64Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(r64Info[i].Count[1], Nx);
                    EXPECT_EQ(r64Info[i].WriterID, inRank);

                    EXPECT_EQ(cr32Info[i].Start[0], 0);
                    EXPECT_EQ(cr32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(cr32Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(cr32Info[i].Count[1], Nx);
                    EXPECT_EQ(cr32Info[i].WriterID, inRank);

                    EXPECT_EQ(cr64Info[i].Start[0], 0);
                    EXPECT_EQ(cr64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(cr64Info[i].Count[0], Ny / 2);
                    EXPECT_EQ(cr64Info[i].Count[1], Nx);
                    EXPECT_EQ(cr64Info[i].WriterID, inRank);

                    i8Min = *std::min_element(currentTestData.I8.begin(),
                                              currentTestData.I8.begin() + Ny * Nx / 2);
                    i8Max = *std::max_element(currentTestData.I8.begin(),
                                              currentTestData.I8.begin() + Ny * Nx / 2);
                    i16Min = *std::min_element(currentTestData.I16.begin(),
                                               currentTestData.I16.begin() + Ny * Nx / 2);
                    i16Max = *std::max_element(currentTestData.I16.begin(),
                                               currentTestData.I16.begin() + Ny * Nx / 2);
                    i32Min = *std::min_element(currentTestData.I32.begin(),
                                               currentTestData.I32.begin() + Ny * Nx / 2);
                    i32Max = *std::max_element(currentTestData.I32.begin(),
                                               currentTestData.I32.begin() + Ny * Nx / 2);
                    i64Min = *std::min_element(currentTestData.I64.begin(),
                                               currentTestData.I64.begin() + Ny * Nx / 2);
                    i64Max = *std::max_element(currentTestData.I64.begin(),
                                               currentTestData.I64.begin() + Ny * Nx / 2);
                    u8Min = *std::min_element(currentTestData.U8.begin(),
                                              currentTestData.U8.begin() + Ny * Nx / 2);
                    u8Max = *std::max_element(currentTestData.U8.begin(),
                                              currentTestData.U8.begin() + Ny * Nx / 2);
                    u16Min = *std::min_element(currentTestData.U16.begin(),
                                               currentTestData.U16.begin() + Ny * Nx / 2);
                    u16Max = *std::max_element(currentTestData.U16.begin(),
                                               currentTestData.U16.begin() + Ny * Nx / 2);
                    u32Min = *std::min_element(currentTestData.U32.begin(),
                                               currentTestData.U32.begin() + Ny * Nx / 2);
                    u32Max = *std::max_element(currentTestData.U32.begin(),
                                               currentTestData.U32.begin() + Ny * Nx / 2);
                    u64Min = *std::min_element(currentTestData.U64.begin(),
                                               currentTestData.U64.begin() + Ny * Nx / 2);
                    u64Max = *std::max_element(currentTestData.U64.begin(),
                                               currentTestData.U64.begin() + Ny * Nx / 2);
                    r32Min = *std::min_element(currentTestData.R32.begin(),
                                               currentTestData.R32.begin() + Ny * Nx / 2);
                    r32Max = *std::max_element(currentTestData.R32.begin(),
                                               currentTestData.R32.begin() + Ny * Nx / 2);
                    r64Min = *std::min_element(currentTestData.R64.begin(),
                                               currentTestData.R64.begin() + Ny * Nx / 2);
                    r64Max = *std::max_element(currentTestData.R64.begin(),
                                               currentTestData.R64.begin() + Ny * Nx / 2);

                    cr32Min = currentTestData.CR32.front();
                    cr32Max = currentTestData.CR32.front();
                    for (auto it = currentTestData.CR32.begin();
                         it != currentTestData.CR32.begin() + Ny * Nx / 2; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr32Min))
                        {
                            cr32Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr32Max))
                        {
                            cr32Max = *it;
                        }
                    }

                    cr64Min = currentTestData.CR64.front();
                    cr64Max = currentTestData.CR64.front();
                    for (auto it = currentTestData.CR64.begin();
                         it != currentTestData.CR64.begin() + Ny * Nx / 2; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr64Min))
                        {
                            cr64Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr64Max))
                        {
                            cr64Max = *it;
                        }
                    }
                }
                else
                {
                    EXPECT_EQ(i8Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(i8Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i8Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(i8Info[i].Count[1], Nx);

                    EXPECT_EQ(i16Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(i16Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i16Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(i16Info[i].Count[1], Nx);

                    EXPECT_EQ(i32Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(i32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i32Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(i32Info[i].Count[1], Nx);

                    EXPECT_EQ(i64Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(i64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(i64Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(i64Info[i].Count[1], Nx);

                    EXPECT_EQ(u8Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(u8Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u8Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(u8Info[i].Count[1], Nx);

                    EXPECT_EQ(u16Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(u16Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u16Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(u16Info[i].Count[1], Nx);

                    EXPECT_EQ(u32Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(u32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u32Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(u32Info[i].Count[1], Nx);

                    EXPECT_EQ(u64Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(u64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(u64Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(u64Info[i].Count[1], Nx);

                    EXPECT_EQ(r32Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(r32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(r32Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(r32Info[i].Count[1], Nx);

                    EXPECT_EQ(r64Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(r64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(r64Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(r64Info[i].Count[1], Nx);

                    EXPECT_EQ(cr32Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(cr32Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(cr32Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(cr32Info[i].Count[1], Nx);

                    EXPECT_EQ(cr64Info[i].Start[0], Ny / 2);
                    EXPECT_EQ(cr64Info[i].Start[1], inRank * Nx);
                    EXPECT_EQ(cr64Info[i].Count[0], Ny - Ny / 2);
                    EXPECT_EQ(cr64Info[i].Count[1], Nx);

                    i8Min = *std::min_element(currentTestData.I8.begin() + Ny * Nx / 2,
                                              currentTestData.I8.begin() + Ny * Nx);
                    i8Max = *std::max_element(currentTestData.I8.begin() + Ny * Nx / 2,
                                              currentTestData.I8.begin() + Ny * Nx);

                    i16Min = *std::min_element(currentTestData.I16.begin() + Ny * Nx / 2,
                                               currentTestData.I16.begin() + Ny * Nx);
                    i16Max = *std::max_element(currentTestData.I16.begin() + Ny * Nx / 2,
                                               currentTestData.I16.begin() + Ny * Nx);

                    i32Min = *std::min_element(currentTestData.I32.begin() + Ny * Nx / 2,
                                               currentTestData.I32.begin() + Ny * Nx);
                    i32Max = *std::max_element(currentTestData.I32.begin() + Ny * Nx / 2,
                                               currentTestData.I32.begin() + Ny * Nx);

                    i64Min = *std::min_element(currentTestData.I64.begin() + Ny * Nx / 2,
                                               currentTestData.I64.begin() + Ny * Nx);
                    i64Max = *std::max_element(currentTestData.I64.begin() + Ny * Nx / 2,
                                               currentTestData.I64.begin() + Ny * Nx);

                    u8Min = *std::min_element(currentTestData.U8.begin() + Ny * Nx / 2,
                                              currentTestData.U8.begin() + Ny * Nx);
                    u8Max = *std::max_element(currentTestData.U8.begin() + Ny * Nx / 2,
                                              currentTestData.U8.begin() + Ny * Nx);

                    u16Min = *std::min_element(currentTestData.U16.begin() + Ny * Nx / 2,
                                               currentTestData.U16.begin() + Ny * Nx);
                    u16Max = *std::max_element(currentTestData.U16.begin() + Ny * Nx / 2,
                                               currentTestData.U16.begin() + Ny * Nx);

                    u32Min = *std::min_element(currentTestData.U32.begin() + Ny * Nx / 2,
                                               currentTestData.U32.begin() + Ny * Nx);
                    u32Max = *std::max_element(currentTestData.U32.begin() + Ny * Nx / 2,
                                               currentTestData.U32.begin() + Ny * Nx);

                    u64Min = *std::min_element(currentTestData.U64.begin() + Ny * Nx / 2,
                                               currentTestData.U64.begin() + Ny * Nx);
                    u64Max = *std::max_element(currentTestData.U64.begin() + Ny * Nx / 2,
                                               currentTestData.U64.begin() + Ny * Nx);

                    r32Min = *std::min_element(currentTestData.R32.begin() + Ny * Nx / 2,
                                               currentTestData.R32.begin() + Ny * Nx);
                    r32Max = *std::max_element(currentTestData.R32.begin() + Ny * Nx / 2,
                                               currentTestData.R32.begin() + Ny * Nx);

                    r64Min = *std::min_element(currentTestData.R64.begin() + Ny * Nx / 2,
                                               currentTestData.R64.begin() + Ny * Nx);
                    r64Max = *std::max_element(currentTestData.R64.begin() + Ny * Nx / 2,
                                               currentTestData.R64.begin() + Ny * Nx);

                    cr32Min = currentTestData.CR32[Ny * Nx / 2];
                    cr32Max = currentTestData.CR32[Ny * Nx / 2];
                    for (auto it = currentTestData.CR32.begin() + Ny * Nx / 2;
                         it != currentTestData.CR32.begin() + Ny * Nx; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr32Min))
                        {
                            cr32Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr32Max))
                        {
                            cr32Max = *it;
                        }
                    }

                    cr64Min = currentTestData.CR64[Ny * Nx / 2];
                    cr64Max = currentTestData.CR64[Ny * Nx / 2];
                    for (auto it = currentTestData.CR64.begin() + Ny * Nx / 2;
                         it != currentTestData.CR64.begin() + Ny * Nx; ++it)
                    {
                        if (std::norm(*it) < std::norm(cr64Min))
                        {
                            cr64Min = *it;
                            continue;
                        }
                        if (std::norm(*it) > std::norm(cr64Max))
                        {
                            cr64Max = *it;
                        }
                    }
                }
                EXPECT_EQ(i8Info[i].Min, i8Min);
                EXPECT_EQ(i8Info[i].Max, i8Max);
                EXPECT_EQ(i16Info[i].Min, i16Min);
                EXPECT_EQ(i16Info[i].Max, i16Max);
                EXPECT_EQ(i32Info[i].Min, i32Min);
                EXPECT_EQ(i32Info[i].Max, i32Max);
                EXPECT_EQ(i64Info[i].Min, i64Min);
                EXPECT_EQ(i64Info[i].Max, i64Max);

                EXPECT_EQ(u8Info[i].Min, u8Min);
                EXPECT_EQ(u8Info[i].Max, u8Max);
                EXPECT_EQ(u16Info[i].Min, u16Min);
                EXPECT_EQ(u16Info[i].Max, u16Max);
                EXPECT_EQ(u32Info[i].Min, u32Min);
                EXPECT_EQ(u32Info[i].Max, u32Max);
                EXPECT_EQ(u64Info[i].Min, u64Min);
                EXPECT_EQ(u64Info[i].Max, u64Max);

                EXPECT_EQ(r32Info[i].Min, r32Min);
                EXPECT_EQ(r32Info[i].Max, r32Max);
                EXPECT_EQ(r64Info[i].Min, r64Min);
                EXPECT_EQ(r64Info[i].Max, r64Max);
            }

            var_i8.SetStepSelection({CurrentStep, 1});
            var_i8.SetSelection(sel1);
            bpReader.Get(var_i8, I8.data());
            var_i8.SetSelection(sel2);
            bpReader.Get(var_i8, I8.data() + Ny * Nx / 2);

            var_i16.SetStepSelection({CurrentStep, 1});
            var_i16.SetSelection(sel1);
            bpReader.Get(var_i16, I16.data());
            var_i16.SetSelection(sel2);
            bpReader.Get(var_i16, I16.data() + Ny * Nx / 2);

            var_i32.SetStepSelection({CurrentStep, 1});
            var_i32.SetSelection(sel1);
            bpReader.Get(var_i32, I32.data());
            var_i32.SetSelection(sel2);
            bpReader.Get(var_i32, I32.data() + Ny * Nx / 2);

            var_i64.SetStepSelection({CurrentStep, 1});
            var_i64.SetSelection(sel1);
            bpReader.Get(var_i64, I64.data());
            var_i64.SetSelection(sel2);
            bpReader.Get(var_i64, I64.data() + Ny * Nx / 2);

            var_u8.SetStepSelection({CurrentStep, 1});
            var_u8.SetSelection(sel1);
            bpReader.Get(var_u8, U8.data());
            var_u8.SetSelection(sel2);
            bpReader.Get(var_u8, U8.data() + Ny * Nx / 2);

            var_u16.SetStepSelection({CurrentStep, 1});
            var_u16.SetSelection(sel1);
            bpReader.Get(var_u16, U16.data());
            var_u16.SetSelection(sel2);
            bpReader.Get(var_u16, U16.data() + Ny * Nx / 2);

            var_u32.SetStepSelection({CurrentStep, 1});
            var_u32.SetSelection(sel1);
            bpReader.Get(var_u32, U32.data());
            var_u32.SetSelection(sel2);
            bpReader.Get(var_u32, U32.data() + Ny * Nx / 2);

            var_u64.SetStepSelection({CurrentStep, 1});
            var_u64.SetSelection(sel1);
            bpReader.Get(var_u64, U64.data());
            var_u64.SetSelection(sel2);
            bpReader.Get(var_u64, U64.data() + Ny * Nx / 2);

            var_r32.SetStepSelection({CurrentStep, 1});
            var_r32.SetSelection(sel1);
            bpReader.Get(var_r32, R32.data());
            var_r32.SetSelection(sel2);
            bpReader.Get(var_r32, R32.data() + Ny * Nx / 2);

            var_r64.SetStepSelection({CurrentStep, 1});
            var_r64.SetSelection(sel1);
            bpReader.Get(var_r64, R64.data());
            var_r64.SetSelection(sel2);
            bpReader.Get(var_r64, R64.data() + Ny * Nx / 2);

            var_cr32.SetStepSelection({CurrentStep, 1});
            var_cr32.SetSelection(sel1);
            bpReader.Get(var_cr32, CR32.data());
            var_cr32.SetSelection(sel2);
            bpReader.Get(var_cr32, CR32.data() + Ny * Nx / 2);

            var_cr64.SetStepSelection({CurrentStep, 1});
            var_cr64.SetSelection(sel1);
            bpReader.Get(var_cr64, CR64.data());
            var_cr64.SetSelection(sel2);
            bpReader.Get(var_cr64, CR64.data() + Ny * Nx / 2);

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            bpReader.PerformGets();

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
                EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
            }

            ++t;
        }
        bpReader.Close();
    }
}

//******************************************************************************
// Test flushing data within the step and that read works properly for all
// blocks across all flushes
//******************************************************************************

TEST_F(BPWriteReadMultiblockTest, MultiblockPerformDataWrite)
{
    if (engineName != "BP5")
    {
        std::cout << "Engine " << engineName << " is not tested for this feature." << std::endl;
        return;
    }
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of elements per blocks (blocksize)
    const size_t Nx = 8;
    // Number of blocks per process (= number of flushes)
    const size_t Nblocks = 3;
    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    const std::string fname("MultiblockPerformDataWrite_MPI.bp");
#else
    adios2::ADIOS adios;
    const std::string fname("MultiblockPerformDataWrite.bp");
#endif
    /* Write output */
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        const adios2::Dims shape{static_cast<size_t>(mpiSize), static_cast<size_t>(Nx * Nblocks)};
        const adios2::Dims start{static_cast<size_t>(mpiRank), 0};
        const adios2::Dims count{1, Nx};

        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("BPFile");
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();

            for (size_t b = 0; b < Nblocks; ++b)
            {
                // Generate test data for each process / block uniquely
                int t = static_cast<int>(step * Nblocks + b);
                SmallTestData currentTestData =
                    generateNewSmallTestData(m_TestData, t, mpiRank, mpiSize);

                const adios2::Box<adios2::Dims> sel({(size_t)mpiRank, b * Nx}, {1, Nx});
                var_i32.SetSelection(sel);
                bpWriter.Put(var_i32, currentTestData.I32.data());

                bpWriter.PerformDataWrite();
            }
            bpWriter.EndStep();
        }

        // Close the file
        bpWriter.Close();
    }

    /* Read back each step, each block and check the unique values.
       Different blocks in each step are coming from different flushes.
    */
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i32.Steps(), NSteps);
        EXPECT_EQ(var_i32.Shape()[0], mpiSize);
        EXPECT_EQ(var_i32.Shape()[1], Nx * Nblocks);

        SmallTestData testData;
        std::array<int32_t, Nx> I32;

        const auto i32AllInfo = bpReader.AllStepsBlocksInfo(var_i32);
        EXPECT_EQ(i32AllInfo.size(), NSteps);

        for (size_t step = 0; step < NSteps; step++)
        {
            var_i32.SetStepSelection({step, 1});
            for (size_t b = 0; b < Nblocks; ++b)
            {
                std::cout << "Read step " << step << " block=" << b << std::endl;
                // Generate test data for each process / block uniquely
                int t = static_cast<int>(step * Nblocks + b);
                SmallTestData currentTestData =
                    generateNewSmallTestData(m_TestData, t, mpiRank, mpiSize);

                const adios2::Box<adios2::Dims> sel({(size_t)mpiRank, b * Nx}, {1, Nx});
                var_i32.SetSelection(sel);
                bpReader.Get(var_i32, I32.data(), adios2::Mode::Sync);

                /* check content of a single block */
                for (size_t i = 0; i < Nx; ++i)
                {
                    std::stringstream ss;
                    ss << "step=" << step << " block=" << b << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();
                    EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                }
            }
        }
        bpReader.Close();
    }
}

//******************************************************************************
// Test reading data where some processes do not contribute to the data
// and some blocks are null
//******************************************************************************

TEST_F(BPWriteReadMultiblockTest, MultiblockNullBlocks)
{
    // Each process would write a 2x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of elements per blocks (blocksize)
    const size_t Nx = 8;
    // Number of blocks per process (= number of flushes)
    const size_t Nblocks = 3;
    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("MultiblockNullBlocks_MPI.bp");
#else
    const std::string fname("MultiblockNullBlocks.bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    /* Write */
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        adios2::Dims shape{static_cast<size_t>(mpiSize), static_cast<size_t>(Nx * (Nblocks - 1))};
        adios2::Dims start{static_cast<size_t>(mpiRank), 0};
        adios2::Dims count{1, Nx};

        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();

            size_t nb = 0;
            for (size_t b = 0; b < Nblocks; ++b)
            {
                // Generate test data for each process / block uniquely
                int t = static_cast<int>(step * Nblocks + b);
                SmallTestData currentTestData =
                    generateNewSmallTestData(m_TestData, t, mpiRank, mpiSize);

                // the first block does not contribute to the variable's data
                if (b == 0)
                {
                    std::array<int32_t, Nx> I32_empty;
                    var_i32.SetSelection(
                        adios2::Box<adios2::Dims>({(size_t)mpiRank, b * Nx}, {0, 0}));
                    bpWriter.Put(var_i32, I32_empty.data());
                }
                else
                {
                    ++nb;
                    start = {static_cast<size_t>(mpiRank), static_cast<size_t>(Nx * (nb - 1))};
                    count = {1, Nx};
                    var_i32.SetSelection({start, count});
                    bpWriter.Put(var_i32, currentTestData.I32.data(), adios2::Mode::Sync);
                }

                bpWriter.PerformDataWrite();
            }
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }
    // Read and check correctness
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        EXPECT_EQ(var_i32.Steps(), NSteps);
        EXPECT_EQ(var_i32.Shape()[0], mpiSize);
        EXPECT_EQ(var_i32.Shape()[1], Nx * (Nblocks - 1));

        SmallTestData testData;
        std::array<int32_t, Nx> I32;

        const auto i32AllInfo = bpReader.AllStepsBlocksInfo(var_i32);
        EXPECT_EQ(i32AllInfo.size(), NSteps);

        for (size_t step = 0; step < NSteps; step++)
        {
            var_i32.SetStepSelection({step, 1});
            for (size_t b = 1; b < Nblocks; ++b)
            {
                std::cout << "Read step " << step << " block=" << b << std::endl;
                // Generate test data for each process / block uniquely
                int t = static_cast<int>(step * Nblocks + b);
                SmallTestData currentTestData =
                    generateNewSmallTestData(m_TestData, t, mpiRank, mpiSize);

                // Block 0 was not written so all blocks are shifted back
                const adios2::Box<adios2::Dims> sel({(size_t)mpiRank, (b - 1) * Nx}, {1, Nx});
                var_i32.SetSelection(sel);
                bpReader.Get(var_i32, I32.data(), adios2::Mode::Sync);

                for (size_t i = 0; i < Nx; ++i)
                {
                    std::stringstream ss;
                    ss << "step=" << step << " block=" << b << " i=" << i << " rank=" << mpiRank;
                    std::string msg = ss.str();
                    EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                }
            }
        }
        bpReader.Close();
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
