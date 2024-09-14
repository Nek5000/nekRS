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

class InlineWriteRead : public ::testing::Test
{
public:
    InlineWriteRead() = default;

    SmallTestData m_TestData;
};

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

// helper
template <class T>
typename adios2::Variable<T>::Info setSelection(adios2::Variable<T> &var_i8, size_t step,
                                                adios2::Engine &inlineReader)
{
    var_i8.SetStepSelection({step, 1});
    auto blocksInfo = inlineReader.BlocksInfo(var_i8, step);
    // ASSERT_EQ(blocksInfo.size(), 1);
    // must be copied, since blocksInfo will go out of scope and be destroyed.
    typename adios2::Variable<T>::Info info = blocksInfo[0];
    var_i8.SetBlockSelection(info.BlockID);
    inlineReader.Get(var_i8, info);
    return info;
}

template <class T>
void testBlocksInfo(adios2::Variable<T> &var, size_t step, adios2::Engine &inlineReader)
{
    var.SetStepSelection({step, 1});
    auto blocksInfo = inlineReader.BlocksInfo(var, step);
    ASSERT_EQ(blocksInfo.size(), 1);
}

TEST_F(InlineWriteRead, InlineWriteRead1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("InlineWriteRead1D8");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
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

        // Create the Engine
        io.SetEngine("Inline");

        // writerID parameter makes sure the reader can find the writer.
        io.SetParameter("verbose", "4");

        adios2::Engine inlineWriter = io.Open(fname + "_write", adios2::Mode::Write);

        adios2::Engine inlineReader = io.Open(fname + "_read", adios2::Mode::Read);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData testData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto wvar_iString = io.InquireVariable<std::string>("iString");
            auto wvar_i8 = io.InquireVariable<int8_t>("i8");
            auto wvar_i16 = io.InquireVariable<int16_t>("i16");
            auto wvar_i32 = io.InquireVariable<int32_t>("i32");
            auto wvar_i64 = io.InquireVariable<int64_t>("i64");
            auto wvar_u8 = io.InquireVariable<uint8_t>("u8");
            auto wvar_u16 = io.InquireVariable<uint16_t>("u16");
            auto wvar_u32 = io.InquireVariable<uint32_t>("u32");
            auto wvar_u64 = io.InquireVariable<uint64_t>("u64");
            auto wvar_r32 = io.InquireVariable<float>("r32");
            auto wvar_r64 = io.InquireVariable<double>("r64");
            auto wvar_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            auto wvar_cr64 = io.InquireVariable<std::complex<double>>("cr64");

            // Make a 1D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});

            EXPECT_THROW(wvar_iString.SetSelection(sel), std::invalid_argument);
            wvar_i8.SetSelection(sel);
            wvar_i16.SetSelection(sel);
            wvar_i32.SetSelection(sel);
            wvar_i64.SetSelection(sel);
            wvar_u8.SetSelection(sel);
            wvar_u16.SetSelection(sel);
            wvar_u32.SetSelection(sel);
            wvar_u64.SetSelection(sel);
            wvar_r32.SetSelection(sel);
            wvar_r64.SetSelection(sel);
            wvar_cr32.SetSelection(sel);
            wvar_cr64.SetSelection(sel);

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            inlineWriter.BeginStep();

            inlineWriter.Put(wvar_iString, testData.S1);
            inlineWriter.Put(wvar_i8, testData.I8.data());
            inlineWriter.Put(wvar_i16, testData.I16.data());
            inlineWriter.Put(wvar_i32, testData.I32.data());
            inlineWriter.Put(wvar_i64, testData.I64.data());
            inlineWriter.Put(wvar_u8, testData.U8.data());
            inlineWriter.Put(wvar_u16, testData.U16.data());
            inlineWriter.Put(wvar_u32, testData.U32.data());
            inlineWriter.Put(wvar_u64, testData.U64.data());
            inlineWriter.Put(wvar_r32, testData.R32.data());
            inlineWriter.Put(wvar_r64, testData.R64.data());
            inlineWriter.Put(wvar_cr32, testData.CR32.data());
            inlineWriter.Put(wvar_cr64, testData.CR64.data());

            inlineWriter.EndStep();

            // start reader section
            auto var_iString = io.InquireVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            ASSERT_EQ(var_iString.Shape().size(), 0);

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

            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            EXPECT_TRUE(var_cr32);
            ASSERT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr32.Shape()[0], mpiSize * Nx);

            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
            EXPECT_TRUE(var_cr64);
            ASSERT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr64.Shape()[0], mpiSize * Nx);

            std::string IString;

            auto writerStep = inlineWriter.CurrentStep();
            auto readerStep = inlineReader.CurrentStep();
            ASSERT_EQ(writerStep, readerStep);

            // Test skipping a step on the read side
            if (step == 1)
            {
                continue;
            }

            inlineReader.BeginStep();
            inlineReader.Get(var_iString, IString);
            auto info_i8 = setSelection<int8_t>(var_i8, step, inlineReader);
            auto info_i16 = setSelection<int16_t>(var_i16, step, inlineReader);
            auto info_i32 = setSelection<int32_t>(var_i32, step, inlineReader);
            auto info_i64 = setSelection<int64_t>(var_i64, step, inlineReader);
            auto info_u8 = setSelection<uint8_t>(var_u8, step, inlineReader);
            auto info_u16 = setSelection<uint16_t>(var_u16, step, inlineReader);
            auto info_u32 = setSelection<uint32_t>(var_u32, step, inlineReader);
            auto info_u64 = setSelection<uint64_t>(var_u64, step, inlineReader);
            auto info_r32 = setSelection<float>(var_r32, step, inlineReader);
            auto info_r64 = setSelection<double>(var_r64, step, inlineReader);
            auto info_cr32 = setSelection<std::complex<float>>(var_cr32, step, inlineReader);
            auto info_cr64 = setSelection<std::complex<double>>(var_cr64, step, inlineReader);

            testBlocksInfo<int8_t>(var_i8, step, inlineReader);
            testBlocksInfo<int16_t>(var_i16, step, inlineReader);
            testBlocksInfo<int32_t>(var_i32, step, inlineReader);
            testBlocksInfo<int64_t>(var_i64, step, inlineReader);
            testBlocksInfo<uint8_t>(var_u8, step, inlineReader);
            testBlocksInfo<uint16_t>(var_u16, step, inlineReader);
            testBlocksInfo<uint32_t>(var_u32, step, inlineReader);
            testBlocksInfo<uint64_t>(var_u64, step, inlineReader);
            testBlocksInfo<float>(var_r32, step, inlineReader);
            testBlocksInfo<double>(var_r64, step, inlineReader);
            testBlocksInfo<std::complex<float>>(var_cr32, step, inlineReader);
            testBlocksInfo<std::complex<double>>(var_cr64, step, inlineReader);

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            inlineReader.EndStep();
            const int8_t *I8 = info_i8.Data();
            const int16_t *I16 = info_i16.Data();
            const int32_t *I32 = info_i32.Data();
            const int64_t *I64 = info_i64.Data();
            const uint8_t *U8 = info_u8.Data();
            const uint16_t *U16 = info_u16.Data();
            const uint32_t *U32 = info_u32.Data();
            const uint64_t *U64 = info_u64.Data();
            const float *R32 = info_r32.Data();
            const double *R64 = info_r64.Data();
            const std::complex<float> *CR32 = info_cr32.Data();
            const std::complex<double> *CR64 = info_cr64.Data();

            EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "step=" << step << " i=" << i << " rank=" << mpiRank;
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
        }
        inlineWriter.Close();
        inlineReader.Close();
    }
}

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

TEST_F(InlineWriteRead, InlineWriteRead2D2x4)
{
    // Each process would write a 2x4 array and all processes would
    // form a 2D 2 * (numberOfProcess*Nx) matrix where Nx is 4 here
    const std::string fname("InlineWriteRead2D2x4");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 4;
    const size_t Ny = 2;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
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

        // Create the Engine
        io.SetEngine("Inline");

        // writerID parameter makes sure the reader can find the writer.
        io.SetParameter("verbose", "4");

        adios2::Engine inlineWriter = io.Open(fname + "_write", adios2::Mode::Write);

        adios2::Engine inlineReader = io.Open(fname + "_read", adios2::Mode::Read);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData testData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto wvar_iString = io.InquireVariable<std::string>("iString");
            auto wvar_i8 = io.InquireVariable<int8_t>("i8");
            auto wvar_i16 = io.InquireVariable<int16_t>("i16");
            auto wvar_i32 = io.InquireVariable<int32_t>("i32");
            auto wvar_i64 = io.InquireVariable<int64_t>("i64");
            auto wvar_u8 = io.InquireVariable<uint8_t>("u8");
            auto wvar_u16 = io.InquireVariable<uint16_t>("u16");
            auto wvar_u32 = io.InquireVariable<uint32_t>("u32");
            auto wvar_u64 = io.InquireVariable<uint64_t>("u64");
            auto wvar_r32 = io.InquireVariable<float>("r32");
            auto wvar_r64 = io.InquireVariable<double>("r64");
            auto wvar_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            auto wvar_cr64 = io.InquireVariable<std::complex<double>>("cr64");

            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({0, static_cast<size_t>(mpiRank * Nx)}, {Ny, Nx});

            EXPECT_THROW(wvar_iString.SetSelection(sel), std::invalid_argument);
            wvar_i8.SetSelection(sel);
            wvar_i16.SetSelection(sel);
            wvar_i32.SetSelection(sel);
            wvar_i64.SetSelection(sel);
            wvar_u8.SetSelection(sel);
            wvar_u16.SetSelection(sel);
            wvar_u32.SetSelection(sel);
            wvar_u64.SetSelection(sel);
            wvar_r32.SetSelection(sel);
            wvar_r64.SetSelection(sel);
            wvar_cr32.SetSelection(sel);
            wvar_cr64.SetSelection(sel);

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            inlineWriter.BeginStep();

            inlineWriter.Put(wvar_iString, testData.S1);
            inlineWriter.Put(wvar_i8, testData.I8.data());
            inlineWriter.Put(wvar_i16, testData.I16.data());
            inlineWriter.Put(wvar_i32, testData.I32.data());
            inlineWriter.Put(wvar_i64, testData.I64.data());
            inlineWriter.Put(wvar_u8, testData.U8.data());
            inlineWriter.Put(wvar_u16, testData.U16.data());
            inlineWriter.Put(wvar_u32, testData.U32.data());
            inlineWriter.Put(wvar_u64, testData.U64.data());
            inlineWriter.Put(wvar_r32, testData.R32.data());
            inlineWriter.Put(wvar_r64, testData.R64.data());
            inlineWriter.Put(wvar_cr32, testData.CR32.data());
            inlineWriter.Put(wvar_cr64, testData.CR64.data());

            inlineWriter.EndStep();

            // start reader section
            auto var_iString = io.InquireVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            ASSERT_EQ(var_iString.Shape().size(), 0);

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i8.Shape()[0], Ny);
            ASSERT_EQ(var_i8.Shape()[1], mpiSize * Nx);

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

            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            EXPECT_TRUE(var_cr32);
            ASSERT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr32.Shape()[0], Ny);
            ASSERT_EQ(var_cr32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
            EXPECT_TRUE(var_cr64);
            ASSERT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr64.Shape()[0], Ny);
            ASSERT_EQ(var_cr64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            std::string IString;

            inlineReader.BeginStep();
            inlineReader.Get(var_iString, IString);
            auto info_i8 = setSelection<int8_t>(var_i8, step, inlineReader);
            auto info_i16 = setSelection<int16_t>(var_i16, step, inlineReader);
            auto info_i32 = setSelection<int32_t>(var_i32, step, inlineReader);
            auto info_i64 = setSelection<int64_t>(var_i64, step, inlineReader);
            auto info_u8 = setSelection<uint8_t>(var_u8, step, inlineReader);
            auto info_u16 = setSelection<uint16_t>(var_u16, step, inlineReader);
            auto info_u32 = setSelection<uint32_t>(var_u32, step, inlineReader);
            auto info_u64 = setSelection<uint64_t>(var_u64, step, inlineReader);
            auto info_r32 = setSelection<float>(var_r32, step, inlineReader);
            auto info_r64 = setSelection<double>(var_r64, step, inlineReader);
            auto info_cr32 = setSelection<std::complex<float>>(var_cr32, step, inlineReader);
            auto info_cr64 = setSelection<std::complex<double>>(var_cr64, step, inlineReader);

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            inlineReader.EndStep();
            const int8_t *I8 = info_i8.Data();
            const int16_t *I16 = info_i16.Data();
            const int32_t *I32 = info_i32.Data();
            const int64_t *I64 = info_i64.Data();
            const uint8_t *U8 = info_u8.Data();
            const uint16_t *U16 = info_u16.Data();
            const uint32_t *U32 = info_u32.Data();
            const uint64_t *U64 = info_u64.Data();
            const float *R32 = info_r32.Data();
            const double *R64 = info_r64.Data();
            const std::complex<float> *CR32 = info_cr32.Data();
            const std::complex<double> *CR64 = info_cr64.Data();

            EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "step=" << step << " i=" << i << " rank=" << mpiRank;
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
        }
        inlineWriter.Close();
        inlineReader.Close();
    }
}

TEST_F(InlineWriteRead, InlineWriteReadContracts)
{
    const std::string fname("InlineWriteReadContracts");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
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

            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
        }

        // Create the Engine
        io.SetEngine("Inline");

        // writerID parameter makes sure the reader can find the writer.
        io.SetParameter("verbose", "4");

        adios2::Engine inlineWriter = io.Open(fname + "_write", adios2::Mode::Write);

        adios2::Engine inlineReader = io.Open(fname + "_read", adios2::Mode::Read);

        // Want to test that the engine correctly fails when using the engine
        // incorrectly
        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData testData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto wvar_iString = io.InquireVariable<std::string>("iString");
            auto wvar_i32 = io.InquireVariable<int32_t>("i32");

            // Make a 1D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});

            EXPECT_THROW(wvar_iString.SetSelection(sel), std::invalid_argument);
            wvar_i32.SetSelection(sel);

            EXPECT_THROW(inlineWriter.EndStep(), std::runtime_error);
            ASSERT_EQ(inlineWriter.BeginStep(), adios2::StepStatus::OK);
            EXPECT_THROW(inlineWriter.BeginStep(), std::runtime_error);
            ASSERT_EQ(inlineReader.BeginStep(), adios2::StepStatus::NotReady);
            EXPECT_THROW(inlineReader.EndStep(), std::runtime_error);

            inlineWriter.Put(wvar_iString, testData.S1, adios2::Mode::Sync);
            EXPECT_THROW(inlineWriter.Put(wvar_i32, testData.I32.data(), adios2::Mode::Sync),
                         std::invalid_argument);
            inlineWriter.Put(wvar_i32, testData.I32.data());
            inlineWriter.EndStep();

            auto var_iString = io.InquireVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            ASSERT_EQ(var_iString.Shape().size(), 0);

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx);

            EXPECT_THROW(inlineReader.EndStep(), std::runtime_error);
            ASSERT_EQ(inlineReader.BeginStep(), adios2::StepStatus::OK);
            EXPECT_THROW(inlineReader.BeginStep(), std::runtime_error);
            ASSERT_EQ(inlineWriter.BeginStep(), adios2::StepStatus::NotReady);
            EXPECT_THROW(inlineWriter.EndStep(), std::runtime_error);

            std::string IString;
            inlineReader.Get(var_iString, IString);
            auto info_i32 = setSelection<int32_t>(var_i32, step, inlineReader);
            testBlocksInfo<int32_t>(var_i32, step, inlineReader);

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            const int32_t *testI32 = info_i32.Data();
            EXPECT_EQ(testI32, nullptr);
            inlineReader.EndStep();

            EXPECT_EQ(IString, currentTestData.S1);
            const int32_t *I32 = info_i32.Data();
            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "step=" << step << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
            }
        }
        inlineWriter.Close();
        inlineReader.Close();
    }
}

TEST_F(InlineWriteRead, InlineWriteReadContracts2)
{
    const std::string fname("InlineWriteReadContracts2");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
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

            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
        }

        // Create the Engine
        io.SetEngine("Inline");

        // writerID parameter makes sure the reader can find the writer.
        io.SetParameter("verbose", "4");

        adios2::Engine inlineWriter = io.Open(fname + "_write", adios2::Mode::Write);

        adios2::Engine inlineReader = io.Open(fname + "_read", adios2::Mode::Read);

        // Want to test that the engine correctly fails when using the engine
        // incorrectly
        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData testData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto wvar_iString = io.InquireVariable<std::string>("iString");
            auto wvar_i32 = io.InquireVariable<int32_t>("i32");

            // Make a 1D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});

            EXPECT_THROW(wvar_iString.SetSelection(sel), std::invalid_argument);
            wvar_i32.SetSelection(sel);

            inlineWriter.Put(wvar_iString, testData.S1);
            EXPECT_THROW(inlineWriter.Put(wvar_i32, testData.I32.data(), adios2::Mode::Sync),
                         std::invalid_argument);
            inlineWriter.Put(wvar_i32, testData.I32.data());
            inlineWriter.PerformPuts();

            auto var_iString = io.InquireVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            ASSERT_EQ(var_iString.Shape().size(), 0);

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx);

            std::string IString;
            inlineReader.Get(var_iString, IString, adios2::Mode::Sync);
            auto info_i32 = setSelection<int32_t>(var_i32, step, inlineReader);
            testBlocksInfo<int32_t>(var_i32, step, inlineReader);

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            const int32_t *testI32 = info_i32.Data();
            EXPECT_EQ(testI32, nullptr);
            inlineReader.PerformGets();

            EXPECT_EQ(IString, currentTestData.S1);
            const int32_t *I32 = info_i32.Data();
            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "step=" << step << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
            }
        }
        inlineWriter.Close();
        inlineReader.Close();
    }
}

TEST_F(InlineWriteRead, IOInvariants)
{
#if ADIOS2_USE_MPI
    int mpiRank = 0, mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");
    io.SetEngine("Inline");

    adios2::Engine inlineWriter = io.Open("writer", adios2::Mode::Write);
    EXPECT_TRUE(inlineWriter);
    // The inline engine does not support multiple writers:
    EXPECT_THROW(io.Open("another_writer", adios2::Mode::Write), std::exception);
    // The inline engine does not support append mode:
    EXPECT_THROW(io.Open("append_mode", adios2::Mode::Append), std::exception);
    adios2::Engine inlineReader = io.Open("reader", adios2::Mode::Read);
    EXPECT_TRUE(inlineReader);
    // The inline engine does not support more than 2 writers or readers:
    EXPECT_THROW(io.Open("reader2", adios2::Mode::Read), std::exception);
}

TEST_F(InlineWriteRead, PointerArithmetic)
{
    int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    adios2::ADIOS adios(MPI_COMM_WORLD);
    std::cout << "Using " << mpiSize << " ranks\n";
#else
    adios2::ADIOS adios;
    std::cout << "Using serial version of inline reader.\n";
#endif
    adios2::IO io = adios.DeclareIO("TestIO");
    io.SetEngine("Inline");

    adios2::Engine writer = io.Open("writer", adios2::Mode::Write);
    adios2::Engine reader = io.Open("reader", adios2::Mode::Read);

    size_t N = 256;
    // Test local array:
    auto local_array = io.DefineVariable<double>("u", {}, {}, {N}, adios2::ConstantDims);
    // Test global array:
    auto global_array =
        io.DefineVariable<double>("v", {mpiSize * N}, {mpiRank * N}, {N}, adios2::ConstantDims);
    for (int64_t timeStep = 0; timeStep < 2; ++timeStep)
    {
        writer.BeginStep();
        std::vector<double> sim_data(N, 3.2);
        writer.Put(local_array, sim_data.data());
        writer.Put(global_array, sim_data.data());
        writer.EndStep();

        reader.BeginStep();
        double *local_data = nullptr;
        reader.Get(local_array, &local_data);
        double *global_data = nullptr;
        reader.Get(global_array, &global_data);
        // The data is valid before Endstep():
        EXPECT_EQ(sim_data.data(), local_data);
        EXPECT_EQ(sim_data.data(), global_data);
        reader.EndStep();

        // And it is valid after EndStep():
        EXPECT_EQ(sim_data.data(), global_data);
        EXPECT_EQ(sim_data.data(), local_data);
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
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
