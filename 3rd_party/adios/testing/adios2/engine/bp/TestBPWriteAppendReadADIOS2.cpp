/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <numeric> //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPWriteAppendReadTestADIOS2 : public ::testing::Test
{
public:
    BPWriteAppendReadTestADIOS2() = default;

    SmallTestData m_TestData;
};

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

// ADIOS2 BP write, native ADIOS1 read
TEST_F(BPWriteAppendReadTestADIOS2, ADIOS2BPWriteAppendRead2D2x4)
{
    // Each process would write a 2x4 array and all processes would
    // form a 2D 2 * (numberOfProcess*Nx) matrix where Nx is 4 here

    const std::string zero = std::to_string(0);
    const std::string s1_Single = std::string("s1_Single_") + zero;
    const std::string s1_Array = std::string("s1_Array_") + zero;
    const std::string i8_Single = std::string("i8_Single_") + zero;
    const std::string i16_Single = std::string("i16_Single_") + zero;
    const std::string i32_Single = std::string("i32_Single_") + zero;
    const std::string i32_Array = std::string("i32_Array_") + zero;
    const std::string i64_Single = std::string("i64_Single_") + zero;
    const std::string u8_Single = std::string("u8_Single_") + zero;
    const std::string u16_Single = std::string("u16_Single_") + zero;
    const std::string u32_Single = std::string("u32_Single_") + zero;
    const std::string u32_Array = std::string("u32_Array_") + zero;
    const std::string u64_Single = std::string("u64_Single_") + zero;
    const std::string r32_Single = std::string("r32_Single_") + zero;
    const std::string r32_Array = std::string("r32_Array_") + zero;
    const std::string r64_Single = std::string("r64_Single_") + zero;
    const std::string cr32_Single = std::string("cr32_Single_") + zero;
    const std::string cr32_Array = std::string("cr32_Array_") + zero;
    const std::string cr64_Single = std::string("cr64_Single_") + zero;
    SmallTestData attributeTestData = generateNewSmallTestData(m_TestData, 0, 0, 0);

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
    const std::string fname("ADIOS2BPWriteAppendRead2D2x4Test_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteAppendRead2D2x4Test.bp");
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("WriteIO");

        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
            const adios2::Dims count{Ny, Nx};

            auto var_iString = io.DefineVariable<std::string>("iString");
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);

            (void)var_iString;
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
        }

        {

            // Declare Single Value Attributes
            io.DefineAttribute<std::string>(s1_Single, attributeTestData.S1);
            io.DefineAttribute<std::string>(s1_Array, attributeTestData.S1array.data(),
                                            attributeTestData.S1array.size());

            io.DefineAttribute<int8_t>(i8_Single, attributeTestData.I8.front());
            io.DefineAttribute<int16_t>(i16_Single, attributeTestData.I16.front());
            io.DefineAttribute<int32_t>(i32_Single, attributeTestData.I32.front());
            io.DefineAttribute<int32_t>(i32_Array, attributeTestData.I32.data(),
                                        attributeTestData.I32.size());
            io.DefineAttribute<int64_t>(i64_Single, attributeTestData.I64.front());

            io.DefineAttribute<uint8_t>(u8_Single, attributeTestData.U8.front());
            io.DefineAttribute<uint16_t>(u16_Single, attributeTestData.U16.front());
            io.DefineAttribute<uint32_t>(u32_Single, attributeTestData.U32.front());
            io.DefineAttribute<uint32_t>(u32_Array, attributeTestData.U32.data(),
                                         attributeTestData.U32.size());
            io.DefineAttribute<uint64_t>(u64_Single, attributeTestData.U64.front());

            io.DefineAttribute<float>(r32_Single, attributeTestData.R32.front());
            io.DefineAttribute<float>(r32_Array, attributeTestData.R32.data(),
                                      attributeTestData.R32.size());
            io.DefineAttribute<double>(r64_Single, attributeTestData.R64.front());

            io.DefineAttribute<std::complex<float>>(cr32_Single, attributeTestData.CR32.front());
            io.DefineAttribute<std::complex<float>>(cr32_Array, attributeTestData.CR32.data(),
                                                    attributeTestData.CR32.size());
            io.DefineAttribute<std::complex<double>>(cr64_Single, attributeTestData.CR64.front());
        }

        io.SetEngine(engineName);
        io.SetParameter("AggregationRatio", "1");

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

            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({0, static_cast<size_t>(mpiRank * Nx)}, {Ny, Nx});
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

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            bpWriter.BeginStep();
            bpWriter.Put(var_iString, currentTestData.S1);
            bpWriter.Put(var_i8, currentTestData.I8.data());
            bpWriter.Put(var_i16, currentTestData.I16.data());
            bpWriter.Put(var_i32, currentTestData.I32.data());
            bpWriter.Put(var_i64, currentTestData.I64.data());
            bpWriter.Put(var_u8, currentTestData.U8.data());
            bpWriter.Put(var_u16, currentTestData.U16.data());
            bpWriter.Put(var_u32, currentTestData.U32.data());
            bpWriter.Put(var_u64, currentTestData.U64.data());
            bpWriter.Put(var_r32, currentTestData.R32.data());
            bpWriter.Put(var_r64, currentTestData.R64.data());

            bpWriter.EndStep();
        }

        // Close the file
        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("AppendIO");

        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
            const adios2::Dims count{Ny, Nx};

            auto var_iString = io.DefineVariable<std::string>("iString");
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);

            (void)var_iString;
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
        }

        {

            // Declare Single Value Attributes
            io.DefineAttribute<std::string>(s1_Single, attributeTestData.S1);
            io.DefineAttribute<std::string>(s1_Array, attributeTestData.S1array.data(),
                                            attributeTestData.S1array.size());

            io.DefineAttribute<int8_t>(i8_Single, attributeTestData.I8.front());
            io.DefineAttribute<int16_t>(i16_Single, attributeTestData.I16.front());
            io.DefineAttribute<int32_t>(i32_Single, attributeTestData.I32.front());
            io.DefineAttribute<int32_t>(i32_Array, attributeTestData.I32.data(),
                                        attributeTestData.I32.size());
            io.DefineAttribute<int64_t>(i64_Single, attributeTestData.I64.front());

            io.DefineAttribute<uint8_t>(u8_Single, attributeTestData.U8.front());
            io.DefineAttribute<uint16_t>(u16_Single, attributeTestData.U16.front());
            io.DefineAttribute<uint32_t>(u32_Single, attributeTestData.U32.front());
            io.DefineAttribute<uint32_t>(u32_Array, attributeTestData.U32.data(),
                                         attributeTestData.U32.size());
            io.DefineAttribute<uint64_t>(u64_Single, attributeTestData.U64.front());

            io.DefineAttribute<float>(r32_Single, attributeTestData.R32.front());
            io.DefineAttribute<float>(r32_Array, attributeTestData.R32.data(),
                                      attributeTestData.R32.size());
            io.DefineAttribute<double>(r64_Single, attributeTestData.R64.front());

            io.DefineAttribute<std::complex<float>>(cr32_Single, attributeTestData.CR32.front());
            io.DefineAttribute<std::complex<float>>(cr32_Array, attributeTestData.CR32.data(),
                                                    attributeTestData.CR32.size());
            io.DefineAttribute<std::complex<double>>(cr64_Single, attributeTestData.CR64.front());
        }

        io.SetEngine(engineName);
        io.AddTransport("file");
        io.SetParameter("AggregationRatio", "1");

        adios2::Engine bpAppender = io.Open(fname, adios2::Mode::Append);

        for (size_t step = NSteps; step < 2 * NSteps; ++step)
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

            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({0, static_cast<size_t>(mpiRank * Nx)}, {Ny, Nx});
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

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            bpAppender.BeginStep();
            bpAppender.Put(var_iString, currentTestData.S1);
            bpAppender.Put(var_i8, currentTestData.I8.data());
            bpAppender.Put(var_i16, currentTestData.I16.data());
            bpAppender.Put(var_i32, currentTestData.I32.data());
            bpAppender.Put(var_i64, currentTestData.I64.data());
            bpAppender.Put(var_u8, currentTestData.U8.data());
            bpAppender.Put(var_u16, currentTestData.U16.data());
            bpAppender.Put(var_u32, currentTestData.U32.data());
            bpAppender.Put(var_u64, currentTestData.U64.data());
            bpAppender.Put(var_r32, currentTestData.R32.data());
            bpAppender.Put(var_r64, currentTestData.R64.data());

            bpAppender.EndStep();
        }

        // Close the file
        bpAppender.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine(engineName);

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto attr_s1 = io.InquireAttribute<std::string>(s1_Single);
        auto attr_s1a = io.InquireAttribute<std::string>(s1_Array);
        auto attr_i8 = io.InquireAttribute<int8_t>(i8_Single);
        auto attr_i16 = io.InquireAttribute<int16_t>(i16_Single);
        auto attr_i32 = io.InquireAttribute<int32_t>(i32_Single);
        auto attr_i32a = io.InquireAttribute<int32_t>(i32_Array);
        auto attr_i64 = io.InquireAttribute<int64_t>(i64_Single);

        auto attr_u8 = io.InquireAttribute<uint8_t>(u8_Single);
        auto attr_u16 = io.InquireAttribute<uint16_t>(u16_Single);
        auto attr_u32 = io.InquireAttribute<uint32_t>(u32_Single);
        auto attr_u32a = io.InquireAttribute<uint32_t>(u32_Array);
        auto attr_u64 = io.InquireAttribute<uint64_t>(u64_Single);

        auto attr_r32 = io.InquireAttribute<float>(r32_Single);
        auto attr_r32a = io.InquireAttribute<float>(r32_Array);
        auto attr_r64 = io.InquireAttribute<double>(r64_Single);

        auto attr_cr32 = io.InquireAttribute<std::complex<float>>(cr32_Single);
        auto attr_cr32a = io.InquireAttribute<std::complex<float>>(cr32_Array);
        auto attr_cr64 = io.InquireAttribute<std::complex<double>>(cr64_Single);

        EXPECT_TRUE(attr_s1);
        ASSERT_EQ(attr_s1.Name(), s1_Single);
        ASSERT_EQ(attr_s1.Data().size() == 1, true);
        ASSERT_EQ(attr_s1.Type(), adios2::GetType<std::string>());
        ASSERT_EQ(attr_s1.Data().front(), attributeTestData.S1);

        EXPECT_TRUE(attr_s1a);
        ASSERT_EQ(attr_s1a.Name(), s1_Array);
        ASSERT_EQ(attr_s1a.Data().size() == 1, true);
        ASSERT_EQ(attr_s1a.Type(), adios2::GetType<std::string>());
        ASSERT_EQ(attr_s1a.Data()[0], attributeTestData.S1array[0]);

        EXPECT_TRUE(attr_i8);
        ASSERT_EQ(attr_i8.Name(), i8_Single);
        ASSERT_EQ(attr_i8.Data().size() == 1, true);
        ASSERT_EQ(attr_i8.Type(), adios2::GetType<int8_t>());
        ASSERT_EQ(attr_i8.Data().front(), attributeTestData.I8.front());

        EXPECT_TRUE(attr_i16);
        ASSERT_EQ(attr_i16.Name(), i16_Single);
        ASSERT_EQ(attr_i16.Data().size() == 1, true);
        ASSERT_EQ(attr_i16.Type(), adios2::GetType<int16_t>());
        ASSERT_EQ(attr_i16.Data().front(), attributeTestData.I16.front());

        EXPECT_TRUE(attr_i32);
        ASSERT_EQ(attr_i32.Name(), i32_Single);
        ASSERT_EQ(attr_i32.Data().size() == 1, true);
        ASSERT_EQ(attr_i32.Type(), adios2::GetType<int32_t>());
        ASSERT_EQ(attr_i32.Data().front(), attributeTestData.I32.front());

        EXPECT_TRUE(attr_i32a);
        ASSERT_EQ(attr_i32a.Name(), i32_Array);
        ASSERT_EQ(attr_i32a.Data().size() == attributeTestData.I32.size(), true);
        ASSERT_EQ(attr_i32a.Type(), adios2::GetType<int32_t>());
        ASSERT_EQ(attr_i32a.Data()[0], attributeTestData.I32[0]);

        EXPECT_TRUE(attr_i64);
        ASSERT_EQ(attr_i64.Name(), i64_Single);
        ASSERT_EQ(attr_i64.Data().size() == 1, true);
        ASSERT_EQ(attr_i64.Type(), adios2::GetType<int64_t>());
        ASSERT_EQ(attr_i64.Data().front(), attributeTestData.I64.front());

        EXPECT_TRUE(attr_u8);
        ASSERT_EQ(attr_u8.Name(), u8_Single);
        ASSERT_EQ(attr_u8.Data().size() == 1, true);
        ASSERT_EQ(attr_u8.Type(), adios2::GetType<uint8_t>());
        ASSERT_EQ(attr_u8.Data().front(), attributeTestData.U8.front());

        EXPECT_TRUE(attr_u16);
        ASSERT_EQ(attr_u16.Name(), u16_Single);
        ASSERT_EQ(attr_u16.Data().size() == 1, true);
        ASSERT_EQ(attr_u16.Type(), adios2::GetType<uint16_t>());
        ASSERT_EQ(attr_u16.Data().front(), attributeTestData.U16.front());

        EXPECT_TRUE(attr_u32);
        ASSERT_EQ(attr_u32.Name(), u32_Single);
        ASSERT_EQ(attr_u32.Data().size() == 1, true);
        ASSERT_EQ(attr_u32.Type(), adios2::GetType<uint32_t>());
        ASSERT_EQ(attr_u32.Data().front(), attributeTestData.U32.front());

        EXPECT_TRUE(attr_u32a);
        ASSERT_EQ(attr_u32a.Name(), u32_Array);
        ASSERT_EQ(attr_u32a.Data().size() == attributeTestData.U32.size(), true);
        ASSERT_EQ(attr_u32a.Type(), adios2::GetType<uint32_t>());
        ASSERT_EQ(attr_u32a.Data()[0], attributeTestData.U32[0]);

        EXPECT_TRUE(attr_u64);
        ASSERT_EQ(attr_u64.Name(), u64_Single);
        ASSERT_EQ(attr_u64.Data().size() == 1, true);
        ASSERT_EQ(attr_u64.Type(), adios2::GetType<uint64_t>());
        ASSERT_EQ(attr_u64.Data().front(), attributeTestData.U64.front());

        EXPECT_TRUE(attr_r32);
        ASSERT_EQ(attr_r32.Name(), r32_Single);
        ASSERT_EQ(attr_r32.Data().size() == 1, true);
        ASSERT_EQ(attr_r32.Type(), adios2::GetType<float>());
        ASSERT_EQ(attr_r32.Data().front(), attributeTestData.R32.front());

        EXPECT_TRUE(attr_r32a);
        ASSERT_EQ(attr_r32a.Name(), r32_Array);
        ASSERT_EQ(attr_r32a.Data().size() == attributeTestData.R32.size(), true);
        ASSERT_EQ(attr_r32a.Type(), adios2::GetType<float>());
        ASSERT_EQ(attr_r32a.Data()[0], attributeTestData.R32[0]);

        EXPECT_TRUE(attr_r64);
        ASSERT_EQ(attr_r64.Name(), r64_Single);
        ASSERT_EQ(attr_r64.Data().size() == 1, true);
        ASSERT_EQ(attr_r64.Type(), adios2::GetType<double>());
        ASSERT_EQ(attr_r64.Data().front(), attributeTestData.R64.front());

        EXPECT_TRUE(attr_cr32);
        ASSERT_EQ(attr_cr32.Name(), cr32_Single);
        ASSERT_EQ(attr_cr32.Data().size() == 1, true);
        ASSERT_EQ(attr_cr32.Type(), adios2::GetType<std::complex<float>>());
        ASSERT_EQ(attr_cr32.Data().front(), attributeTestData.CR32.front());

        EXPECT_TRUE(attr_cr32a);
        ASSERT_EQ(attr_cr32a.Name(), cr32_Array);
        ASSERT_EQ(attr_cr32a.Data().size() == attributeTestData.CR32.size(), true);
        ASSERT_EQ(attr_cr32a.Type(), adios2::GetType<std::complex<float>>());
        ASSERT_EQ(attr_cr32a.Data()[0], attributeTestData.CR32[0]);

        EXPECT_TRUE(attr_cr64);
        ASSERT_EQ(attr_cr64.Name(), cr64_Single);
        ASSERT_EQ(attr_cr64.Data().size() == 1, true);
        ASSERT_EQ(attr_cr64.Type(), adios2::GetType<std::complex<double>>());
        ASSERT_EQ(attr_cr64.Data().front(), attributeTestData.CR64.front());

        EXPECT_EQ(bpReader.Steps(), 2 * NSteps);
        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        ASSERT_EQ(var_iString.Shape().size(), 0);
        ASSERT_EQ(var_iString.Steps(), 2 * NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), 2 * NSteps);
        ASSERT_EQ(var_i8.Shape()[0], Ny);
        ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), 2 * NSteps);
        ASSERT_EQ(var_i16.Shape()[0], Ny);
        ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), 2 * NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), 2 * NSteps);
        ASSERT_EQ(var_i64.Shape()[0], Ny);
        ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), 2 * NSteps);
        ASSERT_EQ(var_u8.Shape()[0], Ny);
        ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), 2 * NSteps);
        ASSERT_EQ(var_u16.Shape()[0], Ny);
        ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), 2 * NSteps);
        ASSERT_EQ(var_u32.Shape()[0], Ny);
        ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), 2 * NSteps);
        ASSERT_EQ(var_u64.Shape()[0], Ny);
        ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), 2 * NSteps);
        ASSERT_EQ(var_r32.Shape()[0], Ny);
        ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), 2 * NSteps);
        ASSERT_EQ(var_r64.Shape()[0], Ny);
        ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

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

        for (size_t t = 0; t < 2 * NSteps; ++t)
        {
            var_i8.SetStepSelection({t, 1});
            var_i16.SetStepSelection({t, 1});
            var_i32.SetStepSelection({t, 1});
            var_i64.SetStepSelection({t, 1});

            var_u8.SetStepSelection({t, 1});
            var_u16.SetStepSelection({t, 1});
            var_u32.SetStepSelection({t, 1});
            var_u64.SetStepSelection({t, 1});

            var_r32.SetStepSelection({t, 1});
            var_r64.SetStepSelection({t, 1});

            bpReader.Get(var_iString, IString);

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
            }
        }
        bpReader.Close();
    }
}

// Write with append combined with aggregation, same aggregation ratio
TEST_F(BPWriteAppendReadTestADIOS2, ADIOS2BPWriteAppendReadAggregate)
{
    int mpiRank = 0, mpiSize = 1;
    const std::size_t Nx = 4;
    const std::size_t Ny = 2;
    const std::size_t NSteps = 2;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteAppendReadAggregate_MPI.bp");
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    const std::string fname("ADIOS2BPWriteAppendReadAggregate.bp");
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("WriteIO");
        io.SetEngine(engineName);
        io.AddTransport("file");
        io.SetParameter("NumAggregators", "2");

        const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);

        /* Write phase I */
        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);
        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);
            bpWriter.BeginStep();
            bpWriter.Put(var_i32, currentTestData.I32.data());
            bpWriter.EndStep();
        }
        bpWriter.Close();

        /* Write phase II: append */
        bpWriter = io.Open(fname, adios2::Mode::Append);
        for (size_t step = NSteps; step < 2 * NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);
            bpWriter.BeginStep();
            bpWriter.Put(var_i32, currentTestData.I32.data());
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine(engineName);

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), 2 * NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::array<int32_t, Nx * Ny> I32;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i32.SetSelection(sel);
        for (size_t t = 0; t < 2 * NSteps; ++t)
        {
            var_i32.SetStepSelection({t, 1});
            bpReader.Get(var_i32, I32.data());
            bpReader.PerformGets();

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
            }
        }
        bpReader.Close();
    }
}

// Write with append combined with aggregation, same aggregation ratio
TEST_F(BPWriteAppendReadTestADIOS2, ADIOS2BPWriteAppendReadVaryingAggregation)
{
    int mpiRank = 0, mpiSize = 1;
    const std::size_t Nx = 4;
    const std::size_t Ny = 2;
    const std::size_t NSteps = 2;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2BPWriteAppendReadVaryingAggregate_MPI.bp");
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    const std::string fname("ADIOS2BPWriteAppendReadVaryingAggregate.bp");
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("WriteIO");
        io.SetEngine(engineName);
        io.AddTransport("file");
        io.SetParameter("NumAggregators", "3");

        const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);

        /* Write phase I */
        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);
        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);
            bpWriter.BeginStep();
            bpWriter.Put(var_i32, currentTestData.I32.data());
            bpWriter.EndStep();
        }
        bpWriter.Close();

        /* Write phase II: append */
        io.SetParameter("NumAggregators", "1");
        bpWriter = io.Open(fname, adios2::Mode::Append);
        for (size_t step = NSteps; step < 2 * NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);
            bpWriter.BeginStep();
            bpWriter.Put(var_i32, currentTestData.I32.data());
            bpWriter.EndStep();
        }
        bpWriter.Close();

        /* Write phase III: append */
        io.SetParameter("NumAggregators", "2");
        bpWriter = io.Open(fname, adios2::Mode::Append);
        for (size_t step = 2 * NSteps; step < 3 * NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);
            bpWriter.BeginStep();
            bpWriter.Put(var_i32, currentTestData.I32.data());
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }

    {
        size_t NumSteps = 3 * NSteps;
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine(engineName);

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NumSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::array<int32_t, Nx * Ny> I32;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i32.SetSelection(sel);

        for (size_t t = 0; t < NumSteps; ++t)
        {
            var_i32.SetStepSelection({t, 1});
            bpReader.Get(var_i32, I32.data());
            bpReader.PerformGets();

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
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
    else
    {
        engineName = "BP4";
    }
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
