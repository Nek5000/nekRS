/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <complex>
#include <cstdint>
#include <string>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPWriteReadAttributes : public ::testing::Test
{
public:
    BPWriteReadAttributes() = default;

    SmallTestData m_TestData;
};

// ADIOS2 write, read for single value attributes
TEST_F(BPWriteReadAttributes, WriteReadSingleTypes)
{

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
    const std::string r128_Single = std::string("r128_Single_") + zero;
    const std::string cr32_Single = std::string("cr32_Single_") + zero;
    const std::string cr64_Single = std::string("cr64_Single_") + zero;

    // When collective meta generation has landed, use
    // generateNewSmallTestData(m_TestData, 0, mpiRank, mpiSize);
    // Generate current testing data
    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, 0, 0);

// Write test data using BP
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "WriteAttributeReadSingleTypes_MPI.bp";
#else
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "WriteAttributeReadSingleTypes.bp";
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare Single Value Attributes
        io.DefineAttribute<std::string>(s1_Single, currentTestData.S1);
        io.DefineAttribute<std::string>(s1_Array, currentTestData.S1array.data(),
                                        currentTestData.S1array.size());

        io.DefineAttribute<int8_t>(i8_Single, currentTestData.I8.front());
        io.DefineAttribute<int16_t>(i16_Single, currentTestData.I16.front());
        io.DefineAttribute<int32_t>(i32_Single, currentTestData.I32.front());
        io.DefineAttribute<int64_t>(i64_Single, currentTestData.I64.front());

        io.DefineAttribute<uint8_t>(u8_Single, currentTestData.U8.front());
        io.DefineAttribute<uint16_t>(u16_Single, currentTestData.U16.front());
        io.DefineAttribute<uint32_t>(u32_Single, currentTestData.U32.front());
        io.DefineAttribute<uint64_t>(u64_Single, currentTestData.U64.front());

        io.DefineAttribute<float>(r32_Single, currentTestData.R32.front());
        io.DefineAttribute<double>(r64_Single, currentTestData.R64.front());
        io.DefineAttribute<long double>(r128_Single, currentTestData.R128.front());

        io.DefineAttribute<std::complex<float>>(cr32_Single, currentTestData.CR32.front());
        io.DefineAttribute<std::complex<double>>(cr64_Single, currentTestData.CR64.front());

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.AddTransport("File");

        adios2::Engine engine = io.Open(fName, adios2::Mode::Write);
        // only attributes are written
        engine.Close();
    }

    {
        adios2::IO ioRead = adios.DeclareIO("ioRead");
        // ioRead.AddTransport("File");
        // ioRead.SetParameter("OpenAsFile", "true");
        if (!engineName.empty())
        {
            ioRead.SetEngine(engineName);
        }

        adios2::Engine bpRead = ioRead.Open(fName, adios2::Mode::ReadRandomAccess);

        auto attr_s1 = ioRead.InquireAttribute<std::string>(s1_Single);
        auto attr_s1a = ioRead.InquireAttribute<std::string>(s1_Array);
        auto attr_i8 = ioRead.InquireAttribute<int8_t>(i8_Single);
        auto attr_i16 = ioRead.InquireAttribute<int16_t>(i16_Single);
        auto attr_i32 = ioRead.InquireAttribute<int32_t>(i32_Single);
        auto attr_i64 = ioRead.InquireAttribute<int64_t>(i64_Single);

        auto attr_u8 = ioRead.InquireAttribute<uint8_t>(u8_Single);
        auto attr_u16 = ioRead.InquireAttribute<uint16_t>(u16_Single);
        auto attr_u32 = ioRead.InquireAttribute<uint32_t>(u32_Single);
        auto attr_u64 = ioRead.InquireAttribute<uint64_t>(u64_Single);

        auto attr_r32 = ioRead.InquireAttribute<float>(r32_Single);
        auto attr_r64 = ioRead.InquireAttribute<double>(r64_Single);
        auto attr_r128 = ioRead.InquireAttribute<long double>(r128_Single);

        auto attr_cr32 = ioRead.InquireAttribute<std::complex<float>>(cr32_Single);
        auto attr_cr64 = ioRead.InquireAttribute<std::complex<double>>(cr64_Single);

        EXPECT_TRUE(attr_s1);
        ASSERT_EQ(attr_s1.Name(), s1_Single);
        ASSERT_EQ(attr_s1.Data().size() == 1, true);
        ASSERT_EQ(attr_s1.Type(), adios2::GetType<std::string>());
        ASSERT_EQ(attr_s1.Data().front(), currentTestData.S1);

        EXPECT_TRUE(attr_s1a);
        ASSERT_EQ(attr_s1a.Name(), s1_Array);
        ASSERT_EQ(attr_s1a.Data().size() == 1, true);
        ASSERT_EQ(attr_s1a.Type(), adios2::GetType<std::string>());
        ASSERT_EQ(attr_s1a.Data()[0], currentTestData.S1array[0]);

        EXPECT_TRUE(attr_i8);
        ASSERT_EQ(attr_i8.Name(), i8_Single);
        ASSERT_EQ(attr_i8.Data().size() == 1, true);
        ASSERT_EQ(attr_i8.Type(), adios2::GetType<int8_t>());
        ASSERT_EQ(attr_i8.Data().front(), currentTestData.I8.front());

        EXPECT_TRUE(attr_i16);
        ASSERT_EQ(attr_i16.Name(), i16_Single);
        ASSERT_EQ(attr_i16.Data().size() == 1, true);
        ASSERT_EQ(attr_i16.Type(), adios2::GetType<int16_t>());
        ASSERT_EQ(attr_i16.Data().front(), currentTestData.I16.front());

        EXPECT_TRUE(attr_i32);
        ASSERT_EQ(attr_i32.Name(), i32_Single);
        ASSERT_EQ(attr_i32.Data().size() == 1, true);
        ASSERT_EQ(attr_i32.Type(), adios2::GetType<int32_t>());
        ASSERT_EQ(attr_i32.Data().front(), currentTestData.I32.front());

        EXPECT_TRUE(attr_i64);
        ASSERT_EQ(attr_i64.Name(), i64_Single);
        ASSERT_EQ(attr_i64.Data().size() == 1, true);
        ASSERT_EQ(attr_i64.Type(), adios2::GetType<int64_t>());
        ASSERT_EQ(attr_i64.Data().front(), currentTestData.I64.front());

        EXPECT_TRUE(attr_u8);
        ASSERT_EQ(attr_u8.Name(), u8_Single);
        ASSERT_EQ(attr_u8.Data().size() == 1, true);
        ASSERT_EQ(attr_u8.Type(), adios2::GetType<uint8_t>());
        ASSERT_EQ(attr_u8.Data().front(), currentTestData.U8.front());

        EXPECT_TRUE(attr_u16);
        ASSERT_EQ(attr_u16.Name(), u16_Single);
        ASSERT_EQ(attr_u16.Data().size() == 1, true);
        ASSERT_EQ(attr_u16.Type(), adios2::GetType<uint16_t>());
        ASSERT_EQ(attr_u16.Data().front(), currentTestData.U16.front());

        EXPECT_TRUE(attr_u32);
        ASSERT_EQ(attr_u32.Name(), u32_Single);
        ASSERT_EQ(attr_u32.Data().size() == 1, true);
        ASSERT_EQ(attr_u32.Type(), adios2::GetType<uint32_t>());
        ASSERT_EQ(attr_u32.Data().front(), currentTestData.U32.front());

        EXPECT_TRUE(attr_u64);
        ASSERT_EQ(attr_u64.Name(), u64_Single);
        ASSERT_EQ(attr_u64.Data().size() == 1, true);
        ASSERT_EQ(attr_u64.Type(), adios2::GetType<uint64_t>());
        ASSERT_EQ(attr_u64.Data().front(), currentTestData.U64.front());

        EXPECT_TRUE(attr_r32);
        ASSERT_EQ(attr_r32.Name(), r32_Single);
        ASSERT_EQ(attr_r32.Data().size() == 1, true);
        ASSERT_EQ(attr_r32.Type(), adios2::GetType<float>());
        ASSERT_EQ(attr_r32.Data().front(), currentTestData.R32.front());

        EXPECT_TRUE(attr_r64);
        ASSERT_EQ(attr_r64.Name(), r64_Single);
        ASSERT_EQ(attr_r64.Data().size() == 1, true);
        ASSERT_EQ(attr_r64.Type(), adios2::GetType<double>());
        ASSERT_EQ(attr_r64.Data().front(), currentTestData.R64.front());

        EXPECT_TRUE(attr_r128);
        ASSERT_EQ(attr_r128.Name(), r128_Single);
        ASSERT_EQ(attr_r128.Data().size() == 1, true);
        ASSERT_EQ(attr_r128.Type(), adios2::GetType<long double>());
        ASSERT_EQ(attr_r128.Data().front(), currentTestData.R128.front());

        EXPECT_TRUE(attr_cr32);
        ASSERT_EQ(attr_cr32.Name(), cr32_Single);
        ASSERT_EQ(attr_cr32.Data().size() == 1, true);
        ASSERT_EQ(attr_cr32.Type(), adios2::GetType<std::complex<float>>());
        ASSERT_EQ(attr_cr32.Data().front(), currentTestData.CR32.front());

        EXPECT_TRUE(attr_cr64);
        ASSERT_EQ(attr_cr64.Name(), cr64_Single);
        ASSERT_EQ(attr_cr64.Data().size() == 1, true);
        ASSERT_EQ(attr_cr64.Type(), adios2::GetType<std::complex<double>>());
        ASSERT_EQ(attr_cr64.Data().front(), currentTestData.CR64.front());

        bpRead.Close();
    }
}

// ADIOS2 write read for array attributes
TEST_F(BPWriteReadAttributes, WriteReadArrayTypes)
{

#if ADIOS2_USE_MPI
    int mpiRank = 0, mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "WriteAttributeReadArrayTypes_MPI.bp";
#else
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "WriteAttributeReadArrayTypes.bp";
#endif

    const std::string zero = std::to_string(0);
    const std::string s1_Array = std::string("s1_Array_") + zero;
    const std::string i8_Array = std::string("i8_Array_") + zero;
    const std::string i16_Array = std::string("i16_Array_") + zero;
    const std::string i32_Array = std::string("i32_Array_") + zero;
    const std::string i64_Array = std::string("i64_Array_") + zero;
    const std::string u8_Array = std::string("u8_Array_") + zero;
    const std::string u16_Array = std::string("u16_Array_") + zero;
    const std::string u32_Array = std::string("u32_Array_") + zero;
    const std::string u64_Array = std::string("u64_Array_") + zero;
    const std::string r32_Array = std::string("r32_Array_") + zero;
    const std::string r64_Array = std::string("r64_Array_") + zero;
    const std::string r128_Array = std::string("r128_Array_") + zero;
    const std::string cr32_Array = std::string("cr32_Array_") + zero;
    const std::string cr64_Array = std::string("cr64_Array_") + zero;

    // When collective meta generation has landed, use
    // generateNewSmallTestData(m_TestData, 0, mpiRank, mpiSize);
    // Generate current testing data
    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, 0, 0);

// Write test data using BP
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare Single Value Attributes
        io.DefineAttribute<std::string>(s1_Array, currentTestData.S3.data(),
                                        currentTestData.S3.size());

        io.DefineAttribute<int8_t>(i8_Array, currentTestData.I8.data(), currentTestData.I8.size());
        io.DefineAttribute<int16_t>(i16_Array, currentTestData.I16.data(),
                                    currentTestData.I16.size());
        io.DefineAttribute<int32_t>(i32_Array, currentTestData.I32.data(),
                                    currentTestData.I32.size());
        io.DefineAttribute<int64_t>(i64_Array, currentTestData.I64.data(),
                                    currentTestData.I64.size());

        io.DefineAttribute<uint8_t>(u8_Array, currentTestData.U8.data(), currentTestData.U8.size());
        io.DefineAttribute<uint16_t>(u16_Array, currentTestData.U16.data(),
                                     currentTestData.U16.size());
        io.DefineAttribute<uint32_t>(u32_Array, currentTestData.U32.data(),
                                     currentTestData.U32.size());
        io.DefineAttribute<uint64_t>(u64_Array, currentTestData.U64.data(),
                                     currentTestData.U64.size());

        io.DefineAttribute<float>(r32_Array, currentTestData.R32.data(),
                                  currentTestData.R32.size());
        io.DefineAttribute<double>(r64_Array, currentTestData.R64.data(),
                                   currentTestData.R64.size());
        io.DefineAttribute<long double>(r128_Array, currentTestData.R128.data(),
                                        currentTestData.R128.size());

        io.DefineAttribute<std::complex<float>>(cr32_Array, currentTestData.CR32.data(),
                                                currentTestData.CR32.size());
        io.DefineAttribute<std::complex<double>>(cr64_Array, currentTestData.CR64.data(),
                                                 currentTestData.CR64.size());

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("File");
        }
        io.AddTransport("file");

        adios2::Engine engine = io.Open(fName, adios2::Mode::Write);
        // only attributes are written
        engine.Close();
    }

    {
        adios2::IO ioRead = adios.DeclareIO("ioRead");

        if (!engineName.empty())
        {
            ioRead.SetEngine(engineName);
        }

        adios2::Engine bpRead = ioRead.Open(fName, adios2::Mode::ReadRandomAccess);

        auto attr_s1 = ioRead.InquireAttribute<std::string>(s1_Array);

        auto attr_i8 = ioRead.InquireAttribute<int8_t>(i8_Array);
        auto attr_i16 = ioRead.InquireAttribute<int16_t>(i16_Array);
        auto attr_i32 = ioRead.InquireAttribute<int32_t>(i32_Array);
        auto attr_i64 = ioRead.InquireAttribute<int64_t>(i64_Array);

        auto attr_u8 = ioRead.InquireAttribute<uint8_t>(u8_Array);
        auto attr_u16 = ioRead.InquireAttribute<uint16_t>(u16_Array);
        auto attr_u32 = ioRead.InquireAttribute<uint32_t>(u32_Array);
        auto attr_u64 = ioRead.InquireAttribute<uint64_t>(u64_Array);

        auto attr_r32 = ioRead.InquireAttribute<float>(r32_Array);
        auto attr_r64 = ioRead.InquireAttribute<double>(r64_Array);
        auto attr_r128 = ioRead.InquireAttribute<long double>(r128_Array);

        auto attr_cr32 = ioRead.InquireAttribute<std::complex<float>>(cr32_Array);
        auto attr_cr64 = ioRead.InquireAttribute<std::complex<double>>(cr64_Array);

        EXPECT_TRUE(attr_s1);
        ASSERT_EQ(attr_s1.Name(), s1_Array);
        ASSERT_EQ(attr_s1.Data().size() == 1, false);
        ASSERT_EQ(attr_s1.Type(), adios2::GetType<std::string>());

        EXPECT_TRUE(attr_i8);
        ASSERT_EQ(attr_i8.Name(), i8_Array);
        ASSERT_EQ(attr_i8.Data().size() == 1, false);
        ASSERT_EQ(attr_i8.Type(), adios2::GetType<int8_t>());

        EXPECT_TRUE(attr_i16);
        ASSERT_EQ(attr_i16.Name(), i16_Array);
        ASSERT_EQ(attr_i16.Data().size() == 1, false);
        ASSERT_EQ(attr_i16.Type(), adios2::GetType<int16_t>());

        EXPECT_TRUE(attr_i32);
        ASSERT_EQ(attr_i32.Name(), i32_Array);
        ASSERT_EQ(attr_i32.Data().size() == 1, false);
        ASSERT_EQ(attr_i32.Type(), adios2::GetType<int32_t>());

        EXPECT_TRUE(attr_i64);
        ASSERT_EQ(attr_i64.Name(), i64_Array);
        ASSERT_EQ(attr_i64.Data().size() == 1, false);
        ASSERT_EQ(attr_i64.Type(), adios2::GetType<int64_t>());

        EXPECT_TRUE(attr_u8);
        ASSERT_EQ(attr_u8.Name(), u8_Array);
        ASSERT_EQ(attr_u8.Data().size() == 1, false);
        ASSERT_EQ(attr_u8.Type(), adios2::GetType<uint8_t>());

        EXPECT_TRUE(attr_u16);
        ASSERT_EQ(attr_u16.Name(), u16_Array);
        ASSERT_EQ(attr_u16.Data().size() == 1, false);
        ASSERT_EQ(attr_u16.Type(), adios2::GetType<uint16_t>());

        EXPECT_TRUE(attr_u32);
        ASSERT_EQ(attr_u32.Name(), u32_Array);
        ASSERT_EQ(attr_u32.Data().size() == 1, false);
        ASSERT_EQ(attr_u32.Type(), adios2::GetType<uint32_t>());

        EXPECT_TRUE(attr_u64);
        ASSERT_EQ(attr_u64.Name(), u64_Array);
        ASSERT_EQ(attr_u64.Data().size() == 1, false);
        ASSERT_EQ(attr_u64.Type(), adios2::GetType<uint64_t>());

        EXPECT_TRUE(attr_r32);
        ASSERT_EQ(attr_r32.Name(), r32_Array);
        ASSERT_EQ(attr_r32.Data().size() == 1, false);
        ASSERT_EQ(attr_r32.Type(), adios2::GetType<float>());

        EXPECT_TRUE(attr_r64);
        ASSERT_EQ(attr_r64.Name(), r64_Array);
        ASSERT_EQ(attr_r64.Data().size() == 1, false);
        ASSERT_EQ(attr_r64.Type(), adios2::GetType<double>());

        EXPECT_TRUE(attr_r128);
        ASSERT_EQ(attr_r128.Name(), r128_Array);
        ASSERT_EQ(attr_r128.Data().size() == 1, false);
        ASSERT_EQ(attr_r128.Type(), adios2::GetType<long double>());

        EXPECT_TRUE(attr_cr32);
        ASSERT_EQ(attr_cr32.Name(), cr32_Array);
        ASSERT_EQ(attr_cr32.Data().size() == 1, false);
        ASSERT_EQ(attr_cr32.Type(), adios2::GetType<std::complex<float>>());

        EXPECT_TRUE(attr_cr64);
        ASSERT_EQ(attr_cr64.Name(), cr64_Array);
        ASSERT_EQ(attr_cr64.Data().size() == 1, false);
        ASSERT_EQ(attr_cr64.Type(), adios2::GetType<std::complex<double>>());

        auto I8 = attr_i8.Data();
        auto I16 = attr_i16.Data();
        auto I32 = attr_i32.Data();
        auto I64 = attr_i64.Data();

        auto U8 = attr_u8.Data();
        auto U16 = attr_u16.Data();
        auto U32 = attr_u32.Data();
        auto U64 = attr_u64.Data();

        const size_t Nx = 10;
        for (size_t i = 0; i < Nx; ++i)
        {
            EXPECT_EQ(I8[i], currentTestData.I8[i]);
            EXPECT_EQ(I16[i], currentTestData.I16[i]);
            EXPECT_EQ(I32[i], currentTestData.I32[i]);
            EXPECT_EQ(I64[i], currentTestData.I64[i]);

            EXPECT_EQ(U8[i], currentTestData.U8[i]);
            EXPECT_EQ(U16[i], currentTestData.U16[i]);
            EXPECT_EQ(U32[i], currentTestData.U32[i]);
            EXPECT_EQ(U64[i], currentTestData.U64[i]);
        }

        bpRead.Close();
    }
}

TEST_F(BPWriteReadAttributes, BPWriteReadSingleTypesVar)
{

    const std::string zero = std::to_string(0);
    const std::string s1_Single = std::string("s1_Single_") + zero;
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
    const std::string r128_Single = std::string("r128_Single_") + zero;
    const std::string cr32_Single = std::string("cr32_Single_") + zero;
    const std::string cr64_Single = std::string("cr64_Single_") + zero;

    // When collective meta generation has landed, use
    // generateNewSmallTestData(m_TestData, 0, mpiRank, mpiSize);
    // Generate current testing data
    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, 0, 0);

    const std::string separator = "/";

// Write test data using BP
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    const std::string fName = "foo" + std::string(&adios2::PathSeparator, 1) +
                              "BPWriteAttributeReadSingleTypesVar_MPI.bp";
#else
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "BPWriteAttributeReadSingleTypesVar.bp";
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        // Declare Single Value Attributes
        auto var = io.DefineVariable<int>("myVar");

        io.DefineAttribute<std::string>(s1_Single, currentTestData.S1, var.Name());
        io.DefineAttribute<int8_t>(i8_Single, currentTestData.I8.front(), var.Name());
        io.DefineAttribute<int16_t>(i16_Single, currentTestData.I16.front(), var.Name());
        io.DefineAttribute<int32_t>(i32_Single, currentTestData.I32.front(), var.Name());
        io.DefineAttribute<int64_t>(i64_Single, currentTestData.I64.front(), var.Name());

        io.DefineAttribute<uint8_t>(u8_Single, currentTestData.U8.front(), var.Name());
        io.DefineAttribute<uint16_t>(u16_Single, currentTestData.U16.front(), var.Name());
        io.DefineAttribute<uint32_t>(u32_Single, currentTestData.U32.front(), var.Name());
        io.DefineAttribute<uint64_t>(u64_Single, currentTestData.U64.front(), var.Name());

        io.DefineAttribute<float>(r32_Single, currentTestData.R32.front(), var.Name());
        io.DefineAttribute<double>(r64_Single, currentTestData.R64.front(), var.Name());
        io.DefineAttribute<long double>(r128_Single, currentTestData.R128.front(), var.Name());

        io.DefineAttribute<std::complex<float>>(cr32_Single, currentTestData.CR32.front(),
                                                var.Name());
        io.DefineAttribute<std::complex<double>>(cr64_Single, currentTestData.CR64.front(),
                                                 var.Name());

        adios2::Engine engine = io.Open(fName, adios2::Mode::Write);
        engine.Put(var, 10);
        engine.Close();
    }

    {
        adios2::IO ioRead = adios.DeclareIO("ioRead");
        if (!engineName.empty())
        {
            ioRead.SetEngine(engineName);
        }

        adios2::Engine bpRead = ioRead.Open(fName, adios2::Mode::ReadRandomAccess);

        auto var = ioRead.InquireVariable<int>("myVar");

        auto attr_s1 = ioRead.InquireAttribute<std::string>(s1_Single, var.Name());
        auto attr_i8 = ioRead.InquireAttribute<int8_t>(i8_Single, var.Name());
        auto attr_i16 = ioRead.InquireAttribute<int16_t>(i16_Single, var.Name());
        auto attr_i32 = ioRead.InquireAttribute<int32_t>(i32_Single, var.Name());
        auto attr_i64 = ioRead.InquireAttribute<int64_t>(i64_Single, var.Name());

        auto attr_u8 = ioRead.InquireAttribute<uint8_t>(u8_Single, var.Name());
        auto attr_u16 = ioRead.InquireAttribute<uint16_t>(u16_Single, var.Name());
        auto attr_u32 = ioRead.InquireAttribute<uint32_t>(u32_Single, var.Name());
        auto attr_u64 = ioRead.InquireAttribute<uint64_t>(u64_Single, var.Name());

        auto attr_r32 = ioRead.InquireAttribute<float>(r32_Single, var.Name());
        auto attr_r64 = ioRead.InquireAttribute<double>(r64_Single, var.Name());
        auto attr_r128 = ioRead.InquireAttribute<long double>(r128_Single, var.Name());

        auto attr_cr32 = ioRead.InquireAttribute<std::complex<float>>(cr32_Single, var.Name());
        auto attr_cr64 = ioRead.InquireAttribute<std::complex<double>>(cr64_Single, var.Name());

        EXPECT_TRUE(attr_s1);
        ASSERT_EQ(attr_s1.Name(), var.Name() + separator + s1_Single);
        ASSERT_EQ(attr_s1.Data().size() == 1, true);
        ASSERT_EQ(attr_s1.Type(), adios2::GetType<std::string>());
        ASSERT_EQ(attr_s1.Data().front(), currentTestData.S1);

        EXPECT_TRUE(attr_i8);
        ASSERT_EQ(attr_i8.Name(), var.Name() + separator + i8_Single);
        ASSERT_EQ(attr_i8.Data().size() == 1, true);
        ASSERT_EQ(attr_i8.Type(), adios2::GetType<int8_t>());
        ASSERT_EQ(attr_i8.Data().front(), currentTestData.I8.front());

        EXPECT_TRUE(attr_i16);
        ASSERT_EQ(attr_i16.Name(), var.Name() + separator + i16_Single);
        ASSERT_EQ(attr_i16.Data().size() == 1, true);
        ASSERT_EQ(attr_i16.Type(), adios2::GetType<int16_t>());
        ASSERT_EQ(attr_i16.Data().front(), currentTestData.I16.front());

        EXPECT_TRUE(attr_i32);
        ASSERT_EQ(attr_i32.Name(), var.Name() + separator + i32_Single);
        ASSERT_EQ(attr_i32.Data().size() == 1, true);
        ASSERT_EQ(attr_i32.Type(), adios2::GetType<int32_t>());
        ASSERT_EQ(attr_i32.Data().front(), currentTestData.I32.front());

        EXPECT_TRUE(attr_i64);
        ASSERT_EQ(attr_i64.Name(), var.Name() + separator + i64_Single);
        ASSERT_EQ(attr_i64.Data().size() == 1, true);
        ASSERT_EQ(attr_i64.Type(), adios2::GetType<int64_t>());
        ASSERT_EQ(attr_i64.Data().front(), currentTestData.I64.front());

        EXPECT_TRUE(attr_u8);
        ASSERT_EQ(attr_u8.Name(), var.Name() + separator + u8_Single);
        ASSERT_EQ(attr_u8.Data().size() == 1, true);
        ASSERT_EQ(attr_u8.Type(), adios2::GetType<uint8_t>());
        ASSERT_EQ(attr_u8.Data().front(), currentTestData.U8.front());

        EXPECT_TRUE(attr_u16);
        ASSERT_EQ(attr_u16.Name(), var.Name() + separator + u16_Single);
        ASSERT_EQ(attr_u16.Data().size() == 1, true);
        ASSERT_EQ(attr_u16.Type(), adios2::GetType<uint16_t>());
        ASSERT_EQ(attr_u16.Data().front(), currentTestData.U16.front());

        EXPECT_TRUE(attr_u32);
        ASSERT_EQ(attr_u32.Name(), var.Name() + separator + u32_Single);
        ASSERT_EQ(attr_u32.Data().size() == 1, true);
        ASSERT_EQ(attr_u32.Type(), adios2::GetType<uint32_t>());
        ASSERT_EQ(attr_u32.Data().front(), currentTestData.U32.front());

        EXPECT_TRUE(attr_u64);
        ASSERT_EQ(attr_u64.Name(), var.Name() + separator + u64_Single);
        ASSERT_EQ(attr_u64.Data().size() == 1, true);
        ASSERT_EQ(attr_u64.Type(), adios2::GetType<uint64_t>());
        ASSERT_EQ(attr_u64.Data().front(), currentTestData.U64.front());

        EXPECT_TRUE(attr_r32);
        ASSERT_EQ(attr_r32.Name(), var.Name() + separator + r32_Single);
        ASSERT_EQ(attr_r32.Data().size() == 1, true);
        ASSERT_EQ(attr_r32.Type(), adios2::GetType<float>());
        ASSERT_EQ(attr_r32.Data().front(), currentTestData.R32.front());

        EXPECT_TRUE(attr_r64);
        ASSERT_EQ(attr_r64.Name(), var.Name() + separator + r64_Single);
        ASSERT_EQ(attr_r64.Data().size() == 1, true);
        ASSERT_EQ(attr_r64.Type(), adios2::GetType<double>());
        ASSERT_EQ(attr_r64.Data().front(), currentTestData.R64.front());

        EXPECT_TRUE(attr_r128);
        ASSERT_EQ(attr_r128.Name(), var.Name() + separator + r128_Single);
        ASSERT_EQ(attr_r128.Data().size() == 1, true);
        ASSERT_EQ(attr_r128.Type(), adios2::GetType<long double>());
        ASSERT_EQ(attr_r128.Data().front(), currentTestData.R128.front());

        EXPECT_TRUE(attr_cr32);
        ASSERT_EQ(attr_cr32.Name(), var.Name() + separator + cr32_Single);
        ASSERT_EQ(attr_cr32.Data().size() == 1, true);
        ASSERT_EQ(attr_cr32.Type(), adios2::GetType<std::complex<float>>());
        ASSERT_EQ(attr_cr32.Data().front(), currentTestData.CR32.front());

        EXPECT_TRUE(attr_cr64);
        ASSERT_EQ(attr_cr64.Name(), var.Name() + separator + cr64_Single);
        ASSERT_EQ(attr_cr64.Data().size() == 1, true);
        ASSERT_EQ(attr_cr64.Type(), adios2::GetType<std::complex<double>>());
        ASSERT_EQ(attr_cr64.Data().front(), currentTestData.CR64.front());

        bpRead.Close();
    }
}

// ADIOS2 write read for array attributes
TEST_F(BPWriteReadAttributes, WriteReadArrayTypesVar)
{

#if ADIOS2_USE_MPI
    int mpiRank = 0, mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "BPWriteAttributeReadArrayTypesVar_MPI.bp";
#else
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "BPWriteAttributeReadArrayTypesVar.bp";
#endif

    const std::string zero = std::to_string(0);
    const std::string s1_Array = std::string("s1_Array_") + zero;
    const std::string i8_Array = std::string("i8_Array_") + zero;
    const std::string i16_Array = std::string("i16_Array_") + zero;
    const std::string i32_Array = std::string("i32_Array_") + zero;
    const std::string i64_Array = std::string("i64_Array_") + zero;
    const std::string u8_Array = std::string("u8_Array_") + zero;
    const std::string u16_Array = std::string("u16_Array_") + zero;
    const std::string u32_Array = std::string("u32_Array_") + zero;
    const std::string u64_Array = std::string("u64_Array_") + zero;
    const std::string r32_Array = std::string("r32_Array_") + zero;
    const std::string r64_Array = std::string("r64_Array_") + zero;
    const std::string r128_Array = std::string("r128_Array_") + zero;
    const std::string cr32_Array = std::string("cr32_Array_") + zero;
    const std::string cr64_Array = std::string("cr64_Array_") + zero;

    const std::string separator = "/";

    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, 0, 0);

// Write test data using BP
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        auto var = io.DefineVariable<int>("myVar");
        // Declare Single Value Attributes
        io.DefineAttribute<std::string>(s1_Array, currentTestData.S3.data(),
                                        currentTestData.S3.size(), var.Name());

        io.DefineAttribute<int8_t>(i8_Array, currentTestData.I8.data(), currentTestData.I8.size(),
                                   var.Name());
        io.DefineAttribute<int16_t>(i16_Array, currentTestData.I16.data(),
                                    currentTestData.I16.size(), var.Name());
        io.DefineAttribute<int32_t>(i32_Array, currentTestData.I32.data(),
                                    currentTestData.I32.size(), var.Name());
        io.DefineAttribute<int64_t>(i64_Array, currentTestData.I64.data(),
                                    currentTestData.I64.size(), var.Name());

        io.DefineAttribute<uint8_t>(u8_Array, currentTestData.U8.data(), currentTestData.U8.size(),
                                    var.Name());
        io.DefineAttribute<uint16_t>(u16_Array, currentTestData.U16.data(),
                                     currentTestData.U16.size(), var.Name());
        io.DefineAttribute<uint32_t>(u32_Array, currentTestData.U32.data(),
                                     currentTestData.U32.size(), var.Name());
        io.DefineAttribute<uint64_t>(u64_Array, currentTestData.U64.data(),
                                     currentTestData.U64.size(), var.Name());

        io.DefineAttribute<float>(r32_Array, currentTestData.R32.data(), currentTestData.R32.size(),
                                  var.Name());
        io.DefineAttribute<double>(r64_Array, currentTestData.R64.data(),
                                   currentTestData.R64.size(), var.Name());
        io.DefineAttribute<long double>(r128_Array, currentTestData.R128.data(),
                                        currentTestData.R128.size(), var.Name());

        io.DefineAttribute<std::complex<float>>(cr32_Array, currentTestData.CR32.data(),
                                                currentTestData.CR32.size(), var.Name());
        io.DefineAttribute<std::complex<double>>(cr64_Array, currentTestData.CR64.data(),
                                                 currentTestData.CR64.size(), var.Name());

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            io.SetEngine("File");
        }

        io.AddTransport("file");

        adios2::Engine engine = io.Open(fName, adios2::Mode::Write);
        engine.Put(var, 10);
        engine.Close();
    }

    {
        adios2::IO ioRead = adios.DeclareIO("ioRead");
        if (!engineName.empty())
        {
            ioRead.SetEngine(engineName);
        }

        adios2::Engine bpRead = ioRead.Open(fName, adios2::Mode::ReadRandomAccess);

        auto var = ioRead.InquireVariable<int>("myVar");

        auto attr_s1 = ioRead.InquireAttribute<std::string>(s1_Array, var.Name());

        auto attr_i8 = ioRead.InquireAttribute<int8_t>(i8_Array, var.Name());
        auto attr_i16 = ioRead.InquireAttribute<int16_t>(i16_Array, var.Name());
        auto attr_i32 = ioRead.InquireAttribute<int32_t>(i32_Array, var.Name());
        auto attr_i64 = ioRead.InquireAttribute<int64_t>(i64_Array, var.Name());

        auto attr_u8 = ioRead.InquireAttribute<uint8_t>(u8_Array, var.Name());
        auto attr_u16 = ioRead.InquireAttribute<uint16_t>(u16_Array, var.Name());
        auto attr_u32 = ioRead.InquireAttribute<uint32_t>(u32_Array, var.Name());
        auto attr_u64 = ioRead.InquireAttribute<uint64_t>(u64_Array, var.Name());

        auto attr_r32 = ioRead.InquireAttribute<float>(r32_Array, var.Name());
        auto attr_r64 = ioRead.InquireAttribute<double>(r64_Array, var.Name());
        auto attr_r128 = ioRead.InquireAttribute<long double>(r128_Array, var.Name());

        auto attr_cr32 = ioRead.InquireAttribute<std::complex<float>>(cr32_Array, var.Name());
        auto attr_cr64 = ioRead.InquireAttribute<std::complex<double>>(cr64_Array, var.Name());

        EXPECT_TRUE(attr_s1);
        ASSERT_EQ(attr_s1.Name(), var.Name() + separator + s1_Array);
        ASSERT_EQ(attr_s1.Data().size() == 1, false);
        ASSERT_EQ(attr_s1.Type(), adios2::GetType<std::string>());

        EXPECT_TRUE(attr_i8);
        ASSERT_EQ(attr_i8.Name(), var.Name() + separator + i8_Array);
        ASSERT_EQ(attr_i8.Data().size() == 1, false);
        ASSERT_EQ(attr_i8.Type(), adios2::GetType<int8_t>());

        EXPECT_TRUE(attr_i16);
        ASSERT_EQ(attr_i16.Name(), var.Name() + separator + i16_Array);
        ASSERT_EQ(attr_i16.Data().size() == 1, false);
        ASSERT_EQ(attr_i16.Type(), adios2::GetType<int16_t>());

        EXPECT_TRUE(attr_i32);
        ASSERT_EQ(attr_i32.Name(), var.Name() + separator + i32_Array);
        ASSERT_EQ(attr_i32.Data().size() == 1, false);
        ASSERT_EQ(attr_i32.Type(), adios2::GetType<int32_t>());

        EXPECT_TRUE(attr_i64);
        ASSERT_EQ(attr_i64.Name(), var.Name() + separator + i64_Array);
        ASSERT_EQ(attr_i64.Data().size() == 1, false);
        ASSERT_EQ(attr_i64.Type(), adios2::GetType<int64_t>());

        EXPECT_TRUE(attr_u8);
        ASSERT_EQ(attr_u8.Name(), var.Name() + separator + u8_Array);
        ASSERT_EQ(attr_u8.Data().size() == 1, false);
        ASSERT_EQ(attr_u8.Type(), adios2::GetType<uint8_t>());

        EXPECT_TRUE(attr_u16);
        ASSERT_EQ(attr_u16.Name(), var.Name() + separator + u16_Array);
        ASSERT_EQ(attr_u16.Data().size() == 1, false);
        ASSERT_EQ(attr_u16.Type(), adios2::GetType<uint16_t>());

        EXPECT_TRUE(attr_u32);
        ASSERT_EQ(attr_u32.Name(), var.Name() + separator + u32_Array);
        ASSERT_EQ(attr_u32.Data().size() == 1, false);
        ASSERT_EQ(attr_u32.Type(), adios2::GetType<uint32_t>());

        EXPECT_TRUE(attr_u64);
        ASSERT_EQ(attr_u64.Name(), var.Name() + separator + u64_Array);
        ASSERT_EQ(attr_u64.Data().size() == 1, false);
        ASSERT_EQ(attr_u64.Type(), adios2::GetType<uint64_t>());

        EXPECT_TRUE(attr_r32);
        ASSERT_EQ(attr_r32.Name(), var.Name() + separator + r32_Array);
        ASSERT_EQ(attr_r32.Data().size() == 1, false);
        ASSERT_EQ(attr_r32.Type(), adios2::GetType<float>());

        EXPECT_TRUE(attr_r64);
        ASSERT_EQ(attr_r64.Name(), var.Name() + separator + r64_Array);
        ASSERT_EQ(attr_r64.Data().size() == 1, false);
        ASSERT_EQ(attr_r64.Type(), adios2::GetType<double>());

        EXPECT_TRUE(attr_r128);
        ASSERT_EQ(attr_r128.Name(), var.Name() + separator + r128_Array);
        ASSERT_EQ(attr_r128.Data().size() == 1, false);
        ASSERT_EQ(attr_r128.Type(), adios2::GetType<long double>());

        EXPECT_TRUE(attr_cr32);
        ASSERT_EQ(attr_cr32.Name(), var.Name() + separator + cr32_Array);
        ASSERT_EQ(attr_cr32.Data().size() == 1, false);
        ASSERT_EQ(attr_cr32.Type(), adios2::GetType<std::complex<float>>());

        EXPECT_TRUE(attr_cr64);
        ASSERT_EQ(attr_cr64.Name(), var.Name() + separator + cr64_Array);
        ASSERT_EQ(attr_cr64.Data().size() == 1, false);
        ASSERT_EQ(attr_cr64.Type(), adios2::GetType<std::complex<double>>());

        auto I8 = attr_i8.Data();
        auto I16 = attr_i16.Data();
        auto I32 = attr_i32.Data();
        auto I64 = attr_i64.Data();

        auto U8 = attr_u8.Data();
        auto U16 = attr_u16.Data();
        auto U32 = attr_u32.Data();
        auto U64 = attr_u64.Data();

        const size_t Nx = 10;
        for (size_t i = 0; i < Nx; ++i)
        {
            EXPECT_EQ(I8[i], currentTestData.I8[i]);
            EXPECT_EQ(I16[i], currentTestData.I16[i]);
            EXPECT_EQ(I32[i], currentTestData.I32[i]);
            EXPECT_EQ(I64[i], currentTestData.I64[i]);

            EXPECT_EQ(U8[i], currentTestData.U8[i]);
            EXPECT_EQ(U16[i], currentTestData.U16[i]);
            EXPECT_EQ(U32[i], currentTestData.U32[i]);
            EXPECT_EQ(U64[i], currentTestData.U64[i]);
        }

        bpRead.Close();
    }
}

TEST_F(BPWriteReadAttributes, WriteReadStreamVarp)
{

    const std::string separator = "\\";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "AttributesWriteReadVar_MPI.bp";
#else
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "AttributesWriteReadVar.bp";
#endif

    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, 0, 0);

// Write test data using BP
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            io.SetEngine("FileStream");
        }

        auto var1 = io.DefineVariable<int32_t>("var1");
        auto var2 = io.DefineVariable<int32_t>("var2", shape, start, count);

        io.DefineAttribute<std::string>("sArray", currentTestData.S3.data(),
                                        currentTestData.S3.size(), var1.Name(), separator);
        io.DefineAttribute<std::string>("sArray", currentTestData.S3.data(),
                                        currentTestData.S3.size(), var2.Name(), separator);

        io.DefineAttribute<uint32_t>("u32Value", 1, var1.Name(), separator);
        io.DefineAttribute<uint32_t>("u32Value", 1, var2.Name(), separator);

#ifndef _WIN32
        io.DefineAttribute<std::string>("smile", "\u263A", var1.Name(), separator);
        io.DefineAttribute<std::string>("smile", "\u263A", var2.Name(), separator);

        io.DefineAttribute<std::string>("utf8", std::string("महसुस"), var1.Name(), separator);
        io.DefineAttribute<std::string>("utf8", std::string("महसुस"), var2.Name(), separator);
#endif
        adios2::Engine bpWriter = io.Open(fName, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            const int32_t step32 = static_cast<int32_t>(step);

            bpWriter.BeginStep();
            bpWriter.Put(var1, step32);
            if (step % 2 == 1)
            {
                bpWriter.Put(var2, currentTestData.I32.data());
            }
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }

    // reader
    {
        auto lf_VerifyAttributes = [](const std::string &variableName, const std::string separator,
                                      adios2::IO &io, const bool fullNames) {
            const std::map<std::string, adios2::Params> attributesInfo =
                io.AvailableAttributes(variableName, separator, fullNames);

            const std::string sArrayName =
                fullNames ? variableName + separator + "sArray" : "sArray";
            auto itSArray = attributesInfo.find(sArrayName);
            EXPECT_NE(itSArray, attributesInfo.end());
            EXPECT_EQ(itSArray->second.at("Type"), "string");
            EXPECT_EQ(itSArray->second.at("Elements"), "3");
            EXPECT_EQ(itSArray->second.at("Value"), R"({ "one", "two", "three" })");

            const std::string u32ValueName =
                fullNames ? variableName + separator + "u32Value" : "u32Value";
            auto itU32Value = attributesInfo.find(u32ValueName);
            EXPECT_NE(itU32Value, attributesInfo.end());
            EXPECT_EQ(itU32Value->second.at("Type"), "uint32_t");
            EXPECT_EQ(itU32Value->second.at("Elements"), "1");
            EXPECT_EQ(itU32Value->second.at("Value"), "1");

#ifndef _WIN32
            const std::string smileName = fullNames ? variableName + separator + "smile" : "smile";
            auto itSmile = attributesInfo.find(smileName);
            EXPECT_NE(itSmile, attributesInfo.end());
            EXPECT_EQ(itSmile->second.at("Type"), "string");
            EXPECT_EQ(itSmile->second.at("Elements"), "1");
            EXPECT_EQ(itSmile->second.at("Value"), std::string("\"\u263A\""));

            const std::string utf8Name = fullNames ? variableName + separator + "utf8" : "utf8";
            auto itUTF8 = attributesInfo.find(utf8Name);
            EXPECT_NE(itUTF8, attributesInfo.end());
            EXPECT_EQ(itUTF8->second.at("Type"), "string");
            EXPECT_EQ(itUTF8->second.at("Elements"), "1");
            EXPECT_EQ(itUTF8->second.at("Value"), "\"" + std::string("महसुस") + "\"");
#endif
        };

        adios2::IO io = adios.DeclareIO("ReaderIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            io.SetEngine("FileStream");
        }
        adios2::Engine bpReader = io.Open(fName, adios2::Mode::Read);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var1 = io.InquireVariable<int32_t>("var1");
            if (var1)
            {
                lf_VerifyAttributes("var1", separator, io, false);
                lf_VerifyAttributes("var1", separator, io, true);
            }

            auto var2 = io.InquireVariable<int32_t>("var2");
            if (var2)
            {
                lf_VerifyAttributes("var2", separator, io, false);
                lf_VerifyAttributes("var2", separator, io, true);
            }

            bpReader.EndStep();
        }
    }
}

TEST_F(BPWriteReadAttributes, WriteReadStreamModifiable)
{

    const std::string separator = "\\";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "AttributesWriteReadModifiable_MPI.bp";
#else
    const std::string fName =
        "foo" + std::string(&adios2::PathSeparator, 1) + "AttributesWriteReadModifiable.bp";
#endif

    const double d3[3] = {-1.1, -1.2, -1.3};
    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, 0, 0);

// Write test data using BP
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        adios2::IO io = adios.DeclareIO("TestIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            io.SetEngine("FileStream");
        }

        auto var1 = io.DefineVariable<int32_t>("var1");
        auto var2 = io.DefineVariable<int32_t>("var2", shape, start, count);

        io.DefineAttribute<double>("dArray", d3, 3, var1.Name(), separator, true);
        io.DefineAttribute<double>("dArray", d3, 3, var2.Name(), separator, true);

        io.DefineAttribute<int32_t>("i32Value", -1, var1.Name(), separator, true);
        io.DefineAttribute<int32_t>("i32Value", -1, var2.Name(), separator, true);

        adios2::Engine bpWriter = io.Open(fName, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            const int32_t step32 = static_cast<int32_t>(step);
            const double stepD = static_cast<double>(step);
            double d[3] = {stepD + 0.1, stepD + 0.2, stepD + 0.3};

            bpWriter.BeginStep();

            io.DefineAttribute<double>("dArray", d, 3, var1.Name(), separator, true);
            io.DefineAttribute<int32_t>("i32Value", step32, var1.Name(), separator, true);
            bpWriter.Put(var1, step32);

            if (step % 2 == 0)
            {
                bpWriter.Put(var2, currentTestData.I32.data());
                io.DefineAttribute<double>("dArray", d, 3, var2.Name(), separator, true);
                io.DefineAttribute<int32_t>("i32Value", step32, var2.Name(), separator, true);
            }

            bpWriter.EndStep();
        }
        bpWriter.Close();
    }

    // reader
    {
        auto lf_VerifyAttributes = [](const int32_t step, const std::string &variableName,
                                      const std::string separator, adios2::IO &io) {
            const std::map<std::string, adios2::Params> attributesInfo =
                io.AvailableAttributes(variableName, separator, false);

            const double stepD = static_cast<double>(step);
            const double d[3] = {stepD + 0.1, stepD + 0.2, stepD + 0.3};

            auto itDArray = attributesInfo.find("dArray");
            EXPECT_NE(itDArray, attributesInfo.end());
            EXPECT_EQ(itDArray->second.at("Type"), "double");
            EXPECT_EQ(itDArray->second.at("Elements"), "3");

            auto a = io.InquireAttribute<double>("dArray", variableName, separator);
            auto adata = a.Data();
            for (int i = 0; i < 3; ++i)
            {
                EXPECT_EQ(adata[i], d[i]);
            }

            const std::string stepS = std::to_string(step);
            auto iti32Value = attributesInfo.find("i32Value");
            EXPECT_NE(iti32Value, attributesInfo.end());
            EXPECT_EQ(iti32Value->second.at("Type"), "int32_t");
            EXPECT_EQ(iti32Value->second.at("Elements"), "1");
            EXPECT_EQ(iti32Value->second.at("Value"), stepS);
        };

        adios2::IO io = adios.DeclareIO("ReaderIO");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
            io.SetParameter("StreamReader", "ON");
        }
        else
        {
            io.SetEngine("FileStream");
        }
        adios2::Engine bpReader = io.Open(fName, adios2::Mode::Read);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            int32_t step = static_cast<int32_t>(bpReader.CurrentStep());
            if (engineName == "BP3")
            {
                // BP3 does not support changing attributes
                step = 0;
            }
            auto var1 = io.InquireVariable<int32_t>("var1");
            if (var1)
            {
                lf_VerifyAttributes(step, "var1", separator, io);
            }

            auto var2 = io.InquireVariable<int32_t>("var2");
            if (var2)
            {
                lf_VerifyAttributes(step, "var2", separator, io);
            }

            bpReader.EndStep();
        }
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
