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

class StreamWriteReadHighLevelAPI : public ::testing::Test
{
public:
    StreamWriteReadHighLevelAPI() = default;

    SmallTestData m_TestData;
};

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

TEST_F(StreamWriteReadHighLevelAPI, ADIOS2BPWriteRead1D8)
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
    const std::string fname("ADIOS2BPWriteRead1D8_hl_MPI.bp");
#else
    const std::string fname("ADIOS2BPWriteRead1D8_hl.bp");
#endif

    // write test data using BP
    {
#if ADIOS2_USE_MPI

        adios2::fstream oStream(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);

#else
        adios2::fstream oStream(fname, adios2::fstream::out, engineName);
#endif

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData stepData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            if (step == 0)
            {
                // write globals
                oStream.write("gi8", m_TestData.I8.front());
                oStream.write("gi16", m_TestData.I16.front());
                oStream.write("gi32", m_TestData.I32.front());
                oStream.write("gi64", m_TestData.I64.front());
                oStream.write("gu8", m_TestData.U8.front());
                oStream.write("gu16", m_TestData.U16.front());
                oStream.write("gu32", m_TestData.U32.front());
                oStream.write("gu64", m_TestData.U64.front());
                oStream.write("gr32", m_TestData.R32.front());
                oStream.write("gr64", m_TestData.R64.front());

                oStream.write_attribute<std::string>("attrStr", "StringValue");
                oStream.write_attribute("attri8", m_TestData.I8.front());
                oStream.write_attribute("attri16", m_TestData.I16.front());
                oStream.write_attribute("attri32", m_TestData.I32.front());
                oStream.write_attribute("attri64", m_TestData.I64.front());
                oStream.write_attribute("attru8", m_TestData.U8.front());
                oStream.write_attribute("attru16", m_TestData.U16.front());
                oStream.write_attribute("attru32", m_TestData.U32.front());
                oStream.write_attribute("attru64", m_TestData.U64.front());
                oStream.write_attribute("attrr32", m_TestData.R32.front());
                oStream.write_attribute("attrr64", m_TestData.R64.front());

                oStream.write_attribute("attrStrarray", m_TestData.S3.data(), m_TestData.S3.size());
                oStream.write_attribute("attri8array", m_TestData.I8.data(), m_TestData.I8.size());
                oStream.write_attribute("attri16array", m_TestData.I16.data(),
                                        m_TestData.I16.size());
                oStream.write_attribute("attri32array", m_TestData.I32.data(),
                                        m_TestData.I32.size());
                oStream.write_attribute("attri64array", m_TestData.I64.data(),
                                        m_TestData.I64.size());
                oStream.write_attribute("attru8array", m_TestData.U8.data(), m_TestData.U8.size());
                oStream.write_attribute("attru16array", m_TestData.U16.data(),
                                        m_TestData.U16.size());
                oStream.write_attribute("attru32array", m_TestData.U32.data(),
                                        m_TestData.U32.size());
                oStream.write_attribute("attru64array", m_TestData.U64.data(),
                                        m_TestData.U64.size());
                oStream.write_attribute("attrr32array", m_TestData.R32.data(),
                                        m_TestData.R32.size());
                oStream.write_attribute("attrr64array", m_TestData.R64.data(),
                                        m_TestData.R64.size());
            }

            oStream.write("iString", stepData.S1);
            oStream.write("i8", stepData.I8.data(), shape, start, count);
            oStream.write("i16", stepData.I16.data(), shape, start, count);
            oStream.write("i32", stepData.I32.data(), shape, start, count);
            oStream.write("i64", stepData.I64.data(), shape, start, count);
            oStream.write("u8", stepData.U8.data(), shape, start, count);
            oStream.write("u16", stepData.U16.data(), shape, start, count);
            oStream.write("u32", stepData.U32.data(), shape, start, count);
            oStream.write("u64", stepData.U64.data(), shape, start, count);
            oStream.write("r32", stepData.R32.data(), shape, start, count);
            oStream.write("r64", stepData.R64.data(), shape, start, count);

            if (step == 0)
            {
                oStream.write_attribute<std::string>("attrStr", m_TestData.S1, "iString");
                oStream.write_attribute("attri8", m_TestData.I8.front(), "i8");
                oStream.write_attribute("attri16", m_TestData.I16.front(), "i16");
                oStream.write_attribute("attri32", m_TestData.I32.front(), "i32");
                oStream.write_attribute("attri64", m_TestData.I64.front(), "i64");
                oStream.write_attribute("attru8", m_TestData.U8.front(), "u8");
                oStream.write_attribute("attru16", m_TestData.U16.front(), "u16");
                oStream.write_attribute("attru32", m_TestData.U32.front(), "u32");
                oStream.write_attribute("attru64", m_TestData.U64.front(), "u64");
                oStream.write_attribute("attrr32", m_TestData.R32.front(), "r32");
                oStream.write_attribute("attrr64", m_TestData.R64.front(), "r64");

                oStream.write_attribute("attrStrarray", m_TestData.S3.data(), m_TestData.S3.size(),
                                        "iString", "::");
                oStream.write_attribute("attri8array", m_TestData.I8.data(), m_TestData.I8.size(),
                                        "i8", "::");
                oStream.write_attribute("attri16array", m_TestData.I16.data(),
                                        m_TestData.I16.size(), "i16", "::");
                oStream.write_attribute("attri32array", m_TestData.I32.data(),
                                        m_TestData.I32.size(), "i32", "::");
                oStream.write_attribute("attri64array", m_TestData.I64.data(),
                                        m_TestData.I64.size(), "i64", "::");
                oStream.write_attribute("attru8array", m_TestData.U8.data(), m_TestData.U8.size(),
                                        "u8", "::");
                oStream.write_attribute("attru16array", m_TestData.U16.data(),
                                        m_TestData.U16.size(), "u16", "::");
                oStream.write_attribute("attru32array", m_TestData.U32.data(),
                                        m_TestData.U32.size(), "u32", "::");
                oStream.write_attribute("attru64array", m_TestData.U64.data(),
                                        m_TestData.U64.size(), "u64", "::");
                oStream.write_attribute("attrr32array", m_TestData.R32.data(),
                                        m_TestData.R32.size(), "r32", "::");
                oStream.write_attribute("attrr64array", m_TestData.R64.data(),
                                        m_TestData.R64.size(), "r64", "::");
            }

            oStream.end_step();
        }
        oStream.close();
    }

    // READ
    {
        adios2::fstream iStream;
        EXPECT_FALSE(iStream);

#if ADIOS2_USE_MPI
        iStream.open(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        iStream.open(fname, adios2::fstream::in, engineName);
#endif

        EXPECT_TRUE(iStream);
        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        int8_t gi8 = 0;
        int16_t gi16 = 0;
        int32_t gi32 = 0;
        int64_t gi64 = 0;
        uint8_t gu8 = 0;
        uint16_t gu16 = 0;
        uint32_t gu32 = 0;
        uint64_t gu64 = 0;
        float gr32 = -1.f;
        double gr64 = -1.f;

        EXPECT_EQ(iStream.steps(), NSteps);

        size_t t = 0;
        for (adios2::fstep iStep; adios2::getstep(iStream, iStep);)
        {
            if (iStep.current_step() == 0)
            {
                iStep.read("gi8", gi8);
                iStep.read("gi16", gi16);
                iStep.read("gi32", gi32);
                iStep.read("gi64", gi64);
                iStep.read("gu8", gu8);
                iStep.read("gu16", gu16);
                iStep.read("gu32", gu32);
                iStep.read("gu64", gu64);
                iStep.read("gr32", gr32);
                iStep.read("gr64", gr64);

                EXPECT_EQ(gi8, m_TestData.I8.front());
                EXPECT_EQ(gi16, m_TestData.I16.front());
                EXPECT_EQ(gi32, m_TestData.I32.front());
                EXPECT_EQ(gi64, m_TestData.I64.front());
                EXPECT_EQ(gu8, m_TestData.U8.front());
                EXPECT_EQ(gu16, m_TestData.U16.front());
                EXPECT_EQ(gu32, m_TestData.U32.front());
                EXPECT_EQ(gu64, m_TestData.U64.front());
                EXPECT_EQ(gr32, m_TestData.R32.front());
                EXPECT_EQ(gr64, m_TestData.R64.front());

                auto vgi8 = iStep.read<int8_t>("gi8");
                auto vgi16 = iStep.read<int16_t>("gi16");
                auto vgi32 = iStep.read<int32_t>("gi32");
                auto vgi64 = iStep.read<int64_t>("gi64");
                auto vgu8 = iStep.read<uint8_t>("gu8");
                auto vgu16 = iStep.read<uint16_t>("gu16");
                auto vgu32 = iStep.read<uint32_t>("gu32");
                auto vgu64 = iStep.read<uint64_t>("gu64");
                auto vgr32 = iStep.read<float>("gr32");
                auto vgr64 = iStep.read<double>("gr64");

                EXPECT_EQ(vgi8.front(), m_TestData.I8.front());
                EXPECT_EQ(vgi16.front(), m_TestData.I16.front());
                EXPECT_EQ(vgi32.front(), m_TestData.I32.front());
                EXPECT_EQ(vgi64.front(), m_TestData.I64.front());
                EXPECT_EQ(vgu8.front(), m_TestData.U8.front());
                EXPECT_EQ(vgu16.front(), m_TestData.U16.front());
                EXPECT_EQ(vgu32.front(), m_TestData.U32.front());
                EXPECT_EQ(vgu64.front(), m_TestData.U64.front());
                EXPECT_EQ(vgr32.front(), m_TestData.R32.front());
                EXPECT_EQ(vgr64.front(), m_TestData.R64.front());

                auto vattrStr = iStep.read_attribute<std::string>("attrStr");
                auto vattri8 = iStep.read_attribute<int8_t>("attri8");
                auto vattri16 = iStep.read_attribute<int16_t>("attri16");
                auto vattri32 = iStep.read_attribute<int32_t>("attri32");
                auto vattri64 = iStep.read_attribute<int64_t>("attri64");
                auto vattru8 = iStep.read_attribute<uint8_t>("attru8");
                auto vattru16 = iStep.read_attribute<uint16_t>("attru16");
                auto vattru32 = iStep.read_attribute<uint32_t>("attru32");
                auto vattru64 = iStep.read_attribute<uint64_t>("attru64");
                auto vattrr32 = iStep.read_attribute<float>("attrr32");
                auto vattrr64 = iStep.read_attribute<double>("attrr64");

                EXPECT_EQ(vattrStr.front(), "StringValue");
                EXPECT_EQ(vattri8.front(), m_TestData.I8.front());
                EXPECT_EQ(vattri16.front(), m_TestData.I16.front());
                EXPECT_EQ(vattri32.front(), m_TestData.I32.front());
                EXPECT_EQ(vattri64.front(), m_TestData.I64.front());
                EXPECT_EQ(vattru8.front(), m_TestData.U8.front());
                EXPECT_EQ(vattru16.front(), m_TestData.U16.front());
                EXPECT_EQ(vattru32.front(), m_TestData.U32.front());
                EXPECT_EQ(vattru64.front(), m_TestData.U64.front());
                EXPECT_EQ(vattrr32.front(), m_TestData.R32.front());
                EXPECT_EQ(vattrr64.front(), m_TestData.R64.front());

                auto vattrStrarray = iStep.read_attribute<std::string>("attrStrarray");
                auto vattri8array = iStep.read_attribute<int8_t>("attri8array");
                auto vattri16array = iStep.read_attribute<int16_t>("attri16array");
                auto vattri32array = iStep.read_attribute<int32_t>("attri32array");
                auto vattri64array = iStep.read_attribute<int64_t>("attri64array");
                auto vattru8array = iStep.read_attribute<uint8_t>("attru8array");
                auto vattru16array = iStep.read_attribute<uint16_t>("attru16array");
                auto vattru32array = iStep.read_attribute<uint32_t>("attru32array");
                auto vattru64array = iStep.read_attribute<uint64_t>("attru64array");
                auto vattrr32array = iStep.read_attribute<float>("attrr32array");
                auto vattrr64array = iStep.read_attribute<double>("attrr64array");

                for (size_t i = 0; i < vattrStrarray.size(); ++i)
                {
                    EXPECT_EQ(vattrStrarray[i], m_TestData.S3[i]);
                }

                for (size_t i = 0; i < vattri8array.size(); ++i)
                {
                    EXPECT_EQ(vattri8array[i], m_TestData.I8[i]);
                    EXPECT_EQ(vattri16array[i], m_TestData.I16[i]);
                    EXPECT_EQ(vattri32array[i], m_TestData.I32[i]);
                    EXPECT_EQ(vattri64array[i], m_TestData.I64[i]);
                    EXPECT_EQ(vattru8array[i], m_TestData.U8[i]);
                    EXPECT_EQ(vattru16array[i], m_TestData.U16[i]);
                    EXPECT_EQ(vattru32array[i], m_TestData.U32[i]);
                    EXPECT_EQ(vattru64array[i], m_TestData.U64[i]);
                    EXPECT_EQ(vattrr32array[i], m_TestData.R32[i]);
                    EXPECT_EQ(vattrr64array[i], m_TestData.R64[i]);
                }

                // var attributes
                auto vvarattrStr = iStep.read_attribute<std::string>("attrStr", "iString");
                auto vvarattri8 = iStep.read_attribute<int8_t>("attri8", "i8");
                auto vvarattri16 = iStep.read_attribute<int16_t>("attri16", "i16");
                auto vvarattri32 = iStep.read_attribute<int32_t>("attri32", "i32");
                auto vvarattri64 = iStep.read_attribute<int64_t>("attri64", "i64");
                auto vvarattru8 = iStep.read_attribute<uint8_t>("attru8", "u8");
                auto vvarattru16 = iStep.read_attribute<uint16_t>("attru16", "u16");
                auto vvarattru32 = iStep.read_attribute<uint32_t>("attru32", "u32");
                auto vvarattru64 = iStep.read_attribute<uint64_t>("attru64", "u64");
                auto vvarattrr32 = iStep.read_attribute<float>("attrr32", "r32");
                auto vvarattrr64 = iStep.read_attribute<double>("attrr64", "r64");

                EXPECT_EQ(vvarattrStr.front(), m_TestData.S1);
                EXPECT_EQ(vvarattri8.front(), m_TestData.I8.front());
                EXPECT_EQ(vvarattri16.front(), m_TestData.I16.front());
                EXPECT_EQ(vvarattri32.front(), m_TestData.I32.front());
                EXPECT_EQ(vvarattri64.front(), m_TestData.I64.front());
                EXPECT_EQ(vvarattru8.front(), m_TestData.U8.front());
                EXPECT_EQ(vvarattru16.front(), m_TestData.U16.front());
                EXPECT_EQ(vvarattru32.front(), m_TestData.U32.front());
                EXPECT_EQ(vvarattru64.front(), m_TestData.U64.front());
                EXPECT_EQ(vvarattrr32.front(), m_TestData.R32.front());
                EXPECT_EQ(vvarattrr64.front(), m_TestData.R64.front());

                auto vvarattrStrarray =
                    iStep.read_attribute<std::string>("attrStrarray", "iString", "::");
                auto vvarattri8array = iStep.read_attribute<int8_t>("attri8array", "i8", "::");
                auto vvarattri16array = iStep.read_attribute<int16_t>("attri16array", "i16", "::");
                auto vvarattri32array = iStep.read_attribute<int32_t>("attri32array", "i32", "::");
                auto vvarattri64array = iStep.read_attribute<int64_t>("attri64array", "i64", "::");
                auto vvarattru8array = iStep.read_attribute<uint8_t>("attru8array", "u8", "::");
                auto vvarattru16array = iStep.read_attribute<uint16_t>("attru16array", "u16", "::");
                auto vvarattru32array = iStep.read_attribute<uint32_t>("attru32array", "u32", "::");
                auto vvarattru64array = iStep.read_attribute<uint64_t>("attru64array", "u64", "::");
                auto vvarattrr32array = iStep.read_attribute<float>("attrr32array", "r32", "::");
                auto vvarattrr64array = iStep.read_attribute<double>("attrr64array", "r64", "::");

                for (size_t i = 0; i < vvarattrStrarray.size(); ++i)
                {
                    EXPECT_EQ(vvarattrStrarray[i], m_TestData.S3[i]);
                }

                for (size_t i = 0; i < vvarattri8array.size(); ++i)
                {
                    EXPECT_EQ(vvarattri8array[i], m_TestData.I8[i]);
                    EXPECT_EQ(vvarattri16array[i], m_TestData.I16[i]);
                    EXPECT_EQ(vvarattri32array[i], m_TestData.I32[i]);
                    EXPECT_EQ(vvarattri64array[i], m_TestData.I64[i]);
                    EXPECT_EQ(vvarattru8array[i], m_TestData.U8[i]);
                    EXPECT_EQ(vvarattru16array[i], m_TestData.U16[i]);
                    EXPECT_EQ(vvarattru32array[i], m_TestData.U32[i]);
                    EXPECT_EQ(vvarattru64array[i], m_TestData.U64[i]);
                    EXPECT_EQ(vvarattrr32array[i], m_TestData.R32[i]);
                    EXPECT_EQ(vvarattrr64array[i], m_TestData.R64[i]);
                }
            }

            auto IString = iStep.read<std::string>("iString");

            std::vector<int8_t> I8(Nx);
            std::vector<int16_t> I16(Nx);
            std::vector<int32_t> I32(Nx);
            std::vector<int64_t> I64(Nx);

            iStep.read<int8_t>("i8", I8.data(), start, count);
            iStep.read<int16_t>("i16", I16.data(), start, count);
            iStep.read<int32_t>("i32", I32.data(), start, count);
            iStep.read<int64_t>("i64", I64.data(), start, count);
            auto U8 = iStep.read<uint8_t>("u8", start, count);
            auto U16 = iStep.read<uint16_t>("u16", start, count);
            auto U32 = iStep.read<uint32_t>("u32", start, count);
            auto U64 = iStep.read<uint64_t>("u64", start, count);
            auto R32 = iStep.read<float>("r32", start, count);
            auto R64 = iStep.read<double>("r64", start, count);

            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString.front(), currentTestData.S1);

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

        iStream.close();
    }
}

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

TEST_F(StreamWriteReadHighLevelAPI, ADIOS2BPwriteRead2D2x4)
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
    const std::string fname("ADIOS2BPwriteRead2D2x4Test_hl_MPI.bp");
#else
    const std::string fname("ADIOS2BPwriteRead2D2x4Test_hl.bp");
#endif

    // write test data using ADIOS2
    {
#if ADIOS2_USE_MPI

        adios2::fstream oStream(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);

#else
        adios2::fstream oStream(fname, adios2::fstream::out, engineName);
#endif

        const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData stepData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            oStream.write("iString", stepData.S1);
            oStream.write("i8", stepData.I8.data(), shape, start, count);
            oStream.write("i16", stepData.I16.data(), shape, start, count);
            oStream.write("i32", stepData.I32.data(), shape, start, count);
            oStream.write("i64", stepData.I64.data(), shape, start, count);
            oStream.write("u8", stepData.U8.data(), shape, start, count);
            oStream.write("u16", stepData.U16.data(), shape, start, count);
            oStream.write("u32", stepData.U32.data(), shape, start, count);
            oStream.write("u64", stepData.U64.data(), shape, start, count);
            oStream.write("r32", stepData.R32.data(), shape, start, count);
            oStream.write("r64", stepData.R64.data(), shape, start, count, adios2::end_step);
        }

        // Close the file
        oStream.close();
    }

    // READ
    {
#if ADIOS2_USE_MPI
        adios2::fstream iStream(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream iStream(fname, adios2::fstream::in, engineName);
#endif
        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};
        const adios2::Box<adios2::Dims> sel(start, count);

        EXPECT_EQ(iStream.steps(), NSteps);

        size_t t = 0;
        adios2::fstep iStep;
        while (adios2::getstep(iStream, iStep))
        {
            auto IString = iStep.read<std::string>("iString");
            auto I8 = iStep.read<int8_t>("i8", start, count);
            auto I16 = iStep.read<int16_t>("i16", start, count);
            auto I32 = iStep.read<int32_t>("i32", start, count);
            auto I64 = iStep.read<int64_t>("i64", start, count);
            auto U8 = iStep.read<uint8_t>("u8", start, count);
            auto U16 = iStep.read<uint16_t>("u16", start, count);
            auto U32 = iStep.read<uint32_t>("u32", start, count);
            auto U64 = iStep.read<uint64_t>("u64", start, count);
            auto R32 = iStep.read<float>("r32", start, count);
            auto R64 = iStep.read<double>("r64", start, count);

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString.front(), currentTestData.S1);

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
        iStream.close();
    }
}

//******************************************************************************
// 2D 4x2 test data
//******************************************************************************

TEST_F(StreamWriteReadHighLevelAPI, ADIOS2BPwriteRead2D4x2)
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
    const std::string fname("ADIOS2BPwriteRead2D4x2Test_hl_MPI.bp");
#else
    const std::string fname("ADIOS2BPwriteRead2D4x2Test_hl.bp");
#endif

    // write test data using ADIOS2
    {
#if ADIOS2_USE_MPI

        adios2::fstream oStream(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);

#else
        adios2::fstream oStream(fname, adios2::fstream::out, engineName);
#endif

        // Declare 2D variables (4 * (NumberOfProcess * Nx))
        // The local process' part (start, count) can be defined now or
        // later
        // before write().
        adios2::Dims shape{Ny, static_cast<std::size_t>(mpiSize * Nx)};
        adios2::Dims start{0, static_cast<std::size_t>(mpiRank * Nx)};
        adios2::Dims count{Ny, Nx};

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData stepData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            oStream.write("iString", stepData.S1);
            oStream.write("i8", stepData.I8.data(), shape, start, count);
            oStream.write("i16", stepData.I16.data(), shape, start, count);
            oStream.write("i32", stepData.I32.data(), shape, start, count);
            oStream.write("i64", stepData.I64.data(), shape, start, count);
            oStream.write("u8", stepData.U8.data(), shape, start, count);
            oStream.write("u16", stepData.U16.data(), shape, start, count);
            oStream.write("u32", stepData.U32.data(), shape, start, count);
            oStream.write("u64", stepData.U64.data(), shape, start, count);
            oStream.write("r32", stepData.R32.data(), shape, start, count);
            oStream.write("r64", stepData.R64.data(), shape, start, count, adios2::end_step);
        }

        EXPECT_THROW(oStream.write<int16_t>("i8", 1), std::invalid_argument);
        oStream.close();
    }

    {
#if ADIOS2_USE_MPI
        adios2::fstream iStream(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream iStream(fname, adios2::fstream::in, engineName);
#endif

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};
        const adios2::Box<adios2::Dims> sel(start, count);

        EXPECT_EQ(iStream.steps(), NSteps);
        size_t t = 0;

        for (adios2::fstep iStep; adios2::getstep(iStream, iStep);)
        {
            auto IString = iStep.read<std::string>("iString");
            auto I8 = iStep.read<int8_t>("i8", start, count);
            auto I16 = iStep.read<int16_t>("i16", start, count);
            auto I32 = iStep.read<int32_t>("i32", start, count);
            auto I64 = iStep.read<int64_t>("i64", start, count);
            auto U8 = iStep.read<uint8_t>("u8", start, count);
            auto U16 = iStep.read<uint16_t>("u16", start, count);
            auto U32 = iStep.read<uint32_t>("u32", start, count);
            auto U64 = iStep.read<uint64_t>("u64", start, count);
            auto R32 = iStep.read<float>("r32", start, count);
            auto R64 = iStep.read<double>("r64", start, count);

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString.front(), currentTestData.S1);

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
        iStream.close();
    }
}

TEST_F(StreamWriteReadHighLevelAPI, DoubleOpenException)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    {
#if ADIOS2_USE_MPI
        const std::string fname("ADIOS2BP_hl_exception_MPI.bp");

        adios2::fstream oStream(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
        EXPECT_THROW(oStream.open("second", adios2::fstream::out, MPI_COMM_WORLD, engineName),
                     std::invalid_argument);

#else
        const std::string fname("ADIOS2BP_hl_exception.bp");

        adios2::fstream oStream(fname, adios2::fstream::out);
        EXPECT_THROW(oStream.open("second", adios2::fstream::out, engineName),
                     std::invalid_argument);
#endif
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
