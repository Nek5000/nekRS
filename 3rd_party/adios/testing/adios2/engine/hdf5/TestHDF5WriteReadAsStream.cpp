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

class HDF5WriteReadAsStreamTestADIOS2 : public ::testing::Test
{
public:
    HDF5WriteReadAsStreamTestADIOS2() = default;

    SmallTestData m_TestData;
};

TEST_F(HDF5WriteReadAsStreamTestADIOS2, ADIOS2HDF5WriteRead1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 5;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2HDF5WriteReadAsStream1D8_MPI.h5");
#else
    const std::string fname("ADIOS2HDF5WriteReadAsStream1D8.h5");
#endif

    // Write test data using HDF5

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetEngine("HDF5");

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
            const adios2::Dims count{Nx};

            io.DefineVariable<std::string>("iString");
            io.DefineVariable<char>("ch", shape, start, count);

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

        adios2::Engine h5Writer = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // EXPECT_EQ(h5Writer.CurrentStep(), step);

            h5Writer.BeginStep();

            if (step == 0)
            {
                h5Writer.Put<std::string>("iString", currentTestData.S1);
            }

            if (step % 2 == 0)
            {
                h5Writer.Put<int8_t>("i8", currentTestData.I8.data());
                h5Writer.Put<int16_t>("i16", currentTestData.I16.data());
                h5Writer.Put<int32_t>("i32", currentTestData.I32.data());
                h5Writer.Put<int64_t>("i64", currentTestData.I64.data());
            }
            if (step % 2 == 1)
            {
                h5Writer.Put<uint8_t>("u8", currentTestData.U8.data());
                h5Writer.Put<uint16_t>("u16", currentTestData.U16.data());
                h5Writer.Put<uint32_t>("u32", currentTestData.U32.data());
                h5Writer.Put<uint64_t>("u64", currentTestData.U64.data());
            }

            if (step < NSteps - 1) // all but the last
            {
                h5Writer.Put<float>("r32", currentTestData.R32.data());
                h5Writer.Put<double>("r64", currentTestData.R64.data());
            }

            h5Writer.Put<char>("ch", currentTestData.CHAR.data());
            h5Writer.Put<std::complex<float>>("cr32", currentTestData.CR32.data());
            h5Writer.Put<std::complex<double>>("cr64", currentTestData.CR64.data());

            h5Writer.EndStep();
        }

        h5Writer.Close();
    }

    // if (false)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine h5Reader = io.Open(fname, adios2::Mode::Read);

        std::string IString;
        std::array<uint8_t, Nx> CHAR;
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

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        size_t t = 0;

        while (h5Reader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = h5Reader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            auto var_iString = io.InquireVariable<std::string>("iString");
            auto var_ch = io.InquireVariable<uint8_t>("ch");
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

            if (currentStep == 0)
            {
                EXPECT_TRUE(var_iString);
                EXPECT_TRUE(var_iString);
                ASSERT_EQ(var_iString.Shape().size(), 0);
                ASSERT_EQ(var_iString.Steps(), 1);
            }
            else
            {
                EXPECT_FALSE(var_iString);
            }

            if (currentStep % 2 == 0)
            {
                EXPECT_TRUE(var_i8);
                EXPECT_TRUE(var_i16);
                EXPECT_TRUE(var_i32);
                EXPECT_TRUE(var_i64);
                EXPECT_FALSE(var_u8);
                EXPECT_FALSE(var_u16);
                EXPECT_FALSE(var_u32);
                EXPECT_FALSE(var_u64);

                ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i8.Steps(), NSteps / 2 + NSteps % 2);
                ASSERT_EQ(var_i8.Shape()[0], static_cast<size_t>(mpiSize * Nx));

                ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i16.Steps(), NSteps / 2 + NSteps % 2);
                ASSERT_EQ(var_i16.Shape()[0], static_cast<size_t>(mpiSize * Nx));

                ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i32.Steps(), NSteps / 2 + NSteps % 2);
                ASSERT_EQ(var_i32.Shape()[0], static_cast<size_t>(mpiSize * Nx));

                ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_i64.Steps(), NSteps / 2 + NSteps % 2);
                ASSERT_EQ(var_i64.Shape()[0], static_cast<size_t>(mpiSize * Nx));

                var_i8.SetSelection(sel);
                var_i16.SetSelection(sel);
                var_i32.SetSelection(sel);
                var_i64.SetSelection(sel);

                h5Reader.Get(var_i8, I8.data());
                h5Reader.Get(var_i16, I16.data());
                h5Reader.Get(var_i32, I32.data());
                h5Reader.Get(var_i64, I64.data());
            }

            if (currentStep % 2 == 1)
            {
                EXPECT_FALSE(var_i8);
                EXPECT_FALSE(var_i16);
                EXPECT_FALSE(var_i32);
                EXPECT_FALSE(var_i64);
                EXPECT_TRUE(var_u8);
                EXPECT_TRUE(var_u16);
                EXPECT_TRUE(var_u32);
                EXPECT_TRUE(var_u64);

                ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u8.Steps(), NSteps / 2);
                ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx);

                ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u16.Steps(), NSteps / 2);
                ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx);

                ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u32.Steps(), NSteps / 2);
                ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx);

                ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_u64.Steps(), NSteps / 2);
                ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx);

                var_u8.SetSelection(sel);
                var_u16.SetSelection(sel);
                var_u32.SetSelection(sel);
                var_u64.SetSelection(sel);

                h5Reader.Get(var_u8, U8.data());
                h5Reader.Get(var_u16, U16.data());
                h5Reader.Get(var_u32, U32.data());
                h5Reader.Get(var_u64, U64.data());
            }

            if (currentStep == NSteps - 1)
            {
                EXPECT_FALSE(var_r32);
                EXPECT_FALSE(var_r64);
            }
            else
            {
                EXPECT_TRUE(var_r32);
                EXPECT_TRUE(var_r64);

                ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r32.Steps(), NSteps - 1);
                ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

                ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
                ASSERT_EQ(var_r64.Steps(), NSteps - 1);
                ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

                var_r32.SetSelection(sel);
                var_r64.SetSelection(sel);

                h5Reader.Get(var_r32, R32.data());
                h5Reader.Get(var_r64, R64.data());
            }

            EXPECT_TRUE(var_ch);
            EXPECT_TRUE(var_cr32);
            EXPECT_TRUE(var_cr64);

            ASSERT_EQ(var_ch.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_ch.Steps(), NSteps);
            ASSERT_EQ(var_ch.Shape()[0], mpiSize * Nx);

            ASSERT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr32.Steps(), NSteps);
            ASSERT_EQ(var_cr32.Shape()[0], mpiSize * Nx);

            ASSERT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr64.Steps(), NSteps);
            ASSERT_EQ(var_cr64.Shape()[0], mpiSize * Nx);

            var_ch.SetSelection(sel);
            var_cr32.SetSelection(sel);
            var_cr64.SetSelection(sel);

            h5Reader.Get(var_ch, CHAR.data());
            h5Reader.Get(var_cr32, CR32.data());
            h5Reader.Get(var_cr64, CR64.data());

            h5Reader.EndStep();

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

                if (var_ch)
                {
                    EXPECT_EQ(static_cast<char>(CHAR[i]), currentTestData.CHAR[i]) << msg;
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
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        h5Reader.Close();
    }
}

TEST_F(HDF5WriteReadAsStreamTestADIOS2, ADIOS2HDF5WriteRead2D2x4)
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

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2HDF5WriteReadAsStream2D2x4Test_MPI.h5");
#else
    const std::string fname("ADIOS2HDF5WriteReadAsStream2D2x4Test.h5");
#endif

    // Write test data using ADIOS2

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Create the HDF5 Engine
        io.SetEngine("HDF5");

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

        adios2::Engine h5Writer = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // EXPECT_EQ(h5Writer.CurrentStep(), step);

            h5Writer.BeginStep();
            h5Writer.Put<int8_t>("i8", currentTestData.I8.data());
            h5Writer.Put<int16_t>("i16", currentTestData.I16.data());
            h5Writer.Put<int32_t>("i32", currentTestData.I32.data());
            h5Writer.Put<int64_t>("i64", currentTestData.I64.data());
            h5Writer.Put<uint8_t>("u8", currentTestData.U8.data());
            h5Writer.Put<uint16_t>("u16", currentTestData.U16.data());
            h5Writer.Put<uint32_t>("u32", currentTestData.U32.data());
            h5Writer.Put<uint64_t>("u64", currentTestData.U64.data());
            h5Writer.Put<float>("r32", currentTestData.R32.data());
            h5Writer.Put<double>("r64", currentTestData.R64.data());
            h5Writer.EndStep();
        }

        h5Writer.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine h5Reader = io.Open(fname, adios2::Mode::Read);

        size_t t = 0;

        while (h5Reader.BeginStep() == adios2::StepStatus::OK)
        {
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

            const size_t currentStep = h5Reader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            h5Reader.Get(var_i8, I8.data());
            h5Reader.Get(var_i16, I16.data());
            h5Reader.Get(var_i32, I32.data());
            h5Reader.Get(var_i64, I64.data());

            h5Reader.Get(var_u8, U8.data());
            h5Reader.Get(var_u16, U16.data());
            h5Reader.Get(var_u32, U32.data());
            h5Reader.Get(var_u64, U64.data());

            h5Reader.Get(var_r32, R32.data());
            h5Reader.Get(var_r64, R64.data());

            h5Reader.PerformGets();
            h5Reader.EndStep();

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

        h5Reader.Close();
    }
}

TEST_F(HDF5WriteReadAsStreamTestADIOS2, ADIOS2HDF5WriteRead2D4x2)
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

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("ADIOS2HDF5WriteReadAsStream2D4x2Test_MPI.h5");
#else
    const std::string fname("ADIOS2HDF5WriteReadAsStream2D4x2Test.h5");
#endif

    // Write test data using ADIOS2

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetEngine("HDF5");

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

        adios2::Engine h5Writer = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            // EXPECT_EQ(h5Writer.CurrentStep(), step);

            h5Writer.BeginStep();
            h5Writer.Put<int8_t>("i8", currentTestData.I8.data());
            h5Writer.Put<int16_t>("i16", currentTestData.I16.data());
            h5Writer.Put<int32_t>("i32", currentTestData.I32.data());
            h5Writer.Put<int64_t>("i64", currentTestData.I64.data());
            h5Writer.Put<uint8_t>("u8", currentTestData.U8.data());
            h5Writer.Put<uint16_t>("u16", currentTestData.U16.data());
            h5Writer.Put<uint32_t>("u32", currentTestData.U32.data());
            h5Writer.Put<uint64_t>("u64", currentTestData.U64.data());
            h5Writer.Put<float>("r32", currentTestData.R32.data());
            h5Writer.Put<double>("r64", currentTestData.R64.data());
            h5Writer.EndStep();
        }

        h5Writer.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine h5Reader = io.Open(fname, adios2::Mode::Read);

        size_t t = 0;

        while (h5Reader.BeginStep() == adios2::StepStatus::OK)
        {

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

            const size_t currentStep = h5Reader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            h5Reader.Get(var_i8, I8.data());
            h5Reader.Get(var_i16, I16.data());
            h5Reader.Get(var_i32, I32.data());
            h5Reader.Get(var_i64, I64.data());

            h5Reader.Get(var_u8, U8.data());
            h5Reader.Get(var_u16, U16.data());
            h5Reader.Get(var_u32, U32.data());
            h5Reader.Get(var_u64, U64.data());

            h5Reader.Get(var_r32, R32.data());
            h5Reader.Get(var_r64, R64.data());

            h5Reader.PerformGets();

            h5Reader.EndStep();

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
        h5Reader.Close();
    }
}

TEST_F(HDF5WriteReadAsStreamTestADIOS2, ReaderWriterDefineVariable)
{
    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 2;
    // Number of cols
    const std::size_t Ny = 4;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fnameFloat("HDF5ReaderWriterDefineVariable_float.h5");
    const std::string fname("HDF5ReaderWriterDefineVariable_all.h5");
#else
    const std::string fnameFloat("HDF5ReaderWriterDefineVariable_float.h5");
    const std::string fname("HDF5ReaderWriterDefineVariable_all.h5");
#endif

    // Write test data using ADIOS2

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
    const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
    const adios2::Dims count{Ny, Nx};
    // simple writer to generate content
    {
        adios2::IO io = adios.DeclareIO("Writer");
        io.SetEngine("HDF5");

        io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);

        adios2::Engine h5Writer = io.Open(fnameFloat, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);
            h5Writer.BeginStep();
            h5Writer.Put<float>("r32", currentTestData.R32.data());
            h5Writer.EndStep();
        }

        h5Writer.Close();
    }
#ifdef TEST_HDF5_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    {
        adios2::IO io = adios.DeclareIO("ReaderWriter");
        io.SetEngine("HDF5");

        adios2::Engine reader = io.Open(fnameFloat, adios2::Mode::Read);
        adios2::Engine writer = io.Open(fname, adios2::Mode::Write);
        for (size_t step = 0; step < NSteps; ++step)
        {
            reader.BeginStep();
            adios2::Variable<float> varR32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(varR32);
            reader.EndStep();

            if (step == 0)
            {
                adios2::Variable<double> varR64 =
                    io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
                EXPECT_TRUE(varR64);
            }
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);
            writer.BeginStep();
            writer.Put<float>("r32", currentTestData.R32.data());
            writer.Put<double>("r64", currentTestData.R64.data());
            writer.EndStep();
        }

        writer.Close();
        reader.Close();
    }
#ifdef TEST_HDF5_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    {
        adios2::IO io = adios.DeclareIO("Reader");
        io.SetEngine("HDF5");
        adios2::Engine reader = io.Open(fname, adios2::Mode::Read);
        while (reader.BeginStep() == adios2::StepStatus::OK)
        {
            adios2::Variable<float> varR32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(varR32);
            adios2::Variable<double> varR64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(varR64);
            reader.EndStep();
        }
        reader.Close();
    }
}

int main(int argc, char **argv)
{
#ifdef TEST_HDF5_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

#ifdef TEST_HDF5_MPI
    MPI_Finalize();
#endif

    return result;
}
