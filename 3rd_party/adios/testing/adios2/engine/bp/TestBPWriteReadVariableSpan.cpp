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

class BPWriteReadSpan : public ::testing::Test
{
public:
    BPWriteReadSpan() = default;

    SmallTestData m_TestData;
};

TEST_F(BPWriteReadSpan, BPWriteRead1D8)
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
    const std::string fname("BPWriteReadSpan1D8_MPI.bp");
#else
    const std::string fname("BPWriteReadSpan1D8.bp");
#endif

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
        else
        {
            io.SetEngine("BP3");
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_Step = io.DefineVariable<size_t>("step");
        /* Why is there no Span for string variable? */
        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);
        auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);
        auto var_u16 =
            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
        auto var_u32 =
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
        auto var_u64 =
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                               adios2::ConstantDims);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                                adios2::ConstantDims);

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

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            bpWriter.BeginStep();

            EXPECT_EQ(bpWriter.CurrentStep(), step);

            bpWriter.Put(var_Step, step);
            // bpWriter.Put<std::string>("iString", currentTestData.S1);

            adios2::Variable<int8_t>::Span i8Span = bpWriter.Put(var_i8);
            adios2::Variable<int16_t>::Span i16Span = bpWriter.Put(var_i16);
            adios2::Variable<int32_t>::Span i32Span = bpWriter.Put(var_i32);
            adios2::Variable<int64_t>::Span i64Span = bpWriter.Put(var_i64);
            adios2::Variable<uint8_t>::Span u8Span = bpWriter.Put(var_u8);
            adios2::Variable<uint16_t>::Span u16Span = bpWriter.Put(var_u16);
            adios2::Variable<uint32_t>::Span u32Span = bpWriter.Put(var_u32);
            adios2::Variable<uint64_t>::Span u64Span = bpWriter.Put(var_u64);
            adios2::Variable<float>::Span r32Span = bpWriter.Put(var_r32);
            adios2::Variable<double>::Span r64Span = bpWriter.Put(var_r64);
            adios2::Variable<std::complex<float>>::Span cr32Span = bpWriter.Put(var_cr32);
            adios2::Variable<std::complex<double>>::Span cr64Span = bpWriter.Put(var_cr64);

            auto ptr = i64Span.data();
            (void)ptr;

            // Testing Data()
            std::copy(currentTestData.I8.begin(), currentTestData.I8.begin() + Nx, i8Span.begin());
            std::copy(currentTestData.I16.begin(), currentTestData.I16.begin() + Nx,
                      i16Span.begin());
            std::copy(currentTestData.I32.begin(), currentTestData.I32.begin() + Nx,
                      i32Span.data());
            std::copy(currentTestData.I64.begin(), currentTestData.I64.begin() + Nx,
                      i64Span.data());
            // Testing operator[] and At
            for (size_t i = 0; i < Nx; ++i)
            {
                u8Span[i] = currentTestData.U8[i];
                u16Span[i] = currentTestData.U16[i];
                u32Span[i] = currentTestData.U32[i];
                u64Span[i] = currentTestData.U64[i];

                r32Span.at(i) = currentTestData.R32[i];
                r64Span.at(i) = currentTestData.R64[i];
                cr32Span.at(i) = currentTestData.CR32[i];
                cr64Span.at(i) = currentTestData.CR64[i];
            }

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
        else
        {
            io.SetEngine("BP3");
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        size_t IStep;
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

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        size_t t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            auto var_iStep = io.InquireVariable<size_t>("step");
            // auto var_iString = io.InquireVariable<std::string>("iString");
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

            EXPECT_EQ(var_iStep.ShapeID(), adios2::ShapeID::GlobalValue);

            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i8.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i16.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i32.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i64.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u8.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u16.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u32.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u64.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r32.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r64.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr32.Shape()[0], static_cast<size_t>(mpiSize * Nx));

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr64.Shape()[0], static_cast<size_t>(mpiSize * Nx));

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
            var_cr32.SetSelection(sel);
            var_cr64.SetSelection(sel);

            bpReader.Get(var_iStep, IStep);
            // bpReader.Get(var_iString, IString);
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

            bpReader.EndStep();

            EXPECT_EQ(IStep, currentStep);
            // EXPECT_EQ(IString, currentTestData.S1);

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

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadSpan, BPWriteRead2D2x4)
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
    const std::string fname("BPWriteReadSpan2D2x4_MPI.bp");
#else
    const std::string fname("BPWriteReadSpan2D2x4.bp");
#endif

    // Write test data using ADIOS2

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
        else
        {
            io.SetEngine("BP3");
        }

        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        const adios2::Dims shape{Ny, static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);
        auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);
        auto var_u16 =
            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
        auto var_u32 =
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
        auto var_u64 =
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                               adios2::ConstantDims);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                                adios2::ConstantDims);

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

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            bpWriter.BeginStep();
            EXPECT_EQ(bpWriter.CurrentStep(), step);
            adios2::Variable<int8_t>::Span i8Span = bpWriter.Put(var_i8);
            adios2::Variable<int16_t>::Span i16Span = bpWriter.Put(var_i16);
            adios2::Variable<int32_t>::Span i32Span = bpWriter.Put(var_i32);
            adios2::Variable<int64_t>::Span i64Span = bpWriter.Put(var_i64);
            adios2::Variable<uint8_t>::Span u8Span = bpWriter.Put(var_u8);
            adios2::Variable<uint16_t>::Span u16Span = bpWriter.Put(var_u16);
            adios2::Variable<uint32_t>::Span u32Span = bpWriter.Put(var_u32);
            adios2::Variable<uint64_t>::Span u64Span = bpWriter.Put(var_u64);
            adios2::Variable<float>::Span r32Span = bpWriter.Put(var_r32);
            adios2::Variable<double>::Span r64Span = bpWriter.Put(var_r64);
            adios2::Variable<std::complex<float>>::Span cr32Span = bpWriter.Put(var_cr32);
            adios2::Variable<std::complex<double>>::Span cr64Span = bpWriter.Put(var_cr64);

            // Testing Data()
            std::copy(currentTestData.I8.begin(), currentTestData.I8.begin() + Nx * Ny,
                      i8Span.data());
            std::copy(currentTestData.I16.begin(), currentTestData.I16.begin() + Nx * Ny,
                      i16Span.data());

            size_t i = 0;
            for (auto &i32 : i32Span)
            {
                i32 = currentTestData.I32[i];
                ++i;
            }
            i = 0;
            for (auto &i64 : i64Span)
            {
                i64 = currentTestData.I64[i];
                ++i;
            }

            // Testing operator[] and At
            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                u8Span[i] = currentTestData.U8[i];
                u16Span[i] = currentTestData.U16[i];
                u32Span[i] = currentTestData.U32[i];
                u64Span[i] = currentTestData.U64[i];

                r32Span.at(i) = currentTestData.R32[i];
                r64Span.at(i) = currentTestData.R64[i];
                cr32Span.at(i) = currentTestData.CR32[i];
                cr64Span.at(i) = currentTestData.CR64[i];
            }

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
        else
        {
            io.SetEngine("BP3");
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        size_t t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i8.Shape()[0], Ny);
            EXPECT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i16.Shape()[0], Ny);
            EXPECT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i32.Shape()[0], Ny);
            EXPECT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i64.Shape()[0], Ny);
            EXPECT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            EXPECT_TRUE(var_u8);
            EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u8.Shape()[0], Ny);
            EXPECT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            EXPECT_TRUE(var_u16);
            EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u16.Shape()[0], Ny);
            EXPECT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            EXPECT_TRUE(var_u32);
            EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u32.Shape()[0], Ny);
            EXPECT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            EXPECT_TRUE(var_u64);
            EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u64.Shape()[0], Ny);
            EXPECT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r32.Shape()[0], Ny);
            EXPECT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r64.Shape()[0], Ny);
            EXPECT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            EXPECT_TRUE(var_cr32);
            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr32.Shape()[0], Ny);
            EXPECT_EQ(var_cr32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
            EXPECT_TRUE(var_cr64);
            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr64.Shape()[0], Ny);
            EXPECT_EQ(var_cr64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

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
            var_cr32.SetSelection(sel);
            var_cr64.SetSelection(sel);

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

            bpReader.Get(var_cr32, CR32.data());
            bpReader.Get(var_cr64, CR64.data());

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
                EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
            }
            ++t;
        }
        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadSpan, BPWriteRead1D8Local)
{

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadSpan1D8Local_MPI.bp");
#else
    const std::string fname("BPWriteReadSpan1D8Local.bp");
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
        else
        {
            io.SetEngine("BP3");
        }

        const adios2::Dims shape{};
        const adios2::Dims start{};
        const adios2::Dims count{Nx};

        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);
        auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);
        auto var_u16 =
            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
        auto var_u32 =
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
        auto var_u64 =
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                               adios2::ConstantDims);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                                adios2::ConstantDims);

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

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            bpWriter.BeginStep();
            EXPECT_EQ(bpWriter.CurrentStep(), step);
            // bpWriter.Put<std::string>("iString", currentTestData.S1);

            adios2::Variable<int8_t>::Span i8Span = bpWriter.Put(var_i8);
            adios2::Variable<int16_t>::Span i16Span = bpWriter.Put(var_i16);
            adios2::Variable<int32_t>::Span i32Span = bpWriter.Put(var_i32);
            adios2::Variable<int64_t>::Span i64Span = bpWriter.Put(var_i64);
            adios2::Variable<uint8_t>::Span u8Span = bpWriter.Put(var_u8);
            adios2::Variable<uint16_t>::Span u16Span = bpWriter.Put(var_u16);
            adios2::Variable<uint32_t>::Span u32Span = bpWriter.Put(var_u32);
            adios2::Variable<uint64_t>::Span u64Span = bpWriter.Put(var_u64);
            adios2::Variable<float>::Span r32Span = bpWriter.Put(var_r32);
            adios2::Variable<double>::Span r64Span = bpWriter.Put(var_r64);
            adios2::Variable<std::complex<float>>::Span cr32Span = bpWriter.Put(var_cr32);
            adios2::Variable<std::complex<double>>::Span cr64Span = bpWriter.Put(var_cr64);

            // Testing Data()
            std::copy(currentTestData.I8.begin(), currentTestData.I8.begin() + Nx, i8Span.begin());
            std::copy(currentTestData.I16.begin(), currentTestData.I16.begin() + Nx,
                      i16Span.begin());
            std::copy(currentTestData.I32.begin(), currentTestData.I32.begin() + Nx,
                      i32Span.data());
            std::copy(currentTestData.I64.begin(), currentTestData.I64.begin() + Nx,
                      i64Span.data());
            // Testing operator[] and At
            for (size_t i = 0; i < Nx; ++i)
            {
                u8Span[i] = currentTestData.U8[i];
                u16Span[i] = currentTestData.U16[i];
                u32Span[i] = currentTestData.U32[i];
                u64Span[i] = currentTestData.U64[i];

                r32Span.at(i) = currentTestData.R32[i];
                r64Span.at(i) = currentTestData.R64[i];
                cr32Span.at(i) = currentTestData.CR32[i];
                cr64Span.at(i) = currentTestData.CR64[i];
            }

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
        else
        {
            io.SetEngine("BP3");
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

            SmallTestData currentTestData = generateNewSmallTestData(
                m_TestData, static_cast<int>(currentStep), mpiRank, mpiSize);

            // auto var_iString = io.InquireVariable<std::string>("iString");
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

            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);
            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);

            const size_t rankBlock = static_cast<size_t>(mpiRank);
            var_i8.SetBlockSelection(rankBlock);
            var_i16.SetBlockSelection(rankBlock);
            var_i32.SetBlockSelection(rankBlock);
            var_i64.SetBlockSelection(rankBlock);
            var_u8.SetBlockSelection(rankBlock);
            var_u16.SetBlockSelection(rankBlock);
            var_u32.SetBlockSelection(rankBlock);
            var_u64.SetBlockSelection(rankBlock);
            var_r32.SetBlockSelection(rankBlock);
            var_r64.SetBlockSelection(rankBlock);
            var_cr32.SetBlockSelection(rankBlock);
            var_cr64.SetBlockSelection(rankBlock);

            // bpReader.Get(var_iString, IString);
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

            bpReader.EndStep();

            // EXPECT_EQ(IString, currentTestData.S1);

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

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadSpan, BPWriteRead2D2x4Local)
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
    const std::string fname("BPWriteReadSpan2D2x4Local_MPI.bp");
#else
    const std::string fname("BPWriteReadSpan2D2x4Local.bp");
#endif

    // Write test data using ADIOS2

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
        else
        {
            io.SetEngine("BP3");
        }

        const adios2::Dims shape{};
        const adios2::Dims start{};
        const adios2::Dims count{Ny, Nx};

        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);
        auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);
        auto var_u16 =
            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
        auto var_u32 =
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
        auto var_u64 =
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                               adios2::ConstantDims);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                                adios2::ConstantDims);

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

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            bpWriter.BeginStep();
            EXPECT_EQ(bpWriter.CurrentStep(), step);
            adios2::Variable<int8_t>::Span i8Span = bpWriter.Put(var_i8);
            adios2::Variable<int16_t>::Span i16Span = bpWriter.Put(var_i16);
            adios2::Variable<int32_t>::Span i32Span = bpWriter.Put(var_i32);
            adios2::Variable<int64_t>::Span i64Span = bpWriter.Put(var_i64);
            adios2::Variable<uint8_t>::Span u8Span = bpWriter.Put(var_u8);
            adios2::Variable<uint16_t>::Span u16Span = bpWriter.Put(var_u16);
            adios2::Variable<uint32_t>::Span u32Span = bpWriter.Put(var_u32);
            adios2::Variable<uint64_t>::Span u64Span = bpWriter.Put(var_u64);
            adios2::Variable<float>::Span r32Span = bpWriter.Put(var_r32);
            adios2::Variable<double>::Span r64Span = bpWriter.Put(var_r64);
            adios2::Variable<std::complex<float>>::Span cr32Span = bpWriter.Put(var_cr32);
            adios2::Variable<std::complex<double>>::Span cr64Span = bpWriter.Put(var_cr64);

            // Testing Data()
            std::copy(currentTestData.I8.begin(), currentTestData.I8.begin() + Nx * Ny,
                      i8Span.data());
            std::copy(currentTestData.I16.begin(), currentTestData.I16.begin() + Nx * Ny,
                      i16Span.data());

            size_t i = 0;
            for (auto &i32 : i32Span)
            {
                i32 = currentTestData.I32[i];
                ++i;
            }
            i = 0;
            for (auto it = i64Span.begin(); it != i64Span.end(); ++it)
            {
                *it = currentTestData.I64[i];
                ++i;
            }

            // Testing operator[] and At
            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                u8Span[i] = currentTestData.U8[i];
                u16Span[i] = currentTestData.U16[i];
                u32Span[i] = currentTestData.U32[i];
                u64Span[i] = currentTestData.U64[i];

                r32Span.at(i) = currentTestData.R32[i];
                r64Span.at(i) = currentTestData.R64[i];
                cr32Span.at(i) = currentTestData.CR32[i];
                cr64Span.at(i) = currentTestData.CR64[i];
            }

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
        else
        {
            io.SetEngine("BP3");
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        size_t t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            EXPECT_TRUE(var_u8);
            EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            EXPECT_TRUE(var_u16);
            EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            EXPECT_TRUE(var_u32);
            EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            EXPECT_TRUE(var_u64);
            EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            EXPECT_TRUE(var_cr32);
            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::LocalArray);

            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
            EXPECT_TRUE(var_cr64);
            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::LocalArray);

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

            const size_t rankBlock = static_cast<size_t>(mpiRank);
            var_i8.SetBlockSelection(rankBlock);
            var_i16.SetBlockSelection(rankBlock);
            var_i32.SetBlockSelection(rankBlock);
            var_i64.SetBlockSelection(rankBlock);
            var_u8.SetBlockSelection(rankBlock);
            var_u16.SetBlockSelection(rankBlock);
            var_u32.SetBlockSelection(rankBlock);
            var_u64.SetBlockSelection(rankBlock);
            var_r32.SetBlockSelection(rankBlock);
            var_r64.SetBlockSelection(rankBlock);
            var_cr32.SetBlockSelection(rankBlock);
            var_cr64.SetBlockSelection(rankBlock);

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
            bpReader.Get(var_cr32, CR32.data());
            bpReader.Get(var_cr64, CR64.data());

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
                EXPECT_EQ(CR32[i], currentTestData.CR32[i]) << msg;
                EXPECT_EQ(CR64[i], currentTestData.CR64[i]) << msg;
            }
            ++t;
        }
        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

TEST_F(BPWriteReadSpan, BPWriteRead1D8FillValue)
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
    const std::string fname("BPWriteReadSpan1D8FillValue_MPI.bp");
#else
    const std::string fname("BPWriteReadSpan1D8FillValue.bp");
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
        else
        {
            io.SetEngine("BP3");
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count, adios2::ConstantDims);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count, adios2::ConstantDims);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count, adios2::ConstantDims);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count, adios2::ConstantDims);
        auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count, adios2::ConstantDims);
        auto var_u16 =
            io.DefineVariable<uint16_t>("u16", shape, start, count, adios2::ConstantDims);
        auto var_u32 =
            io.DefineVariable<uint32_t>("u32", shape, start, count, adios2::ConstantDims);
        auto var_u64 =
            io.DefineVariable<uint64_t>("u64", shape, start, count, adios2::ConstantDims);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count,
                                                               adios2::ConstantDims);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count,
                                                                adios2::ConstantDims);

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

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            EXPECT_EQ(bpWriter.CurrentStep(), step);
            // bpWriter.Put<std::string>("iString", std::to_string(step));

            adios2::Variable<int8_t>::Span i8Span =
                bpWriter.Put(var_i8, true, static_cast<int8_t>(step));
            adios2::Variable<int16_t>::Span i16Span =
                bpWriter.Put(var_i16, true, static_cast<int16_t>(step));
            adios2::Variable<int32_t>::Span i32Span =
                bpWriter.Put(var_i32, true, static_cast<int32_t>(step));
            adios2::Variable<int64_t>::Span i64Span =
                bpWriter.Put(var_i64, true, static_cast<int64_t>(step));

            adios2::Variable<uint8_t>::Span u8Span =
                bpWriter.Put(var_u8, true, static_cast<uint8_t>(step));
            adios2::Variable<uint16_t>::Span u16Span =
                bpWriter.Put(var_u16, true, static_cast<uint16_t>(step));
            adios2::Variable<uint32_t>::Span u32Span =
                bpWriter.Put(var_u32, true, static_cast<uint32_t>(step));
            adios2::Variable<uint64_t>::Span u64Span =
                bpWriter.Put(var_u64, true, static_cast<uint64_t>(step));

            adios2::Variable<float>::Span r32Span =
                bpWriter.Put(var_r32, true, static_cast<float>(step));

            adios2::Variable<double>::Span r64Span =
                bpWriter.Put(var_r64, true, static_cast<double>(step));

            adios2::Variable<std::complex<float>>::Span cr32Span =
                bpWriter.Put(var_cr32, true, {static_cast<float>(step), static_cast<float>(step)});
            adios2::Variable<std::complex<double>>::Span cr64Span = bpWriter.Put(
                var_cr64, true, {static_cast<double>(step), static_cast<double>(step)});

            (void)i8Span;
            (void)i16Span;
            (void)i32Span;
            (void)i64Span;
            (void)u8Span;
            (void)u16Span;
            (void)u32Span;
            (void)u64Span;
            (void)r32Span;
            (void)r64Span;
            (void)cr32Span;
            (void)cr64Span;

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
        else
        {
            io.SetEngine("BP3");
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

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        size_t t = 0;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t currentStep = bpReader.CurrentStep();
            EXPECT_EQ(currentStep, static_cast<size_t>(t));

            // auto var_iString = io.InquireVariable<std::string>("iString");
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

            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i8.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_i8.Min(), static_cast<int8_t>(currentStep));
            // EXPECT_EQ(var_i8.Max(), static_cast<int8_t>(currentStep));

            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i16.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_i16.Min(), static_cast<int16_t>(currentStep));
            // EXPECT_EQ(var_i16.Max(), static_cast<int16_t>(currentStep));

            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i32.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_i32.Min(), static_cast<int32_t>(currentStep));
            // EXPECT_EQ(var_i32.Max(), static_cast<int32_t>(currentStep));

            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i64.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_i64.Min(), static_cast<int64_t>(currentStep));
            // EXPECT_EQ(var_i64.Max(), static_cast<int64_t>(currentStep));

            EXPECT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u8.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_u8.Min(), static_cast<uint8_t>(currentStep));
            // EXPECT_EQ(var_u8.Max(), static_cast<uint8_t>(currentStep));

            EXPECT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u16.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_u16.Min(), static_cast<uint16_t>(currentStep));
            // EXPECT_EQ(var_u16.Max(), static_cast<uint16_t>(currentStep));

            EXPECT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u32.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_u32.Min(), static_cast<uint32_t>(currentStep));
            // EXPECT_EQ(var_u32.Max(), static_cast<uint32_t>(currentStep));

            EXPECT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_u64.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_u64.Min(), static_cast<uint64_t>(currentStep));
            // EXPECT_EQ(var_u64.Max(), static_cast<uint64_t>(currentStep));

            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r32.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_r32.Min(), static_cast<float>(currentStep));
            // EXPECT_EQ(var_r32.Max(), static_cast<float>(currentStep));

            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r64.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_r64.Min(), static_cast<double>(currentStep));
            // EXPECT_EQ(var_r64.Max(), static_cast<double>(currentStep));

            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr32.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_cr32.Min(),
            //          std::complex<float>(static_cast<float>(currentStep),
            //                              static_cast<float>(currentStep)));
            // EXPECT_EQ(var_cr32.Max(),
            //          std::complex<float>(static_cast<float>(currentStep),
            //                              static_cast<float>(currentStep)));

            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr64.Shape()[0], static_cast<size_t>(mpiSize * Nx));
            // EXPECT_EQ(var_cr64.Min(),
            //          std::complex<double>(static_cast<double>(currentStep),
            //                               static_cast<double>(currentStep)));
            // EXPECT_EQ(var_cr64.Max(),
            //          std::complex<double>(static_cast<double>(currentStep),
            //                               static_cast<double>(currentStep)));

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
            var_cr32.SetSelection(sel);
            var_cr64.SetSelection(sel);

            // bpReader.Get(var_iString, IString);
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

            bpReader.EndStep();

            // EXPECT_EQ(IString, std::to_string(currentStep));

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], static_cast<int8_t>(currentStep)) << msg;
                EXPECT_EQ(I16[i], static_cast<int16_t>(currentStep)) << msg;
                EXPECT_EQ(I32[i], static_cast<int32_t>(currentStep)) << msg;
                EXPECT_EQ(I64[i], static_cast<int64_t>(currentStep)) << msg;
                EXPECT_EQ(U8[i], static_cast<uint8_t>(currentStep)) << msg;
                EXPECT_EQ(U16[i], static_cast<uint16_t>(currentStep)) << msg;
                EXPECT_EQ(U32[i], static_cast<uint32_t>(currentStep)) << msg;
                EXPECT_EQ(U64[i], static_cast<uint64_t>(currentStep)) << msg;
                EXPECT_EQ(R32[i], static_cast<float>(currentStep)) << msg;
                EXPECT_EQ(R64[i], static_cast<double>(currentStep)) << msg;
                EXPECT_EQ(CR32[i], std::complex<float>(static_cast<float>(currentStep),
                                                       static_cast<float>(currentStep)))
                    << msg;
                EXPECT_EQ(CR64[i], std::complex<double>(static_cast<double>(currentStep),
                                                        static_cast<double>(currentStep)))
                    << msg;
            }

            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

#ifdef ADIOS2_HAVE_BZIP2
TEST_F(BPWriteReadSpan, BPWriteSpanOperatorException)
{

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteSpanOperatorException_MPI.bp");
#else
    const std::string fname("BPWriteSpanOperatorException.bp");
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
        else
        {
            io.SetEngine("BP3");
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        adios2::Operator BZIP2Op =
            adios.DefineOperator("BZIP2Compressor", adios2::ops::LosslessBZIP2);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();

            // using bzip2, it could have been any operator to generate an
            // exception
            var_r32.AddOperation(BZIP2Op, {{adios2::ops::bzip2::key::blockSize100k, "1e-4"}});
            var_r64.AddOperation(BZIP2Op, {{adios2::ops::bzip2::key::blockSize100k, "1e-4"}});

            EXPECT_THROW(bpWriter.Put(var_r32, true, static_cast<float>(step)),
                         std::invalid_argument);

            EXPECT_THROW(bpWriter.Put(var_r64, true, static_cast<double>(step)),
                         std::invalid_argument);

            bpWriter.EndStep();
        }

        bpWriter.Close();
    }
}
#endif

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
