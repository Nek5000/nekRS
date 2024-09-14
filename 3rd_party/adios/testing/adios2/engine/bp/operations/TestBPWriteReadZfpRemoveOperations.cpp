/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <algorithm>
#include <iostream>
#include <numeric> //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

/// Test suit that only compresses on odd steps to test
/// Variable<T>::RemoveOperations functionality

void ZFPRate1D(const std::string rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx);
    std::vector<double> r64s(Nx);

    // range 0 to 999
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRZFPOdd1D_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWRZFPOdd1D_" + rate + ".bp");
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

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        // add operations
        adios2::Operator ZFPOp = adios.DefineOperator("ZFPCompressor", adios2::ops::LossyZFP);

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();

            if (step % 2 == 0)
            {
                // remove current operations to avoid compression in even steps
                var_r32.RemoveOperations();
                var_r64.RemoveOperations();

                EXPECT_TRUE(var_r32.Operations().empty());
                EXPECT_TRUE(var_r64.Operations().empty());
            }
            else
            {
                var_r32.AddOperation(ZFPOp, {{adios2::ops::zfp::key::rate, rate}});
                var_r64.AddOperation(
                    ZFPOp, {{adios2::ops::zfp::key::rate, std::to_string(2 * std::stod(rate))}});

                EXPECT_FALSE(var_r32.Operations().empty());
                EXPECT_FALSE(var_r64.Operations().empty());
            }

            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            auto mmR32 = std::minmax_element(r32s.begin(), r32s.end());
            EXPECT_EQ(var_r32.Min(), *mmR32.first);
            EXPECT_EQ(var_r32.Max(), *mmR32.second);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            auto mmR64 = std::minmax_element(r64s.begin(), r64s.end());
            EXPECT_EQ(var_r64.Min(), *mmR64.first);
            EXPECT_EQ(var_r64.Max(), *mmR64.second);

            const adios2::Dims start{mpiRank * Nx};
            const adios2::Dims count{Nx};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZFPRate2D(const std::string rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;
    const size_t Ny = 50;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx * Ny);
    std::vector<double> r64s(Nx * Ny);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRZFPOdd2D_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWRZFPOdd2D_" + rate + ".bp");
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

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
        const adios2::Dims count{Nx, Ny};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        // add operations
        adios2::Operator zfpOp = adios.DefineOperator("ZFPCompressor", adios2::ops::LossyZFP);

        var_r32.AddOperation(zfpOp, {{"Accuracy", rate}});
        var_r64.AddOperation(zfpOp, {{"accuracy", rate}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();

            if (step % 2 == 0)
            {
                // remove current operations to avoid compression in even steps
                var_r32.RemoveOperations();
                var_r64.RemoveOperations();

                EXPECT_TRUE(var_r32.Operations().empty());
                EXPECT_TRUE(var_r64.Operations().empty());
            }
            else
            {
                var_r32.AddOperation(zfpOp, {{adios2::ops::zfp::key::rate, rate}});
                var_r64.AddOperation(
                    zfpOp, {{adios2::ops::zfp::key::rate, std::to_string(2 * std::stod(rate))}});

                EXPECT_FALSE(var_r32.Operations().empty());
                EXPECT_FALSE(var_r64.Operations().empty());
            }

            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
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

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
        ASSERT_EQ(var_r32.Shape()[1], Ny);
        auto mmR32 = std::minmax_element(r32s.begin(), r32s.end());
        EXPECT_EQ(var_r32.Min(), *mmR32.first);
        EXPECT_EQ(var_r32.Max(), *mmR32.second);

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
        ASSERT_EQ(var_r64.Shape()[1], Ny);
        auto mmR64 = std::minmax_element(r64s.begin(), r64s.end());
        EXPECT_EQ(var_r64.Min(), *mmR64.first);
        EXPECT_EQ(var_r64.Max(), *mmR64.second);

        const adios2::Dims start{mpiRank * Nx, 0};
        const adios2::Dims count{Nx, Ny};
        const adios2::Box<adios2::Dims> sel(start, count);
        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZFPRate3D(const std::string rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 10;
    const size_t Ny = 20;
    const size_t Nz = 15;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx * Ny * Nz);
    std::vector<double> r64s(Nx * Ny * Nz);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRZFPOdd3D_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWRZFPOdd3D_" + rate + ".bp");
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

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny, Nz};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0, 0};
        const adios2::Dims count{Nx, Ny, Nz};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        // add operations
        adios2::Operator zfpOp = adios.DefineOperator("ZFPCompressor", adios2::ops::LossyZFP);

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();

            if (step % 2 == 0)
            {
                // remove current operations to avoid compression in even steps
                var_r32.RemoveOperations();
                var_r64.RemoveOperations();

                EXPECT_TRUE(var_r32.Operations().empty());
                EXPECT_TRUE(var_r64.Operations().empty());
            }
            else
            {
                var_r32.AddOperation(zfpOp, {{adios2::ops::zfp::key::rate, rate}});
                var_r64.AddOperation(
                    zfpOp, {{adios2::ops::zfp::key::rate, std::to_string(2 * std::stod(rate))}});

                EXPECT_FALSE(var_r32.Operations().empty());
                EXPECT_FALSE(var_r64.Operations().empty());
            }

            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
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

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
        ASSERT_EQ(var_r32.Shape()[1], Ny);
        ASSERT_EQ(var_r32.Shape()[2], Nz);
        auto mmR32 = std::minmax_element(r32s.begin(), r32s.end());
        EXPECT_EQ(var_r32.Min(), *mmR32.first);
        EXPECT_EQ(var_r32.Max(), *mmR32.second);

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
        ASSERT_EQ(var_r64.Shape()[1], Ny);
        ASSERT_EQ(var_r64.Shape()[2], Nz);
        auto mmR64 = std::minmax_element(r64s.begin(), r64s.end());
        EXPECT_EQ(var_r64.Min(), *mmR64.first);
        EXPECT_EQ(var_r64.Max(), *mmR64.second);

        const adios2::Dims start{mpiRank * Nx, 0, 0};
        const adios2::Dims count{Nx, Ny, Nz};
        const adios2::Box<adios2::Dims> sel(start, count);
        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny * Nz; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

class BPWRZFPODD : public ::testing::TestWithParam<std::string>
{
public:
    BPWRZFPODD() = default;

    virtual void SetUp() {}
    virtual void TearDown() {}
};

TEST_P(BPWRZFPODD, ADIOS2BPWRZFPODD1D) { ZFPRate1D(GetParam()); }
TEST_P(BPWRZFPODD, ADIOS2BPWRZFPODD2D) { ZFPRate2D(GetParam()); }
TEST_P(BPWRZFPODD, ADIOS2BPWRZFPODD3D) { ZFPRate3D(GetParam()); }

INSTANTIATE_TEST_SUITE_P(ZFPRate, BPWRZFPODD, ::testing::Values("8", "9", "10"));

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
