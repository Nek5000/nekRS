/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <cstdint>
#include <cstring>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include <adios2.h>
#include <cuda_runtime.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

void MGARDAccuracy2D(const std::string tolerance)
{

    int mpiRank = 0, mpiSize = 1;
    const size_t Nx = 100;
    const size_t Ny = 50;
    const size_t NSteps = 1;

    std::vector<double> r64s(Nx * Ny);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRMGARDCU2D_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARDCU2D_" + tolerance + ".bp");
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

        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        // add MGARD operations
        adios2::Operator mgardOp = adios.DefineOperator("mgardCompressor", adios2::ops::LossyMGARD);
        var_r64.AddOperation(mgardOp, {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        double *gpu64s = nullptr;
        cudaMalloc(&gpu64s, Nx * Ny * sizeof(double));
        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            cudaMemcpy(gpu64s, r64s.data(), Nx * Ny * sizeof(double), cudaMemcpyHostToDevice);
            var_r64.SetMemorySpace(adios2::MemorySpace::GPU);
            bpWriter.Put<double>("r64", gpu64s);
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

        unsigned int t = 0;
        std::vector<double> decompressedR64s(Nx * Ny);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);

            const adios2::Dims start{mpiRank * Nx, 0};
            const adios2::Dims count{Nx, Ny};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r64.SetSelection(sel);

            double *gpu64s = nullptr;
            cudaMalloc(&gpu64s, Nx * Ny * sizeof(double));
            bpReader.Get(var_r64, gpu64s);
            bpReader.EndStep();
            cudaMemcpy(decompressedR64s.data(), gpu64s, Nx * Ny * sizeof(double),
                       cudaMemcpyDeviceToHost);

            double maxDiff = 0;
            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                double diff = std::abs(r64s[i] - decompressedR64s[i]);
                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }
            ++t;

            auto itMax = std::max_element(r64s.begin(), r64s.end());
            const double relativeMaxDiff = maxDiff / *itMax;
            ASSERT_LT(relativeMaxDiff, std::stod(tolerance));
            std::cout << "Relative Max Diff " << relativeMaxDiff << " tolerance " << tolerance
                      << "\n";
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void MGARDAccuracySmall(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx);
    std::iota(r32s.begin(), r32s.end(), 0.f);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRMGARD1D_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD1D_" + tolerance + ".bp");
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

        // add operations
        adios2::Operator mgardOp = adios.DefineOperator("mgardCompressor", adios2::ops::LossyMGARD);

        var_r32.AddOperation(mgardOp, {{adios2::ops::mgard::key::tolerance, tolerance},
                                       {adios2::ops::mgard::key::s, "inf"}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        float *gpu32s = nullptr;
        cudaMalloc(&gpu32s, Nx * sizeof(float));
        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            cudaMemcpy(gpu32s, r32s.data(), Nx * sizeof(float), cudaMemcpyHostToDevice);
            bpWriter.Put<float>("r32", gpu32s);
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

        unsigned int t = 0;
        std::vector<float> decompressedR32s(Nx);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

            const adios2::Dims start{mpiRank * Nx};
            const adios2::Dims count{Nx};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);

            float *gpu32s = nullptr;
            cudaMalloc(&gpu32s, Nx * sizeof(float));
            bpReader.Get(var_r32, gpu32s);
            bpReader.EndStep();
            cudaMemcpy(decompressedR32s.data(), gpu32s, Nx * sizeof(float), cudaMemcpyDeviceToHost);

            double maxDiff = 0, relativeMaxDiff = 0;

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                double diff = std::abs(r32s[i] - decompressedR32s[i]);

                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }
            ++t;

            auto r32s_Max = std::max_element(r32s.begin(), r32s.end());
            relativeMaxDiff = maxDiff / *r32s_Max;

            ASSERT_LT(relativeMaxDiff, std::stod(tolerance));
            std::cout << "Relative Max Diff " << relativeMaxDiff << " tolerance " << tolerance
                      << "\n";
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

class BPWriteReadMGARD : public ::testing::TestWithParam<std::string>
{
public:
    BPWriteReadMGARD() = default;
    virtual void SetUp(){};
    virtual void TearDown(){};
};

TEST_P(BPWriteReadMGARD, BPWRMGARDCU2D) { MGARDAccuracy2D(GetParam()); }
TEST_P(BPWriteReadMGARD, BPWRMGARDCU1D) { MGARDAccuracySmall(GetParam()); }

INSTANTIATE_TEST_SUITE_P(MGARDAccuracy, BPWriteReadMGARD,
                         ::testing::Values("0.01", "0.001", "0.0001", "0.00001"));

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
