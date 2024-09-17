/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "CudaRoutines.h"

#include <adios2.h>

#include <algorithm>
#include <array>
#include <gtest/gtest.h>
#include <iostream>
#include <numeric> //std::iota

std::string engineName; // comes from command line

const float EPSILON = std::numeric_limits<float>::epsilon();
const float INCREMENT = 10.0f;

void ZFPRateCUDA(const std::string rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    // Number of rows
    const size_t Nx = 100;

    // Number of steps
    const size_t NSteps = 1;

    int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRZFP1D_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWRZFP1D_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    const size_t NxTotal = Nx * mpiSize;

    // Initialize the simulation data
    std::vector<float> r32s(NxTotal, .0f);
    std::iota(r32s.begin(), r32s.end(), .0f);

    float *gpuSimData = nullptr;
    cudaMalloc(&gpuSimData, Nx * sizeof(float));
    cudaMemcpy(gpuSimData, ((float *)&r32s[0] + (Nx * mpiRank)), Nx * sizeof(float),
               cudaMemcpyHostToDevice);

    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(NxTotal)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);

        // add operations
        adios2::Operator ZFPOp = adios.DefineOperator("ZFPCompressor", adios2::ops::LossyZFP);

        var_r32.AddOperation(ZFPOp, {{adios2::ops::zfp::key::rate, rate}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Update values in the simulation data
            cuda_increment(Nx, 1, 0, gpuSimData, INCREMENT);

            bpWriter.BeginStep();
            var_r32.SetMemorySpace(adios2::MemorySpace::GPU);
            bpWriter.Put(var_r32, gpuSimData);
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

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], NxTotal);

        auto mmR32 = std::minmax_element(r32s.begin(), r32s.end());
        EXPECT_EQ(var_r32.Min() - INCREMENT, *mmR32.first);
        EXPECT_EQ(var_r32.Max() - INCREMENT, *mmR32.second);

        unsigned int t = 0;
        for (; bpReader.BeginStep() == adios2::StepStatus::OK; ++t)
        {
            std::vector<float> r32o(NxTotal);
            bpReader.Get(var_r32, r32o);
            bpReader.EndStep();

            // Remove INCREMENT from each element
            std::transform(r32o.begin(), r32o.end(), r32o.begin(),
                           std::bind(std::minus<float>(), std::placeholders::_1, INCREMENT));

            for (size_t i = 0; i < NxTotal; i++)
            {
                char msg[1 << 8] = {0};
                snprintf(msg, sizeof(msg), "t=%d i=%zu rank=%d r32o=%f r32s=%f", t, i, mpiRank,
                         r32o[i], r32s[i]);
                ASSERT_LT(std::abs(r32o[i] - r32s[i]), EPSILON) << msg;
            }
        }
        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

class BPWRZFPCUDA : public ::testing::TestWithParam<std::string>
{
public:
    BPWRZFPCUDA() = default;

    virtual void SetUp() {}
    virtual void TearDown() {}
};

TEST_P(BPWRZFPCUDA, ADIOS2BPWRZFPCUDA) { ZFPRateCUDA(GetParam()); }

INSTANTIATE_TEST_SUITE_P(ZFPRate, BPWRZFPCUDA, ::testing::Values("16", "32", "64"));

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
