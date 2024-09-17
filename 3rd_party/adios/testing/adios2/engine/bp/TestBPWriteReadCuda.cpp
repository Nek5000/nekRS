/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "operations/CudaRoutines.h"

#include <adios2.h>

#include <algorithm>
#include <array>
#include <gtest/gtest.h>
#include <iostream>
#include <numeric> //std::iota

std::string engineName; // comes from command line

const float EPSILON = std::numeric_limits<float>::epsilon();
const float INCREMENT = 10.0f;

void CUDAWrongMemSpace()
{
    const std::string fname("BPWRCUFail.bp");
    const size_t Nx = 5;

    adios2::ADIOS adios;
    std::vector<float> r32s(Nx, .0f);
    std::iota(r32s.begin(), r32s.end(), .0f);
    { // write
        adios2::IO io = adios.DeclareIO("TestIO");
        const adios2::Dims shape{Nx};
        const adios2::Dims start{0};
        const adios2::Dims count{Nx};
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
        auto var_r32_cpu = io.DefineVariable<float>("r32cpu", shape, start, count);

        float *gpuSimData = nullptr;
        cudaMalloc(&gpuSimData, Nx * sizeof(float));
        cudaMemcpy(gpuSimData, (float *)&r32s[0], Nx * sizeof(float), cudaMemcpyHostToDevice);

        io.SetEngine("BP5");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        bpWriter.BeginStep();
        var_r32.SetMemorySpace(adios2::MemorySpace::Host);
        EXPECT_DEATH(bpWriter.Put(var_r32, gpuSimData), "");
        var_r32_cpu.SetMemorySpace(adios2::MemorySpace::GPU);
        bpWriter.Put(var_r32_cpu, r32s.data());
        bpWriter.EndStep();

        bpWriter.Close();
    }
    { // read
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("BP5");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        bpReader.BeginStep();
        auto var_r32 = io.InquireVariable<float>("r32cpu");
        EXPECT_TRUE(var_r32);

        std::vector<float> r32o(Nx);
        float *gpuSimData;
        cudaMalloc(&gpuSimData, Nx * sizeof(float));
        var_r32.SetMemorySpace(adios2::MemorySpace::Host);
        EXPECT_THROW(bpReader.Get(var_r32, gpuSimData, adios2::Mode::Sync), std::ios_base::failure);
        var_r32.SetMemorySpace(adios2::MemorySpace::GPU);
        EXPECT_THROW(bpReader.Get(var_r32, r32o.data(), adios2::Mode::Sync),
                     std::ios_base::failure);
        // bpReader.EndStep();
        // bpReader.Close();
    }
}

void CUDADetectMemSpace(const std::string mode)
{
    const std::string fname("BPWRCUDetect" + mode + ".bp");
    adios2::Mode ioMode = adios2::Mode::Deferred;
    if (mode == "Sync")
        ioMode = adios2::Mode::Sync;

    // Number of rows
    const size_t Nx = 5;
    // Number of columns
    const size_t Ny = 2;
    const size_t NTotal = Nx * Ny;
    // Number of steps
    const size_t NSteps = 10;

    adios2::ADIOS adios;
    // simulation data
    std::vector<float> r32s(NTotal, .0f);
    std::iota(r32s.begin(), r32s.end(), .0f);

    { // write
        adios2::IO io = adios.DeclareIO("TestIO");
        const adios2::Dims shape{Ny, Nx};
        const adios2::Dims start{0, 0};
        const adios2::Dims count{Ny, Nx};
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
        EXPECT_TRUE(var_r32);

        float *gpuSimData = nullptr;
        cudaMalloc(&gpuSimData, NTotal * sizeof(float));
        cudaMemcpy(gpuSimData, (float *)&r32s[0], NTotal * sizeof(float), cudaMemcpyHostToDevice);

        io.SetEngine("BP5");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Update values in the simulation data
            cuda_increment(NTotal, 1, 0, gpuSimData, INCREMENT);
            std::transform(r32s.begin(), r32s.end(), r32s.begin(),
                           std::bind(std::plus<float>(), std::placeholders::_1, INCREMENT));

            bpWriter.BeginStep();
            if (step % 2 == 0)
                bpWriter.Put(var_r32, gpuSimData, ioMode);
            else
                bpWriter.Put(var_r32, r32s.data(), ioMode);
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }
    // reset the initial data
    std::iota(r32s.begin(), r32s.end(), .0f);
    { // read
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("BP5");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        for (; bpReader.BeginStep() == adios2::StepStatus::OK; ++t)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], Ny);
            ASSERT_EQ(var_r32.Shape()[1], Nx);

            std::vector<float> r32o(NTotal);
            float *gpuSimData;
            cudaMalloc(&gpuSimData, NTotal * sizeof(float));
            if (t % 2 == 0)
            {
                bpReader.Get(var_r32, r32o.data(), ioMode);
                bpReader.EndStep();
            }
            else
            {
                bpReader.Get(var_r32, gpuSimData, ioMode);
                bpReader.EndStep();
                cudaMemcpy(r32o.data(), gpuSimData, NTotal * sizeof(float), cudaMemcpyDeviceToHost);
            }
            // Remove INCREMENT from each element
            std::transform(
                r32o.begin(), r32o.end(), r32o.begin(),
                std::bind(std::minus<float>(), std::placeholders::_1, (t + 1) * INCREMENT));

            for (size_t i = 0; i < NTotal; i++)
            {
                char msg[1 << 8] = {0};
                snprintf(msg, sizeof(msg), "t=%d i=%zu r32o=%f r32s=%f", t, i, r32o[i], r32s[i]);
                ASSERT_LT(std::abs(r32o[i] - r32s[i]), EPSILON) << msg;
            }
        }
        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void CUDAWriteReadMemorySelection()
{
    const std::string fname("BPWRCUSel1D.bp");
    const size_t Nx = 10;
    const size_t NSteps = 2;
    const size_t ghostCells = 1;
    std::vector<float> r32s(Nx + 2 * ghostCells);
    std::iota(r32s.begin(), r32s.end(), .0f);

    adios2::ADIOS adios;
    {
        // cuda simulation buffer
        float *gpuSimData = nullptr;
        cudaMalloc(&gpuSimData, (Nx + 2 * ghostCells) * sizeof(float));
        cudaMemcpy(gpuSimData, r32s.data(), (Nx + 2 * ghostCells) * sizeof(float),
                   cudaMemcpyHostToDevice);

        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetEngine("BP5");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx)};
        const adios2::Dims start{static_cast<size_t>(0)};
        const adios2::Dims count{Nx};
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);

        const adios2::Dims memoryStart = {ghostCells};
        const adios2::Dims memoryCount = {Nx + 2 * ghostCells};
        var_r32.SetMemorySelection({memoryStart, memoryCount});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            cuda_increment(Nx + 2 * ghostCells, 1, 0, gpuSimData, INCREMENT);

            bpWriter.BeginStep();
            var_r32.SetMemorySpace(adios2::MemorySpace::GPU);
            bpWriter.Put(var_r32, gpuSimData);
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }
    {
        // remove ghost cells from the input vector when checking correctness
        r32s.erase(r32s.begin(), r32s.begin() + ghostCells);
        r32s.erase(r32s.end() - ghostCells, r32s.end());

        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("BP5");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        for (; bpReader.BeginStep() == adios2::StepStatus::OK; ++t)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], Nx);

            auto mmR32 = std::minmax_element(r32s.begin(), r32s.end());
            EXPECT_EQ(var_r32.Min() - (t + 1) * INCREMENT, *mmR32.first);
            EXPECT_EQ(var_r32.Max() - (t + 1) * INCREMENT, *mmR32.second);

            std::vector<float> r32o(Nx);
            float *gpuSimData;
            cudaMalloc(&gpuSimData, Nx * sizeof(float));
            var_r32.SetMemorySpace(adios2::MemorySpace::GPU);
            bpReader.Get(var_r32, gpuSimData);
            bpReader.EndStep();
            cudaMemcpy(r32o.data(), gpuSimData, Nx * sizeof(float), cudaMemcpyDeviceToHost);

            // Remove INCREMENT from each element
            std::transform(
                r32o.begin(), r32o.end(), r32o.begin(),
                std::bind(std::minus<float>(), std::placeholders::_1, (t + 1) * INCREMENT));
            for (size_t i = 0; i < Nx; i++)
            {
                char msg[1 << 8] = {0};
                snprintf(msg, sizeof(msg), "t=%d i=%zu r32o=%f r32s=%f", t, i, r32o[i], r32s[i]);
                ASSERT_LT(std::abs(r32o[i] - r32s[i]), EPSILON) << msg;
            }
        }
        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void CUDAWriteReadMPI1D(const std::string mode)
{
    adios2::Mode ioMode = adios2::Mode::Deferred;
    if (mode == "Sync")
        ioMode = adios2::Mode::Sync;

    // Number of rows
    const size_t Nx = 100;
    // Number of steps
    const size_t NSteps = 10;

    int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRCU1D_" + mode + "_MPI.bp");
#else
    const std::string fname("BPWRCU1D_" + mode + ".bp");
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

    {
        // cuda simulation buffer
        float *gpuSimData = nullptr;
        cudaMalloc(&gpuSimData, Nx * sizeof(float));
        cudaMemcpy(gpuSimData, ((float *)&r32s[0] + (Nx * mpiRank)), Nx * sizeof(float),
                   cudaMemcpyHostToDevice);
        // host simulation buffer
        std::vector<float> simData(r32s.begin() + (Nx * mpiRank),
                                   r32s.begin() + (Nx * (mpiRank + 1)));

        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetEngine("BP5");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(NxTotal)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
        auto var_r32_host = io.DefineVariable<float>("r32host", shape, start, count);

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Update values in the simulation data
            cuda_increment(Nx, 1, 0, gpuSimData, INCREMENT);
            std::transform(simData.begin(), simData.end(), simData.begin(),
                           std::bind(std::plus<float>(), std::placeholders::_1, INCREMENT));

            bpWriter.BeginStep();
            var_r32.SetMemorySpace(adios2::MemorySpace::GPU);
            bpWriter.Put(var_r32, gpuSimData, ioMode);
            var_r32_host.SetMemorySpace(adios2::MemorySpace::Host);
            bpWriter.Put(var_r32_host, simData.data(), ioMode);
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    {
        adios2::IO io = adios.DeclareIO("ReadIO");
        io.SetEngine("BP5");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        for (; bpReader.BeginStep() == adios2::StepStatus::OK; ++t)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], NxTotal);

            auto mmR32 = std::minmax_element(r32s.begin(), r32s.end());
            EXPECT_EQ(var_r32.Min() - (t + 1) * INCREMENT, *mmR32.first);
            EXPECT_EQ(var_r32.Max() - (t + 1) * INCREMENT, *mmR32.second);

            std::vector<float> r32o(NxTotal);
            float *gpuSimData;
            cudaMalloc(&gpuSimData, NxTotal * sizeof(float));
            var_r32.SetMemorySpace(adios2::MemorySpace::GPU);
            bpReader.Get(var_r32, gpuSimData, ioMode);
            bpReader.EndStep();
            cudaMemcpy(r32o.data(), gpuSimData, NxTotal * sizeof(float), cudaMemcpyDeviceToHost);

            // Remove INCREMENT from each element
            std::transform(
                r32o.begin(), r32o.end(), r32o.begin(),
                std::bind(std::minus<float>(), std::placeholders::_1, (t + 1) * INCREMENT));

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

class BPWRCUDA : public ::testing::TestWithParam<std::string>
{
public:
    BPWRCUDA() = default;

    virtual void SetUp() {}
    virtual void TearDown() {}
};

TEST_P(BPWRCUDA, ADIOS2BPWRCUDA1D) { CUDAWriteReadMPI1D(GetParam()); }
TEST_P(BPWRCUDA, ADIOS2BPCUDADetect) { CUDADetectMemSpace(GetParam()); }
TEST_P(BPWRCUDA, ADIOS2BPCUDAWrong) { CUDAWrongMemSpace(); }
TEST_P(BPWRCUDA, ADIOS2BPCUDAMemSel) { CUDAWriteReadMemorySelection(); }

INSTANTIATE_TEST_SUITE_P(CudaRW, BPWRCUDA, ::testing::Values("deferred", "sync"));

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
