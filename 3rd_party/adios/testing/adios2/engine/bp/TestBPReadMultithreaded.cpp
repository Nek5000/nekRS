/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Test "SelectSteps" parameter for reading a BP file
 *
 *  Created on: Feb 4, 2021
 *      Author: Norbert Podhorszki
 */

#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line
constexpr std::size_t NSteps = 4;
const std::size_t Nx = 10;
using DataArray = std::array<int32_t, Nx>;

class BPReadMultithreadedTest : public ::testing::Test
{
public:
    BPReadMultithreadedTest() = default;

    DataArray GenerateData(int step, int rank, int size)
    {
        DataArray d;
        int j = rank + 1 + step * size;
        for (size_t i = 0; i < d.size(); ++i)
        {
            d[i] = j;
        }
        return d;
    }

    std::string ArrayToString(int32_t *data, size_t nelems)
    {
        std::stringstream ss;
        ss << "[";
        for (size_t i = 0; i < nelems; ++i)
        {
            ss << data[i];
            if (i < nelems - 1)
            {
                ss << " ";
            }
        }
        ss << "]";
        return ss.str();
    }

    bool OutputWritten = false;

    void CreateOutput()
    {
        // This is not really a test but to create a dataset for all read tests
        if (OutputWritten)
        {
            return;
        }
        int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        std::string filename = "BPReadMultithreaded" + std::to_string(mpiSize) + ".bp";
        adios2::IO ioWrite = adios.DeclareIO("TestIOWrite");
        ioWrite.SetEngine(engineName);
        adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
        // Number of elements per process
        const std::size_t Nx = 10;
        adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
        adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        auto v1 = ioWrite.DefineVariable<int32_t>("v1", shape, start, count);
        auto v2 = ioWrite.DefineVariable<int32_t>("v2", shape, start, count);
        auto v3 = ioWrite.DefineVariable<int32_t>("v3", shape, start, count);
        auto v4 = ioWrite.DefineVariable<int32_t>("v4", shape, start, count);
        for (size_t step = 0; step < NSteps; ++step)
        {
            int s = static_cast<int>(step);
            auto d = GenerateData(s, mpiRank, mpiSize);
            engine.BeginStep();
            engine.Put(v1, d.data());
            engine.Put(v2, d.data());
            engine.Put(v3, d.data());
            engine.Put(v4, d.data());
            engine.EndStep();
        }
        engine.Close();
        OutputWritten = true;
    }
};

class BPReadMultithreadedTestP : public BPReadMultithreadedTest,
                                 public ::testing::WithParamInterface<int>
{
protected:
    int GetThreads() { return GetParam(); };
};

TEST_P(BPReadMultithreadedTestP, ReadFile)
{
    int mpiRank = 0, mpiSize = 1;
    int nThreads = GetThreads();
    std::cout << "---- Test Multithreaded ReadRandomAccess with " << nThreads << " threads ----"
              << std::endl;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    CreateOutput();
    std::string filename = "BPReadMultithreaded" + std::to_string(mpiSize) + ".bp";
    adios2::IO ioRead = adios.DeclareIO("TestIORead");
    ioRead.SetEngine(engineName);
    ioRead.SetParameter("Threads", std::to_string(nThreads));
    adios2::Engine reader = ioRead.Open(filename, adios2::Mode::ReadRandomAccess);
    EXPECT_TRUE(reader);

    const size_t nsteps = reader.Steps();
    EXPECT_EQ(nsteps, NSteps);

    adios2::Variable<int> v1 = ioRead.InquireVariable<int32_t>("v1");
    adios2::Variable<int> v3 = ioRead.InquireVariable<int32_t>("v2");
    adios2::Variable<int> v2 = ioRead.InquireVariable<int32_t>("v3");
    adios2::Variable<int> v4 = ioRead.InquireVariable<int32_t>("v4");

    std::vector<int> res1, res2, res3, res4;
    res1.resize(NSteps * Nx);
    res2.resize(NSteps * Nx);
    res3.resize(NSteps * Nx);
    res4.resize(NSteps * Nx);
    v1.SetSelection({{Nx * mpiRank}, {Nx}});
    v2.SetSelection({{Nx * mpiRank}, {Nx}});
    v3.SetSelection({{Nx * mpiRank}, {Nx}});
    v4.SetSelection({{Nx * mpiRank}, {Nx}});
    v1.SetStepSelection(adios2::Box<size_t>(0, nsteps));
    v2.SetStepSelection(adios2::Box<size_t>(0, nsteps));
    v3.SetStepSelection(adios2::Box<size_t>(0, nsteps));
    v4.SetStepSelection(adios2::Box<size_t>(0, nsteps));
    reader.Get<int>(v1, res1, adios2::Mode::Deferred);
    reader.Get<int>(v2, res2, adios2::Mode::Deferred);
    reader.Get<int>(v3, res3, adios2::Mode::Deferred);
    reader.Get<int>(v4, res4, adios2::Mode::Deferred);
    reader.PerformGets();

    for (size_t step = 0; step < nsteps; step++)
    {
        int s = static_cast<int>(step);
        auto d = GenerateData(s, mpiRank, mpiSize);
        EXPECT_EQ(res1[step * Nx], d[0]);
        EXPECT_EQ(res2[step * Nx], d[0]);
        EXPECT_EQ(res3[step * Nx], d[0]);
        EXPECT_EQ(res4[step * Nx], d[0]);
    }

    reader.Close();
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

TEST_P(BPReadMultithreadedTestP, ReadStream)
{
    int mpiRank = 0, mpiSize = 1;
    int nThreads = GetThreads();
    std::cout << "---- Test Multithreaded stream Read with " << nThreads << " threads ----"
              << std::endl;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    CreateOutput();
    std::string filename = "BPReadMultithreaded" + std::to_string(mpiSize) + ".bp";
    adios2::IO ioRead = adios.DeclareIO("TestIORead");
    ioRead.SetEngine(engineName);
    ioRead.SetParameter("Threads", std::to_string(nThreads));
    adios2::Engine reader = ioRead.Open(filename, adios2::Mode::Read);
    EXPECT_TRUE(reader);

    // Number of elements per process
    const std::size_t Nx = 10;
    adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
    adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
    adios2::Dims count{static_cast<unsigned int>(Nx)};

    for (size_t step = 0; step < NSteps; ++step)
    {
        reader.BeginStep();
        adios2::Variable<int> v1 = ioRead.InquireVariable<int32_t>("v1");
        adios2::Variable<int> v3 = ioRead.InquireVariable<int32_t>("v2");
        adios2::Variable<int> v2 = ioRead.InquireVariable<int32_t>("v3");
        adios2::Variable<int> v4 = ioRead.InquireVariable<int32_t>("v4");

        std::vector<int> res1, res2, res3, res4;
        res1.resize(NSteps * Nx);
        res2.resize(NSteps * Nx);
        res3.resize(NSteps * Nx);
        res4.resize(NSteps * Nx);
        v1.SetSelection({{Nx * mpiRank}, {Nx}});
        v2.SetSelection({{Nx * mpiRank}, {Nx}});
        v3.SetSelection({{Nx * mpiRank}, {Nx}});
        v4.SetSelection({{Nx * mpiRank}, {Nx}});
        reader.Get<int>(v1, res1, adios2::Mode::Deferred);
        reader.Get<int>(v2, res2, adios2::Mode::Deferred);
        reader.Get<int>(v3, res3, adios2::Mode::Deferred);
        reader.Get<int>(v4, res4, adios2::Mode::Deferred);
        reader.EndStep();

        int s = static_cast<int>(step);
        auto d = GenerateData(s, mpiRank, mpiSize);
        EXPECT_EQ(res1[0], d[0]);
        EXPECT_EQ(res2[0], d[0]);
        EXPECT_EQ(res3[0], d[0]);
        EXPECT_EQ(res4[0], d[0]);
    }

    auto status = reader.BeginStep(adios2::StepMode::Read, 1.0f);
    EXPECT_EQ(status, adios2::StepStatus::EndOfStream);
    reader.Close();
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

INSTANTIATE_TEST_SUITE_P(BPReadMultithreadedTest, BPReadMultithreadedTestP,
                         ::testing::Values(1, 2, 3, 0));

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    MPI_Init(nullptr, nullptr);
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
