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
const std::size_t Nx = 10;
using DataArray = std::array<int32_t, Nx>;

class BPAppendAfterSteps : public ::testing::Test
{
public:
    BPAppendAfterSteps() = default;

    DataArray GenerateData(size_t step, int rank, int size)
    {
        DataArray d;
        int j = rank + 1 + (int)step * size;
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
};

class BPAppendAfterStepsP : public BPAppendAfterSteps,
                            public ::testing::WithParamInterface<std::tuple<size_t, int>>
{
protected:
    size_t GetSteps() { return std::get<0>(GetParam()); };
    int GetAppendAfterSteps() { return std::get<1>(GetParam()); };
};

TEST_P(BPAppendAfterStepsP, Test)
{
    int mpiRank = 0, mpiSize = 1;

    size_t nSteps = GetSteps();
    int nAppendAfterSteps = GetAppendAfterSteps();

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    std::cout << "Test AppendAfterSteps parameter with writing " << nSteps
              << " steps, then appending " << nSteps << " steps again with parameter "
              << nAppendAfterSteps << std::endl;

    std::string filename = "AppendAfterSteps_N" + std::to_string(mpiSize) + "_Steps" +
                           std::to_string(nSteps) + "_Append_" + std::to_string(nAppendAfterSteps) +
                           ".bp";
    size_t totalNSteps = 0;

    {
        /* Write nSteps steps */
        adios2::IO ioWrite = adios.DeclareIO("TestIOWrite");
        ioWrite.SetEngine(engineName);
        adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
        adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
        adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        auto var0 = ioWrite.DefineVariable<int32_t>("var", shape, start, count);
        auto var1 = ioWrite.DefineVariable<size_t>("step");
        for (size_t step = 0; step < nSteps; ++step)
        {
            auto d = GenerateData(step, mpiRank, mpiSize);
            engine.BeginStep();
            std::string aname = "a" + std::to_string(step);
            ioWrite.DefineAttribute<size_t>(aname, step);
            engine.Put(var0, d.data());
            engine.Put(var1, step);
            engine.EndStep();
        }
        engine.Close();

        /* Write again nSteps steps but append to file at nAppendAfterSteps*/
        ioWrite.SetParameter("AppendAfterSteps", std::to_string(nAppendAfterSteps));
        engine = ioWrite.Open(filename, adios2::Mode::Append);

        size_t beginStep = 0;
        if (nAppendAfterSteps > (int)nSteps)
        {
            beginStep = nSteps;
        }
        else if (nAppendAfterSteps < 0)
        {
            int t = (int)nSteps + 1 + nAppendAfterSteps;
            if (t < 0)
            {
                t = 0;
            }
            beginStep = (size_t)t;
        }
        else
        {
            beginStep = nAppendAfterSteps;
        }
        totalNSteps = beginStep + nSteps;

        for (size_t step = beginStep; step < beginStep + nSteps; ++step)
        {
            auto d = GenerateData(step, mpiRank, mpiSize);
            engine.BeginStep();
            std::string aname = "a" + std::to_string(step);
            ioWrite.DefineAttribute<size_t>(aname, step);
            engine.Put(var0, d.data());
            engine.Put(var1, step);
            engine.EndStep();
        }
        engine.Close();
    }

    {
        adios2::IO ioRead = adios.DeclareIO("TestIORead");
        ioRead.SetEngine(engineName);
        adios2::Engine engine_s = ioRead.Open(filename, adios2::Mode::ReadRandomAccess);
        EXPECT_TRUE(engine_s);

        const size_t nsteps = engine_s.Steps();
        EXPECT_EQ(nsteps, totalNSteps);

        adios2::Variable<int> var = ioRead.InquireVariable<int32_t>("var");
        adios2::Variable<size_t> varStep = ioRead.InquireVariable<size_t>("step");

        for (size_t readStep = 0; readStep < totalNSteps; readStep++)
        {
            var.SetStepSelection(adios2::Box<size_t>(readStep, 1));
            std::vector<int> res;
            var.SetSelection({{Nx * mpiRank}, {Nx}});
            engine_s.Get<int>(var, res, adios2::Mode::Sync);
            auto d = GenerateData(readStep, mpiRank, mpiSize);
            EXPECT_EQ(res[0], d[0]);

            varStep.SetStepSelection(adios2::Box<size_t>(readStep, 1));
            size_t stepInFile;
            engine_s.Get<size_t>(varStep, stepInFile);
            EXPECT_EQ(stepInFile, readStep);

            std::string aname = "a" + std::to_string(readStep);
            adios2::Attribute<size_t> a = ioRead.InquireAttribute<size_t>(aname);
            size_t stepInAttribute = a.Data()[0];
            EXPECT_EQ(stepInAttribute, readStep);
        }

        engine_s.Close();
    }
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

INSTANTIATE_TEST_SUITE_P(BPAppendAfterSteps, BPAppendAfterStepsP,
                         ::testing::Values(std::make_tuple(1, 0), std::make_tuple(1, 1),
                                           std::make_tuple(2, 0), std::make_tuple(2, 1),
                                           std::make_tuple(2, 2), std::make_tuple(2, 3),
                                           std::make_tuple(2, 1), std::make_tuple(2, -1),
                                           std::make_tuple(2, -2), std::make_tuple(2, -3)));

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
