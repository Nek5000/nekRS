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
constexpr std::size_t NSteps = 10;
const std::size_t Nx = 10;
using DataArray = std::array<int32_t, Nx>;

class BPParameterSelectSteps : public ::testing::Test
{
public:
    BPParameterSelectSteps() = default;

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
        std::string filename = "ParameterSelectSteps" + std::to_string(mpiSize) + ".bp";
        adios2::IO ioWrite = adios.DeclareIO("TestIOWrite");
        ioWrite.SetEngine(engineName);
        adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
        // Number of elements per process
        const std::size_t Nx = 10;
        adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
        adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        auto var0 = ioWrite.DefineVariable<int32_t>("var", shape, start, count);
        auto var1 = ioWrite.DefineVariable<size_t>("step");
        for (size_t step = 0; step < NSteps; ++step)
        {
            int s = static_cast<int>(step);
            auto d = GenerateData(s, mpiRank, mpiSize);
            engine.BeginStep();
            std::string aname = "a" + std::to_string(step);
            ioWrite.DefineAttribute<size_t>(aname, step);
            engine.Put(var0, d.data());
            engine.Put(var1, step);
            engine.EndStep();
        }
        engine.Close();
        OutputWritten = true;
    }
};

class BPParameterSelectStepsP
: public BPParameterSelectSteps,
  public ::testing::WithParamInterface<std::tuple<std::string, std::vector<size_t>>>
{
protected:
    std::string GetSelectionString() { return std::get<0>(GetParam()); };
    std::vector<size_t> GetSteps() { return std::get<1>(GetParam()); };
};

TEST_P(BPParameterSelectStepsP, Read)
{
    int mpiRank = 0, mpiSize = 1;
    std::string selection = GetSelectionString();
    // with this selection these are the original (absolute) steps
    // that this reader will now see as steps 0,1,2,...
    std::vector<size_t> absoluteSteps = GetSteps();

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
    std::string filename = "ParameterSelectSteps" + std::to_string(mpiSize) + ".bp";
    adios2::IO ioRead = adios.DeclareIO("TestIORead");
    ioRead.SetEngine(engineName);
    ioRead.SetParameter("SelectSteps", selection);
    adios2::Engine engine_s = ioRead.Open(filename, adios2::Mode::ReadRandomAccess);
    EXPECT_TRUE(engine_s);

    const size_t nsteps = engine_s.Steps();
    EXPECT_EQ(nsteps, absoluteSteps.size());

    adios2::Variable<int> var = ioRead.InquireVariable<int32_t>("var");
    adios2::Variable<size_t> varStep = ioRead.InquireVariable<size_t>("step");

    for (size_t readStep = 0; readStep < nsteps; readStep++)
    {
        var.SetStepSelection(adios2::Box<size_t>(readStep, 1));
        std::vector<int> res;
        var.SetSelection({{Nx * mpiRank}, {Nx}});
        engine_s.Get<int>(var, res, adios2::Mode::Sync);
        int s = static_cast<int>(absoluteSteps[readStep]);
        auto d = GenerateData(s, mpiRank, mpiSize);
        EXPECT_EQ(res[0], d[0]);

        varStep.SetStepSelection(adios2::Box<size_t>(readStep, 1));
        size_t stepInFile;
        engine_s.Get<size_t>(varStep, stepInFile);
        EXPECT_EQ(stepInFile, absoluteSteps[readStep]);

        std::string aname = "a" + std::to_string(absoluteSteps[readStep]);
        adios2::Attribute<size_t> a = ioRead.InquireAttribute<size_t>(aname);
        size_t stepInAttribute = a.Data()[0];
        EXPECT_EQ(stepInAttribute, absoluteSteps[readStep]);
    }

    engine_s.Close();
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

TEST_P(BPParameterSelectStepsP, Stream)
{
    int mpiRank = 0, mpiSize = 1;
    std::string selection = GetSelectionString();
    // with this selection these are the original (absolute) steps
    // that this reader will now see as steps 0,1,2,...
    std::vector<size_t> absoluteSteps = GetSteps();

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    std::string filename = "ParameterSelectStepsStream" + std::to_string(mpiSize) + ".bp";
    adios2::IO ioWrite = adios.DeclareIO("TestIOWrite");
    ioWrite.SetEngine(engineName);
    adios2::Engine writer = ioWrite.Open(filename, adios2::Mode::Write);

    adios2::IO ioRead = adios.DeclareIO("TestIORead");
    ioRead.SetEngine(engineName);
    ioRead.SetParameter("SelectSteps", selection);
    ioRead.SetParameter("OpenTimeoutSecs", "30.0");
    adios2::Engine reader = ioRead.Open(filename, adios2::Mode::Read);
    EXPECT_TRUE(reader);

    // Number of elements per process
    const std::size_t Nx = 10;
    adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
    adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
    adios2::Dims count{static_cast<unsigned int>(Nx)};

    auto var0 = ioWrite.DefineVariable<int32_t>("var", shape, start, count);
    auto var1 = ioWrite.DefineVariable<size_t>("step");

    for (size_t step = 0; step < NSteps; ++step)
    {
        int s = static_cast<int>(step);
        auto d = GenerateData(s, mpiRank, mpiSize);
        writer.BeginStep();
        std::string aname = "a" + std::to_string(step);
        ioWrite.DefineAttribute<size_t>(aname, step);
        writer.Put(var0, d.data());
        writer.Put(var1, step);
        writer.EndStep();

        if (!mpiRank)
        {
            std::cout << "Writer done step " << step << std::endl;
        }
        auto status = reader.BeginStep(adios2::StepMode::Read, 0.0f);
        if (!mpiRank)
        {
            std::cout << "Reader BeginStep() step " << step << " status " << status << std::endl;
        }
        if (status == adios2::StepStatus::OK)
        {
            size_t readStep = reader.CurrentStep();
            if (!mpiRank)
            {
                std::cout << "Reader got read step " << readStep
                          << ". Check if it equals to writer step " << absoluteSteps[readStep]
                          << std::endl;
            }
            std::vector<int32_t> res;
            adios2::Variable<int32_t> var = ioRead.InquireVariable<int32_t>("var");
            var.SetSelection({{Nx * mpiRank}, {Nx}});
            reader.Get<int32_t>(var, res, adios2::Mode::Sync);
            int s = static_cast<int>(absoluteSteps[readStep]);
            auto d = GenerateData(s, mpiRank, mpiSize);
            EXPECT_EQ(res[0], d[0]);

            adios2::Variable<size_t> varStep = ioRead.InquireVariable<size_t>("step");
            size_t stepInFile;
            reader.Get<size_t>(varStep, stepInFile);
            EXPECT_EQ(stepInFile, absoluteSteps[readStep]);

            std::string aname = "a" + std::to_string(absoluteSteps[readStep]);
            adios2::Attribute<size_t> a = ioRead.InquireAttribute<size_t>(aname);
            size_t stepInAttribute = a.Data()[0];
            EXPECT_EQ(stepInAttribute, absoluteSteps[readStep]);

            reader.EndStep();
        }
        else
        {
            EXPECT_EQ(status, adios2::StepStatus::NotReady);
        }
    }
    writer.Close();
    auto status = reader.BeginStep(adios2::StepMode::Read, 1.0f);
    EXPECT_EQ(status, adios2::StepStatus::EndOfStream);
    reader.Close();
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

const std::vector<size_t> s_0n1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
const std::vector<size_t> s_152 = {1, 3, 5};
const std::vector<size_t> s_3n3 = {3, 6, 9};
const std::vector<size_t> s_1n2_0n2 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

INSTANTIATE_TEST_SUITE_P(BPParameterSelectSteps, BPParameterSelectStepsP,
                         ::testing::Values(std::make_tuple("0:n:1", s_0n1),
                                           std::make_tuple("1:5:2", s_152),
                                           std::make_tuple("3:n:3", s_3n3),
                                           std::make_tuple("1:n:2 0:n:2", s_1n2_0n2)));

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
