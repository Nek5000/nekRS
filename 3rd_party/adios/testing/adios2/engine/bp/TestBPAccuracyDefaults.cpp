/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestBPAccuracyDefaults.cpp :
 *
 *  Created on: Dec 5, 2023
 *      Author: Norbert Podhorszki
 */

#include <cstdint>
#include <cstring>

#include <iostream>
#include <numeric> //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

std::string engineName;       // comes from command line
std::string engineParameters; // comes from command line

class AccuracyTests : public ::testing::Test
{
public:
    AccuracyTests() = default;
};

// Check if SetAccuracy/GetAccuracy default behavior works
TEST_F(AccuracyTests, DefaultAccuracy)
{

    int mpiRank = 0, mpiSize = 1;

    const std::size_t Nx = 10;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("DefaultAccuracy_MPI.bp");
#else
    const std::string fname("DefaultAccuracy.bp");
#endif

    std::vector<int64_t> localData(Nx);
    std::iota(localData.begin(), localData.end(), mpiRank * Nx);

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("DefaultAccuracyWrite");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        if (!engineParameters.empty())
        {
            io.SetParameters(engineParameters);
        }

        adios2::Variable<int64_t> var =
            io.DefineVariable<int64_t>("range", {static_cast<std::size_t>(Nx * mpiSize)},
                                       {static_cast<std::size_t>(Nx * mpiRank)}, {Nx});

        adios2::Accuracy accuracyRequested = {0.001, std::numeric_limits<double>::infinity(), true};
        var.SetAccuracy(accuracyRequested);
        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        bpWriter.Put<int64_t>("range", localData.data());
        bpWriter.Close();

        auto accuracyGot = var.GetAccuracy();
        EXPECT_EQ(accuracyGot.error, 0.0); // no error whatsoever
        EXPECT_EQ(accuracyGot.norm, accuracyRequested.norm);
        EXPECT_EQ(accuracyGot.relative, accuracyRequested.relative);
    }
    // Reader
    {
        adios2::IO io = adios.DeclareIO("DefaultAccuracyRead");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        if (!engineParameters.empty())
        {
            io.SetParameters(engineParameters);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);
        adios2::Variable<int64_t> varRange = io.InquireVariable<int64_t>("range");
        adios2::Accuracy accuracyRequested = {0.001, 1.17, false};
        varRange.SetAccuracy(accuracyRequested);

        const std::size_t gNx = static_cast<std::size_t>(Nx * mpiSize);
        std::vector<int64_t> globalData(gNx);
        bpReader.Get(varRange, globalData);
        bpReader.PerformGets();

        auto accuracyGot = varRange.GetAccuracy();
        EXPECT_EQ(accuracyGot.error, 0.0); // no error whatsoever
        EXPECT_EQ(accuracyGot.error, 0.0); // no error whatsoever
        EXPECT_EQ(accuracyGot.norm, accuracyRequested.norm);
        EXPECT_EQ(accuracyGot.relative, accuracyRequested.relative);

        std::vector<int64_t> iStartEndData;
        iStartEndData.reserve(gNx); // maximum possible

        for (size_t i = 1; i < gNx; ++i)
        {
            varRange.SetSelection({{i}, {gNx - i}});

            bpReader.Get("range", iStartEndData);
            bpReader.PerformGets();

            for (size_t j = i; j < gNx; ++j)
            {
                EXPECT_EQ(globalData[j], iStartEndData[j - i]);
            }
        }
        bpReader.Close();
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
    if (argc > 2)
    {
        engineParameters = std::string(argv[2]);
    }
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
