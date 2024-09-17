/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestBPLargeMetadata.cpp :
 *
 *  Created on: Aug 1, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPLargeMetadata : public ::testing::Test
{
public:
    BPLargeMetadata() = default;

    SmallTestData m_TestData;
};

TEST_F(BPLargeMetadata, BPWrite1D_LargeMetadata)
{
    // Each process would write a 4x2 array and all processes would
    // form a 2D 4 * (NumberOfProcess * Nx) matrix where Nx is 2 here

    int mpiRank = 0, mpiSize = 1;

    const std::size_t Nx = 10;
    const std::size_t NSteps = 1;
    const std::size_t NVars = 1000;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWrite1D_LargeMetadata_MPI.bp");
#else
    const std::string fname("BPWrite1D_LargeMetadata.bp");
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");
        // io.SetEngine("ADIOS1");
        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Dims shape{static_cast<size_t>(mpiSize * Nx)};
        adios2::Dims start{static_cast<size_t>(mpiRank * Nx)};
        adios2::Dims count{Nx};

        std::vector<adios2::Variable<float>> varsR32(NVars);
        std::vector<adios2::Variable<double>> varsR64(NVars);

        for (size_t i = 0; i < NVars; ++i)
        {
            varsR32[i] = io.DefineVariable<float>("varR32_" + std::to_string(i), shape, start,
                                                  count, adios2::ConstantDims);
            varsR64[i] = io.DefineVariable<double>("varR64_" + std::to_string(i), shape, start,
                                                   count, adios2::ConstantDims);
        }

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            bpWriter.BeginStep();
            for (size_t i = 0; i < NVars; ++i)
            {
                bpWriter.Put(varsR32[i], currentTestData.R32.data());
                bpWriter.Put(varsR64[i], currentTestData.R64.data());
            }
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }
}

TEST_F(BPLargeMetadata, ManyLongStrings)
{
    const std::string longString = "test_string "
                                   "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                                   "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                                   "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                                   "aaaaaaaaaaaaaaaaaaaaaaaaa";

    const std::size_t NVars = 100;

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    const std::string fname("BPWrite1D_LargeMetadataStrings_MPI.bp");
#else
    adios2::ADIOS adios;
    const std::string fname("BPWrite1D_LargeMetadataStrings.bp");
#endif

    adios2::IO io = adios.DeclareIO("myIO");
    if (!engineName.empty())
    {
        io.SetEngine(engineName);
    }
    {
        adios2::Engine writer = io.Open(fname, adios2::Mode::Write);

        for (size_t i = 0; i < NVars; ++i)
        {
            const std::string variableName = "string" + std::to_string(i);
            auto variableString = io.DefineVariable<std::string>(variableName);
            writer.Put(variableString, longString);
        }

        writer.Close();
    }
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // reader
    {
        io.RemoveAllVariables();
        adios2::Engine reader = io.Open(fname, adios2::Mode::ReadRandomAccess);

        std::string inString;

        for (size_t i = 0; i < NVars; ++i)
        {
            const std::string variableName = "string" + std::to_string(i);
            auto variableString = io.InquireVariable<std::string>(variableName);
            EXPECT_TRUE(variableString);
            reader.Get(variableString, inString);
            EXPECT_EQ(inString, longString);
        }

        reader.Close();
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
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
