/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>
#include <ctime>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

class CommonWriteTest : public ::testing::Test
{
public:
    CommonWriteTest() = default;
};

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

// ADIOS2 COMMON write
TEST_F(CommonWriteTest, ADIOS2CommonWrite)
{
    // form a mpiSize * Nx 1D array
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
    MPI_Comm_size(testComm, &mpiSize);
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    std::string varname = "r64";

    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    unsigned int myStart = (int)Nx * mpiRank;
    unsigned int myCount = (int)Nx;
    {
        adios2::Dims shape{static_cast<unsigned int>(Nx * mpiSize)};
        adios2::Dims start{static_cast<unsigned int>(myStart)};
        adios2::Dims count{static_cast<unsigned int>(myCount)};

        (void)io.DefineVariable<double>(varname, shape, start, count);
    }

    // Create the Engine
    io.SetEngine(engine);
    io.SetParameters(engineParams);
    /* ensure min/max are done */
    io.SetParameters({{"StatsLevel", "1"}});

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    for (size_t step = 0; step < 10; step++)
    {
        // Generate test data for each process uniquely
        std::vector<double> data_forward;

        engine.BeginStep();
        auto var = io.InquireVariable<double>(varname);

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces
        adios2::Box<adios2::Dims> sel({myStart}, {myCount});
        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count

        var.SetSelection(sel);

        adios2::Variable<double>::Span VarSpan = engine.Put(var);
        generateSimpleForwardData(VarSpan.data(), (int)step, myStart, myCount, (int)Nx * mpiSize);

        engine.EndStep();
    }

    // Close the file
    engine.Close();
}

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);

    ParseArgs(argc, argv);

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engine == "SST" || engine == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;
    MPI_Init_thread(nullptr, nullptr, thread_support_level, &provided);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &testComm);
#endif

    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
#ifdef CRAY_MPICH_VERSION
    MPI_Barrier(MPI_COMM_WORLD);
#else
    MPI_Finalize();
#endif
#endif

    return result;
}
