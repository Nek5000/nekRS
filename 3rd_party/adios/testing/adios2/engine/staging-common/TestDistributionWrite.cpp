/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <chrono>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <thread>

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
    adios2::IO io1 = adios.DeclareIO("TestIO");

    std::string varname1 = "r64";
    std::string varname2 = "r64_2";

    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    unsigned int myStart1 = (int)Nx * mpiRank;
    unsigned int myCount1 = (int)Nx;
    adios2::Dims shape1{static_cast<unsigned int>(Nx * mpiSize)};
    adios2::Dims start1{static_cast<unsigned int>(myStart1)};
    adios2::Dims count1{static_cast<unsigned int>(myCount1)};

    {
        //        auto var1 =
        (void)io1.DefineVariable<double>(varname1, shape1, start1, count1);
        (void)io1.DefineVariable<double>(varname2, shape1, start1, count1);
        (void)io1.DefineVariable<size_t>("Step");
    }

    // Create the Engine
    io1.SetEngine(engine);
    if (RoundRobin)
    {
        engineParams["StepDistributionMode"] = "RoundRobin";
    }
    else if (OnDemand)
    {
        engineParams["StepDistributionMode"] = "OnDemand";
        std::cout << "Running this test in OnDemand mode may produce "
                     "unpredictable results.  Not suitable for CI."
                  << std::endl;
    }
    else
    {
        // default
    }
    io1.SetParameters(engineParams);

    adios2::Engine engine1 = io1.Open(fname, adios2::Mode::Write);

    for (size_t step = 0; step < (size_t)NSteps; step++)
    {

        // Generate test data for each process uniquely
        std::vector<double> data_forward;

        generateSimpleForwardData(data_forward, (int)step, myStart1, myCount1, (int)Nx * mpiSize);

        engine1.BeginStep();
        auto var1 = io1.InquireVariable<double>(varname1);
        auto step_var = io1.InquireVariable<size_t>("Step");

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces
        adios2::Box<adios2::Dims> sel1({myStart1}, {myCount1});

        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count
        const adios2::Mode sync = GlobalWriteMode;

        var1.SetSelection(sel1);
        engine1.Put(var1, data_forward.data(), sync);
        if (step > 5)
        {
            auto var2 = io1.InquireVariable<double>(varname2);
            var2.SetSelection(sel1);
            engine1.Put(var2, data_forward.data(), sync);
        }
        engine1.Put(step_var, step);
        engine1.EndStep();
    }
    // Close the file
    engine1.Close();
}

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    MPI_Init(nullptr, nullptr);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &testComm);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);

    ParseArgs(argc, argv);

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
