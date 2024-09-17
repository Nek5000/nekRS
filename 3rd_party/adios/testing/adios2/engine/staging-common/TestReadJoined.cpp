/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <thread>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

class CommonReadTest : public ::testing::Test
{
public:
    CommonReadTest() = default;
};

typedef std::chrono::duration<double> Seconds;

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

const int nsteps = 3;
const size_t Ncols = 4;
const std::vector<int> nblocksPerProcess = {2, 3, 2, 1, 3, 2};
int nMyTotalRows[nsteps];
int nTotalRows[nsteps];

// ADIOS2 Common read
TEST_F(CommonReadTest, ADIOS2CommonRead1D8)
{
    // Each process would write a 1x8 array and all processes would
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
    adios2::IO inIO = adios.DeclareIO("Input");

    inIO.SetEngine(engine);
    inIO.SetParameters(engineParams);

    adios2::Engine reader = inIO.Open(fname, adios2::Mode::Read);

    if (!mpiRank)
    {
        std::cout << "Reading as stream with BeginStep/EndStep on " << mpiSize
                  << "processes :" << std::endl;
    }

    int step = 0;
    while (true)
    {
        adios2::StepStatus status = reader.BeginStep(adios2::StepMode::Read);

        if (status != adios2::StepStatus::OK)
        {
            break;
        }

        auto rows_var = inIO.InquireVariable<int>("totalrows");
        auto var = inIO.InquireVariable<double>("table");
        EXPECT_TRUE(var);

        if (!mpiRank)
        {
            std::cout << "Step " << step << " table shape (" << var.Shape()[0] << ", "
                      << var.Shape()[1] << ")" << std::endl;
        }

        int Nrows;
        reader.Get(rows_var, Nrows);
        std::cout << "Reader expecting " << Nrows << std::endl;
        EXPECT_EQ(var.Shape()[0], Nrows);
        EXPECT_EQ(var.Shape()[1], Ncols);

        var.SetSelection({{0, 0}, {(size_t)Nrows, Ncols}});

        // Check data on rank 0
        if (!mpiRank)
        {
            std::vector<double> data(Nrows * Ncols);
            reader.Get(var, data.data());
            reader.PerformGets();
            for (size_t i = 0; i < (size_t)Nrows; ++i)
            {
                for (size_t j = 0; j < Ncols; ++j)
                {
                    EXPECT_GE(data[i * Ncols + j], (step + 1) * 1.0);
                    EXPECT_LT(data[i * Ncols + j], (nsteps + 1) * 1.0 + 0.9999);
                }
            }
        }

        reader.EndStep();
        ++step;
    }
    reader.Close();
    EXPECT_EQ(step, nsteps);
}

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    ParseArgs(argc, argv);

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engine == "SST" || engine == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, thread_support_level, &provided);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 2;
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
