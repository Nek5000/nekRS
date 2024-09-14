/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>
#include <ctime>

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <thread>

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

const int nsteps = 3;
const size_t Ncols = 4;
const std::vector<int> nblocksPerProcess = {2, 3, 2, 1, 3, 2};
int nMyTotalRows[nsteps];
int nTotalRows[nsteps];

// ADIOS2 COMMON write
TEST_F(CommonWriteTest, ADIOS2CommonWrite)
{
    // form a mpiSize * Nx 1D array
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
    MPI_Comm_size(testComm, &mpiSize);
    const int nblocks =
        (mpiRank < static_cast<int>(nblocksPerProcess.size()) ? nblocksPerProcess[mpiRank] : 1);
#else
    const int nblocks = nblocksPerProcess[0];
#endif

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO outIO = adios.DeclareIO("Output");

    EXPECT_LE(mpiSize, nblocksPerProcess.size());
    outIO.SetEngine(engine);
    outIO.SetParameters(engineParams);

    adios2::Engine writer = outIO.Open(fname, adios2::Mode::Write);
    auto bpp_var = outIO.DefineVariable<int>("blocksperprocess", {nblocksPerProcess.size()}, {0},
                                             {nblocksPerProcess.size()});

    auto rows_var = outIO.DefineVariable<int>("totalrows");

    auto var = outIO.DefineVariable<double>("table", {adios2::JoinedDim, Ncols}, {}, {1, Ncols});

    if (!mpiRank)
    {
        std::cout << "Writing to " << fname << std::endl;
    }

    for (int step = 0; step < nsteps; step++)
    {
        // Application variables for output random size per process, 5..10
        // each
        std::vector<size_t> Nrows;
        nMyTotalRows[step] = 0;
        for (int i = 0; i < nblocks; ++i)
        {
            int n = rand() % 6 + 5;
            Nrows.push_back(static_cast<size_t>(n));
            nMyTotalRows[step] += n;
        }

        nTotalRows[step] = nMyTotalRows[step];
#if ADIOS2_USE_MPI
        MPI_Allreduce(&(nMyTotalRows[step]), &(nTotalRows[step]), 1, MPI_INT, MPI_SUM, testComm);
#endif

        if (!mpiRank)
        {
            std::cout << "Writing " << nTotalRows[step] << " rows in step " << step << std::endl;
        }

        writer.BeginStep();
        if ((step == 0) && (mpiRank == 0))
        {
            writer.Put(bpp_var, nblocksPerProcess.data());
        }
        if (mpiRank == 0)
        {
            std::cout << "Writer Generating " << nTotalRows[step] << " in total" << std::endl;
            writer.Put(rows_var, nTotalRows[step]);
        }
        for (int block = 0; block < nblocks; ++block)
        {
            std::vector<double> mytable(Nrows[block] * Ncols);
            for (size_t row = 0; row < Nrows[block]; row++)
            {
                for (size_t col = 0; col < Ncols; col++)
                {
                    mytable[row * Ncols + col] =
                        static_cast<double>((step + 1) * 1.0 + mpiRank * 0.1 + block * 0.01 +
                                            row * 0.001 + col * 0.0001);
                }
            }

            var.SetSelection({{}, {Nrows[block], Ncols}});

            std::cout << "Step " << step << " rank " << mpiRank << " block " << block << " count ("
                      << var.Count()[0] << ", " << var.Count()[1] << ")" << std::endl;

            writer.Put(var, mytable.data(), adios2::Mode::Sync);
        }
        writer.EndStep();
    }
    writer.Close();
}

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);

    DelayMS = 0; // zero for common writer

    ParseArgs(argc, argv);

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engine == "SST" || engine == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
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
