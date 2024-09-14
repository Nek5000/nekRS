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
    adios2::IO io1 = adios.DeclareIO("TestIO1");
    adios2::IO tmp_io = adios.DeclareIO("TestIO2");

    adios2::IO *io2;

    if (SharedIO)
    {
        io2 = &io1;
    }
    else
    {
        io2 = &tmp_io;
    }
    std::string varname1 = "r64";
    std::string varname2 = "r64_2";

    if (SharedVar)
    {
        varname2 = "r64";
    }
    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    unsigned int myStart1 = (int)Nx * mpiRank, myStart2 = (int)Nx * mpiRank;
    unsigned int myCount1 = (int)Nx, myCount2 = (int)Nx;
    if (mpiRank == 0)
    {
        /* first guy gets twice allotment var 1 */
        myCount1 = 2 * (int)Nx;
    }
    else
    {
        /* everyone else shifts up */
        myStart1 += (int)Nx;
    }
    if (mpiRank == (mpiSize - 1))
    {
        /* last guy  gets twice allotment var 2 */
        myCount2 = 2 * (int)Nx;
    }
    {
        adios2::Dims shape1{static_cast<unsigned int>(Nx * (mpiSize + 1))};
        adios2::Dims start1{static_cast<unsigned int>(myStart1)};
        adios2::Dims count1{static_cast<unsigned int>(myCount1)};
        adios2::Dims shape2{static_cast<unsigned int>(Nx * (mpiSize + 1))};
        adios2::Dims start2{static_cast<unsigned int>(myStart2)};
        adios2::Dims count2{static_cast<unsigned int>(myCount2)};

        //        auto var1 =
        (void)io1.DefineVariable<double>(varname1, shape1, start1, count1);
        if (!SharedVar)
        {
            //
            (void)io2->DefineVariable<double>(varname2, shape2, start2, count2);
        }
    }

    // Create the Engine
    io1.SetEngine(engine);
    io1.SetParameters(engineParams);
    io2->SetEngine(engine);
    io2->SetParameters(engineParams);

    std::string fname1 = fname + "1";
    std::string fname2 = fname + "2";
    adios2::Engine engine1 = io1.Open(fname1, adios2::Mode::Write);
    adios2::Engine engine2 = io2->Open(fname2, adios2::Mode::Write);

    {
        size_t step = 0;
        // Generate test data for each process uniquely
        std::vector<double> data_forward;
        std::vector<double> data_reverse;

        generateSimpleForwardData(data_forward, (int)step, myStart1, myCount1,
                                  (int)Nx * (mpiSize + 1));
        generateSimpleReverseData(data_reverse, (int)step, myStart2, myCount2,
                                  (int)Nx * (mpiSize + 1));

        engine1.BeginStep();
        engine2.BeginStep();
        auto var1 = io1.InquireVariable<double>(varname1);
        auto var2 = io2->InquireVariable<double>(varname2);

        if (SharedVar)
        {
            var2 = var1;
        }

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces
        adios2::Box<adios2::Dims> sel1({myStart1}, {myCount1});
        adios2::Box<adios2::Dims> sel2({myStart2}, {myCount2});
        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count
        const adios2::Mode sync = GlobalWriteMode;

        /* var1 and var2 may be the same variable (--shared_var), so nothing
         * should come between the SetSelection/Put pair for each. */
        var1.SetSelection(sel1);
        engine1.Put(var1, data_forward.data(), sync);

        var2.SetSelection(sel2);
        engine2.Put(var2, data_reverse.data(), sync);

        engine1.EndStep();
        engine2.EndStep();
    }

    // Close the file
    engine1.Close();
    engine2.Close();
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
