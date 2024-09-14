#include <adios2.h>
#include <chrono>
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <thread>
#include <time.h>
#include <vector>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

void ReadVariable(const std::string &name, adios2::IO &io, adios2::Engine &reader, size_t step)
{
    adios2::Variable<double> variable = io.InquireVariable<double>(name);

    if (variable)
    {
        auto blocksInfo = reader.BlocksInfo(variable, step);

        std::cout << "    " << name << " has " << blocksInfo.size() << " blocks in step " << step
                  << std::endl;

        // create a data vector for each block
        std::vector<double> dataSet;
        dataSet.resize(blocksInfo.size());

        // schedule a read operation for each block separately
        for (auto &info : blocksInfo)
        {
            variable.SetBlockSelection(info.BlockID);
            int t0 = clock();
            reader.Get(variable, dataSet, adios2::Mode::Sync);
            printf("done in %f\n", 1.0 * (clock() - t0) / CLOCKS_PER_SEC);
        }
    }
    else
    {
        std::cout << "    Variable " << name << " not found in step " << step << std::endl;
    }
}

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

class CommonWriteTest : public ::testing::Test
{
public:
    CommonWriteTest() = default;
};

// ADIOS2 COMMON write
TEST_F(CommonWriteTest, ADIOS2CommonWrite)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
#endif

    // Write test data using ADIOS2

    // v0 has the same size on every process at every step
    const size_t Nglobal = 40;
    std::vector<double> v0(Nglobal);

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Variable<double> varV0 = io.DefineVariable<double>("v0", {}, {}, {Nglobal});

    adios2::Engine writer = io.Open(fname, adios2::Mode::Write);

    const int NSTEPS = 5;

    for (int step = 0; step < NSTEPS; step++)
    {
        writer.BeginStep();
        for (int block = 0; block < 500; block++)
        {
            // v0
            for (size_t i = 0; i < Nglobal; i++)
            {
                v0[i] = mpiRank * 1.0 + step * 0.1;
            }
            writer.Put<double>(varV0, v0.data(), adios2::Mode::Sync);
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
