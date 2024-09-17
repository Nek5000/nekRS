/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>
#include <ctime>

#include <chrono>
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <thread>
#include <time.h>
#include <vector>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

static void ReadVariable(const std::string &name, adios2::IO &io, adios2::Engine &reader,
                         size_t step)
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

class LocalReadTest : public ::testing::Test
{
public:
    LocalReadTest() = default;
};

// ADIOS2 COMMON write
TEST_F(LocalReadTest, ADIOS2LocalRead)
{

    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Engine reader = io.Open(fname, adios2::Mode::Read);

    while (true)
    {

        // Begin step
        adios2::StepStatus read_status = reader.BeginStep(adios2::StepMode::Read, 10.0f);
        if (read_status == adios2::StepStatus::NotReady)
        {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (read_status != adios2::StepStatus::OK)
        {
            break;
        }

        size_t step = reader.CurrentStep();
        std::cout << "Process step " << step << ": " << std::endl;
        if (step == 0)
        {
            ReadVariable("v0", io, reader, step);
        }
        reader.EndStep();
    }

    reader.Close();
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
