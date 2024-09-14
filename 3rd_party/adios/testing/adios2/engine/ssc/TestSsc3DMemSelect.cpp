/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestSsc3DMemSelect.cpp
 *
 *  Created on: Aug 12, 2022
 *      Author: Jason Wang
 */

#include <adios2.h>
#include <gtest/gtest.h>
#include <numeric>
#include <thread>

using namespace adios2;
int mpiRank = 0;
int mpiSize = 1;
MPI_Comm mpiComm;

size_t steps = 100;

Dims shape = {4, 4, 4};
std::vector<int> global_data = {
    0,   1,   2,   3,   10,  11,  12,  13,  20,  21,  22,  23,  30,  31,  32,  33,

    100, 101, 102, 103, 110, 111, 112, 113, 120, 121, 122, 123, 130, 131, 132, 133,

    200, 201, 202, 203, 210, 211, 212, 213, 220, 221, 222, 223, 230, 231, 232, 233,

    300, 301, 302, 303, 310, 311, 312, 313, 320, 321, 322, 323, 330, 331, 332, 333};

Dims start = {1, 2, 1};
Dims count = {2, 1, 2};
std::vector<int> writer_data = {121, 122, 221, 222};

Dims memStart = {0, 1, 1};
Dims memCount = {3, 3, 3};
std::vector<int> reader_data = {11,  12,  13,  21,  22,  23,  31,  32,  33,

                                111, 112, 113, 121, 122, 123, 131, 132, 133,

                                211, 212, 213, 221, 222, 223, 231, 232, 233};

class SscEngineTest : public ::testing::Test
{
public:
    SscEngineTest() = default;
};

template <class T>
void PrintData(const T *data, const size_t step, const Dims &start, const Dims &count)
{
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    std::cout << "Step: " << step << " Size:" << size << "\n";
    size_t printsize = 128;

    if (size < printsize)
    {
        printsize = size;
    }
    int s = 0;
    for (size_t i = 0; i < printsize; ++i)
    {
        ++s;
        std::cout << data[i] << " ";
        if (s == count[1])
        {
            std::cout << std::endl;
            s = 0;
        }
    }

    std::cout << "]" << std::endl;
}

void VerifyData(const int *data, size_t step, const Dims &start, const Dims &count,
                const Dims &shape)
{
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    bool compressed = false;
    for (size_t i = 0; i < size; ++i)
    {
        if (!compressed)
        {
            ASSERT_EQ(data[i], reader_data[i]);
        }
    }
}

void Writer(const Params &engineParams)
{
    adios2::ADIOS adios(mpiComm);
    adios2::IO io = adios.DeclareIO("io");
    io.SetEngine("ssc");
    io.SetParameters(engineParams);
    auto varInts = io.DefineVariable<int>("varInts", shape, start, count);
    adios2::Engine engine = io.Open("stream", adios2::Mode::Write);
    for (size_t i = 0; i < steps; ++i)
    {
        engine.BeginStep();
        engine.Put(varInts, writer_data.data(), adios2::Mode::Sync);
        engine.EndStep();
    }
    engine.Close();
}

void Reader(const Params &engineParams)
{
    adios2::ADIOS adios(mpiComm);
    adios2::IO io = adios.DeclareIO("io");
    io.SetEngine("ssc");
    io.SetParameters(engineParams);
    adios2::Engine engine = io.Open("stream", adios2::Mode::Read);
    std::vector<int> myInts = reader_data;
    while (true)
    {
        adios2::StepStatus status = engine.BeginStep(StepMode::Read, 5);
        if (status == adios2::StepStatus::OK)
        {
            const auto &vars = io.AvailableVariables();
            ASSERT_EQ(vars.size(), 1);
            size_t currentStep = engine.CurrentStep();
            adios2::Variable<int> varInts = io.InquireVariable<int>("varInts");
            varInts.SetSelection({start, count});
            varInts.SetMemorySelection({memStart, memCount});
            engine.Get(varInts, myInts.data(), adios2::Mode::Sync);
            VerifyData(myInts.data(), currentStep, memStart, memCount, shape);
            engine.EndStep();
        }
        else if (status == adios2::StepStatus::EndOfStream)
        {
            break;
        }
        else if (status == adios2::StepStatus::NotReady)
        {
            continue;
        }
    }
    engine.Close();
}

TEST_F(SscEngineTest, TestSsc3DMemSelect)
{
    {
        adios2::Params engineParams = {{"Verbose", "0"}};
        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        int mpiGroup = worldRank / (worldSize / 2);
        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);
        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);
        if (mpiGroup == 0)
        {
            Writer(engineParams);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (mpiGroup == 1)
        {
            Reader(engineParams);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    {
        adios2::Params engineParams = {{"Verbose", "0"}, {"EngineMode", "naive"}};
        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        int mpiGroup = worldRank / (worldSize / 2);
        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);
        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);
        if (mpiGroup == 0)
        {
            Writer(engineParams);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (mpiGroup == 1)
        {
            Reader(engineParams);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int worldRank, worldSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
