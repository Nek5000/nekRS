/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestDataMan3DMemSelect.cpp
 *
 *  Created on: Nov 24, 2018
 *      Author: Jason Wang
 */

#include <adios2.h>
#include <gtest/gtest.h>
#include <numeric>
#include <thread>

using namespace adios2;
size_t print_lines = 0;

Dims shape = {4, 4, 4};
std::vector<int> global_data = {
    0,   1,   2,   3,   10,  11,  12,  13,  20,  21,  22,  23,  30,  31,  32,  33,

    100, 101, 102, 103, 110, 111, 112, 113, 120, 121, 122, 123, 130, 131, 132, 133,

    200, 201, 202, 203, 210, 211, 212, 213, 220, 221, 222, 223, 230, 231, 232, 233,

    300, 301, 302, 303, 310, 311, 312, 313, 320, 321, 322, 323, 330, 331, 332, 333};

Dims start = {1, 2, 1};
Dims count = {2, 1, 2};
std::vector<int> writer_data = {121, 122, 221, 222};

Dims memstart = {0, 1, 1};
Dims memcount = {3, 3, 3};
std::vector<int> reader_data = {11,  12,  13,  21,  22,  23,  31,  32,  33,

                                111, 112, 113, 121, 122, 123, 131, 132, 133,

                                211, 212, 213, 221, 222, 223, 231, 232, 233};

class DataManEngineTest : public ::testing::Test
{
public:
    DataManEngineTest() = default;
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

void DataManWriterP2PMemSelect(const Dims &shape, const Dims &start, const Dims &count,
                               const size_t steps, const adios2::Params &engineParams)
{
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("WAN");
    io.SetEngine("DataMan");
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

void DataManReaderP2PMemSelect(const Dims &shape, const Dims &start, const Dims &count,
                               const Dims &memStart, const Dims &memCount, const size_t steps,
                               const adios2::Params &engineParams)
{
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("WAN");
    io.SetEngine("DataMan");
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
    print_lines = 0;
}

#ifdef ADIOS2_HAVE_ZEROMQ
TEST_F(DataManEngineTest, 3D_MemSelect)
{

    size_t steps = 3500;
    adios2::Params engineParams = {{"IPAddress", "127.0.0.1"}, {"Port", "12350"}};
    // run workflow

    auto r = std::thread(DataManReaderP2PMemSelect, shape, start, count, memstart, memcount, steps,
                         engineParams);

    auto w = std::thread(DataManWriterP2PMemSelect, shape, start, count, steps, engineParams);

    w.join();

    r.join();
}

#endif // ZEROMQ

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    return result;
}
