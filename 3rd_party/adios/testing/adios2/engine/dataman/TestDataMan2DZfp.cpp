/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestDataMan2DZfp.cpp
 *
 *  Created on: Nov 24, 2020
 *      Author: Jason Wang
 */

#include <adios2.h>
#include <cmath>
#include <gtest/gtest.h>
#include <numeric>
#include <thread>

using namespace adios2;
size_t print_lines = 0;

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

template <class T>
void GenData(std::vector<T> &data, const size_t step, const Dims &start, const Dims &count,
             const Dims &shape)
{
    if (start.size() == 2)
    {
        for (size_t i = 0; i < count[0]; ++i)
        {
            for (size_t j = 0; j < count[1]; ++j)
            {
                data[i * count[1] + j] = (i + start[1]) * shape[1] + j + start[0] + 0.01;
            }
        }
    }
}

template <class T>
void VerifyData(const std::complex<T> *data, size_t step, const Dims &start, const Dims &count,
                const Dims &shape)
{
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    std::vector<std::complex<T>> tmpdata(size);
    GenData(tmpdata, step, start, count, shape);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_EQ(abs(data[i].real() - tmpdata[i].real()) < 0.01, true);
        ASSERT_EQ(abs(data[i].imag() - tmpdata[i].imag()) < 0.01, true);
    }
}

template <class T>
void VerifyData(const T *data, size_t step, const Dims &start, const Dims &count, const Dims &shape)
{
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    std::vector<T> tmpdata(size);
    GenData(tmpdata, step, start, count, shape);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_EQ(abs((double)(data[i] - tmpdata[i])) < 0.01, true);
    }
}

void DataManWriterP2PMemSelect(const Dims &shape, const Dims &start, const Dims &count,
                               const size_t steps, const adios2::Params &engineParams)
{
    size_t datasize = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("WAN");
    io.SetEngine("DataMan");
    io.SetParameters(engineParams);
    io.AddOperation("varFloats", "zfp", {{"accuracy", "0.01"}});
    io.AddOperation("varDoubles", "zfp", {{"accuracy", "0.1"}});
    io.AddOperation("varComplexes", "zfp", {{"accuracy", "0.1"}});
    io.AddOperation("varDComplexes", "zfp", {{"accuracy", "0.1"}});
    std::vector<char> myChars(datasize);
    std::vector<unsigned char> myUChars(datasize);
    std::vector<short> myShorts(datasize);
    std::vector<unsigned short> myUShorts(datasize);
    std::vector<int> myInts(datasize);
    std::vector<unsigned int> myUInts(datasize);
    std::vector<float> myFloats(datasize);
    std::vector<double> myDoubles(datasize);
    std::vector<std::complex<float>> myComplexes(datasize);
    std::vector<std::complex<double>> myDComplexes(datasize);
    auto varChars = io.DefineVariable<char>("varChars", shape, start, count);
    auto varUChars = io.DefineVariable<unsigned char>("varUChars", shape, start, count);
    auto varShorts = io.DefineVariable<short>("varShorts", shape, start, count);
    auto varUShorts = io.DefineVariable<unsigned short>("varUShorts", shape, start, count);
    auto varInts = io.DefineVariable<int>("varInts", shape, start, count);
    auto varUInts = io.DefineVariable<unsigned int>("varUInts", shape, start, count);
    auto varFloats = io.DefineVariable<float>("varFloats", shape, start, count);
    auto varDoubles = io.DefineVariable<double>("varDoubles", shape, start, count);
    auto varComplexes = io.DefineVariable<std::complex<float>>("varComplexes", shape, start, count);
    auto varDComplexes =
        io.DefineVariable<std::complex<double>>("varDComplexes", shape, start, count);
    io.DefineAttribute<int>("AttInt", 110);
    adios2::Engine engine = io.Open("stream", adios2::Mode::Write);
    for (size_t i = 0; i < steps; ++i)
    {
        engine.BeginStep();
        GenData(myChars, i, start, count, shape);
        GenData(myUChars, i, start, count, shape);
        GenData(myShorts, i, start, count, shape);
        GenData(myUShorts, i, start, count, shape);
        GenData(myInts, i, start, count, shape);
        GenData(myUInts, i, start, count, shape);
        GenData(myFloats, i, start, count, shape);
        GenData(myDoubles, i, start, count, shape);
        GenData(myComplexes, i, start, count, shape);
        GenData(myDComplexes, i, start, count, shape);
        engine.Put(varChars, myChars.data(), adios2::Mode::Sync);
        engine.Put(varUChars, myUChars.data(), adios2::Mode::Sync);
        engine.Put(varShorts, myShorts.data(), adios2::Mode::Sync);
        engine.Put(varUShorts, myUShorts.data(), adios2::Mode::Sync);
        engine.Put(varInts, myInts.data(), adios2::Mode::Sync);
        engine.Put(varUInts, myUInts.data(), adios2::Mode::Sync);
        engine.Put(varFloats, myFloats.data(), adios2::Mode::Sync);
        engine.Put(varDoubles, myDoubles.data(), adios2::Mode::Sync);
        engine.Put(varComplexes, myComplexes.data(), adios2::Mode::Sync);
        engine.Put(varDComplexes, myDComplexes.data(), adios2::Mode::Sync);
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

    size_t datasize =
        std::accumulate(memCount.begin(), memCount.end(), 1, std::multiplies<size_t>());
    std::vector<char> myChars(datasize);
    std::vector<unsigned char> myUChars(datasize);
    std::vector<short> myShorts(datasize);
    std::vector<unsigned short> myUShorts(datasize);
    std::vector<int> myInts(datasize);
    std::vector<unsigned int> myUInts(datasize);
    std::vector<float> myFloats(datasize);
    std::vector<double> myDoubles(datasize);
    std::vector<std::complex<float>> myComplexes(datasize);
    std::vector<std::complex<double>> myDComplexes(datasize);
    bool received_steps = false;
    size_t currentStep;
    while (true)
    {
        adios2::StepStatus status = engine.BeginStep();
        if (status == adios2::StepStatus::OK)
        {
            received_steps = true;
            const auto &vars = io.AvailableVariables();
            ASSERT_EQ(vars.size(), 10);
            currentStep = engine.CurrentStep();
            GenData(myChars, currentStep, memStart, memCount, shape);
            GenData(myUChars, currentStep, memStart, memCount, shape);
            GenData(myShorts, currentStep, memStart, memCount, shape);
            GenData(myUShorts, currentStep, memStart, memCount, shape);
            GenData(myInts, currentStep, memStart, memCount, shape);
            GenData(myUInts, currentStep, memStart, memCount, shape);
            GenData(myFloats, currentStep, memStart, memCount, shape);
            GenData(myDoubles, currentStep, memStart, memCount, shape);
            GenData(myComplexes, currentStep, memStart, memCount, shape);
            GenData(myDComplexes, currentStep, memStart, memCount, shape);
            adios2::Variable<char> varChars = io.InquireVariable<char>("varChars");
            adios2::Variable<unsigned char> varUChars =
                io.InquireVariable<unsigned char>("varUChars");
            adios2::Variable<short> varShorts = io.InquireVariable<short>("varShorts");
            adios2::Variable<unsigned short> varUShorts =
                io.InquireVariable<unsigned short>("varUShorts");
            adios2::Variable<int> varInts = io.InquireVariable<int>("varInts");
            adios2::Variable<unsigned int> varUInts = io.InquireVariable<unsigned int>("varUInts");
            adios2::Variable<float> varFloats = io.InquireVariable<float>("varFloats");
            adios2::Variable<double> varDoubles = io.InquireVariable<double>("varDoubles");
            adios2::Variable<std::complex<float>> varComplexes =
                io.InquireVariable<std::complex<float>>("varComplexes");
            adios2::Variable<std::complex<double>> varDComplexes =
                io.InquireVariable<std::complex<double>>("varDComplexes");
            auto charsBlocksInfo = engine.AllStepsBlocksInfo(varChars);

            varChars.SetSelection({start, count});
            varUChars.SetSelection({start, count});
            varShorts.SetSelection({start, count});
            varUShorts.SetSelection({start, count});
            varInts.SetSelection({start, count});
            varUInts.SetSelection({start, count});
            varFloats.SetSelection({start, count});
            varDoubles.SetSelection({start, count});
            varComplexes.SetSelection({start, count});
            varDComplexes.SetSelection({start, count});

            varChars.SetMemorySelection({memStart, memCount});
            varUChars.SetMemorySelection({memStart, memCount});
            varShorts.SetMemorySelection({memStart, memCount});
            varUShorts.SetMemorySelection({memStart, memCount});
            varInts.SetMemorySelection({memStart, memCount});
            varUInts.SetMemorySelection({memStart, memCount});
            varFloats.SetMemorySelection({memStart, memCount});
            varDoubles.SetMemorySelection({memStart, memCount});
            varComplexes.SetMemorySelection({memStart, memCount});
            varDComplexes.SetMemorySelection({memStart, memCount});

            engine.Get(varChars, myChars.data(), adios2::Mode::Sync);
            engine.Get(varUChars, myUChars.data(), adios2::Mode::Sync);
            engine.Get(varShorts, myShorts.data(), adios2::Mode::Sync);
            engine.Get(varUShorts, myUShorts.data(), adios2::Mode::Sync);
            engine.Get(varInts, myInts.data(), adios2::Mode::Sync);
            engine.Get(varUInts, myUInts.data(), adios2::Mode::Sync);
            engine.Get(varFloats, myFloats.data(), adios2::Mode::Sync);
            engine.Get(varDoubles, myDoubles.data(), adios2::Mode::Sync);
            engine.Get(varComplexes, myComplexes.data(), adios2::Mode::Sync);
            engine.Get(varDComplexes, myDComplexes.data(), adios2::Mode::Sync);
            VerifyData(myChars.data(), currentStep, memStart, memCount, shape);
            VerifyData(myUChars.data(), currentStep, memStart, memCount, shape);
            VerifyData(myShorts.data(), currentStep, memStart, memCount, shape);
            VerifyData(myUShorts.data(), currentStep, memStart, memCount, shape);
            VerifyData(myInts.data(), currentStep, memStart, memCount, shape);
            VerifyData(myUInts.data(), currentStep, memStart, memCount, shape);
            VerifyData(myFloats.data(), currentStep, memStart, memCount, shape);
            VerifyData(myDoubles.data(), currentStep, memStart, memCount, shape);
            VerifyData(myComplexes.data(), currentStep, memStart, memCount, shape);
            VerifyData(myDComplexes.data(), currentStep, memStart, memCount, shape);
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
    if (received_steps)
    {
        auto attInt = io.InquireAttribute<int>("AttInt");
        ASSERT_EQ(110, attInt.Data()[0]);
        ASSERT_NE(111, attInt.Data()[0]);
    }
    engine.Close();
    print_lines = 0;
}

#ifdef ADIOS2_HAVE_ZEROMQ
TEST_F(DataManEngineTest, 2D_Zfp)
{
    // set parameters
    Dims shape = {10, 10};
    Dims start = {2, 2};
    Dims count = {5, 5};
    Dims memstart = start;
    Dims memcount = count;
    memstart = {1, 1};
    memcount = {7, 9};

    size_t steps = 5000;
    adios2::Params engineParams = {{"IPAddress", "127.0.0.1"}, {"Port", "12340"}, {"Verbose", "0"}};

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
