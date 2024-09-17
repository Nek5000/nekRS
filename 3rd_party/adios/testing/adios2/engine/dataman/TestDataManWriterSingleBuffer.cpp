/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <numeric>
#include <thread>

#include <adios2.h>
#include <adios2/common/ADIOSMacros.h>
#include <adios2/helper/adiosFunctions.h>
#include <gtest/gtest.h>

using namespace adios2;

int mpiRank = 0;
int mpiSize = 1;

size_t print_lines = 0;
size_t to_print_lines = 10;

template <class T>
void GenData(std::vector<std::complex<T>> &data, const size_t step)
{
    for (size_t i = 0; i < data.size(); ++i)
    {
        data[i] = {static_cast<T>(i + mpiRank * 10000 + step * 100),
                   static_cast<T>(i + mpiRank * 10000)};
    }
}

template <class T>
void GenData(std::vector<T> &data, const size_t step)
{
    for (size_t i = 0; i < data.size(); ++i)
    {
        data[i] = i + mpiRank * 10000 + step * 100;
    }
}

template <class T>
void PrintData(const T *data, const size_t size, const size_t step)
{
    std::cout << "Rank: " << mpiRank << " Step: " << step << " [";
    size_t printsize = 32;
    if (size < printsize)
    {
        printsize = size;
    }
    for (size_t i = 0; i < printsize; ++i)
    {
        std::cout << data[i] << " ";
    }
    std::cout << "]" << std::endl;
}

template <class T>
void VerifyData(const std::complex<T> *data, const size_t size, size_t step)
{
    std::vector<std::complex<T>> tmpdata(size);
    GenData(tmpdata, step);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_EQ(data[i], tmpdata[i]);
    }
}

template <class T>
void VerifyData(const T *data, const size_t size, size_t step)
{
    std::vector<T> tmpdata(size);
    GenData(tmpdata, step);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_EQ(data[i], tmpdata[i]);
    }
}

template <class T>
void VerifyData(const std::vector<T> &data, const size_t step)
{
    VerifyData(data.data(), data.size(), step);
}

void DataManWriter(const Dims &shape, const Dims &start, const Dims &count, const size_t steps,
                   const adios2::Params &engineParams)
{
    size_t datasize = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("WAN");
    io.SetEngine("DataMan");
    io.SetParameters(engineParams);
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
    auto varUInt64s = io.DefineVariable<uint64_t>("varUInt64s");
    io.DefineAttribute<int>("AttInt", 110);
    adios2::Engine engine = io.Open("stream", adios2::Mode::Write);
    for (uint64_t i = 0; i < steps; ++i)
    {
        engine.BeginStep();
        GenData(myChars, i);
        GenData(myUChars, i);
        GenData(myShorts, i);
        GenData(myUShorts, i);
        GenData(myInts, i);
        GenData(myUInts, i);
        GenData(myFloats, i);
        GenData(myDoubles, i);
        GenData(myComplexes, i);
        GenData(myDComplexes, i);
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
        engine.Put(varUInt64s, i);
        engine.EndStep();
    }
    engine.Close();
}

void DataManReader(const Dims &shape, const Dims &start, const Dims &count, const size_t steps,
                   const adios2::Params &engineParams)
{
    size_t datasize = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("WAN");
    io.SetEngine("DataMan");
    io.SetParameters(engineParams);
    adios2::Engine engine = io.Open("stream", adios2::Mode::Read);

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
        adios2::StepStatus status = engine.BeginStep(StepMode::Read, 5);
        if (status == adios2::StepStatus::OK)
        {
            received_steps = true;
            const auto &vars = io.AvailableVariables();
            ASSERT_EQ(vars.size(), 11);
            currentStep = engine.CurrentStep();
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
            adios2::Variable<uint64_t> varUInt64s = io.InquireVariable<uint64_t>("varUInt64s");
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
            uint64_t stepValue;
            engine.Get(varUInt64s, &stepValue, adios2::Mode::Sync);
            ASSERT_EQ(currentStep, stepValue);
            VerifyData(myChars, currentStep);
            VerifyData(myUChars, currentStep);
            VerifyData(myShorts, currentStep);
            VerifyData(myUShorts, currentStep);
            VerifyData(myInts, currentStep);
            VerifyData(myUInts, currentStep);
            VerifyData(myFloats, currentStep);
            VerifyData(myDoubles, currentStep);
            VerifyData(myComplexes, currentStep);
            VerifyData(myDComplexes, currentStep);
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
}

class DataManEngineTest : public ::testing::Test
{
public:
    DataManEngineTest() = default;
};

#ifdef ADIOS2_HAVE_ZEROMQ
TEST_F(DataManEngineTest, WriterSingleBuffer)
{
    // set parameters
    Dims shape = {10};
    Dims start = {0};
    Dims count = {10};
    size_t steps = 5000;

    // run workflow
    adios2::Params readerEngineParams = {{"IPAddress", "127.0.0.1"}, {"Port", "12400"}};
    auto r = std::thread(DataManReader, shape, start, count, steps, readerEngineParams);
    adios2::Params writerEngineParams = {
        {"IPAddress", "127.0.0.1"}, {"Port", "12400"}, {"DoubleBuffer", "false"}};
    auto w = std::thread(DataManWriter, shape, start, count, steps, writerEngineParams);
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
