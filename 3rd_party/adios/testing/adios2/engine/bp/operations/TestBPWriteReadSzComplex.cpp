/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestBPWriteReadSzComplex.cpp
 *
 *  Created on: Jun 20, 2021
 *      Author: Jason Wang
 */

#include <adios2.h>
#include <cmath>
#include <gtest/gtest.h>
#include <numeric>
#include <thread>

using namespace adios2;

std::string ioName = "TestIO";
std::string fileName = "TestBPWriteReadSzComplex";
std::string accuracy = "0.01";

class BPEngineTest : public ::testing::Test
{
public:
    BPEngineTest() = default;
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
void GenData(std::vector<std::complex<T>> &data, const size_t step, const Dims &start,
             const Dims &count, const Dims &shape)
{
    if (start.size() == 2)
    {
        for (size_t i = 0; i < count[0]; ++i)
        {
            for (size_t j = 0; j < count[1]; ++j)
            {
                data[i * count[1] + j] = {static_cast<T>((i + start[1]) * shape[1] + j + start[0] +
                                                         std::stof(accuracy) * 0.01 * (T)step),
                                          static_cast<T>((i + start[1]) * shape[1] + j + start[0] +
                                                         std::stof(accuracy) * 0.01 * (T)step) +
                                              1};
            }
        }
    }
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
                data[i * count[1] + j] = (i + start[1]) * shape[1] + j + start[0] +
                                         std::stof(accuracy) * 0.00001 * (T)step;
            }
        }
    }
}

template <class T>
void VerifyData(const std::complex<T> *data, size_t step, const Dims &start, const Dims &count,
                const Dims &shape, bool &compressed)
{
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    std::vector<std::complex<T>> tmpdata(size);
    GenData(tmpdata, step, start, count, shape);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_LT(std::abs((double)data[i].real() - (double)tmpdata[i].real()),
                  std::stof(accuracy));
        ASSERT_LT(std::abs((double)data[i].imag() - (double)tmpdata[i].imag()),
                  std::stof(accuracy));
        if (data[i].real() != tmpdata[i].real() || data[i].imag() != tmpdata[i].imag())
        {
            compressed = true;
        }
    }
}

template <class T>
void VerifyData(const T *data, size_t step, const Dims &start, const Dims &count, const Dims &shape,
                bool &compressed)
{
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    std::vector<T> tmpdata(size);
    GenData(tmpdata, step, start, count, shape);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_LT(std::abs((double)(data[i] - tmpdata[i])), std::stof(accuracy));
        if (data[i] != tmpdata[i])
        {
            compressed = true;
        }
    }
}

void Writer(const Dims &shape, const Dims &start, const Dims &count, const size_t steps)
{
    size_t datasize = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    fileName = "TestBPWriteReadSzComplex_MPI";
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO(ioName);
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
    auto bpChars = io.DefineVariable<char>("bpChars", shape, start, count);
    auto bpUChars = io.DefineVariable<unsigned char>("bpUChars", shape, start, count);
    auto bpShorts = io.DefineVariable<short>("bpShorts", shape, start, count);
    auto bpUShorts = io.DefineVariable<unsigned short>("bpUShorts", shape, start, count);
    auto bpInts = io.DefineVariable<int>("bpInts", shape, start, count);
    auto bpUInts = io.DefineVariable<unsigned int>("bpUInts", shape, start, count);
    auto bpFloats = io.DefineVariable<float>("bpFloats", shape, start, count);
    adios2::Operator compressor = adios.DefineOperator("szCompressor", adios2::ops::LossySZ);
    bpFloats.AddOperation(compressor, {{"accuracy", accuracy}});
    auto bpDoubles = io.DefineVariable<double>("bpDoubles", shape, start, count);
    bpDoubles.AddOperation(compressor, {{"accuracy", accuracy}});
    auto bpComplexes = io.DefineVariable<std::complex<float>>("bpComplexes", shape, start, count);
    bpComplexes.AddOperation(compressor, {{"accuracy", accuracy}});
    auto bpDComplexes =
        io.DefineVariable<std::complex<double>>("bpDComplexes", shape, start, count);
    bpDComplexes.AddOperation(compressor, {{"accuracy", accuracy}});
    io.DefineAttribute<int>("AttInt", 110);
    adios2::Engine writerEngine = io.Open(fileName, adios2::Mode::Write);
    for (int i = 0; i < static_cast<int>(steps); ++i)
    {
        writerEngine.BeginStep();
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
        writerEngine.Put(bpChars, myChars.data(), adios2::Mode::Sync);
        writerEngine.Put(bpUChars, myUChars.data(), adios2::Mode::Sync);
        writerEngine.Put(bpShorts, myShorts.data(), adios2::Mode::Sync);
        writerEngine.Put(bpUShorts, myUShorts.data(), adios2::Mode::Sync);
        writerEngine.Put(bpInts, myInts.data(), adios2::Mode::Sync);
        writerEngine.Put(bpUInts, myUInts.data(), adios2::Mode::Sync);
        writerEngine.Put(bpFloats, myFloats.data(), adios2::Mode::Sync);
        writerEngine.Put(bpDoubles, myDoubles.data(), adios2::Mode::Sync);
        writerEngine.Put(bpComplexes, myComplexes.data(), adios2::Mode::Sync);
        writerEngine.Put(bpDComplexes, myDComplexes.data(), adios2::Mode::Sync);
        writerEngine.EndStep();
    }
    writerEngine.Close();
}

void Reader(const Dims &shape, const Dims &start, const Dims &count, const size_t steps)
{
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO(ioName);
    adios2::Engine readerEngine = io.Open(fileName, adios2::Mode::Read);

    size_t datasize = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
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
    size_t currentStep;
    bool floatCompressed = false;
    bool doubleCompressed = false;
    bool complexCompressed = false;
    bool dcomplexCompressed = false;
    bool otherCompressed = false;
    while (true)
    {
        adios2::StepStatus status = readerEngine.BeginStep();
        if (status == adios2::StepStatus::OK)
        {
            const auto &vars = io.AvailableVariables();
            ASSERT_EQ(vars.size(), 10);
            currentStep = readerEngine.CurrentStep();
            GenData(myChars, currentStep, start, count, shape);
            GenData(myUChars, currentStep, start, count, shape);
            GenData(myShorts, currentStep, start, count, shape);
            GenData(myUShorts, currentStep, start, count, shape);
            GenData(myInts, currentStep, start, count, shape);
            GenData(myUInts, currentStep, start, count, shape);
            GenData(myFloats, currentStep, start, count, shape);
            GenData(myDoubles, currentStep, start, count, shape);
            GenData(myComplexes, currentStep, start, count, shape);
            GenData(myDComplexes, currentStep, start, count, shape);
            adios2::Variable<char> bpChars = io.InquireVariable<char>("bpChars");
            adios2::Variable<unsigned char> bpUChars =
                io.InquireVariable<unsigned char>("bpUChars");
            adios2::Variable<short> bpShorts = io.InquireVariable<short>("bpShorts");
            adios2::Variable<unsigned short> bpUShorts =
                io.InquireVariable<unsigned short>("bpUShorts");
            adios2::Variable<int> bpInts = io.InquireVariable<int>("bpInts");
            adios2::Variable<unsigned int> bpUInts = io.InquireVariable<unsigned int>("bpUInts");
            adios2::Variable<float> bpFloats = io.InquireVariable<float>("bpFloats");
            adios2::Variable<double> bpDoubles = io.InquireVariable<double>("bpDoubles");
            adios2::Variable<std::complex<float>> bpComplexes =
                io.InquireVariable<std::complex<float>>("bpComplexes");
            adios2::Variable<std::complex<double>> bpDComplexes =
                io.InquireVariable<std::complex<double>>("bpDComplexes");
            auto charsBlocksInfo = readerEngine.AllStepsBlocksInfo(bpChars);

            bpChars.SetSelection({start, count});
            bpUChars.SetSelection({start, count});
            bpShorts.SetSelection({start, count});
            bpUShorts.SetSelection({start, count});
            bpInts.SetSelection({start, count});
            bpUInts.SetSelection({start, count});
            bpFloats.SetSelection({start, count});
            bpDoubles.SetSelection({start, count});
            bpComplexes.SetSelection({start, count});
            bpDComplexes.SetSelection({start, count});

            readerEngine.Get(bpChars, myChars.data(), adios2::Mode::Sync);
            readerEngine.Get(bpUChars, myUChars.data(), adios2::Mode::Sync);
            readerEngine.Get(bpShorts, myShorts.data(), adios2::Mode::Sync);
            readerEngine.Get(bpUShorts, myUShorts.data(), adios2::Mode::Sync);
            readerEngine.Get(bpInts, myInts.data(), adios2::Mode::Sync);
            readerEngine.Get(bpUInts, myUInts.data(), adios2::Mode::Sync);
            readerEngine.Get(bpFloats, myFloats.data(), adios2::Mode::Sync);
            readerEngine.Get(bpDoubles, myDoubles.data(), adios2::Mode::Sync);
            readerEngine.Get(bpComplexes, myComplexes.data(), adios2::Mode::Sync);
            readerEngine.Get(bpDComplexes, myDComplexes.data(), adios2::Mode::Sync);

            VerifyData(myChars.data(), currentStep, start, count, shape, otherCompressed);
            VerifyData(myUChars.data(), currentStep, start, count, shape, otherCompressed);
            VerifyData(myShorts.data(), currentStep, start, count, shape, otherCompressed);
            VerifyData(myUShorts.data(), currentStep, start, count, shape, otherCompressed);
            VerifyData(myInts.data(), currentStep, start, count, shape, otherCompressed);
            VerifyData(myUInts.data(), currentStep, start, count, shape, otherCompressed);
            VerifyData(myFloats.data(), currentStep, start, count, shape, floatCompressed);
            VerifyData(myDoubles.data(), currentStep, start, count, shape, doubleCompressed);
            VerifyData(myComplexes.data(), currentStep, start, count, shape, complexCompressed);
            VerifyData(myDComplexes.data(), currentStep, start, count, shape, dcomplexCompressed);
            readerEngine.EndStep();
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

    auto attInt = io.InquireAttribute<int>("AttInt");
    ASSERT_EQ(110, attInt.Data()[0]);
    ASSERT_EQ(otherCompressed, false);
    ASSERT_EQ(floatCompressed, true);
    ASSERT_EQ(doubleCompressed, true);
    ASSERT_EQ(complexCompressed, true);
    ASSERT_EQ(dcomplexCompressed, true);
    readerEngine.Close();
}

TEST_F(BPEngineTest, SzComplex)
{

    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    Dims shape = {(size_t)mpiSize, 100};
    Dims start = {(size_t)mpiRank, 0};
    Dims count = {1, 80};
    size_t steps = 500;

    Writer(shape, start, count, steps);
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    Reader(shape, start, count, steps);
}

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif
    return result;
}
