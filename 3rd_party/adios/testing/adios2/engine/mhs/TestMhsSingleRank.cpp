/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <adios2.h>
#include <gtest/gtest.h>
#include <numeric>
#include <thread>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

using namespace adios2;

char runMode;

class MhsEngineTest : public ::testing::Test
{
public:
    MhsEngineTest() = default;
};

template <class T>
void GenData(std::complex<T> *data, const size_t row, const Dims &count)
{
    for (size_t i = 0; i < count[1]; ++i)
    {
        for (size_t j = 0; j < count[2]; ++j)
        {
            data[i * count[2] + j] = {static_cast<T>(i * count[2] + j + row),
                                      static_cast<T>(i * count[2])};
        }
    }
}

template <class T>
void GenData(T *data, const size_t row, const Dims &count)
{
    for (size_t i = 0; i < count[1]; ++i)
    {
        for (size_t j = 0; j < count[2]; ++j)
        {
            data[i * count[2] + j] = static_cast<T>(i * count[2] + j + row);
        }
    }
}

template <class T>
void GenData(std::vector<T> &data, const size_t row, const Dims &count)
{
    GenData(data.data(), row, count);
}

template <class T>
void VerifyData(const T *data, const size_t rows, const Dims &shape)
{
    size_t columnSize = 1;
    for (const auto &i : shape)
    {
        columnSize *= i;
    }
    size_t rowSize = 1;
    for (size_t i = 1; i < shape.size(); ++i)
    {
        rowSize *= shape[i];
    }
    std::vector<T> tmpdata(columnSize);
    size_t position = 0;
    for (size_t i = 0; i < rows; ++i)
    {
        GenData(tmpdata.data() + position, i, shape);
        position += rowSize;
    }
    for (size_t i = 0; i < columnSize; ++i)
    {
        ASSERT_EQ(data[i], tmpdata[i]);
    }
}

void Reader(const Dims &shape, const Dims &start, const Dims &count, const size_t rows,
            const adios2::Params &engineParams, const std::string &name)
{
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("ms");
    io.SetEngine("mhs");
    io.SetParameters(engineParams);
    adios2::Engine readerEngine = io.Open(name, adios2::Mode::Read);
    size_t datasize = 1;
    for (const auto &i : shape)
    {
        datasize *= i;
    }
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

    readerEngine.BeginStep();
    const auto &vars = io.AvailableVariables();
    std::cout << "All available variables : ";
    for (const auto &var : vars)
    {
        std::cout << var.first << ", ";
    }
    std::cout << std::endl;
    ASSERT_EQ(vars.size(), 10);

    adios2::Variable<char> bpChars = io.InquireVariable<char>("bpChars");
    adios2::Variable<unsigned char> bpUChars = io.InquireVariable<unsigned char>("bpUChars");
    adios2::Variable<short> bpShorts = io.InquireVariable<short>("bpShorts");
    adios2::Variable<unsigned short> bpUShorts = io.InquireVariable<unsigned short>("bpUShorts");
    adios2::Variable<int> bpInts = io.InquireVariable<int>("bpInts");
    adios2::Variable<unsigned int> bpUInts = io.InquireVariable<unsigned int>("bpUInts");
    adios2::Variable<float> bpFloats = io.InquireVariable<float>("bpFloats");
    adios2::Variable<double> bpDoubles = io.InquireVariable<double>("bpDoubles");
    adios2::Variable<std::complex<float>> bpComplexes =
        io.InquireVariable<std::complex<float>>("bpComplexes");
    adios2::Variable<std::complex<double>> bpDComplexes =
        io.InquireVariable<std::complex<double>>("bpDComplexes");

    bpChars.SetSelection({start, shape});
    bpUChars.SetSelection({start, shape});
    bpShorts.SetSelection({start, shape});
    bpUShorts.SetSelection({start, shape});
    bpInts.SetSelection({start, shape});
    bpUInts.SetSelection({start, shape});
    bpFloats.SetSelection({start, shape});
    bpDoubles.SetSelection({start, shape});
    bpComplexes.SetSelection({start, shape});
    bpDComplexes.SetSelection({start, shape});

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

    VerifyData(myChars.data(), rows, shape);
    VerifyData(myUChars.data(), rows, shape);
    VerifyData(myShorts.data(), rows, shape);
    VerifyData(myUShorts.data(), rows, shape);
    VerifyData(myInts.data(), rows, shape);
    VerifyData(myUInts.data(), rows, shape);
    VerifyData(myFloats.data(), rows, shape);
    VerifyData(myDoubles.data(), rows, shape);
    VerifyData(myComplexes.data(), rows, shape);
    VerifyData(myDComplexes.data(), rows, shape);
    readerEngine.EndStep();
    readerEngine.Close();
}

void Writer(const Dims &shape, const Dims &start, const Dims &count, const size_t rows,
            const adios2::Params &engineParams, const std::string &name)
{
    size_t datasize = 1;
    for (const auto &i : count)
    {
        datasize *= i;
    }
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("ms");
    io.SetEngine("mhs");
    io.SetParameters(engineParams);
    io.AddTransport("sirius", {{"variable", "bpFloats"}});
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
    auto bpDoubles = io.DefineVariable<double>("bpDoubles", shape, start, count);
    auto bpComplexes = io.DefineVariable<std::complex<float>>("bpComplexes", shape, start, count);
    auto bpDComplexes =
        io.DefineVariable<std::complex<double>>("bpDComplexes", shape, start, count);
    adios2::Engine writerEngine = io.Open(name, adios2::Mode::Write);
    writerEngine.BeginStep();
    for (size_t i = 0; i < rows; ++i)
    {
        Dims startRow = start;
        startRow[0] = i;
        bpChars.SetSelection({startRow, count});
        bpUChars.SetSelection({startRow, count});
        bpShorts.SetSelection({startRow, count});
        bpUShorts.SetSelection({startRow, count});
        bpInts.SetSelection({startRow, count});
        bpUInts.SetSelection({startRow, count});
        bpFloats.SetSelection({startRow, count});
        bpDoubles.SetSelection({startRow, count});
        bpComplexes.SetSelection({startRow, count});
        bpDComplexes.SetSelection({startRow, count});
        GenData(myChars, i, count);
        GenData(myUChars, i, count);
        GenData(myShorts, i, count);
        GenData(myUShorts, i, count);
        GenData(myInts, i, count);
        GenData(myUInts, i, count);
        GenData(myFloats, i, count);
        GenData(myDoubles, i, count);
        GenData(myComplexes, i, count);
        GenData(myDComplexes, i, count);
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
    }
    writerEngine.EndStep();
    writerEngine.Close();
}

TEST_F(MhsEngineTest, TestMhsSingleRank)
{
    std::string filename = "TestMhsSingleRank";
    adios2::Params engineParams = {{"Verbose", "0"}, {"Tiers", "1"}};

    size_t rows = 100;
    Dims shape = {rows, 1, 128};
    Dims start = {0, 0, 0};
    Dims count = {1, 1, 128};

    Writer(shape, start, count, rows, engineParams, filename);

    Reader(shape, start, count, rows, engineParams, filename);

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif
    return result;
}
