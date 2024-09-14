/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "TestMhsCommon.h"
#include <adios2.h>
#include <gtest/gtest.h>
#include <mpi.h>
#include <numeric>
#include <thread>

using namespace adios2;
int mpiRank = 0;
int mpiSize = 1;

char runMode;

class MhsEngineTest : public ::testing::Test
{
public:
    MhsEngineTest() = default;
};

void Reader(const Dims &shape, const Dims &start, const Dims &count, const size_t rows,
            const adios2::Params &engineParams, const std::string &name)
{

    if (mpiRank != 0)
    {
        return;
    }

    adios2::ADIOS adios(MPI_COMM_SELF);
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

    for (size_t step = 0; step < 10; step++)
    {
        readerEngine.BeginStep();

        const auto &vars = io.AvailableVariables();
        std::cout << "All available variables : ";
        for (const auto &var : vars)
        {
            std::cout << var.first << ", ";
        }
        std::cout << std::endl;

        adios2::Variable<char> bpChars = io.InquireVariable<char>("bpChars");
        adios2::Variable<unsigned char> bpUChars = io.InquireVariable<unsigned char>("bpUChars");
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
        VerifyData(myChars.data(), step, {0, 0, 0}, shape, shape, "bpChars");

        readerEngine.Get(bpUChars, myUChars.data(), adios2::Mode::Sync);
        VerifyData(myUChars.data(), step, {0, 0, 0}, shape, shape, "bpUChars");

        readerEngine.Get(bpShorts, myShorts.data(), adios2::Mode::Sync);
        VerifyData(myShorts.data(), step, {0, 0, 0}, shape, shape, "bpShorts");

        readerEngine.Get(bpUShorts, myUShorts.data(), adios2::Mode::Sync);
        VerifyData(myUShorts.data(), step, {0, 0, 0}, shape, shape, "bpUShorts");

        readerEngine.Get(bpInts, myInts.data(), adios2::Mode::Sync);
        VerifyData(myInts.data(), step, {0, 0, 0}, shape, shape, "bpInts");

        readerEngine.Get(bpUInts, myUInts.data(), adios2::Mode::Sync);
        VerifyData(myUInts.data(), step, {0, 0, 0}, shape, shape, "bpUInts");

        readerEngine.Get(bpFloats, myFloats.data(), adios2::Mode::Sync);
        VerifyData(myFloats.data(), step, {0, 0, 0}, shape, shape, "bpFloats");

        readerEngine.Get(bpDoubles, myDoubles.data(), adios2::Mode::Sync);
        VerifyData(myDoubles.data(), step, {0, 0, 0}, shape, shape, "bpDoubles");

        readerEngine.Get(bpComplexes, myComplexes.data(), adios2::Mode::Sync);
        VerifyData(myComplexes.data(), step, {0, 0, 0}, shape, shape, "bpComplexes");

        readerEngine.Get(bpDComplexes, myDComplexes.data(), adios2::Mode::Sync);
        VerifyData(myDComplexes.data(), step, {0, 0, 0}, shape, shape, "bpDComplexes");

        readerEngine.EndStep();
    }
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
    adios2::ADIOS adios(MPI_COMM_WORLD);
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

    for (size_t step = 0; step < 10; step++)
    {
        writerEngine.BeginStep();
        for (int i = mpiRank; i < static_cast<int>(rows); i += mpiSize)
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
            GenData(myChars, step, startRow, count, shape);
            GenData(myUChars, step, startRow, count, shape);
            GenData(myShorts, step, startRow, count, shape);
            GenData(myUShorts, step, startRow, count, shape);
            GenData(myInts, step, startRow, count, shape);
            GenData(myUInts, step, startRow, count, shape);
            GenData(myFloats, step, startRow, count, shape);
            GenData(myDoubles, step, startRow, count, shape);
            GenData(myComplexes, step, startRow, count, shape);
            GenData(myDComplexes, step, startRow, count, shape);
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
    }
    writerEngine.Close();
}

TEST_F(MhsEngineTest, TestMhsMultiRank)
{
    std::string filename = "TestMhsMultiRank";
    adios2::Params engineParams = {{"Verbose", "0"}, {"Tiers", "4"}};

    size_t rows = 32;
    Dims shape = {rows, 8, 16};
    Dims start = {0, 0, 0};
    Dims count = {1, 8, 16};

    Writer(shape, start, count, rows, engineParams, filename);
    MPI_Barrier(MPI_COMM_WORLD);

    Reader(shape, start, count, rows, engineParams, filename);
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
