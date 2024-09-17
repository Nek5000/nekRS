/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "TestSscCommon.h"
#include <adios2.h>
#include <gtest/gtest.h>
#include <mpi.h>
#include <numeric>
#include <thread>

using namespace adios2;
int mpiRank = 0;
int mpiSize = 1;
MPI_Comm mpiComm;

class SscEngineTest : public ::testing::Test
{
public:
    SscEngineTest() = default;
};

void Writer(const Dims &shape, const Dims &start, const Dims &count, const size_t steps,
            const adios2::Params &engineParams, const std::string &name)
{
    size_t datasize = std::accumulate(shape.begin(), shape.end(), static_cast<size_t>(1),
                                      std::multiplies<size_t>());
    adios2::ADIOS adios(mpiComm);
    adios2::IO io = adios.DeclareIO("Test");
    io.SetEngine("ssc");
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
    auto varChars = io.DefineVariable<char>("varChars", shape, start, shape);
    auto varUChars = io.DefineVariable<unsigned char>("varUChars", shape, start, shape);
    auto varShorts = io.DefineVariable<short>("varShorts", shape, start, shape);
    auto varUShorts = io.DefineVariable<unsigned short>("varUShorts", shape, start, shape);
    auto varInts = io.DefineVariable<int>("varInts", shape, start, shape);
    auto varUInts = io.DefineVariable<unsigned int>("varUInts", shape, start, shape);
    auto varFloats = io.DefineVariable<float>("varFloats", shape, start, shape);
    auto varDoubles = io.DefineVariable<double>("varDoubles", shape, start, shape);
    auto varComplexes = io.DefineVariable<std::complex<float>>("varComplexes", shape, start, shape);
    auto varDComplexes =
        io.DefineVariable<std::complex<double>>("varDComplexes", shape, start, shape);
    auto varIntScalar = io.DefineVariable<int>("varIntScalar");
    auto varString = io.DefineVariable<std::string>("varString");
    io.DefineAttribute<int>("AttInt", 110);
    adios2::Engine engine = io.Open(name, adios2::Mode::Write);
    engine.LockWriterDefinitions();
    for (size_t i = 0; i < steps; ++i)
    {
        engine.BeginStep();
        if (mpiRank == 0)
        {
            GenData(myChars, i, start, shape, shape);
            GenData(myUChars, i, start, shape, shape);
            GenData(myShorts, i, start, shape, shape);
            GenData(myUShorts, i, start, shape, shape);
            GenData(myInts, i, start, shape, shape);
            GenData(myUInts, i, start, shape, shape);
            GenData(myFloats, i, start, shape, shape);
            GenData(myDoubles, i, start, shape, shape);
            GenData(myComplexes, i, start, shape, shape);
            GenData(myDComplexes, i, start, shape, shape);
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
            engine.Put(varIntScalar, static_cast<int>(i));
            std::string s = "sample string sample string sample string";
            engine.Put(varString, s);
        }
        engine.EndStep();
    }
    engine.Close();
}

void Reader(const Dims &shape, const Dims &start, const Dims &count, const size_t steps,
            const adios2::Params &engineParams, const std::string &name)
{
    adios2::ADIOS adios(mpiComm);
    adios2::IO io = adios.DeclareIO("Test");
    io.SetEngine("ssc");
    io.SetParameters(engineParams);
    adios2::Engine engine = io.Open(name, adios2::Mode::Read);

    size_t datasize = std::accumulate(count.begin(), count.end(), static_cast<size_t>(1),
                                      std::multiplies<size_t>());
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

    engine.LockReaderSelections();

    while (true)
    {
        adios2::StepStatus status = engine.BeginStep(StepMode::Read, 5);
        if (status == adios2::StepStatus::OK)
        {
            auto varIntScalar = io.InquireVariable<int>("varIntScalar");
            auto blocksInfo = engine.BlocksInfo(varIntScalar, engine.CurrentStep());

            for (const auto &bi : blocksInfo)
            {
                ASSERT_EQ(bi.IsValue, true);
                ASSERT_EQ(bi.Value, engine.CurrentStep());
                ASSERT_EQ(varIntScalar.Min(), engine.CurrentStep());
                ASSERT_EQ(varIntScalar.Max(), engine.CurrentStep());
            }

            const auto &vars = io.AvailableVariables();
            ASSERT_EQ(vars.size(), 12);
            size_t currentStep = engine.CurrentStep();
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
            adios2::Variable<std::string> varString = io.InquireVariable<std::string>("varString");

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

            VerifyData(myChars.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myUChars.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myShorts.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myUShorts.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myInts.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myUInts.data(), currentStep, start, count, shape, mpiRank);

            engine.Get(varFloats, myFloats.data(), adios2::Mode::Deferred);
            engine.Get(varDoubles, myDoubles.data(), adios2::Mode::Deferred);
            engine.Get(varComplexes, myComplexes.data(), adios2::Mode::Deferred);
            engine.Get(varDComplexes, myDComplexes.data(), adios2::Mode::Deferred);
            engine.PerformGets();

            VerifyData(myFloats.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myDoubles.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myComplexes.data(), currentStep, start, count, shape, mpiRank);
            VerifyData(myDComplexes.data(), currentStep, start, count, shape, mpiRank);

            std::string s;
            engine.Get(varString, s);
            engine.PerformGets();
            ASSERT_EQ(s, "sample string sample string sample string");
            ASSERT_EQ(varString.Min(), "sample string sample string sample string");
            ASSERT_EQ(varString.Max(), "sample string sample string sample string");

            int i;
            engine.Get(varIntScalar, &i);
            engine.PerformGets();
            ASSERT_EQ(i, currentStep);
            engine.EndStep();
        }
        else if (status == adios2::StepStatus::EndOfStream)
        {
            std::cout << "[Rank " + std::to_string(mpiRank) + "] SscTest reader end of stream!"
                      << std::endl;
            break;
        }
    }
    auto attInt = io.InquireAttribute<int>("AttInt");
    ASSERT_EQ(110, attInt.Data()[0]);
    engine.Close();
}

TEST_F(SscEngineTest, TestSscZeroBlock)
{
    {
        std::string filename = "TestSscZeroBlock";
        adios2::Params engineParams = {{"Verbose", "0"}};
        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        int mpiGroup = worldRank / (worldSize / 2);
        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);
        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);
        Dims shape = {10, (size_t)mpiSize * 2};
        Dims start = {2, (size_t)mpiRank * 2};
        Dims count = {5, 2};
        size_t steps = 100;
        if (mpiGroup == 0)
        {
            Writer(shape, start, count, steps, engineParams, filename);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (mpiGroup == 1)
        {
            Reader(shape, start, count, steps, engineParams, filename);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    {
        std::string filename = "TestSscZeroBlock";
        adios2::Params engineParams = {{"Verbose", "0"}, {"EngineMode", "naive"}};
        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        int mpiGroup = worldRank / (worldSize / 2);
        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);
        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);
        Dims shape = {10, (size_t)mpiSize * 2};
        Dims start = {2, (size_t)mpiRank * 2};
        Dims count = {5, 2};
        size_t steps = 100;
        if (mpiGroup == 0)
        {
            Writer(shape, start, count, steps, engineParams, filename);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (mpiGroup == 1)
        {
            Reader(shape, start, count, steps, engineParams, filename);
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
