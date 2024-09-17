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
    size_t datasize = std::accumulate(count.begin(), count.end(), static_cast<size_t>(1),
                                      std::multiplies<size_t>());
    adios2::ADIOS adios(mpiComm);
    adios2::IO io = adios.DeclareIO("WAN");
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
    auto varIntScalar = io.DefineVariable<int>("varIntScalar");
    auto varString = io.DefineVariable<std::string>("varString");
    io.DefineAttribute<int>("AttInt", 110);
    adios2::Engine engine = io.Open(name, adios2::Mode::Write);
    //    engine.LockWriterDefinitions();
    for (int i = 0; i < static_cast<int>(steps); ++i)
    {
        engine.BeginStep();

        Dims vstart = {0, (size_t)mpiRank * 2};
        //        Dims vcount = {(size_t)i+1, 2};
        //        Dims vshape = {(size_t)i+1, (size_t)mpiSize * 2};
        Dims vcount = {10, 2};
        Dims vshape = {10, (size_t)mpiSize * 2};

        varChars.SetShape(vshape);
        varChars.SetSelection({vstart, vcount});
        varUChars.SetShape(vshape);
        varUChars.SetSelection({vstart, vcount});
        varShorts.SetShape(vshape);
        varShorts.SetSelection({vstart, vcount});
        varUShorts.SetShape(vshape);
        varUShorts.SetSelection({vstart, vcount});
        varInts.SetShape(vshape);
        varInts.SetSelection({vstart, vcount});
        varUInts.SetShape(vshape);
        varUInts.SetSelection({vstart, vcount});
        varFloats.SetShape(vshape);
        varFloats.SetSelection({vstart, vcount});
        varDoubles.SetShape(vshape);
        varDoubles.SetSelection({vstart, vcount});
        varComplexes.SetShape(vshape);
        varComplexes.SetSelection({vstart, vcount});
        varDComplexes.SetShape(vshape);
        varDComplexes.SetSelection({vstart, vcount});

        GenData(myChars, i, vstart, vcount, vshape);
        GenData(myUChars, i, vstart, vcount, vshape);
        GenData(myShorts, i, vstart, vcount, vshape);
        GenData(myUShorts, i, vstart, vcount, vshape);
        GenData(myInts, i, vstart, vcount, vshape);
        GenData(myUInts, i, vstart, vcount, vshape);
        GenData(myFloats, i, vstart, vcount, vshape);
        GenData(myDoubles, i, vstart, vcount, vshape);
        GenData(myComplexes, i, vstart, vcount, vshape);
        GenData(myDComplexes, i, vstart, vcount, vshape);

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
        engine.Put(varIntScalar, i);
        std::string s = "sample string sample string sample string";
        engine.Put(varString, s);
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
    //    engine.LockReaderSelections();

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

            auto vshape = varChars.Shape();
            myChars.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                           std::multiplies<size_t>()));
            myUChars.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                            std::multiplies<size_t>()));
            myShorts.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                            std::multiplies<size_t>()));
            myUShorts.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                             std::multiplies<size_t>()));
            myInts.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                          std::multiplies<size_t>()));
            myUInts.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                           std::multiplies<size_t>()));
            myFloats.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                            std::multiplies<size_t>()));
            myDoubles.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                             std::multiplies<size_t>()));
            myComplexes.resize(std::accumulate(vshape.begin(), vshape.end(), static_cast<size_t>(1),
                                               std::multiplies<size_t>()));
            myDComplexes.resize(std::accumulate(vshape.begin(), vshape.end(),
                                                static_cast<size_t>(1), std::multiplies<size_t>()));

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
            std::string s;
            engine.Get(varString, s, adios2::Mode::Sync);
            ASSERT_EQ(s, "sample string sample string sample string");
            ASSERT_EQ(varString.Min(), "sample string sample string sample string");
            ASSERT_EQ(varString.Max(), "sample string sample string sample string");

            int i;
            engine.Get(varIntScalar, &i, adios2::Mode::Sync);
            ASSERT_EQ(i, currentStep);

            VerifyData(myChars.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myUChars.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myShorts.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myUShorts.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myInts.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myUInts.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myFloats.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myDoubles.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myComplexes.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
            VerifyData(myDComplexes.data(), currentStep, {0, 0}, vshape, vshape, mpiRank);
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

TEST_F(SscEngineTest, TestSscVaryingSteps)
{
    {
        std::string filename = "TestSscVaryingSteps";
        adios2::Params engineParams = {};

        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        int mpiGroup = worldRank / (worldSize / 2);
        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);

        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);

        Dims shape = {1, (size_t)mpiSize * 2};
        Dims start = {0, (size_t)mpiRank * 2};
        Dims count = {1, 2};
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
        std::string filename = "TestSscVaryingStepsNaive";
        adios2::Params engineParams = {{"Verbose", "0"}, {"EngineMode", "naive"}};

        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        int mpiGroup = worldRank / (worldSize / 2);
        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);

        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);

        Dims shape = {1, (size_t)mpiSize * 2};
        Dims start = {0, (size_t)mpiRank * 2};
        Dims count = {1, 2};
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
