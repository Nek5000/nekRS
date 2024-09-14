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
int worldRank = 0;
int worldSize = 2;
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

    struct particle
    {
        int8_t a;
        int b[4];
    };
    auto particleDef = io.DefineStruct("particle", sizeof(particle));
    particleDef.AddField("a", offsetof(struct particle, a),
                         adios2::DataType::Int8); // simple field
    particleDef.AddField("b", offsetof(struct particle, b), adios2::DataType::Int32,
                         sizeof(int)); // static array
    auto varStruct = io.DefineStructVariable("particles", particleDef, shape, start, count);
    EXPECT_THROW(particleDef.AddField("c", 12, adios2::DataType::Int32), std::runtime_error);
    std::vector<particle> myParticles(datasize);
    for (size_t i = 0; i < datasize; ++i)
    {
        myParticles[i].a = 5;
        myParticles[i].b[0] = 0;
        myParticles[i].b[1] = 10;
        myParticles[i].b[2] = 20;
        myParticles[i].b[3] = 30;
    }

    io.DefineAttribute<int>("AttInt", 110);
    adios2::Engine engine = io.Open(name, adios2::Mode::Write);
    engine.LockWriterDefinitions();
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
        engine.Put(varStruct, myParticles.data());
        engine.Put(varIntScalar, static_cast<int>(i));
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
    struct particle
    {
        int8_t a;
        int b[4];
    };
    std::vector<particle> myParticles(datasize);

    engine.LockReaderSelections();

    auto particleDef1 = io.DefineStruct("particle1", sizeof(particle));
    particleDef1.AddField("a", offsetof(struct particle, a), adios2::DataType::Int8);
    particleDef1.AddField("b", offsetof(struct particle, b), adios2::DataType::Int32, 4);

    auto particleDef2 = io.DefineStruct("particle2", sizeof(particle) + 4);
    particleDef2.AddField("a", offsetof(struct particle, a), adios2::DataType::Int8, 1);
    particleDef2.AddField("b", offsetof(struct particle, b), adios2::DataType::Int32, 4);
    particleDef2.AddField("c", sizeof(particle) /* OK offset */, adios2::DataType::Int32);

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

            {
                auto varStruct = io.InquireVariable("particles");
                varStruct.SetSelection({start, count});
                ASSERT_TRUE(varStruct);
                engine.Get(varStruct, myParticles.data(), adios2::Mode::Sync);
                for (size_t i = 0; i < datasize; ++i)
                {
                    ASSERT_EQ(myParticles[i].a, 5);
                    ASSERT_EQ(myParticles[i].b[0], 0);
                    ASSERT_EQ(myParticles[i].b[1], 10);
                    ASSERT_EQ(myParticles[i].b[2], 20);
                    ASSERT_EQ(myParticles[i].b[3], 30);
                }
                auto structBlocks = engine.BlocksInfo(varStruct, currentStep);
                ASSERT_EQ(structBlocks.size(), worldSize - mpiSize);
                for (size_t i = 0; i < structBlocks.size(); ++i)
                {
                    ASSERT_FALSE(structBlocks[i].IsValue);
                    ASSERT_FALSE(structBlocks[i].IsReverseDims);
                    ASSERT_EQ(structBlocks[i].Step, currentStep);
                }
            }

            {
                auto varStruct = io.InquireStructVariable("particles", particleDef1);
                varStruct.SetSelection({start, count});
                ASSERT_TRUE(varStruct);
                engine.Get(varStruct, myParticles.data(), adios2::Mode::Sync);
                for (size_t i = 0; i < datasize; ++i)
                {
                    ASSERT_EQ(myParticles[i].a, 5);
                    ASSERT_EQ(myParticles[i].b[0], 0);
                    ASSERT_EQ(myParticles[i].b[1], 10);
                    ASSERT_EQ(myParticles[i].b[2], 20);
                    ASSERT_EQ(myParticles[i].b[3], 30);
                }
                auto structBlocks = engine.BlocksInfo(varStruct, currentStep);
                ASSERT_EQ(structBlocks.size(), worldSize - mpiSize);
                for (size_t i = 0; i < structBlocks.size(); ++i)
                {
                    ASSERT_FALSE(structBlocks[i].IsValue);
                    ASSERT_FALSE(structBlocks[i].IsReverseDims);
                    ASSERT_EQ(structBlocks[i].Step, currentStep);
                }
            }
            {
                auto varStruct = io.InquireStructVariable("particles", particleDef2);
                ASSERT_FALSE(varStruct);
            }

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

TEST_F(SscEngineTest, TestSscStruct)
{
    {
        std::string filename = "TestSscStruct";
        adios2::Params engineParams = {{"Verbose", "0"}};
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
        std::string filename = "TestSscStructNaive";
        adios2::Params engineParams = {{"Verbose", "0"}, {"EngineMode", "naive"}};
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
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
