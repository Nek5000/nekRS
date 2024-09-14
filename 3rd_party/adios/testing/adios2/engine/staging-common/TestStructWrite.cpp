/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>
#include <ctime>

#include <chrono>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <thread>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData2.h"

#include "ParseArgs.h"

using namespace adios2;

class StructWriteTest : public ::testing::Test
{
public:
    StructWriteTest() = default;
};

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif
int mpiRank = 0;
int mpiSize = 1;

// ADIOS2 COMMON write
TEST_F(StructWriteTest, ADIOS2CommonWrite)
{

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
    MPI_Comm_size(testComm, &mpiSize);

    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Dims shape = {10, (size_t)mpiSize * 2};
    adios2::Dims start = {2, (size_t)mpiRank * 2};
    adios2::Dims count = {5, 2};
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
    particleDef.AddField("a", offsetof(struct particle, a), adios2::DataType::Int8);
    particleDef.AddField("b", offsetof(struct particle, b), adios2::DataType::Int32, 4);
    auto varStruct = io.DefineStructVariable("particles", particleDef, shape, start, count);
    EXPECT_THROW(particleDef.AddField("c", 4, adios2::DataType::Int32), std::runtime_error);
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
    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);
    engine.LockWriterDefinitions();
    for (size_t i = 0; i < (size_t)NSteps; ++i)
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

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);

    DelayMS = 0; // zero for common writer

    ParseArgs(argc, argv);

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engine == "SST" || engine == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, thread_support_level, &provided);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &testComm);
#endif

    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
#ifdef CRAY_MPICH_VERSION
    MPI_Barrier(MPI_COMM_WORLD);
#else
    MPI_Finalize();
#endif
#endif

    return result;
}
