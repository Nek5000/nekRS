/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

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

class StructReadTest : public ::testing::Test
{
public:
    StructReadTest() = default;
};

typedef std::chrono::duration<double> Seconds;

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif
int mpiRank = 0;
int mpiSize = 1;

// ADIOS2 Common read
TEST_F(StructReadTest, ADIOS2StructRead)
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
    adios2::Engine engine = io.Open(fname, adios2::Mode::Read);

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
    struct particle
    {
        int8_t a;
        int b[4];
    };
    std::vector<particle> myParticles(datasize);

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

            int i;
            engine.Get(varIntScalar, &i);
            engine.PerformGets();

            {
                auto varStruct = io.InquireStructVariable("particles");
                StructDefinition ReadStruct = varStruct.GetReadStructDef();
                if (varStruct && !ReadStruct)
                {
                    // can't do a Get without setting the read structure
                    EXPECT_THROW(engine.Get(varStruct, myParticles.data(), adios2::Mode::Sync),
                                 std::logic_error);

                    StructDefinition WriteStruct = varStruct.GetWriteStructDef();
                    ASSERT_TRUE(WriteStruct);
                    std::cout << std::endl
                              << "Writer side structure was named \"" << WriteStruct.StructName()
                              << "\" and has size " << WriteStruct.StructSize() << std::endl;
                    for (size_t i = 0; i < WriteStruct.Fields(); i++)
                    {
                        std::cout << "\tField " << i << " - Name: \"" << WriteStruct.Name(i)
                                  << "\", Offset: " << WriteStruct.Offset(i)
                                  << ", Type: " << WriteStruct.Type(i)
                                  << ", ElementCount : " << WriteStruct.ElementCount(i)
                                  << std::endl;
                    }
                    std::cout << std::endl;
                    auto particleDef1 = io.DefineStruct("particle", sizeof(particle));
                    particleDef1.AddField("a", offsetof(particle, a), adios2::DataType::Int8);
                    particleDef1.AddField("b", offsetof(particle, b), adios2::DataType::Int32, 4);
                    varStruct.SetReadStructDef(particleDef1);
                }
                else if (varStruct)
                {
                    // set this already, but try something else
                    static bool first = true;
                    if (first)
                    {
                        first = false;
                        StructDefinition SaveDef = varStruct.GetReadStructDef();
                        ASSERT_TRUE(SaveDef);
                        auto particleDef2 = io.DefineStruct("particle", sizeof(particle));
                        particleDef2.AddField("a", offsetof(particle, a), adios2::DataType::Int8);
                        particleDef2.AddField("c", offsetof(particle, b), adios2::DataType::Int32,
                                              4);
                        varStruct.SetReadStructDef(particleDef2);

                        // can't change the read structure
                        EXPECT_THROW(engine.Get(varStruct, myParticles.data(), adios2::Mode::Sync),
                                     std::logic_error);
                        // restore the old one so we can succeed below
                        varStruct.SetReadStructDef(SaveDef);
                    }
                }
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

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    ParseArgs(argc, argv);

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engine == "SST" || engine == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, thread_support_level, &provided);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = 2;
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
