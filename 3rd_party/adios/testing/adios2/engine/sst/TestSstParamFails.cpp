/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

class SstParamFailTest : public ::testing::Test
{
public:
    SstParamFailTest() = default;
};

int TimeGapExpected = 0;
std::string fname = "ADIOS2Sst";

static std::string Trim(std::string &str)
{
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

/*
 * Engine parameters spec is a poor-man's JSON.  name:value pairs are separated
 * by commas.  White space is trimmed off front and back.  No quotes or anything
 * fancy allowed.
 */
static adios2::Params ParseEngineParams(std::string Input)
{
    std::istringstream ss(Input);
    std::string Param;
    adios2::Params Ret = {};

    while (std::getline(ss, Param, ','))
    {
        std::istringstream ss2(Param);
        std::string ParamName;
        std::string ParamValue;
        std::getline(ss2, ParamName, ':');
        if (!std::getline(ss2, ParamValue, ':'))
        {
            throw std::invalid_argument("Engine parameter \"" + Param + "\" missing value");
        }
        Ret[Trim(ParamName)] = Trim(ParamValue);
    }
    return Ret;
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, ParamNoValue)
{
    EXPECT_THROW(adios2::Params engineParams = ParseEngineParams("MarshalMethod:"),
                 std::invalid_argument);
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, ParamNoValue2)
{
    EXPECT_THROW(adios2::Params engineParams = ParseEngineParams("MarshalMethod"),
                 std::invalid_argument);
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, MarshalMethodUnknown)
{
    adios2::Params engineParams = ParseEngineParams("MarshalMethod:unknown");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    adios2::IO io = adios.DeclareIO("TestIO");

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    adios2::Engine engine;
    EXPECT_THROW(engine = io.Open(fname, adios2::Mode::Write), std::invalid_argument);
    // Close the file
    EXPECT_FALSE(engine);
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, CompressionMethodUnknown)
{
    adios2::Params engineParams = ParseEngineParams("CompressionMethod:unknown");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    adios2::IO io = adios.DeclareIO("TestIO");

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    adios2::Engine engine;
    EXPECT_THROW(engine = io.Open(fname, adios2::Mode::Write), std::invalid_argument);
    // Close the file
    EXPECT_FALSE(engine);
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, QueueFullPolicyUnknown)
{
    adios2::Params engineParams = ParseEngineParams("QueueFullPolicy:unknown");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    adios2::IO io = adios.DeclareIO("TestIO");

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    adios2::Engine engine;
    EXPECT_THROW(engine = io.Open(fname, adios2::Mode::Write), std::invalid_argument);
    // Close the file
    EXPECT_FALSE(engine);
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, RendezvousReaderCountBad)
{
    adios2::Params engineParams = ParseEngineParams("RendezvousReaderCount:unknown");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    adios2::IO io = adios.DeclareIO("TestIO");

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    adios2::Engine engine;
    EXPECT_THROW(engine = io.Open(fname, adios2::Mode::Write), std::invalid_argument);
    // Close the file
    EXPECT_FALSE(engine);
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, QueueLimitBad)
{
    adios2::Params engineParams = ParseEngineParams("QueueLimit:unknown");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    adios2::IO io = adios.DeclareIO("TestIO");

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    adios2::Engine engine;
    EXPECT_THROW(engine = io.Open(fname, adios2::Mode::Write), std::invalid_argument);
    // Close the file
    EXPECT_FALSE(engine);
}

// ADIOS2 Sst param fail
TEST_F(SstParamFailTest, PreciousFirstTimestepParamTest)
{
    adios2::Params engineParams = ParseEngineParams("FirstTimestepPrecious:unknown");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    adios2::IO io = adios.DeclareIO("TestIO");
    adios2::Engine engine;

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    EXPECT_THROW(engine = io.Open(fname, adios2::Mode::Write), std::invalid_argument);
    // Close the file
    EXPECT_FALSE(engine);
}

//******************************************************************************
// main
//******************************************************************************

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
