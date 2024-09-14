/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>
#include <ctime>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

class SstWriteFails : public ::testing::Test
{
public:
    SstWriteFails() = default;
};

adios2::Params engineParams = {}; // parsed from command line
std::string fname = "ADIOS2SstFAIL";

int CompressSz = 0;
int CompressZfp = 0;

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

// ADIOS2 SST write
TEST_F(SstWriteFails, InvalidPut)
{
    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    auto scalar_r64 = io.DefineVariable<double>("scalar_r64");

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    adios2::Params NoReaders = {{"RendezvousReaderCount", "0"}};
    io.SetParameters(NoReaders);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    double data_scalar_R64 = 0.0;

    EXPECT_THROW(engine.Put(scalar_r64, data_scalar_R64), std::logic_error);

    // Close the file
    engine.Close();
}

// ADIOS2 SST write
TEST_F(SstWriteFails, InvalidBeginStep)
{
    // Write test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    auto scalar_r64 = io.DefineVariable<double>("scalar_r64");

    // Create the Engine
    io.SetEngine("Sst");
    io.SetParameters(engineParams);

    adios2::Params NoReaders = {{"RendezvousReaderCount", "0"}};
    io.SetParameters(NoReaders);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    double data_scalar_R64 = 0.0;

    engine.BeginStep();
    engine.Put(scalar_r64, data_scalar_R64);
    EXPECT_THROW(engine.BeginStep(), std::logic_error);

    // Close the file
    engine.Close();
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

    while ((argc > 1) && (argv[1][0] == '-'))
    {
        if (std::string(argv[1]) == "--expect_time_gap")
        {
            //  TimeGapExpected++;   Nothing on write side
        }
        else if (std::string(argv[1]) == "--compress_sz")
        {
            CompressSz++;
        }
        else if (std::string(argv[1]) == "--compress_zfp")
        {
            CompressZfp++;
        }
        else if (std::string(argv[1]) == "--filename")
        {
            fname = std::string(argv[2]);
            argv++;
            argc--;
        }
        else
        {
            throw std::invalid_argument("Unknown argument \"" + std::string(argv[1]) + "\"");
        }
        argv++;
        argc--;
    }
    if (argc > 1)
    {
        engineParams = ParseEngineParams(argv[1]);
    }

    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
