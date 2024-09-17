/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 *  Created on: Jan 2018
 *      Author: Norbert Podhorszki
 */

#include <chrono>
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <numeric>
#include <sstream>
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <adios2.h>

#include "../SmallTestData.h"

int numprocs, wrank;
std::string engineName;           // comes from command line
bool serializeWriterReader;       // comes from command line
adios2::Params engineParams = {}; // parsed from command line

std::string Trim(std::string &str)
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
adios2::Params ParseEngineParams(std::string Input)
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

static unsigned int ParseUintParam(const std::string &optionName, char *arg)
{
    char *end;
    long retval = std::strtol(arg, &end, 10);
    if (end[0] || errno == ERANGE)
    {
        throw std::invalid_argument("Invalid value given for " + optionName + ": " +
                                    std::string(arg));
    }
    if (retval < 0)
    {
        throw std::invalid_argument("Negative value given for " + optionName + ": " +
                                    std::string(arg));
    }
    return static_cast<unsigned int>(retval);
}

class Common : public ::testing::Test
{
public:
    Common() = default;
    SmallTestData m_TestData;
};

// ADIOS2 write, read  attributes
TEST_F(Common, NewAttributeEveryStep)
{
    std::string streamName = "attributes_" + engineName;
    if (engineName == "HDF5")
    {
        streamName += ".h5";
    }
    else
    {
        streamName += ".bp";
    }
    const int steps = 3;
    const int nwriters = 1;
    const int nreaders = 1;
    if (!wrank)
    {
        std::cout << "test " << nwriters << " writers  with " << nreaders << " readers "
                  << (serializeWriterReader ? "sequentially" : "concurrently") << std::endl;
    }

    if (nwriters + nreaders > numprocs)
    {
        if (!wrank)
        {
            std::cout << "skip test: writers+readers > available processors " << std::endl;
        }
        return;
    }

    unsigned int color;
    if (wrank < nwriters)
    {
        color = 0; // writers
    }
    else if (wrank < nwriters + nreaders)
    {
        color = 1; // readers
    }
    else
    {
        color = 2; // not participating in test
    }
    int rank, nproc;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    if (color == 0)
    {
        std::cout << "Process wrank " << wrank << " plays Writer rank " << rank << std::endl;
        if (!rank)
        {
            std::cout << "There are " << nproc << " Writers" << std::endl;
        }

        adios2::ADIOS adios(comm);
        adios2::IO io = adios.DeclareIO("writer");
        io.SetEngine(engineName);
        io.SetParameters(engineParams);

        const size_t shape = static_cast<size_t>(nproc);
        const size_t start = static_cast<size_t>(rank);
        adios2::Variable<double> var =
            io.DefineVariable<double>("v", {shape}, {start}, {1}, adios2::ConstantDims);
        io.DefineAttribute<std::string>("v/unit", "km/s");

        adios2::Engine writer = io.Open(streamName, adios2::Mode::Write, comm);

        for (size_t step = 0; step < steps; ++step)
        {
            writer.BeginStep(adios2::StepMode::Append);
            const double d = static_cast<double>(step + 1);
            if (rank == 0)
                writer.Put<double>(var, &d);
            const std::string aname = "a" + std::to_string(step);
            io.DefineAttribute<uint64_t>(aname, step + 1);
            writer.EndStep();
        }
        writer.Close();
        if (serializeWriterReader)
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    else if (color == 1)
    {
        if (serializeWriterReader)
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }
        std::cout << "Process wrank " << wrank << " plays Reader rank " << rank << std::endl;
        int rank, nproc;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &nproc);
        if (!rank)
        {
            std::cout << "There are " << nproc << " Readers" << std::endl;
        }

        adios2::ADIOS adios(comm);
        adios2::IO io = adios.DeclareIO("reader");
        io.SetEngine(engineName);
        io.SetParameters(engineParams);
        adios2::Engine reader = io.Open(streamName, adios2::Mode::Read, comm);

        for (size_t step = 0; step < steps; ++step)
        {
            adios2::StepStatus status = reader.BeginStep(adios2::StepMode::Read);
            if (status != adios2::StepStatus::OK)
            {
                throw std::runtime_error("Expected step " + std::to_string(step) +
                                         " to be available in Reader but "
                                         "BeginStep() returned not-OK");
            }

            auto var = io.InquireVariable<double>("v");
            if (!var)
            {
                throw std::ios_base::failure("Missing 'v' variable in step " +
                                             std::to_string(step));
            }

            auto aUnit = io.InquireAttribute<std::string>("v/unit");
            if (!aUnit)
            {
                throw std::ios_base::failure("Missing 'v/unit' attribute in step " +
                                             std::to_string(step));
            }

            const std::string aname = "a" + std::to_string(step);
            auto aStep = io.InquireAttribute<uint64_t>(aname);
            if (!aStep)
            {
                throw std::ios_base::failure("Missing '" + aname + "' attribute in step " +
                                             std::to_string(step));
            }
            uint64_t expectedAttributeValue = step + 1;
            EXPECT_TRUE(aStep);
            ASSERT_EQ(aStep.Data().size() == 1, true);
            ASSERT_EQ(aStep.Type(), adios2::GetType<uint64_t>());
            ASSERT_EQ(aStep.Data().front(), expectedAttributeValue);

            if (!rank)
            {
                std::cout << "In step " << step << " Readers got attribute 'v/unit' with value "
                          << aUnit.Data().front() << std::endl;
                std::cout << "In step " << step << " Readers got attribute '" << aname
                          << "' with value " << aStep.Data().front() << std::endl;
            }

            double d[nwriters];
            reader.Get(var, d);
            reader.EndStep();

            double expectedScalarValue = static_cast<double>(step + 1);

            EXPECT_EQ(d[0], expectedScalarValue)
                << "Error in read, did not receive the expected values for "
                   "'v':"
                << " rank " << rank << ", step " << step << ", expected " << expectedScalarValue
                << ", received " << d[0];
        }
        reader.Close();
    }
    else if (color == 2)
    {
        if (serializeWriterReader)
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    // Separate each individual test with a big gap in time
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void threadTimeoutRun(size_t t)
{
    std::this_thread::sleep_for(std::chrono::seconds(t));
    throw std::runtime_error("Timeout reached");
}

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
    std::thread threadTimeout(threadTimeoutRun, 300);
    threadTimeout.detach();

    int provided;
    ::testing::InitGoogleTest(&argc, argv);

    engineName = std::string(argv[1]);

    int threadSupportLevel = MPI_THREAD_SINGLE;
    if (engineName == "SST")
    {
        threadSupportLevel = MPI_THREAD_MULTIPLE;
    }

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, threadSupportLevel, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    unsigned int p = ParseUintParam("serializeWriterReader", argv[2]);
    serializeWriterReader = (p != 0);

    if (argc > 3)
    {
        engineParams = ParseEngineParams(argv[3]);
    }

    if (!wrank)
    {
        std::cout << "Test " << engineName << " engine with " << numprocs << " processes "
                  << std::endl;
    }

    int result;
    result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
