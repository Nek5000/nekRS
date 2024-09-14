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

#include "ParseArgs.h"

static int numprocs, wrank;
std::string engineName; // comes from command line

struct RunParams
{
    /* 2D decomposition of processes
     * npx_w x npy_w Writers
     * npx_r x npy_r Readers
     */
    size_t npx_w;
    size_t npy_w;
    size_t npx_r;
    size_t npy_r;
    RunParams(size_t xw, size_t yw, size_t xr, size_t yr)
    : npx_w{xw}, npy_w{yw}, npx_r{xr}, npy_r{yr} {};
};

/* This function is executed by INSTANTIATE_TEST_SUITE_P
   before main() and MPI_Init()!!! */
std::vector<RunParams> CreateRunParams()
{
    std::vector<RunParams> params;
    // 2 process test
    params.push_back(RunParams(1, 1, 1, 1));
    // 3 process tests
    params.push_back(RunParams(2, 1, 1, 1));
    params.push_back(RunParams(1, 2, 1, 1));
    params.push_back(RunParams(1, 1, 2, 1));
    params.push_back(RunParams(1, 1, 1, 2));
    // 4 process tests
    params.push_back(RunParams(1, 3, 1, 1));
    params.push_back(RunParams(2, 1, 2, 1));
    params.push_back(RunParams(2, 1, 1, 2));
    // 8 process tests
    params.push_back(RunParams(1, 1, 1, 7));
    params.push_back(RunParams(1, 7, 1, 1));
    params.push_back(RunParams(2, 2, 2, 2));
    // 16 process tests
    params.push_back(RunParams(3, 5, 1, 1));
    params.push_back(RunParams(1, 1, 5, 3));
    params.push_back(RunParams(4, 2, 2, 4));
    return params;
}

class TestStagingMPMD : public ::testing::TestWithParam<RunParams>
{
public:
    TestStagingMPMD() = default;
    const std::string streamName = "TestStream";

    void MainWriters(MPI_Comm comm, size_t npx, size_t npy, int steps, unsigned int sleeptime)
    {
        int rank, nproc;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &nproc);
        if (!rank)
        {
            std::cout << "There are " << nproc << " Writers" << std::endl;
        }
        size_t ndx = 50;
        size_t ndy = 60;
        size_t gndx = npx * ndx;
        size_t gndy = npy * ndy;
        size_t posx = rank % npx;
        size_t posy = rank / npx;
        size_t offsx = posx * ndx;
        size_t offsy = posy * ndy;

        std::vector<float> myArray(ndx * ndy);

        adios2::ADIOS adios(comm);
        adios2::IO io = adios.DeclareIO("writer");
        io.SetEngine(engineName);
        io.SetParameters(engineParams);

        adios2::Variable<float> varArray = io.DefineVariable<float>(
            "myArray", {gndx, gndy}, {offsx, offsy}, {ndx, ndy}, adios2::ConstantDims);

        adios2::Variable<double> varScalar = io.DefineVariable<double>("myScalar");

        adios2::Engine writer = io.Open(streamName, adios2::Mode::Write, comm);

        for (int step = 0; step < steps; ++step)
        {
            size_t idx = 0;
            for (size_t x = 0; x < ndx; ++x)
            {
                for (size_t y = 0; y < ndy; ++y)
                {
                    myArray[idx] = GetValue(offsx + x, offsy + y, step);
                    ++idx;
                }
            }
            writer.BeginStep(adios2::StepMode::Append);
            writer.Put<float>(varArray, myArray.data());
            if (wrank == 0)
                writer.Put<double>(varScalar, 1.5 * (step + 1));
            writer.EndStep();
            std::this_thread::sleep_for(std::chrono::milliseconds(sleeptime));
        }
        writer.Close();
    }

    /* encode position x,y as y.x floating point value in each array cell
     * add 1000*step to each value
     *
     * A 2 by 2 example with 3x3 local arrays
     *   0.0    0.001  0.002 |  0.003  0.004  0.005
     *   1.0    1.001  1.002 |  1.003  1.004  1.005
     *   2.0    2.001  2.002 |  2.003  2.004  2.005
     *   --------------------+---------------------
     *   3.0    3.001  3.002 |  3.003  3.004  3.005
     *   4.0    4.001  4.002 |  4.003  4.004  4.005
     *   5.0    5.001  5.002 |  5.003  5.004  5.005
     *
     * Next step each value is bigger by 1000
     */
    float GetValue(size_t offsx, size_t offsy, size_t step)
    {
        return 1000.0f * step + offsy + offsx / 1000.0f;
    }

    void CheckData(const std::vector<float> &array, size_t gndx, size_t gndy, size_t offsx,
                   size_t offsy, size_t ndx, size_t ndy, size_t step, int rank)
    {
        size_t idx = 0;
        for (size_t x = 0; x < ndx; ++x)
        {
            for (size_t y = 0; y < ndy; ++y)
            {
                float expectedValue = GetValue(offsx + x, offsy + y, step);
                EXPECT_EQ(array[idx], expectedValue)
                    << "Error in read, did not receive the expected value:"
                    << " rank " << rank << ", step " << step << ", gdim {" << gndx << ',' << gndy
                    << '}' << ", offs {" << offsx << ',' << offsy << '}' << ", ldim {" << ndx << ','
                    << ndy << '}' << ", lpos {" << x << ',' << y << '}';
                ++idx;
            }
        }
    }

    void MainReaders(MPI_Comm comm, size_t npx, size_t npy, unsigned int sleeptime,
                     float reader_timeout)
    {
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

        size_t posx = rank % npx;
        size_t posy = rank / npx;
        size_t step = 0;
        adios2::Variable<float> vMyArray;
        std::vector<float> myArray;
        adios2::Variable<double> vMyScalar;
        double myIncomingScalar;

        while (true)
        {
            adios2::StepStatus status = reader.BeginStep(adios2::StepMode::Read, reader_timeout);
            if (status != adios2::StepStatus::OK)
            {
                break;
            }

            vMyArray = io.InquireVariable<float>("myArray");
            if (!vMyArray)
            {
                throw std::ios_base::failure("Missing 'myArray' variable.");
            }

            vMyScalar = io.InquireVariable<double>("myScalar");
            if (!vMyScalar)
            {
                throw std::ios_base::failure("Missing 'myScalar' variable.");
            }

            // 2D decomposition of global array reading
            size_t gndx = vMyArray.Shape()[0];
            size_t gndy = vMyArray.Shape()[1];
            size_t ndx = gndx / npx;
            size_t ndy = gndy / npy;
            size_t offsx = ndx * posx;
            size_t offsy = ndy * posy;
            if (posx == npx - 1)
            {
                // right-most processes need to read all the rest of rows
                ndx = gndx - ndx * (npx - 1);
            }

            if (posy == npy - 1)
            {
                // bottom processes need to read all the rest of columns
                ndy = gndy - ndy * (npy - 1);
            }

            adios2::Dims count({ndx, ndy});
            adios2::Dims start({offsx, offsy});

            vMyArray.SetSelection({start, count});
            size_t elementsSize = count[0] * count[1];
            myArray.resize(elementsSize);

            reader.Get(vMyArray, myArray.data());
            reader.Get(vMyScalar, myIncomingScalar);
            reader.EndStep();
            float expectedScalarValue = 1.5f * (step + 1);
            EXPECT_EQ(myIncomingScalar, expectedScalarValue)
                << "Error in read, did not receive the expected value:"
                << " rank " << rank << ", step " << step;
            CheckData(myArray, gndx, gndy, offsx, offsy, ndx, ndy, step, rank);
            std::this_thread::sleep_for(std::chrono::milliseconds(sleeptime));
            ++step;
        }
        reader.Close();
    }

    void TestCommon(RunParams p, int steps, unsigned int writer_sleeptime,
                    unsigned int reader_sleeptime, float reader_timeout)
    {
        std::cout << "test " << p.npx_w << "x" << p.npy_w << " writers " << p.npx_r << "x"
                  << p.npy_r << " readers " << std::endl;

        size_t nwriters = p.npx_w * p.npy_w;
        size_t nreaders = p.npx_r * p.npy_r;
        if (nwriters + nreaders > static_cast<size_t>(numprocs))
        {
            if (!wrank)
            {
                std::cout << "skip test: writers+readers > available processors " << std::endl;
            }
            return;
        }

        int rank;
        MPI_Comm comm;

        unsigned int color;
        if (static_cast<size_t>(wrank) < nwriters)
        {
            color = 0; // writers
        }
        else if (static_cast<size_t>(wrank) < nwriters + nreaders)
        {
            color = 1; // readers
        }
        else
        {
            color = 2; // not participating in test
        }
        MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);
        MPI_Comm_rank(comm, &rank);

        if (color == 0)
        {
            std::cout << "Process wrank " << wrank << " rank " << rank << " calls MainWriters "
                      << std::endl;
            MainWriters(comm, p.npx_w, p.npy_w, steps, writer_sleeptime);
        }
        else if (color == 1)
        {
            std::cout << "Process wrank " << wrank << " rank " << rank << " calls MainReaders "
                      << std::endl;
            MainReaders(comm, p.npx_r, p.npy_r, reader_sleeptime, reader_timeout);
        }
        std::cout << "Process wrank " << wrank << " rank " << rank << " enters MPI barrier..."
                  << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        // Separate each individual test with a big gap in time
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
};

TEST_P(TestStagingMPMD, SingleStep)
{
    RunParams p = GetParam();
    TestCommon(p, 1, 0, 0, -1.0);
}

TEST_P(TestStagingMPMD, MultipleSteps)
{
    RunParams p = GetParam();
    TestCommon(p, 10, 0, 0, -1.0);
}

TEST_P(TestStagingMPMD, MultipleStepsWithTimeout)
{
    RunParams p = GetParam();
    // timeout in reader's BeginStep instead of blocking wait
    TestCommon(p, 10, 0, 0, 2.0);
}

TEST_P(TestStagingMPMD, SlowWriter)
{
    RunParams p = GetParam();
    TestCommon(p, 4, 100, 0, -1.0);
}

TEST_P(TestStagingMPMD, SlowReader)
{
    RunParams p = GetParam();
    TestCommon(p, 4, 0, 100, -1.0);
}

INSTANTIATE_TEST_SUITE_P(NxM, TestStagingMPMD, ::testing::ValuesIn(CreateRunParams()));

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

    if (argc > 2)
    {
        engineParams = ParseEngineParams(argv[2]);
    }

    if (!wrank)
    {
        std::cout << "Test " << engineName << " engine with " << numprocs << " processes "
                  << std::endl;
    }

    int result;
    result = RUN_ALL_TESTS();

#ifdef CRAY_MPICH_VERSION
    MPI_Barrier(MPI_COMM_WORLD);
#else
    MPI_Finalize();
#endif

    return result;
}
