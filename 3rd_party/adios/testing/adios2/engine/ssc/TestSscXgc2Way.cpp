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
int mpiGroup;
MPI_Comm mpiComm;

class SscEngineTest : public ::testing::Test
{
public:
    SscEngineTest() = default;
};

void xgc(const Dims &shape, const Dims &start, const Dims &count, const size_t steps,
         const adios2::Params &engineParams)
{
    size_t datasize = std::accumulate(count.begin(), count.end(), static_cast<size_t>(1),
                                      std::multiplies<size_t>());

    std::vector<float> x_to_g_data(datasize);
    std::vector<float> g_to_x_data;

    adios2::ADIOS adios(mpiComm);

    adios2::IO x_to_g_io = adios.DeclareIO("x_to_g");
    x_to_g_io.SetEngine("ssc");
    x_to_g_io.SetParameters(engineParams);

    adios2::IO g_to_x_io = adios.DeclareIO("g_to_x");
    g_to_x_io.SetEngine("ssc");
    g_to_x_io.SetParameters(engineParams);

    auto x_to_g_var = x_to_g_io.DefineVariable<float>("x_to_g", shape, start, count);

    adios2::Engine x_to_g_engine = x_to_g_io.Open("x_to_g", adios2::Mode::Write);
    adios2::Engine g_to_x_engine = g_to_x_io.Open("g_to_x", adios2::Mode::Read);
    x_to_g_engine.LockWriterDefinitions();
    g_to_x_engine.LockReaderSelections();

    for (size_t i = 0; i < steps; ++i)
    {
        x_to_g_engine.BeginStep();
        GenData(x_to_g_data, i, start, count, shape);
        x_to_g_engine.Put(x_to_g_var, x_to_g_data.data(), adios2::Mode::Sync);
        x_to_g_engine.EndStep();

        g_to_x_engine.BeginStep();
        auto g_to_x_var = g_to_x_io.InquireVariable<float>("g_to_x");
        auto readShape = g_to_x_var.Shape();
        g_to_x_data.resize(std::accumulate(readShape.begin(), readShape.end(),
                                           static_cast<size_t>(1), std::multiplies<size_t>()));
        g_to_x_engine.Get(g_to_x_var, g_to_x_data.data(), adios2::Mode::Sync);
        VerifyData(g_to_x_data.data(), i, Dims(readShape.size(), 0), readShape, readShape, mpiRank);
        g_to_x_engine.EndStep();
    }

    x_to_g_engine.Close();
    g_to_x_engine.BeginStep();
    g_to_x_engine.Close();
}

void gene(const Dims &shape, const Dims &start, const Dims &count, const size_t steps,
          const adios2::Params &engineParams)
{
    adios2::ADIOS adios(mpiComm);

    adios2::IO x_to_g_io = adios.DeclareIO("x_to_g");
    x_to_g_io.SetEngine("ssc");
    x_to_g_io.SetParameters(engineParams);

    adios2::IO g_to_x_io = adios.DeclareIO("g_to_x");
    g_to_x_io.SetEngine("ssc");
    g_to_x_io.SetParameters(engineParams);

    auto g_to_x_var = g_to_x_io.DefineVariable<float>("g_to_x", shape, start, count);

    adios2::Engine x_to_g_engine = x_to_g_io.Open("x_to_g", adios2::Mode::Read);
    adios2::Engine g_to_x_engine = g_to_x_io.Open("g_to_x", adios2::Mode::Write);
    g_to_x_engine.LockWriterDefinitions();
    x_to_g_engine.LockReaderSelections();

    size_t datasize = std::accumulate(shape.begin(), shape.end(), static_cast<size_t>(1),
                                      std::multiplies<size_t>());
    std::vector<float> x_to_g_data;
    std::vector<float> g_to_x_data(datasize);

    for (size_t i = 0; i < steps; ++i)
    {
        x_to_g_engine.BeginStep(StepMode::Read, 5);
        auto x_to_g_var = x_to_g_io.InquireVariable<float>("x_to_g");
        auto readShape = x_to_g_var.Shape();
        x_to_g_data.resize(std::accumulate(readShape.begin(), readShape.end(),
                                           static_cast<size_t>(1), std::multiplies<size_t>()));
        x_to_g_engine.Get(x_to_g_var, x_to_g_data.data(), adios2::Mode::Sync);
        VerifyData(x_to_g_data.data(), i, Dims(readShape.size(), 0), readShape, readShape, mpiRank);
        x_to_g_engine.EndStep();

        g_to_x_engine.BeginStep();
        GenData(g_to_x_data, i, start, count, shape);
        g_to_x_engine.Put(g_to_x_var, g_to_x_data.data(), adios2::Mode::Sync);
        g_to_x_engine.EndStep();
    }

    x_to_g_engine.BeginStep();
    x_to_g_engine.Close();
    g_to_x_engine.Close();
}

TEST_F(SscEngineTest, TestSscXgc2Way)
{

    {
        Dims start, count, shape;
        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

        if (worldSize < 6)
        {
            return;
        }

        if (worldRank < 2)
        {
            mpiGroup = 0;
        }
        else
        {
            mpiGroup = 1;
        }

        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);

        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);

        size_t steps = 20;

        if (mpiGroup == 0)
        {
            shape = {(size_t)mpiSize, 10};
            start = {(size_t)mpiRank, 0};
            count = {1, 10};
            adios2::Params engineParams = {{"RendezvousAppCount", "2"},
                                           {"MaxStreamsPerApp", "2"},
                                           {"OpenTimeoutSecs", "3"},
                                           {"Verbose", "0"}};
            xgc(shape, start, count, steps, engineParams);
        }

        if (mpiGroup == 1)
        {
            shape = {(size_t)mpiSize, 10};
            start = {(size_t)mpiRank, 0};
            count = {1, 10};
            adios2::Params engineParams = {{"RendezvousAppCount", "2"},
                                           {"MaxStreamsPerApp", "2"},
                                           {"OpenTimeoutSecs", "3"},
                                           {"Verbose", "0"}};
            gene(shape, start, shape, steps, engineParams);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    {
        Dims start, count, shape;
        int worldRank, worldSize;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

        if (worldSize < 6)
        {
            return;
        }

        if (worldRank < 2)
        {
            mpiGroup = 0;
        }
        else
        {
            mpiGroup = 1;
        }

        MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, worldRank, &mpiComm);

        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);

        size_t steps = 20;

        if (mpiGroup == 0)
        {
            shape = {(size_t)mpiSize, 10};
            start = {(size_t)mpiRank, 0};
            count = {1, 10};
            adios2::Params engineParams = {{"RendezvousAppCount", "2"},
                                           {"MaxStreamsPerApp", "2"},
                                           {"OpenTimeoutSecs", "3"},
                                           {"EngineMode", "naive"},
                                           {"Verbose", "0"}};
            xgc(shape, start, count, steps, engineParams);
        }

        if (mpiGroup == 1)
        {
            shape = {(size_t)mpiSize, 10};
            start = {(size_t)mpiRank, 0};
            count = {1, 10};
            adios2::Params engineParams = {{"RendezvousAppCount", "2"},
                                           {"MaxStreamsPerApp", "2"},
                                           {"OpenTimeoutSecs", "3"},
                                           {"EngineMode", "naive"},
                                           {"Verbose", "0"}};
            gene(shape, start, shape, steps, engineParams);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
