#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

std::vector<float> create_random_data(int n)
{
    std::vector<float> v(n);

    std::generate(begin(v), end(v), []() { return ((float)(rand() % 100)); });
    return v;
}

class SstOnDemandWriteTest : public ::testing::Test
{
public:
    SstOnDemandWriteTest() = default;
};

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

// ADIOS2 Sst Write
TEST_F(SstOnDemandWriteTest, ADIOS2SstOnDemandWrite)
{
    int rank = 0;
    int size = 1;

    int variablesSize = 100;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &rank);
    MPI_Comm_size(testComm, &size);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif

    auto myFloats = create_random_data((int)Nx * variablesSize);

    try
    {
        adios2::IO sstIO = adios.DeclareIO("sstOnDemand");

        sstIO.SetEngine(engine);
        engineParams["StepDistributionMode"] = "OnDemand";
        sstIO.SetParameters(engineParams);
        std::vector<adios2::Variable<float>> sstFloats(variablesSize);
        for (int v = 0; v < variablesSize; ++v)
        {
            std::string namev("sstFloats");
            namev += std::to_string(v);
            sstFloats[v] = sstIO.DefineVariable<float>(namev, {size * Nx}, {rank * Nx}, {Nx});
        }
        auto stepVar = sstIO.DefineVariable<int>("Step");

        adios2::Engine sstWriter = sstIO.Open(fname, adios2::Mode::Write);
        double put_time = 0;
        auto start_step = std::chrono::steady_clock::now();
        for (int timeStep = 0; timeStep < NSteps; ++timeStep)
        {
            sstWriter.BeginStep();
            for (int v = 0; v < variablesSize; ++v)
            {
                myFloats[v * Nx] = (float)v + timeStep * variablesSize;
                auto start_put = std::chrono::steady_clock::now();
                sstWriter.Put<float>(sstFloats[v], myFloats.data() + v * Nx);
                auto end_put = std::chrono::steady_clock::now();
                put_time += (end_put - start_put).count() / 1000;
                // std::this_thread::sleep_for (std::chrono::seconds(10));
#ifdef DEBUG
                std::cout << fname << ": Put step " << timeStep << " variable" << v << " "
                          << myFloats[v * Nx] << std::endl;
#endif
            }
            sstWriter.Put<int>(stepVar, timeStep);
            sstWriter.EndStep();
        }
        auto end_step = std::chrono::steady_clock::now();
        double total_time = ((double)(end_step - start_step).count()) / (size * 1000.0);

        double global_put_sum;
        double global_sum;
#if ADIOS2_USE_MPI
        MPI_Reduce(&put_time, &global_put_sum, 1, MPI_DOUBLE, MPI_SUM, 0, testComm);
        MPI_Reduce(&total_time, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, testComm);
#else
        global_sum = total_time;
        global_put_sum = put_time;
#endif
        // Time in microseconds
        if (rank == 0)
            std::cout << "SST,Write," << size << "," << Nx << "," << variablesSize << "," << NSteps
                      << "," << global_put_sum / size << "," << global_sum / size << std::endl;
        sstWriter.Close();
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM from rank " << rank
                  << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
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
