#include <chrono>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

#include <adios2.h>
#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

class SstOnDemandReadTest : public ::testing::Test
{
public:
    SstOnDemandReadTest() = default;
};

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

// ADIOS2 Sst read
TEST_F(SstOnDemandReadTest, ADIOS2SstOnDemandRead)
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
    try
    {
        adios2::IO sstIO = adios.DeclareIO("sstOnDemand");

        sstIO.SetEngine(engine);
        adios2::Engine sstReader = sstIO.Open(fname, adios2::Mode::Read);

        double get_time = 0;
        const std::size_t my_start = Nx * rank;
        const adios2::Dims pos_start{my_start};
        const adios2::Dims count{Nx};
        const adios2::Box<adios2::Dims> sel(pos_start, count);

        auto start_step = std::chrono::steady_clock::now();
        int steps = 0;
        std::vector<float> myFloats(variablesSize * Nx);
        while (sstReader.BeginStep() == adios2::StepStatus::OK)
        {
            adios2::Variable<int> stepVar = sstIO.InquireVariable<int>("Step");
            int writerStep;

            for (int v = 0; v < variablesSize; ++v)
            {
                std::string namev("sstFloats");
                namev += std::to_string(v);
                adios2::Variable<float> sstFloats = sstIO.InquireVariable<float>(namev);

                sstFloats.SetSelection(sel);
                auto start_get = std::chrono::steady_clock::now();
                sstReader.Get(sstFloats, myFloats.data() + (v * Nx));
                auto end_get = std::chrono::steady_clock::now();
                get_time += (end_get - start_get).count() / 1000;
                // std::this_thread::sleep_for (std::chrono::seconds(1));
            }
            sstReader.Get(stepVar, writerStep);
            sstReader.EndStep();
            steps += 1;
#ifdef DEBUG
            size_t currentStep = sstReader.CurrentStep();
            for (unsigned int v = 0; v < variablesSize; ++v)
            {
                std::cout << name << ": Get step " << currentStep << " variable" << v << " "
                          << myFloats[v * Nx] << std::endl;
            }
#endif
        }
        auto end_step = std::chrono::steady_clock::now();
        double total_time = ((double)(end_step - start_step).count()) / (size * 1000.0);
        get_time /= size;

        double global_get_sum = 0;
        double global_sum = 0;
#if ADIOS2_USE_MPI
        MPI_Reduce(&get_time, &global_get_sum, 1, MPI_DOUBLE, MPI_SUM, 0, testComm);
        MPI_Reduce(&total_time, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, testComm);
#else
        global_sum = total_time;
        global_get_sum = get_time;
#endif

        // Time in microseconds
        if (rank == 0)
        {
            std::cout << "SST,Read," << size << "," << Nx << "," << variablesSize << "," << steps
                      << "," << global_get_sum << "," << global_sum << std::endl;
        }
        EXPECT_EQ(NSteps, steps);
        sstReader.Close();
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM "
                     "from rank "
                  << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
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
