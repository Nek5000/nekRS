/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <fstream>
#include <iostream>
#include <stdexcept>

#include <adios2.h>
#include <gtest/gtest.h>

std::string engineName; // comes from command line

class BPBufferSizeTest : public ::testing::Test
{
public:
    BPBufferSizeTest() = default;
};

std::ifstream vm_input;
void PrintVMStatus()
{
    std::string m_Name = "/proc/self/status";

    if (!vm_input.is_open())
    {
        vm_input.open(m_Name.c_str());
    }

    if (vm_input.is_open())
    {
        vm_input.seekg(0);
        for (std::string line; getline(vm_input, line);)
        {
            if (line.find("VmRSS") == 0)
                std::cout << line << " ";
            if (line.find("VmSize") == 0)
                std::cout << line << " ";
        }
        vm_input.close();
    }
}

void VMStatusClose()
{
    if (vm_input.is_open())
    {
        vm_input.close();
    }
}

size_t GetAndPrintBufferSize(adios2::Engine &engine, const std::string &info,
                             size_t step = 999999999)
{

    size_t s = engine.DebugGetDataBufferSize();
    std::cout << std::left << std::setw(35) << info;
    if (step < 999999999)
    {
        std::cout << " step " << std::setw(4) << std::to_string(step);
    }
    else
    {
        std::cout << std::setw(10) << " ";
    }
    std::cout << " buffer size = " << std::setw(12) << std::to_string(s) << " ";
    PrintVMStatus();
    std::cout << std::endl;
    return s;
}

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

// Put(Sync) and Put(Deferred) should have the same buffer consumption
TEST_F(BPBufferSizeTest, SyncDeferredIdenticalUsage)
{
    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 10485760; // 10M elements, 80MB variable
    std::vector<double> data;
    data.resize(Nx);

    // Number of steps
    const std::size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    std::string fnameSync = "ADIOS2BPBufferSizeSync_MPI.bp";
    std::string fnameDeferred = "ADIOS2BPBufferSizeDeferred_MPI.bp";
    std::string fnameDeferredPP = "ADIOS2BPBufferSizeDeferredPP_MPI.bp";
#else
    std::string fnameSync = "ADIOS2BPBufferSizeSync.bp";
    std::string fnameDeferred = "ADIOS2BPBufferSizeDeferred.bp";
    std::string fnameDeferredPP = "ADIOS2BPBufferSizeDeferredPP.bp";

#endif

    // Write test data using BP
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io = adios.DeclareIO("TestIO");
        adios2::Dims shape{static_cast<unsigned int>(Nx * mpiSize)};
        adios2::Dims start{static_cast<unsigned int>(Nx * mpiRank)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        auto var1 = io.DefineVariable<double>("r64_1", shape, start, count);
        auto var2 = io.DefineVariable<double>("r64_2", shape, start, count);
        auto var3 = io.DefineVariable<double>("r64_3", shape, start, count);

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("BPFile");
        }
        io.AddTransport("NULL");
        io.SetParameter("Profile", "OFF");

        // Write with Sync
        size_t bufsize_sync_beginstep[NSteps];
        size_t bufsize_sync_v1[NSteps];
        size_t bufsize_sync_v2[NSteps];
        size_t bufsize_sync_v3[NSteps];
        size_t bufsize_sync_endstep[NSteps];
        {
            adios2::Engine engine = io.Open(fnameSync, adios2::Mode::Write);

            GetAndPrintBufferSize(engine, "After Open():");
            for (size_t step = 0; step < NSteps; ++step)
            {
                engine.BeginStep();
                bufsize_sync_beginstep[step] = GetAndPrintBufferSize(engine, "After BeginStep():");
                engine.Put(var1, data.data(), adios2::Mode::Sync);
                bufsize_sync_v1[step] = GetAndPrintBufferSize(engine, "After Put(v1, Sync):", step);
                engine.Put(var2, data.data(), adios2::Mode::Sync);
                bufsize_sync_v2[step] = GetAndPrintBufferSize(engine, "After Put(v2, Sync):", step);
                engine.Put(var3, data.data(), adios2::Mode::Sync);
                bufsize_sync_v3[step] = GetAndPrintBufferSize(engine, "After Put(v3, Sync):", step);
                engine.EndStep();
                bufsize_sync_endstep[step] =
                    GetAndPrintBufferSize(engine, "After EndStep():", step);
            }

            engine.Close();
        }

        // Write with Deferred
        size_t bufsize_deferred_beginstep[NSteps];
        size_t bufsize_deferred_v1[NSteps];
        size_t bufsize_deferred_v2[NSteps];
        size_t bufsize_deferred_v3[NSteps];
        size_t bufsize_deferred_endstep[NSteps];
        {
            adios2::Engine engine = io.Open(fnameDeferred, adios2::Mode::Write);

            GetAndPrintBufferSize(engine, "After Open():");
            for (size_t step = 0; step < NSteps; ++step)
            {
                engine.BeginStep();
                bufsize_deferred_beginstep[step] =
                    GetAndPrintBufferSize(engine, "After BeginStep():");
                engine.Put(var1, data.data(), adios2::Mode::Deferred);
                bufsize_deferred_v1[step] =
                    GetAndPrintBufferSize(engine, "After Put(v1, Deferred):", step);
                engine.Put(var2, data.data(), adios2::Mode::Deferred);
                bufsize_deferred_v2[step] =
                    GetAndPrintBufferSize(engine, "After Put(v2, Deferred):", step);
                engine.Put(var3, data.data(), adios2::Mode::Deferred);
                bufsize_deferred_v3[step] =
                    GetAndPrintBufferSize(engine, "After Put(v3, Deferred):", step);
                engine.EndStep();
                bufsize_deferred_endstep[step] =
                    GetAndPrintBufferSize(engine, "After EndStep():", step);
            }

            engine.Close();
        }

        // Write with Deferred+PerformPuts
        size_t bufsize_deferred_pp_beginstep[NSteps];
        size_t bufsize_deferred_pp_v1[NSteps];
        size_t bufsize_deferred_pp_v2[NSteps];
        size_t bufsize_deferred_pp_v3[NSteps];
        size_t bufsize_deferred_pp_endstep[NSteps];
        {
            adios2::Engine engine = io.Open(fnameDeferredPP, adios2::Mode::Write);

            GetAndPrintBufferSize(engine, "After Open():");
            for (size_t step = 0; step < NSteps; ++step)
            {
                engine.BeginStep();
                bufsize_deferred_pp_beginstep[step] =
                    GetAndPrintBufferSize(engine, "After BeginStep():");
                engine.Put(var1, data.data(), adios2::Mode::Deferred);
                engine.PerformPuts();
                bufsize_deferred_pp_v1[step] =
                    GetAndPrintBufferSize(engine, "After Put(v1, Def)+PerformPuts():", step);
                engine.Put(var2, data.data(), adios2::Mode::Deferred);
                engine.PerformPuts();
                bufsize_deferred_pp_v2[step] =
                    GetAndPrintBufferSize(engine, "After Put(v2, Def)+PerformPuts():", step);
                engine.Put(var3, data.data(), adios2::Mode::Deferred);
                engine.PerformPuts();
                bufsize_deferred_pp_v3[step] =
                    GetAndPrintBufferSize(engine, "After Put(v3, Def)+PerformPuts():", step);
                engine.EndStep();
                bufsize_deferred_pp_endstep[step] =
                    GetAndPrintBufferSize(engine, "After EndStep():", step);
            }

            engine.Close();
        }

        VMStatusClose();

        /* Compare buffer size usage.
         * Buffer size should never be substantially more than the total data
         * size
         * */
        const size_t TotalDataSize = Nx * sizeof(double) * 3;
        const size_t MaxExtra = 18 * 1024 * 1024; /* 18MB extra allowed in buffer */
        for (size_t step = 0; step < NSteps; ++step)
        {
            EXPECT_LT(bufsize_sync_beginstep[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_sync_v1[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_sync_v2[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_sync_v3[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_sync_endstep[step], TotalDataSize + MaxExtra);

            EXPECT_LT(bufsize_deferred_beginstep[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_v1[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_v2[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_v3[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_endstep[step], TotalDataSize + MaxExtra);

            EXPECT_LT(bufsize_deferred_pp_beginstep[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_pp_v1[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_pp_v2[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_pp_v3[step], TotalDataSize + MaxExtra);
            EXPECT_LT(bufsize_deferred_pp_endstep[step], TotalDataSize + MaxExtra);
        }
    }
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

    if (argc > 1)
    {
        engineName = std::string(argv[1]);
    }
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
