/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>

#include <chrono>
#include <iostream>
#include <signal.h>
#include <sstream>
#include <stdexcept>
#include <thread>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"

class CommonServerTest : public ::testing::Test
{
public:
    CommonServerTest() = default;
};

inline bool file_exists(const std::string &name)
{
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

// ADIOS2 COMMON write
TEST_F(CommonServerTest, ADIOS2CommonServer)
{
    int mpiRank = 0, mpiSize = 1;
    int GlobalCloseNow = 0;

    std::remove(shutdown_name.c_str());
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Server test data using ADIOS2

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    {
        adios2::Dims shape{static_cast<unsigned int>(Nx * mpiSize)};
        adios2::Dims start{static_cast<unsigned int>(Nx * mpiRank)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};
        adios2::Dims shape2{static_cast<unsigned int>(Nx * mpiSize), 2};
        adios2::Dims start2{static_cast<unsigned int>(Nx * mpiRank), 0};
        adios2::Dims count2{static_cast<unsigned int>(Nx), 2};
        adios2::Dims shape3{2, static_cast<unsigned int>(Nx * mpiSize)};
        adios2::Dims start3{0, static_cast<unsigned int>(Nx * mpiRank)};
        adios2::Dims count3{2, static_cast<unsigned int>(Nx)};
        adios2::Dims time_shape{static_cast<unsigned int>(mpiSize)};
        adios2::Dims time_start{static_cast<unsigned int>(mpiRank)};
        adios2::Dims time_count{1};
        io.DefineVariable<double>("scalar_r64");
        io.DefineVariable<int8_t>("i8", shape, start, count);
        io.DefineVariable<int16_t>("i16", shape, start, count);
        io.DefineVariable<int32_t>("i32", shape, start, count);
        io.DefineVariable<int64_t>("i64", shape, start, count);
        io.DefineVariable<float>("r32", shape, start, count);
        io.DefineVariable<double>("r64", shape, start, count);
        io.DefineVariable<std::complex<float>>("c32", shape, start, count);
        io.DefineVariable<std::complex<double>>("c64", shape, start, count);
        io.DefineVariable<double>("r64_2d", shape2, start2, count2);
        io.DefineVariable<double>("r64_2d_rev", shape3, start3, count3);
        io.DefineVariable<int64_t>("time", time_shape, time_start, time_count);
    }

    // Create the Engine
    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    // Advance to the next time step
    std::time_t EndTime = std::time(NULL) + DurationSeconds;
    size_t step = 0;

    while ((std::time(NULL) < EndTime) && !GlobalCloseNow)
    {
        // Generate test data for each process uniquely
        generateCommonTestData((int)step, mpiRank, mpiSize, (int)Nx, (int)Nx);

        engine.BeginStep();
        // Retrieve the variables that previously went out of scope
        auto scalar_r64 = io.InquireVariable<double>("scalar_r64");
        auto var_i8 = io.InquireVariable<int8_t>("i8");
        auto var_i16 = io.InquireVariable<int16_t>("i16");
        auto var_i32 = io.InquireVariable<int32_t>("i32");
        auto var_i64 = io.InquireVariable<int64_t>("i64");
        auto var_r32 = io.InquireVariable<float>("r32");
        auto var_r64 = io.InquireVariable<double>("r64");
        auto var_c32 = io.InquireVariable<std::complex<float>>("c32");
        auto var_c64 = io.InquireVariable<std::complex<double>>("c64");
        auto var_r64_2d = io.InquireVariable<double>("r64_2d");
        auto var_r64_2d_rev = io.InquireVariable<double>("r64_2d_rev");
        auto var_time = io.InquireVariable<int64_t>("time");

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces
        adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});
        adios2::Box<adios2::Dims> sel2({mpiRank * Nx, 0}, {Nx, 2});
        adios2::Box<adios2::Dims> sel3({0, mpiRank * Nx}, {2, Nx});
        adios2::Box<adios2::Dims> sel_time({static_cast<unsigned long>(mpiRank)}, {1});
        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);
        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);
        var_c32.SetSelection(sel);
        var_c64.SetSelection(sel);
        var_r64_2d.SetSelection(sel2);
        var_r64_2d_rev.SetSelection(sel3);
        var_time.SetSelection(sel_time);

        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count
        const adios2::Mode sync = adios2::Mode::Sync;

        engine.Put(scalar_r64, data_scalar_R64, sync);
        engine.Put(var_i8, data_I8.data(), sync);
        engine.Put(var_i16, data_I16.data(), sync);
        engine.Put(var_i32, data_I32.data(), sync);
        engine.Put(var_i64, data_I64.data(), sync);
        engine.Put(var_r32, data_R32.data(), sync);
        engine.Put(var_r64, data_R64.data(), sync);
        engine.Put(var_c32, data_C32.data(), sync);
        engine.Put(var_c64, data_C64.data(), sync);
        engine.Put(var_r64_2d, &data_R64_2d[0], sync);
        engine.Put(var_r64_2d_rev, &data_R64_2d_rev[0], sync);
        // Advance to the next time step
        std::time_t localtime = std::time(NULL);
        engine.Put(var_time, (int64_t *)&localtime);
        if (LockGeometry)
        {
            // we'll never change our data decomposition
            engine.LockWriterDefinitions();
        }
        if (AdvancingAttrs)
        {
            const std::string r64_Single = std::string("r64_PerStep_") + std::to_string(step);
            io.DefineAttribute<double>(r64_Single, (double)(step * 10.0));
        }
        engine.EndStep();
        std::cout << "Writer finished step " << step << std::endl;
        std::this_thread::sleep_for(
            std::chrono::milliseconds(DelayMS)); /* sleep for DelayMS milliseconds */
        step++;
        {
            int MyCloseNow = 0;
            if (file_exists(shutdown_name))
            {
                MyCloseNow = 1;
            }
#if ADIOS2_USE_MPI
            MPI_Allreduce(&MyCloseNow, &GlobalCloseNow, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
#else
            GlobalCloseNow = MyCloseNow;
#endif
        }
        if (GlobalCloseNow)
        {
            std::cout << "Writer closing stream because file \"" << shutdown_name
                      << "\" was noticed" << std::endl;
        }
    }
    std::cout << "Writer closing stream normally" << std::endl;
    // Close the file
    engine.Close();
}

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
