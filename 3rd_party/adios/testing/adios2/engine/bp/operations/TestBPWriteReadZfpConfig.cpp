/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <numeric> //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#define str_helper(X) #X
#define mystr(X) str_helper(X)

void ZfpRate1D(const std::string configFile)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const auto begin = configFile.find("_rate") + 5;
    const auto end = configFile.find_last_of(".") - begin;
    const std::string rate = configFile.substr(begin, end);

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx);
    std::vector<double> r64s(Nx);

    // range 0 to 999
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadZfpConfig1D_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWriteReadZfpConfig1D_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile,
                        MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile);
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

            const adios2::Dims start{mpiRank * Nx};
            const adios2::Dims count{Nx};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZfpRate2D(const std::string configFile)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const auto begin = configFile.find("_rate") + 5;
    const auto end = configFile.find_last_of(".") - begin;
    const std::string rate = configFile.substr(begin, end);

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;
    const size_t Ny = 50;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx * Ny);
    std::vector<double> r64s(Nx * Ny);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadZfpConfig2D_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWriteReadZfpConfig2D_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile,
                        MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile);
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
        const adios2::Dims count{Nx, Ny};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r32.Shape()[1], Ny);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);

            const adios2::Dims start{mpiRank * Nx, 0};
            const adios2::Dims count{Nx, Ny};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZfpRate3D(const std::string configFile)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const auto begin = configFile.find("_rate") + 5;
    const auto end = configFile.find_last_of(".") - begin;
    const std::string rate = configFile.substr(begin, end);

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 10;
    const size_t Ny = 20;
    const size_t Nz = 15;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx * Ny * Nz);
    std::vector<double> r64s(Nx * Ny * Nz);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadZfpConfig3D_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWriteReadZfpConfig3D_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile,
                        MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile);
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny, Nz};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0, 0};
        const adios2::Dims count{Nx, Ny, Nz};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r32.Shape()[1], Ny);
            ASSERT_EQ(var_r32.Shape()[2], Nz);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);
            ASSERT_EQ(var_r64.Shape()[2], Nz);

            const adios2::Dims start{mpiRank * Nx, 0, 0};
            const adios2::Dims count{Nx, Ny, Nz};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx * Ny * Nz; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZfpRate1DSel(const std::string configFile)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const auto begin = configFile.find("_rate") + 5;
    const auto end = configFile.find_last_of(".") - begin;
    const std::string rate = configFile.substr(begin, end);

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx);
    std::vector<double> r64s(Nx);

    // range 0 to 999
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadZfpConfig1DSel_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWriteReadZfpConfig1DSel_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile,
                        MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile);
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        const adios2::Dims shape{static_cast<std::size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<std::size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], static_cast<std::size_t>(mpiSize) * Nx);

            const adios2::Dims start{static_cast<std::size_t>(mpiRank) * Nx + Nx / 2};
            const adios2::Dims count{Nx / 2};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx / 2; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[Nx / 2 + i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[Nx / 2 + i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZfpRate2DSel(const std::string configFile)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const auto begin = configFile.find("_rate") + 5;
    const auto end = configFile.find_last_of(".") - begin;
    const std::string rate = configFile.substr(begin, end);

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;
    const size_t Ny = 50;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx * Ny);
    std::vector<double> r64s(Nx * Ny);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadZfpConfig2DSel_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWriteReadZfpConfig2DSel_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile,
                        MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile);
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
        const adios2::Dims count{Nx, Ny};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r32.Shape()[1], Ny);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);

            const adios2::Dims start{mpiRank * Nx + Nx / 2, 0};
            const adios2::Dims count{Nx / 2, Ny};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx / 2 * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[Nx / 2 * Ny + i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[Nx / 2 * Ny + i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZfpRate3DSel(const std::string configFile)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const auto begin = configFile.find("_rate") + 5;
    const auto end = configFile.find_last_of(".") - begin;
    const std::string rate = configFile.substr(begin, end);

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 10;
    const size_t Ny = 20;
    const size_t Nz = 15;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx * Ny * Nz);
    std::vector<double> r64s(Nx * Ny * Nz);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadZfpConfig3DSel_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWriteReadZfpConfig3DSel_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile,
                        MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile);
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny, Nz};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0, 0};
        const adios2::Dims count{Nx, Ny, Nz};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r32.Shape()[1], Ny);
            ASSERT_EQ(var_r32.Shape()[2], Nz);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);
            ASSERT_EQ(var_r64.Shape()[2], Nz);

            const adios2::Dims start{static_cast<std::size_t>(mpiRank) * Nx + Nx / 2, 0, 0};
            const adios2::Dims count{Nx / 2, Ny, Nz};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx / 2 * Ny * Nz; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                ASSERT_LT(std::abs(decompressedR32s[i] - r32s[Nx / 2 * Ny * Nz + i]), 1E-4) << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[Nx / 2 * Ny * Nz + i]), 1E-4) << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void ZfpRate2DSmallSel(const std::string configFile)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const auto begin = configFile.find("_rate") + 5;
    const auto end = configFile.find_last_of(".") - begin;
    const std::string rate = configFile.substr(begin, end);

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 5;
    const size_t Ny = 5;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                               0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
                               0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24};
    std::vector<double> r64s = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                                0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
                                0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24};
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadZfpConfig2DSmallSel_" + rate + "_MPI.bp");
#else
    const std::string fname("BPWriteReadZfpConfig2DSmallSel_" + rate + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile,
                        MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(std::string(mystr(XML_CONFIG_DIR)) + adios2::PathSeparator + configFile);
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
        const adios2::Dims count{Nx, Ny};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        (void)var_r32;
        (void)var_r64;

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r32.Shape()[1], Ny);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);

            const adios2::Dims start{static_cast<std::size_t>(mpiRank) * Nx + 1, 1};
            const adios2::Dims count{2, 2};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            ASSERT_LT(std::abs(decompressedR32s[0] - 0.06), 0.01);
            ASSERT_LT(std::abs(decompressedR64s[0] - 0.06), 0.01);

            ASSERT_LT(std::abs(decompressedR32s[1] - 0.07), 0.01);
            ASSERT_LT(std::abs(decompressedR64s[1] - 0.07), 0.01);

            ASSERT_LT(std::abs(decompressedR32s[2] - 0.11), 0.01);
            ASSERT_LT(std::abs(decompressedR64s[2] - 0.11), 0.01);

            ASSERT_LT(std::abs(decompressedR32s[3] - 0.12), 0.01);
            ASSERT_LT(std::abs(decompressedR64s[3] - 0.12), 0.01);

            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

class BPWriteReadZfpConfig : public ::testing::TestWithParam<std::string>
{
public:
    BPWriteReadZfpConfig() = default;

    virtual void SetUp() {}
    virtual void TearDown() {}
};

TEST_P(BPWriteReadZfpConfig, ADIOS2BPWriteReadZfp1D) { ZfpRate1D(GetParam()); }
TEST_P(BPWriteReadZfpConfig, ADIOS2BPWriteReadZfp2D) { ZfpRate2D(GetParam()); }
TEST_P(BPWriteReadZfpConfig, ADIOS2BPWriteReadZfp3D) { ZfpRate3D(GetParam()); }
TEST_P(BPWriteReadZfpConfig, ADIOS2BPWriteReadZfp1DSel) { ZfpRate1DSel(GetParam()); }
TEST_P(BPWriteReadZfpConfig, ADIOS2BPWriteReadZfp2DSel) { ZfpRate2DSel(GetParam()); }
TEST_P(BPWriteReadZfpConfig, ADIOS2BPWriteReadZfp3DSel) { ZfpRate3DSel(GetParam()); }
TEST_P(BPWriteReadZfpConfig, ADIOS2BPWriteReadZfp2DSmallSel) { ZfpRate2DSmallSel(GetParam()); }

INSTANTIATE_TEST_SUITE_P(ZfpConfigFile, BPWriteReadZfpConfig,
                         ::testing::Values("configZfp_rate8.xml", "configZfp_rate8Simple.xml",
                                           "configZfp_rate9.xml", "configZfp_rate9Simple.xml",
                                           "configZfp_rate10.xml", "configZfp_rate10Simple.xml",
                                           "configZfp_rate8.yaml", "configZfp_rate9.yaml",
                                           "configZfp_rate10.yaml"));

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
