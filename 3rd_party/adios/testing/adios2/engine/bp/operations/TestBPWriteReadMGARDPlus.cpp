/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <algorithm>
#include <iostream> //std::cout
#include <numeric>  //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

void MGARDAccuracy1D(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;

    // Number of steps
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx);
    std::vector<double> r64s(Nx);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRMGARD1D_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD1D_" + tolerance + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        var_r32.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});
        var_r64.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            // bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            //        auto var_r32 = io.InquireVariable<float>("r32");
            //        EXPECT_TRUE(var_r32);
            //        ASSERT_EQ(var_r32.ShapeID(),
            //        adios2::ShapeID::GlobalArray); ASSERT_EQ(var_r32.Steps(),
            //        NSteps); ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            //        ASSERT_EQ(var_r32.Shape()[1], Ny);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

            const adios2::Dims start{mpiRank * Nx};
            const adios2::Dims count{Nx};
            const adios2::Box<adios2::Dims> sel(start, count);
            // var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            double maxDiff = 0;

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                double diff = std::abs(r64s[i] - decompressedR64s[i]);

                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }
            ++t;

            auto itMax = std::max_element(r64s.begin(), r64s.end());

            const double relativeMaxDiff = maxDiff / *itMax;
            ASSERT_LT(relativeMaxDiff, std::stod(tolerance));
            std::cout << "Relative Max Diff " << relativeMaxDiff << " tolerance " << tolerance
                      << "\n";
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void MGARDAccuracy2D(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

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
    const std::string fname("BPWRMGARD2D_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD2D_" + tolerance + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
        const adios2::Dims count{Nx, Ny};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        var_r32.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});
        var_r64.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            // bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            //        auto var_r32 = io.InquireVariable<float>("r32");
            //        EXPECT_TRUE(var_r32);
            //        ASSERT_EQ(var_r32.ShapeID(),
            //        adios2::ShapeID::GlobalArray); ASSERT_EQ(var_r32.Steps(),
            //        NSteps); ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            //        ASSERT_EQ(var_r32.Shape()[1], Ny);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);

            const adios2::Dims start{mpiRank * Nx, 0};
            const adios2::Dims count{Nx, Ny};
            const adios2::Box<adios2::Dims> sel(start, count);
            // var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            double maxDiff = 0;

            for (size_t i = 0; i < Nx * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                double diff = std::abs(r64s[i] - decompressedR64s[i]);

                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }
            ++t;

            auto itMax = std::max_element(r64s.begin(), r64s.end());

            const double relativeMaxDiff = maxDiff / *itMax;
            ASSERT_LT(relativeMaxDiff, std::stod(tolerance));
            std::cout << "Relative Max Diff " << relativeMaxDiff << " tolerance " << tolerance
                      << "\n";
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void MGARDAccuracy3D(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;
    const size_t Ny = 50;
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
    const std::string fname("BPWRMGARD3D_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD3D_" + tolerance + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny, Nz};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0, 0};
        const adios2::Dims count{Nx, Ny, Nz};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        var_r32.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});
        var_r64.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            // bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            //        auto var_r32 = io.InquireVariable<float>("r32");
            //        EXPECT_TRUE(var_r32);
            //        ASSERT_EQ(var_r32.ShapeID(),
            //        adios2::ShapeID::GlobalArray); ASSERT_EQ(var_r32.Steps(),
            //        NSteps); ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            //        ASSERT_EQ(var_r32.Shape()[1], Ny);
            //        ASSERT_EQ(var_r32.Shape()[2], Nz);

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
            // var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            double maxDiff = 0;

            for (size_t i = 0; i < Nx * Ny * Nz; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                double diff = std::abs(r64s[i] - decompressedR64s[i]);

                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }
            ++t;

            auto itMax = std::max_element(r64s.begin(), r64s.end());

            const double relativeMaxDiff = maxDiff / *itMax;
            ASSERT_LT(relativeMaxDiff, std::stod(tolerance));
            std::cout << "Relative Max Diff " << relativeMaxDiff << " tolerance " << tolerance
                      << "\n";
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void MGARDAccuracy1DSel(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 1000;

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
    const std::string fname("BPWRMGARD1DSel_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD1DSel_" + tolerance + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        var_r32.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});
        var_r64.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            // bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            //        auto var_r32 = io.InquireVariable<float>("r32");
            //        EXPECT_TRUE(var_r32);
            //        ASSERT_EQ(var_r32.ShapeID(),
            //        adios2::ShapeID::GlobalArray); ASSERT_EQ(var_r32.Steps(),
            //        NSteps); ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

            const adios2::Dims start{mpiRank * Nx + Nx / 2};
            const adios2::Dims count{Nx / 2};
            const adios2::Box<adios2::Dims> sel(start, count);
            // var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx / 2; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                //                ASSERT_LT(std::abs(decompressedR32s[i] -
                //                r32s[Nx / 2 + i]),
                //                          tolerance)
                //                    << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[Nx / 2 + i]), std::stod(tolerance))
                    << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void MGARDAccuracy2DSel(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

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
    const std::string fname("BPWRMGARD2DSel_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD2DSel_" + tolerance + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
        const adios2::Dims count{Nx, Ny};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        var_r32.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});
        var_r64.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            // bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            //        auto var_r32 = io.InquireVariable<float>("r32");
            //        EXPECT_TRUE(var_r32);
            //        ASSERT_EQ(var_r32.ShapeID(),
            //        adios2::ShapeID::GlobalArray); ASSERT_EQ(var_r32.Steps(),
            //        NSteps); ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            //        ASSERT_EQ(var_r32.Shape()[1], Ny);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);

            const adios2::Dims start{mpiRank * Nx + Nx / 2, 0};
            const adios2::Dims count{Nx / 2, Ny};
            const adios2::Box<adios2::Dims> sel(start, count);
            // var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            for (size_t i = 0; i < Nx / 2 * Ny; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                //                ASSERT_LT(std::abs(decompressedR32s[i] -
                //                r32s[Nx / 2 * Ny + i]),
                //                          tolerance)
                //                    << msg;
                ASSERT_LT(std::abs(decompressedR64s[i] - r64s[Nx / 2 * Ny + i]),
                          std::stod(tolerance))
                    << msg;
            }
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void MGARDAccuracy3DSel(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

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
    const std::string fname("BPWRMGARD3DSel_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD3DSel_" + tolerance + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny, Nz};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0, 0};
        const adios2::Dims count{Nx, Ny, Nz};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        var_r32.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});
        var_r64.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            // bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);
        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {

            //        auto var_r32 = io.InquireVariable<float>("r32");
            //        EXPECT_TRUE(var_r32);
            //        ASSERT_EQ(var_r32.ShapeID(),
            //        adios2::ShapeID::GlobalArray); ASSERT_EQ(var_r32.Steps(),
            //        NSteps); ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            //        ASSERT_EQ(var_r32.Shape()[1], Ny);
            //        ASSERT_EQ(var_r32.Shape()[2], Nz);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);
            ASSERT_EQ(var_r64.Shape()[2], Nz);

            const adios2::Dims start{mpiRank * Nx + Nx / 2, 0, 0};
            const adios2::Dims count{Nx / 2, Ny, Nz};
            const adios2::Box<adios2::Dims> sel(start, count);
            // var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            double maxDiff = 0;

            for (size_t i = 0; i < Nx / 2 * Ny * Nz; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                //                ASSERT_LT(
                //                    std::abs(decompressedR32s[i] - r32s[Nx / 2
                //                    * Ny * Nz + i]),
                //                    tolerance)
                //                    << msg;
                double diff = std::abs(r64s[Nx / 2 * Ny * Nz + i] - decompressedR64s[i]);

                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }

            auto itMax = std::max_element(r64s.begin(), r64s.end());
            ASSERT_LT(maxDiff / *itMax, std::stod(tolerance));

            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

void MGARDAccuracy2DSmallSel(const std::string tolerance)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

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
    const std::string fname("BPWRMGARD2DSmallSel_" + tolerance + "_MPI.bp");
#else
    const std::string fname("BPWRMGARD2DSmallSel_" + tolerance + ".bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
        const adios2::Dims count{Nx, Ny};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);

        var_r32.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});
        var_r64.AddOperation("MGARDPlus", {{adios2::ops::mgard::key::tolerance, tolerance}});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            // bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        unsigned int t = 0;
        std::vector<float> decompressedR32s;
        std::vector<double> decompressedR64s;

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            //        auto var_r32 = io.InquireVariable<float>("r32");
            //        EXPECT_TRUE(var_r32);
            //        ASSERT_EQ(var_r32.ShapeID(),
            //        adios2::ShapeID::GlobalArray); ASSERT_EQ(var_r32.Steps(),
            //        NSteps); ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);
            //        ASSERT_EQ(var_r32.Shape()[1], Ny);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_r64.Shape()[1], Ny);

            const adios2::Dims start{static_cast<std::size_t>(mpiRank) * Nx + 1, 1};
            const adios2::Dims count{2, 2};
            const adios2::Box<adios2::Dims> sel(start, count);
            // var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // bpReader.Get(var_r32, decompressedR32s);
            bpReader.Get(var_r64, decompressedR64s);
            bpReader.EndStep();

            // ASSERT_LT(std::abs(decompressedR32s[0] - 0.06), tolerance);
            ASSERT_LT(std::abs(decompressedR64s[0] - 0.06), std::stod(tolerance));

            // ASSERT_LT(std::abs(decompressedR32s[1] - 0.07), tolerance);
            ASSERT_LT(std::abs(decompressedR64s[1] - 0.07), std::stod(tolerance));

            // ASSERT_LT(std::abs(decompressedR32s[2] - 0.11), tolerance);
            ASSERT_LT(std::abs(decompressedR64s[2] - 0.11), std::stod(tolerance));

            // ASSERT_LT(std::abs(decompressedR32s[3] - 0.12), tolerance);
            ASSERT_LT(std::abs(decompressedR64s[3] - 0.12), std::stod(tolerance));

            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
    }
}

class BPWriteReadMGARDPlus : public ::testing::TestWithParam<std::string>
{
public:
    BPWriteReadMGARDPlus() = default;
    virtual void SetUp(){};
    virtual void TearDown(){};
};

TEST_P(BPWriteReadMGARDPlus, BPWRMGARD2D) { MGARDAccuracy2D(GetParam()); }

TEST_P(BPWriteReadMGARDPlus, BPWRMGARD3D) { MGARDAccuracy3D(GetParam()); }

INSTANTIATE_TEST_SUITE_P(MGARDAccuracy, BPWriteReadMGARDPlus,
                         ::testing::Values("0.01", "0.001", "0.0001", "0.00001"));

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    MPI_Init(nullptr, nullptr);
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
