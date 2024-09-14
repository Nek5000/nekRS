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

std::string engineName; // comes from command line

void ZfpRate1D(const double rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("BPWriteReadZfp1D_" + std::to_string(rate) + "_hl.bp");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;

    // Number of steps
    const size_t NSteps = 3;

    std::vector<float> r32s(Nx);
    std::vector<double> r64s(Nx);

    // range 0 to 999
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
    const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
    const adios2::Dims count{Nx};

    // Writer
    {
#if ADIOS2_USE_MPI
        adios2::fstream fw(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fw(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            fw.write("r32", r32s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(rate)}}}});
            fw.write("r64", r64s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(2 * rate)}}}}, adios2::end_step);
        }
        fw.close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // reader
    {
#if ADIOS2_USE_MPI
        adios2::fstream fr(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fr(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        adios2::fstep fs;
        while (adios2::getstep(fr, fs))
        {
            std::vector<float> decompressedR32s = fs.read<float>("r32", start, count);

            std::vector<double> decompressedR64s = fs.read<double>("r64", start, count);

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

        fr.close();
    }
}

void ZfpRate2D(const double rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("BPWriteReadZfp2D_" + std::to_string(rate) + "_hl.bp");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 100;
    const size_t Ny = 50;

    // Number of steps
    const size_t NSteps = 3;

    std::vector<float> r32s(Nx * Ny);
    std::vector<double> r64s(Nx * Ny);

    // range 0 to 100*50
    std::iota(r32s.begin(), r32s.end(), 0.f);
    std::iota(r64s.begin(), r64s.end(), 0.);

    const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
    const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
    const adios2::Dims count{Nx, Ny};

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // writer
    {
#if ADIOS2_USE_MPI
        adios2::fstream fw(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fw(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            fw.write("r32", r32s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(rate)}}}});
            fw.write("r64", r64s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(rate)}}}}, adios2::end_step);
        }
        fw.close();
    }

    // reader
    {
#if ADIOS2_USE_MPI
        adios2::fstream fr(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fr(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        adios2::fstep fs;
        while (adios2::getstep(fr, fs))
        {
            std::vector<float> decompressedR32s = fs.read<float>("r32", start, count);

            std::vector<double> decompressedR64s = fs.read<double>("r64", start, count);

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

        fr.close();
    }
}

void ZfpRate3D(const double rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("BPWriteReadZfp3D_" + std::to_string(rate) + "_hl.bp");

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

    const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny, Nz};
    const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0, 0};
    const adios2::Dims count{Nx, Ny, Nz};

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // writer
    {
#if ADIOS2_USE_MPI
        adios2::fstream fw(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fw(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            fw.write("r32", r32s.data(), shape, start, count,
                     {{"zfp", {{"accuracy", std::to_string(rate)}}}});
            fw.write("r64", r64s.data(), shape, start, count,
                     {{"zfp", {{"accuracy", std::to_string(rate)}}}}, adios2::end_step);
        }
        fw.close();
    }

    // reader
    {
#if ADIOS2_USE_MPI
        adios2::fstream fr(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fr(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        adios2::fstep fs;
        while (adios2::getstep(fr, fs))
        {
            std::vector<float> decompressedR32s = fs.read<float>("r32", start, count);

            std::vector<double> decompressedR64s = fs.read<double>("r64", start, count);

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

        fr.close();
    }
}

void ZfpRate1DSel(const double rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("BPWriteReadZfp1DSel_" + std::to_string(rate) + "_hl.bp");

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

    const adios2::Dims shape{static_cast<std::size_t>(Nx * mpiSize)};
    const adios2::Dims start{static_cast<std::size_t>(Nx * mpiRank)};
    const adios2::Dims count{Nx};

    const adios2::Dims startSel{static_cast<std::size_t>(mpiRank) * Nx + Nx / 2};
    const adios2::Dims countSel{Nx / 2};

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Writer
    {
#if ADIOS2_USE_MPI
        adios2::fstream fw(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fw(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            fw.write("r32", r32s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(rate)}}}});
            fw.write("r64", r64s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(2 * rate)}}}}, adios2::end_step);
        }
        fw.close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // reader
    {
#if ADIOS2_USE_MPI
        adios2::fstream fr(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fr(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        adios2::fstep fs;
        while (adios2::getstep(fr, fs))
        {
            std::vector<float> decompressedR32s = fs.read<float>("r32", startSel, countSel);

            std::vector<double> decompressedR64s = fs.read<double>("r64", startSel, countSel);

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

        fr.close();
    }
}

void ZfpRate2DSel(const double rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("BPWriteReadZfp2DSel_" + std::to_string(rate) + "_hl.bp");

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

    const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
    const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
    const adios2::Dims count{Nx, Ny};

    const adios2::Dims startSel{mpiRank * Nx + Nx / 2, 0};
    const adios2::Dims countSel{Nx / 2, Ny};

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Writer
    {
#if ADIOS2_USE_MPI
        adios2::fstream fw(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fw(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            fw.write("r32", r32s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(rate)}}}});
            fw.write("r64", r64s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(2 * rate)}}}}, adios2::end_step);
        }
        fw.close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // reader
    {
#if ADIOS2_USE_MPI
        adios2::fstream fr(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fr(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        adios2::fstep fs;
        while (adios2::getstep(fr, fs))
        {
            std::vector<float> decompressedR32s = fs.read<float>("r32", startSel, countSel);

            std::vector<double> decompressedR64s = fs.read<double>("r64", startSel, countSel);

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

        fr.close();
    }
}

void ZfpRate3DSel(const double rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("BPWriteReadZfp3DSel_" + std::to_string(rate) + "_hl.bp");

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

    const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny, Nz};
    const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0, 0};
    const adios2::Dims count{Nx, Ny, Nz};

    const adios2::Dims startSel{static_cast<std::size_t>(mpiRank) * Nx + Nx / 2, 0, 0};
    const adios2::Dims countSel{Nx / 2, Ny, Nz};

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Writer
    {
#if ADIOS2_USE_MPI
        adios2::fstream fw(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fw(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            fw.write("r32", r32s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(rate)}}}});
            fw.write("r64", r64s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(2 * rate)}}}}, adios2::end_step);
        }
        fw.close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // reader
    {
#if ADIOS2_USE_MPI
        adios2::fstream fr(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fr(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        adios2::fstep fs;
        while (adios2::getstep(fr, fs))
        {
            std::vector<float> decompressedR32s = fs.read<float>("r32", startSel, countSel);

            std::vector<double> decompressedR64s = fs.read<double>("r64", startSel, countSel);

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

        fr.close();
    }
}

void ZfpRate2DSmallSel(const double rate)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname("BPWriteReadZfp2DSmallSel_" + std::to_string(rate) + "_hl.bp");

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 5;
    const size_t Ny = 5;

    // Number of steps
    const size_t NSteps = 1;

    const std::vector<float> r32s = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                                     0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
                                     0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24};
    const std::vector<double> r64s = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                                      0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
                                      0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24};

    const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), Ny};
    const adios2::Dims start{static_cast<size_t>(Nx * mpiRank), 0};
    const adios2::Dims count{Nx, Ny};

    const adios2::Dims startSel{static_cast<std::size_t>(mpiRank) * Nx + 1, 1};
    const adios2::Dims countSel{2, 2};

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Writer
    {
#if ADIOS2_USE_MPI
        adios2::fstream fw(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fw(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            fw.write("r32", r32s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(rate)}}}});
            fw.write("r64", r64s.data(), shape, start, count,
                     {{"zfp", {{"rate", std::to_string(2 * rate)}}}}, adios2::end_step);
        }
        fw.close();
    }

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // reader
    {
#if ADIOS2_USE_MPI
        adios2::fstream fr(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream fr(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        adios2::fstep fs;
        while (adios2::getstep(fr, fs))
        {
            std::vector<float> decompressedR32s = fs.read<float>("r32", startSel, countSel);

            std::vector<double> decompressedR64s = fs.read<double>("r64", startSel, countSel);

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

        fr.close();
    }
}

class BPWriteReadZfpHighLevelAPI : public ::testing::TestWithParam<double>
{
public:
    BPWriteReadZfpHighLevelAPI() = default;

    virtual void SetUp() {}
    virtual void TearDown() {}
};

TEST_P(BPWriteReadZfpHighLevelAPI, ADIOS2BPWriteReadZfp1D) { ZfpRate1D(GetParam()); }
TEST_P(BPWriteReadZfpHighLevelAPI, ADIOS2BPWriteReadZfp2D) { ZfpRate2D(GetParam()); }
TEST_P(BPWriteReadZfpHighLevelAPI, ADIOS2BPWriteReadZfp3D) { ZfpRate3D(GetParam()); }
TEST_P(BPWriteReadZfpHighLevelAPI, ADIOS2BPWriteReadZfp1DSel) { ZfpRate1DSel(GetParam()); }
TEST_P(BPWriteReadZfpHighLevelAPI, ADIOS2BPWriteReadZfp2DSel) { ZfpRate2DSel(GetParam()); }
TEST_P(BPWriteReadZfpHighLevelAPI, ADIOS2BPWriteReadZfp3DSel) { ZfpRate3DSel(GetParam()); }
TEST_P(BPWriteReadZfpHighLevelAPI, ADIOS2BPWriteReadZfp2DSmallSel)
{
    ZfpRate2DSmallSel(GetParam());
}

INSTANTIATE_TEST_SUITE_P(ZfpRate, BPWriteReadZfpHighLevelAPI, ::testing::Values(8., 9., 10));

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
