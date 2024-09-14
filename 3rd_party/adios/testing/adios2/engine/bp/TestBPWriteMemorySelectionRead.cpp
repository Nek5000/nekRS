/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line
bool DoWrite = true;
bool DoRead = true;

namespace
{

template <class T>
inline void AssignStep1D(const size_t step, std::vector<T> &vector, const size_t ghostCells = 0)
{
    std::for_each(vector.begin() + ghostCells, vector.end() - ghostCells,
                  [step](T &value) { value = static_cast<T>(step); });
}

template <>
void AssignStep1D(const size_t step, std::vector<std::complex<float>> &vector,
                  const size_t ghostCells)
{
    std::for_each(
        vector.begin() + ghostCells, vector.end() - ghostCells, [step](std::complex<float> &value) {
            value = std::complex<float>(static_cast<float>(step), static_cast<float>(step));
        });
}

template <>
void AssignStep1D(const size_t step, std::vector<std::complex<double>> &vector,
                  const size_t ghostCells)
{
    std::for_each(vector.begin() + ghostCells, vector.end() - ghostCells,
                  [step](std::complex<double> &value) {
                      value = std::complex<double>(static_cast<double>(step),
                                                   static_cast<double>(step));
                  });
}

template <class T>
inline void AssignStep2D(const size_t step, std::vector<T> &vector, const size_t Nx,
                         const size_t Ny, const size_t ghostCellsX, const size_t ghostCellsY)
{
    for (size_t j = ghostCellsY; j < Ny + ghostCellsY; ++j)
    {
        const size_t indexJ = j * (Nx + 2 * ghostCellsX);

        for (size_t i = ghostCellsX; i < Nx + ghostCellsX; ++i)
        {
            const size_t index = indexJ + i;
            vector[index] = static_cast<T>(step);
        }
    }
}

template <>
void AssignStep2D(const size_t step, std::vector<std::complex<float>> &vector, const size_t Nx,
                  const size_t Ny, const size_t ghostCellsX, const size_t ghostCellsY)
{
    for (size_t j = ghostCellsY; j < Ny + ghostCellsY; ++j)
    {
        const size_t indexJ = j * (Nx + 2 * ghostCellsX);

        for (size_t i = ghostCellsX; i < Nx + ghostCellsX; ++i)
        {
            const size_t index = indexJ + i;
            vector[index] = std::complex<float>(static_cast<float>(step), static_cast<float>(step));
        }
    }
}

template <>
void AssignStep2D(const size_t step, std::vector<std::complex<double>> &vector, const size_t Nx,
                  const size_t Ny, const size_t ghostCellsX, const size_t ghostCellsY)
{
    for (size_t j = ghostCellsY; j < Ny + ghostCellsY; ++j)
    {
        const size_t indexJ = j * (Nx + 2 * ghostCellsX);

        for (size_t i = ghostCellsX; i < Nx + ghostCellsX; ++i)
        {
            const size_t index = indexJ + i;
            vector[index] =
                std::complex<double>(static_cast<double>(step), static_cast<double>(step));
        }
    }
}

template <class T>
inline void AssignStep3D(const size_t step, std::vector<T> &vector, const size_t Nx,
                         const size_t Ny, const size_t Nz, const size_t ghostCellsX,
                         const size_t ghostCellsY, const size_t ghostCellsZ)
{
    for (size_t k = ghostCellsZ; k < Nz + ghostCellsZ; ++k)
    {
        const size_t indexK = k * (Ny + 2 * ghostCellsY) * (Nx + 2 * ghostCellsX);

        for (size_t j = ghostCellsY; j < Ny + ghostCellsY; ++j)
        {
            const size_t indexJ = j * (Nx + 2 * ghostCellsX);

            for (size_t i = ghostCellsX; i < Nx + ghostCellsX; ++i)
            {
                const size_t index = indexK + indexJ + i;
                vector[index] = static_cast<T>(step);
            }
        }
    }
}

template <>
void AssignStep3D(const size_t step, std::vector<std::complex<float>> &vector, const size_t Nx,
                  const size_t Ny, const size_t Nz, const size_t ghostCellsX,
                  const size_t ghostCellsY, const size_t ghostCellsZ)
{
    for (size_t k = ghostCellsZ; k < Nz + ghostCellsZ; ++k)
    {
        const size_t indexK = k * (Ny + 2 * ghostCellsY) * (Nx + 2 * ghostCellsX);

        for (size_t j = ghostCellsY; j < Ny + ghostCellsY; ++j)
        {
            const size_t indexJ = j * (Nx + 2 * ghostCellsX);

            for (size_t i = ghostCellsX; i < Nx + ghostCellsX; ++i)
            {
                const size_t index = indexK + indexJ + i;
                vector[index] =
                    std::complex<float>(static_cast<float>(step), static_cast<float>(step));
            }
        }
    }
}

template <>
void AssignStep3D(const size_t step, std::vector<std::complex<double>> &vector, const size_t Nx,
                  const size_t Ny, const size_t Nz, const size_t ghostCellsX,
                  const size_t ghostCellsY, const size_t ghostCellsZ)
{
    for (size_t k = ghostCellsZ; k < Nz + ghostCellsZ; ++k)
    {
        const size_t indexK = k * (Ny + 2 * ghostCellsY) * (Nx + 2 * ghostCellsX);

        for (size_t j = ghostCellsY; j < Ny + ghostCellsY; ++j)
        {
            const size_t indexJ = j * (Nx + 2 * ghostCellsX);

            for (size_t i = ghostCellsX; i < Nx + ghostCellsX; ++i)
            {
                const size_t index = indexK + indexJ + i;
                vector[index] =
                    std::complex<double>(static_cast<double>(step), static_cast<double>(step));
            }
        }
    }
}

} // end anonymous namespace

#if ADIOS2_USE_MPI
MPI_Comm testComm;
#endif

void BPSteps1D(const size_t ghostCells)
{

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 10;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
    MPI_Comm_size(testComm, &mpiSize);
    const std::string fname("BPSteps1D_" + std::to_string(ghostCells) + "_MPI");
#else
    const std::string fname("BPSteps1D_" + std::to_string(ghostCells));
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    if (DoWrite)
    {
        adios2::IO io = adios.DeclareIO("WriteIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        io.SetParameters("StatsLevel=1");
        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count);

        const adios2::Dims memoryStart = {ghostCells};
        const adios2::Dims memoryCount = {Nx + 2 * ghostCells};

        var_i8.SetMemorySelection({memoryStart, memoryCount});
        var_i16.SetMemorySelection({memoryStart, memoryCount});
        var_i32.SetMemorySelection({memoryStart, memoryCount});
        var_i64.SetMemorySelection({memoryStart, memoryCount});
        var_r32.SetMemorySelection({memoryStart, memoryCount});
        var_r64.SetMemorySelection({memoryStart, memoryCount});
        var_cr32.SetMemorySelection({memoryStart, memoryCount});
        var_cr64.SetMemorySelection({memoryStart, memoryCount});

        std::vector<int8_t> dataI8(Nx + 2 * ghostCells, -1);
        std::vector<int16_t> dataI16(Nx + 2 * ghostCells, -1);
        std::vector<int32_t> dataI32(Nx + 2 * ghostCells, -1);
        std::vector<int64_t> dataI64(Nx + 2 * ghostCells, -1);
        std::vector<float> dataR32(Nx + 2 * ghostCells, -1.f);
        std::vector<double> dataR64(Nx + 2 * ghostCells, -1.);
        std::vector<std::complex<float>> dataCR32(Nx + 2 * ghostCells, {-1.f, -1.f});
        std::vector<std::complex<double>> dataCR64(Nx + 2 * ghostCells, {-1., -1.});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t i = 0; i < NSteps; ++i)
        {
            AssignStep1D(i, dataI8, ghostCells);
            AssignStep1D(i, dataI16, ghostCells);
            AssignStep1D(i, dataI32, ghostCells);
            AssignStep1D(i, dataI64, ghostCells);
            AssignStep1D(i, dataR32, ghostCells);
            AssignStep1D(i, dataR64, ghostCells);
            AssignStep1D(i, dataCR32, ghostCells);
            AssignStep1D(i, dataCR64, ghostCells);

            bpWriter.BeginStep();
            bpWriter.Put(var_i8, dataI8.data());
            bpWriter.Put(var_i16, dataI16.data());
            bpWriter.Put(var_i32, dataI32.data());
            bpWriter.Put(var_i64, dataI64.data());
            bpWriter.Put(var_r32, dataR32.data());
            bpWriter.Put(var_r64, dataR64.data());
            bpWriter.Put(var_cr32, dataCR32.data());
            bpWriter.Put(var_cr64, dataCR64.data());
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }
#if ADIOS2_USE_MPI
    MPI_Barrier(testComm);
#endif
    // Reader
    if (DoRead)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);
        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            const size_t step = bpReader.CurrentStep();

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx);
            ASSERT_EQ(var_i8.Min(), static_cast<int8_t>(step));

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx);

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx);

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx);

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            EXPECT_TRUE(var_cr32);
            ASSERT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr32.Shape()[0], mpiSize * Nx);

            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
            EXPECT_TRUE(var_cr64);
            ASSERT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr64.Shape()[0], mpiSize * Nx);
            std::vector<int8_t> I8;
            std::vector<int16_t> I16;
            std::vector<int32_t> I32;
            std::vector<int64_t> I64;
            std::vector<float> R32;
            std::vector<double> R64;
            std::vector<std::complex<float>> CR32;
            std::vector<std::complex<double>> CR64;

            bpReader.Get(var_i8, I8);
            bpReader.Get(var_i16, I16);
            bpReader.Get(var_i32, I32);
            bpReader.Get(var_i64, I64);
            bpReader.Get(var_r32, R32);
            bpReader.Get(var_r64, R64);
            bpReader.Get(var_cr32, CR32);
            bpReader.Get(var_cr64, CR64);
            bpReader.EndStep();

            EXPECT_EQ(I8.size(), mpiSize * Nx);
            EXPECT_EQ(I16.size(), mpiSize * Nx);
            EXPECT_EQ(I32.size(), mpiSize * Nx);
            EXPECT_EQ(I64.size(), mpiSize * Nx);
            EXPECT_EQ(R32.size(), mpiSize * Nx);
            EXPECT_EQ(R64.size(), mpiSize * Nx);
            EXPECT_EQ(CR32.size(), mpiSize * Nx);
            EXPECT_EQ(CR64.size(), mpiSize * Nx);

            EXPECT_EQ(I8.front(), static_cast<int8_t>(step));
            EXPECT_EQ(I16.front(), static_cast<int16_t>(step));
            EXPECT_EQ(I32.front(), static_cast<int32_t>(step));
            EXPECT_EQ(I64.front(), static_cast<int64_t>(step));
            EXPECT_EQ(R32.front(), static_cast<float>(step));
            EXPECT_EQ(R64.front(), static_cast<double>(step));
            EXPECT_EQ(CR32.front(),
                      std::complex<float>(static_cast<float>(step), static_cast<float>(step)));
            EXPECT_EQ(CR64.front(),
                      std::complex<double>(static_cast<double>(step), static_cast<double>(step)));

            EXPECT_EQ(std::adjacent_find(I8.begin(), I8.end(), std::not_equal_to<int8_t>()),
                      I8.end());
            EXPECT_EQ(std::adjacent_find(I16.begin(), I16.end(), std::not_equal_to<int16_t>()),
                      I16.end());
            EXPECT_EQ(std::adjacent_find(I32.begin(), I32.end(), std::not_equal_to<int32_t>()),
                      I32.end());
            EXPECT_EQ(std::adjacent_find(I64.begin(), I64.end(), std::not_equal_to<int64_t>()),
                      I64.end());
            EXPECT_EQ(std::adjacent_find(R32.begin(), R32.end(), std::not_equal_to<float>()),
                      R32.end());
            EXPECT_EQ(std::adjacent_find(R64.begin(), R64.end(), std::not_equal_to<double>()),
                      R64.end());
            EXPECT_EQ(std::adjacent_find(CR32.begin(), CR32.end(),
                                         std::not_equal_to<std::complex<float>>()),
                      CR32.end());
            EXPECT_EQ(std::adjacent_find(CR64.begin(), CR64.end(),
                                         std::not_equal_to<std::complex<double>>()),
                      CR64.end());
        }

        bpReader.Close();
    }
}

void BPSteps2D4x2(const size_t ghostCells)
{

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 2;
    const size_t Ny = 4;

    const size_t ghostCellsX = ghostCells;
    const size_t ghostCellsY = 2 * ghostCells;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
    MPI_Comm_size(testComm, &mpiSize);
    const std::string fname("BPSteps2D4x2_" + std::to_string(ghostCells) + "_MPI");
#else
    const std::string fname("BPSteps2D4x2_" + std::to_string(ghostCells));
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    if (DoWrite)
    {
        adios2::IO io = adios.DeclareIO("WriteIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Ny * mpiSize), Nx};
        const adios2::Dims start{static_cast<size_t>(Ny * mpiRank), 0};
        const adios2::Dims count{Ny, Nx};

        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count);

        const adios2::Dims memoryStart = {ghostCellsY, ghostCellsX};
        const adios2::Dims memoryCount = {Ny + 2 * ghostCellsY, Nx + 2 * ghostCellsX};

        var_i8.SetMemorySelection({memoryStart, memoryCount});
        var_i16.SetMemorySelection({memoryStart, memoryCount});
        var_i32.SetMemorySelection({memoryStart, memoryCount});
        var_i64.SetMemorySelection({memoryStart, memoryCount});
        var_r32.SetMemorySelection({memoryStart, memoryCount});
        var_r64.SetMemorySelection({memoryStart, memoryCount});
        var_cr32.SetMemorySelection({memoryStart, memoryCount});
        var_cr64.SetMemorySelection({memoryStart, memoryCount});

        const size_t dataSize = (Ny + 2 * ghostCellsY) * (Nx + 2 * ghostCellsX);
        std::vector<int8_t> dataI8(dataSize, -1);
        std::vector<int16_t> dataI16(dataSize, -1);
        std::vector<int32_t> dataI32(dataSize, -1);
        std::vector<int64_t> dataI64(dataSize, -1);
        std::vector<float> dataR32(dataSize, -1.f);
        std::vector<double> dataR64(dataSize, -1.);
        std::vector<std::complex<float>> dataCR32(dataSize, {-1.f, -1.f});
        std::vector<std::complex<double>> dataCR64(dataSize, {-1., -1.});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t i = 0; i < NSteps; ++i)
        {
            AssignStep2D(i, dataI8, Nx, Ny, ghostCellsX, ghostCellsY);
            AssignStep2D(i, dataI16, Nx, Ny, ghostCellsX, ghostCellsY);
            AssignStep2D(i, dataI32, Nx, Ny, ghostCellsX, ghostCellsY);
            AssignStep2D(i, dataI64, Nx, Ny, ghostCellsX, ghostCellsY);
            AssignStep2D(i, dataR32, Nx, Ny, ghostCellsX, ghostCellsY);
            AssignStep2D(i, dataR64, Nx, Ny, ghostCellsX, ghostCellsY);
            AssignStep2D(i, dataCR32, Nx, Ny, ghostCellsX, ghostCellsY);
            AssignStep2D(i, dataCR64, Nx, Ny, ghostCellsX, ghostCellsY);

            bpWriter.BeginStep();
            bpWriter.Put(var_i8, dataI8.data());
            bpWriter.Put(var_i16, dataI16.data());
            bpWriter.Put(var_i32, dataI32.data());
            bpWriter.Put(var_i64, dataI64.data());
            bpWriter.Put(var_r32, dataR32.data());
            bpWriter.Put(var_r64, dataR64.data());
            bpWriter.Put(var_cr32, dataCR32.data(), adios2::Mode::Sync);
            bpWriter.Put(var_cr64, dataCR64.data(), adios2::Mode::Sync);
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }
#if ADIOS2_USE_MPI
    MPI_Barrier(testComm);
#endif
    // Reader
    if (DoRead)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i8.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_i8.Shape()[1], Nx);

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i16.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_i16.Shape()[1], Nx);

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i32.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_i32.Shape()[1], Nx);

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_i64.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_i64.Shape()[1], Nx);

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_r32.Shape()[1], Nx);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_r64.Shape()[1], Nx);

            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            EXPECT_TRUE(var_cr32);
            ASSERT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr32.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_cr32.Shape()[1], Nx);

            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
            EXPECT_TRUE(var_cr64);
            ASSERT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_cr64.Shape()[0], mpiSize * Ny);
            ASSERT_EQ(var_cr64.Shape()[1], Nx);

            std::vector<int8_t> I8;
            std::vector<int16_t> I16;
            std::vector<int32_t> I32;
            std::vector<int64_t> I64;
            std::vector<float> R32;
            std::vector<double> R64;
            std::vector<std::complex<float>> CR32;
            std::vector<std::complex<double>> CR64;

            bpReader.Get(var_i8, I8);
            bpReader.Get(var_i16, I16);
            bpReader.Get(var_i32, I32);
            bpReader.Get(var_i64, I64);
            bpReader.Get(var_r32, R32);
            bpReader.Get(var_r64, R64);
            bpReader.Get(var_cr32, CR32);
            bpReader.Get(var_cr64, CR64);
            bpReader.EndStep();

            const size_t dataSize = mpiSize * Ny * Nx;
            EXPECT_EQ(I8.size(), dataSize);
            EXPECT_EQ(I16.size(), dataSize);
            EXPECT_EQ(I32.size(), dataSize);
            EXPECT_EQ(I64.size(), dataSize);
            EXPECT_EQ(R32.size(), dataSize);
            EXPECT_EQ(R64.size(), dataSize);
            EXPECT_EQ(CR32.size(), dataSize);
            EXPECT_EQ(CR64.size(), dataSize);

            const size_t step = bpReader.CurrentStep();
            EXPECT_EQ(I8.front(), static_cast<int8_t>(step));
            EXPECT_EQ(I16.front(), static_cast<int16_t>(step));
            EXPECT_EQ(I32.front(), static_cast<int32_t>(step));
            EXPECT_EQ(I64.front(), static_cast<int64_t>(step));
            EXPECT_EQ(R32.front(), static_cast<float>(step));
            EXPECT_EQ(R64.front(), static_cast<double>(step));
            EXPECT_EQ(CR32.front(),
                      std::complex<float>(static_cast<float>(step), static_cast<float>(step)));
            EXPECT_EQ(CR64.front(),
                      std::complex<double>(static_cast<double>(step), static_cast<double>(step)));

            EXPECT_EQ(std::adjacent_find(I8.begin(), I8.end(), std::not_equal_to<int8_t>()),
                      I8.end());
            EXPECT_EQ(std::adjacent_find(I16.begin(), I16.end(), std::not_equal_to<int16_t>()),
                      I16.end());
            EXPECT_EQ(std::adjacent_find(I32.begin(), I32.end(), std::not_equal_to<int32_t>()),
                      I32.end());
            EXPECT_EQ(std::adjacent_find(I64.begin(), I64.end(), std::not_equal_to<int64_t>()),
                      I64.end());
            EXPECT_EQ(std::adjacent_find(R32.begin(), R32.end(), std::not_equal_to<float>()),
                      R32.end());
            EXPECT_EQ(std::adjacent_find(R64.begin(), R64.end(), std::not_equal_to<double>()),
                      R64.end());
            EXPECT_EQ(std::adjacent_find(CR32.begin(), CR32.end(),
                                         std::not_equal_to<std::complex<float>>()),
                      CR32.end());
            EXPECT_EQ(std::adjacent_find(CR64.begin(), CR64.end(),
                                         std::not_equal_to<std::complex<double>>()),
                      CR64.end());
        }

        bpReader.Close();
    }
}

void BPSteps3D8x2x4(const size_t ghostCells)
{

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 4;
    const size_t Ny = 2;
    const size_t Nz = 8;

    const size_t ghostCellsX = 2 * ghostCells;
    const size_t ghostCellsY = 4 * ghostCells;
    const size_t ghostCellsZ = ghostCells;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(testComm, &mpiRank);
    MPI_Comm_size(testComm, &mpiSize);
    const std::string fname("BPSteps3D8x2x4_" + std::to_string(ghostCells) + "_MPI");
#else
    const std::string fname("BPSteps3D8x2x4_" + std::to_string(ghostCells));
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(testComm);
#else
    adios2::ADIOS adios;
#endif
    if (DoWrite)
    {
        adios2::IO io = adios.DeclareIO("WriteIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nz * mpiSize), Ny, Nx};
        const adios2::Dims start{static_cast<size_t>(Nz * mpiRank), 0, 0};
        const adios2::Dims count{Nz, Ny, Nx};

        auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
        auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
        auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
        auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
        auto var_cr32 = io.DefineVariable<std::complex<float>>("cr32", shape, start, count);
        auto var_cr64 = io.DefineVariable<std::complex<double>>("cr64", shape, start, count);

        const adios2::Dims memoryStart = {ghostCellsZ, ghostCellsY, ghostCellsX};
        const adios2::Dims memoryCount = {Nz + 2 * ghostCellsZ, Ny + 2 * ghostCellsY,
                                          Nx + 2 * ghostCellsX};

        var_i8.SetMemorySelection({memoryStart, memoryCount});
        var_i16.SetMemorySelection({memoryStart, memoryCount});
        var_i32.SetMemorySelection({memoryStart, memoryCount});
        var_i64.SetMemorySelection({memoryStart, memoryCount});
        var_r32.SetMemorySelection({memoryStart, memoryCount});
        var_r64.SetMemorySelection({memoryStart, memoryCount});
        var_cr32.SetMemorySelection({memoryStart, memoryCount});
        var_cr64.SetMemorySelection({memoryStart, memoryCount});

        const size_t dataSize =
            (Nz + 2 * ghostCellsZ) * (Ny + 2 * ghostCellsY) * (Nx + 2 * ghostCellsX);
        std::vector<int8_t> dataI8(dataSize, -1);
        std::vector<int16_t> dataI16(dataSize, -1);
        std::vector<int32_t> dataI32(dataSize, -1);
        std::vector<int64_t> dataI64(dataSize, -1);
        std::vector<float> dataR32(dataSize, -1.f);
        std::vector<double> dataR64(dataSize, -1.);
        std::vector<std::complex<float>> dataCR32(dataSize, {-1.f, -1.f});
        std::vector<std::complex<double>> dataCR64(dataSize, {-1., -1.});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t i = 0; i < NSteps; ++i)
        {
            AssignStep3D(i, dataI8, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);
            AssignStep3D(i, dataI16, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);
            AssignStep3D(i, dataI32, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);
            AssignStep3D(i, dataI64, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);
            AssignStep3D(i, dataR32, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);
            AssignStep3D(i, dataR64, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);
            AssignStep3D(i, dataCR32, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);
            AssignStep3D(i, dataCR64, Nx, Ny, Nz, ghostCellsX, ghostCellsY, ghostCellsZ);

            bpWriter.BeginStep();
            bpWriter.Put(var_i8, dataI8.data());
            bpWriter.Put(var_i16, dataI16.data());
            bpWriter.Put(var_i32, dataI32.data());
            bpWriter.Put(var_i64, dataI64.data());
            bpWriter.Put(var_r32, dataR32.data());
            bpWriter.Put(var_r64, dataR64.data());
            bpWriter.Put(var_cr32, dataCR32.data());
            bpWriter.Put(var_cr64, dataCR64.data());
            bpWriter.EndStep();
        }
        bpWriter.Close();
    }
#if ADIOS2_USE_MPI
    MPI_Barrier(testComm);
#endif
    // Reader
    if (DoRead)
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            // const size_t currentStep = bpReader.CurrentStep();

            auto var_i8 = io.InquireVariable<int8_t>("i8");
            EXPECT_TRUE(var_i8);
            EXPECT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i8.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_i8.Shape()[1], Ny);
            EXPECT_EQ(var_i8.Shape()[2], Nx);
            // EXPECT_EQ(var_i8.Min(), static_cast<int8_t>(currentStep));
            // EXPECT_EQ(var_i8.Max(), static_cast<int8_t>(currentStep));

            auto var_i16 = io.InquireVariable<int16_t>("i16");
            EXPECT_TRUE(var_i16);
            EXPECT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i16.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_i16.Shape()[1], Ny);
            EXPECT_EQ(var_i16.Shape()[2], Nx);
            // EXPECT_EQ(var_i16.Min(), static_cast<int16_t>(currentStep));
            // EXPECT_EQ(var_i16.Max(), static_cast<int16_t>(currentStep));

            auto var_i32 = io.InquireVariable<int32_t>("i32");
            EXPECT_TRUE(var_i32);
            EXPECT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i32.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_i32.Shape()[1], Ny);
            EXPECT_EQ(var_i32.Shape()[2], Nx);
            // EXPECT_EQ(var_i32.Min(), static_cast<int32_t>(currentStep));
            // EXPECT_EQ(var_i32.Max(), static_cast<int32_t>(currentStep));

            auto var_i64 = io.InquireVariable<int64_t>("i64");
            EXPECT_TRUE(var_i64);
            EXPECT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_i64.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_i64.Shape()[1], Ny);
            EXPECT_EQ(var_i64.Shape()[2], Nx);
            // EXPECT_EQ(var_i64.Min(), static_cast<int64_t>(currentStep));
            // EXPECT_EQ(var_i64.Max(), static_cast<int64_t>(currentStep));

            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            EXPECT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r32.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_r32.Shape()[1], Ny);
            EXPECT_EQ(var_r32.Shape()[2], Nx);
            // EXPECT_EQ(var_r32.Min(), static_cast<float>(currentStep));
            // EXPECT_EQ(var_r32.Max(), static_cast<float>(currentStep));

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            EXPECT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_r64.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_r64.Shape()[1], Ny);
            EXPECT_EQ(var_r64.Shape()[2], Nx);
            // EXPECT_EQ(var_r64.Min(), static_cast<double>(currentStep));
            // EXPECT_EQ(var_r64.Max(), static_cast<double>(currentStep));

            auto var_cr32 = io.InquireVariable<std::complex<float>>("cr32");
            EXPECT_TRUE(var_cr32);
            EXPECT_EQ(var_cr32.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr32.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_cr32.Shape()[1], Ny);
            EXPECT_EQ(var_cr32.Shape()[2], Nx);
            // EXPECT_EQ(var_cr32.Min(),
            //          std::complex<float>(static_cast<float>(currentStep),
            //                              static_cast<float>(currentStep)));
            // EXPECT_EQ(var_cr32.Max(),
            //          std::complex<float>(static_cast<float>(currentStep),
            //                              static_cast<float>(currentStep)));

            auto var_cr64 = io.InquireVariable<std::complex<double>>("cr64");
            EXPECT_TRUE(var_cr64);
            EXPECT_EQ(var_cr64.ShapeID(), adios2::ShapeID::GlobalArray);
            EXPECT_EQ(var_cr64.Shape()[0], mpiSize * Nz);
            EXPECT_EQ(var_cr64.Shape()[1], Ny);
            EXPECT_EQ(var_cr64.Shape()[2], Nx);
            // EXPECT_EQ(var_cr64.Min(),
            //          std::complex<double>(static_cast<double>(currentStep),
            //                               static_cast<double>(currentStep)));
            // EXPECT_EQ(var_cr64.Max(),
            //          std::complex<double>(static_cast<double>(currentStep),
            //                               static_cast<double>(currentStep)));

            std::vector<int8_t> I8;
            std::vector<int16_t> I16;
            std::vector<int32_t> I32;
            std::vector<int64_t> I64;
            std::vector<float> R32;
            std::vector<double> R64;
            std::vector<std::complex<float>> CR32;
            std::vector<std::complex<double>> CR64;

            bpReader.Get(var_i8, I8);
            bpReader.Get(var_i16, I16);
            bpReader.Get(var_i32, I32);
            bpReader.Get(var_i64, I64);
            bpReader.Get(var_r32, R32);
            bpReader.Get(var_r64, R64);
            bpReader.Get(var_cr32, CR32);
            bpReader.Get(var_cr64, CR64);
            bpReader.EndStep();

            const size_t dataSize = mpiSize * Nz * Ny * Nx;
            EXPECT_EQ(I8.size(), dataSize);
            EXPECT_EQ(I16.size(), dataSize);
            EXPECT_EQ(I32.size(), dataSize);
            EXPECT_EQ(I64.size(), dataSize);
            EXPECT_EQ(R32.size(), dataSize);
            EXPECT_EQ(R64.size(), dataSize);
            EXPECT_EQ(CR32.size(), dataSize);
            EXPECT_EQ(CR64.size(), dataSize);

            const size_t step = bpReader.CurrentStep();
            EXPECT_EQ(I8.front(), static_cast<int8_t>(step));
            EXPECT_EQ(I16.front(), static_cast<int16_t>(step));
            EXPECT_EQ(I32.front(), static_cast<int32_t>(step));
            EXPECT_EQ(I64.front(), static_cast<int64_t>(step));
            EXPECT_EQ(R32.front(), static_cast<float>(step));
            EXPECT_EQ(R64.front(), static_cast<double>(step));
            EXPECT_EQ(CR32.front(),
                      std::complex<float>(static_cast<float>(step), static_cast<float>(step)));
            EXPECT_EQ(CR64.front(),
                      std::complex<double>(static_cast<double>(step), static_cast<double>(step)));

            EXPECT_EQ(std::adjacent_find(I8.begin(), I8.end(), std::not_equal_to<int8_t>()),
                      I8.end());
            EXPECT_EQ(std::adjacent_find(I16.begin(), I16.end(), std::not_equal_to<int16_t>()),
                      I16.end());
            EXPECT_EQ(std::adjacent_find(I32.begin(), I32.end(), std::not_equal_to<int32_t>()),
                      I32.end());
            EXPECT_EQ(std::adjacent_find(I64.begin(), I64.end(), std::not_equal_to<int64_t>()),
                      I64.end());
            EXPECT_EQ(std::adjacent_find(R32.begin(), R32.end(), std::not_equal_to<float>()),
                      R32.end());
            EXPECT_EQ(std::adjacent_find(R64.begin(), R64.end(), std::not_equal_to<double>()),
                      R64.end());
            EXPECT_EQ(std::adjacent_find(CR32.begin(), CR32.end(),
                                         std::not_equal_to<std::complex<float>>()),
                      CR32.end());
            EXPECT_EQ(std::adjacent_find(CR64.begin(), CR64.end(),
                                         std::not_equal_to<std::complex<double>>()),
                      CR64.end());
        }

        bpReader.Close();
    }
}

class BPWriteMemSelReadVector : public ::testing::TestWithParam<size_t>
{
public:
    BPWriteMemSelReadVector() = default;
    virtual void SetUp() {}
    virtual void TearDown() {}
};

TEST_P(BPWriteMemSelReadVector, BPMemorySelectionSteps1D) { BPSteps1D(GetParam()); }

TEST_P(BPWriteMemSelReadVector, BPMemorySelectionSteps2D4x2) { BPSteps2D4x2(GetParam()); }

TEST_P(BPWriteMemSelReadVector, BPMemorySelectionSteps3D4x2x8) { BPSteps3D8x2x4(GetParam()); }

INSTANTIATE_TEST_SUITE_P(ghostCells, BPWriteMemSelReadVector, ::testing::Values(1));

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    int bare_arg = 0;

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-do_write") == 0)
        {
            DoWrite = true;
            DoRead = false;
        }
        else if (strcmp(argv[i], "-do_read") == 0)
        {
            DoWrite = false;
            DoRead = true;
        }
        else if (argv[i][0] == '-')
        {
            std::cerr << "Unknown argument: " << argv[i] << std::endl;
            exit(1);
        }
        else
        {
            std::string fname;
            std::string engineParams;
            if (bare_arg == 0)
            {
                /* first arg without -- is engine */
                engineName = std::string(argv[1]);
                bare_arg++;
            }
            else if (bare_arg == 1)
            {
                /* second arg without -- is filename */
                //                fname = std::string(argv[1]);
                bare_arg++;
            }
            else if (bare_arg == 2)
            {
                //                engineParams = ParseEngineParams(argv[1]);
                bare_arg++;
            }
            else
            {

                throw std::invalid_argument("Unknown argument \"" + std::string(argv[1]) + "\"");
            }
        }
    }

#if ADIOS2_USE_MPI
    int provided;
    int thread_support_level =
        (engineName == "SST" || engineName == "sst") ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, thread_support_level, &provided);

    int key;
    MPI_Comm_rank(MPI_COMM_WORLD, &key);

    const unsigned int color = (DoRead & !DoWrite) ? 1 : 0;

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
