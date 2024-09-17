/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * datamanWriterKokkos.cpp  Simple example of writing multiple steps of a 2D float Kokkos::View
 * through ADIOS2 DataMan
 */
#include <adios2.h>
#include <adios2/cxx11/KokkosView.h>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <thread>
#include <vector>

#include <Kokkos_Core.hpp>

size_t Nx = 10;
size_t Ny = 10;
size_t steps = 2;
adios2::Dims shape;
adios2::Dims start;
adios2::Dims count;

int mpiRank, mpiSize;

template <class T, class MemSpace>
void PrintData(Kokkos::View<T **, MemSpace> &gpuData, const size_t step)
{
    auto data = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gpuData);
    std::cout << "Rank: " << mpiRank << " Step: " << step << " [";
    for (int i = 0; i < data.extent_int(0); ++i)
        for (int j = 0; j < data.extent_int(1); ++j)
            std::cout << data(i, j) << " ";
    std::cout << "]" << std::endl;
}

template <class T, class MemSpace, class ExecSpace>
Kokkos::View<T **, MemSpace> GenerateData(const size_t step, const size_t Ny, const size_t mpiRank)
{
    Kokkos::View<T **, MemSpace> gpuSimData("simBuffer", Nx, Ny);
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemSpace>::accessible, "");
    Kokkos::parallel_for(
        "initBuffer", Kokkos::RangePolicy<ExecSpace>(0, Nx), KOKKOS_LAMBDA(int i) {
            for (int j = 0; j < Ny; j++)
                gpuSimData(i, j) = static_cast<float>(i * Ny + j) + mpiRank * 10000 + step;
        });
    Kokkos::fence();
    ExecSpace exe_space;
    std::cout << "Create data for step " << step << " on memory space: " << exe_space.name()
              << std::endl;
    return gpuSimData;
}

int main(int argc, char *argv[])
{
    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    // initialize data dimensions
    count = {Nx, Ny};
    start = {mpiRank * Nx, 0};
    shape = {mpiSize * Nx, Ny};

    // initialize adios2
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO dataManIO = adios.DeclareIO("whatever");
    dataManIO.SetEngine("DataMan");
    dataManIO.SetParameters({{"IPAddress", "127.0.0.1"},
                             {"Port", "12306"},
                             {"Timeout", "5"},
                             {"RendezvousReaderCount", "1"}});

    // open stream
    adios2::Engine dataManWriter = dataManIO.Open("HelloDataMan", adios2::Mode::Write);

    // define variable
    auto floatArrayVar = dataManIO.DefineVariable<float>("FloatArray", shape, start, count);

    // write data
    for (size_t i = 0; i < steps; ++i)
    {
        auto floatVector = GenerateData<float, Kokkos::DefaultExecutionSpace::memory_space,
                                        Kokkos::DefaultExecutionSpace>(i, Ny, mpiRank);
        dataManWriter.BeginStep();
        dataManWriter.Put(floatArrayVar, floatVector, adios2::Mode::Sync);
        PrintData(floatVector, dataManWriter.CurrentStep());
        dataManWriter.EndStep();
    }

    dataManWriter.Close();
    MPI_Finalize();

    return 0;
}
