/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpWriteReadKokkosView.cpp  Simple example of writing and reading Kokkos Views through ADIOS2 BP
 * engine using different combinations of Layouts for the write/read buffers
 */
#include <ios>
#include <iostream>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <adios2.h>
#include <adios2/cxx11/KokkosView.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

template <class Layout>
int BPWrite(const std::string fname, const size_t Nx, const size_t Ny, const size_t nSteps,
            int rank, int size)
{
    Kokkos::View<float **, Layout> gpuSimData("simBuffer", Nx, Ny);
    Kokkos::parallel_for(
        "initBuffer", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {Nx, Ny}),
        KOKKOS_LAMBDA(int i, int j) { gpuSimData(i, j) = static_cast<float>(1 + i); });
    Kokkos::fence();

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("WriteIO");

    const adios2::Dims shape{static_cast<size_t>(size * Nx), Ny};
    const adios2::Dims start{static_cast<size_t>(rank * Nx), 0};
    const adios2::Dims count{Nx, Ny};
    auto data = io.DefineVariable<float>("bpFloats", shape, start, count);

    adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

    for (unsigned int step = 0; step < nSteps; ++step)
    {
        bpWriter.BeginStep();
        bpWriter.Put(data, gpuSimData);
        bpWriter.EndStep();

        Kokkos::parallel_for(
            "updateBuffer", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {Nx, Ny}),
            KOKKOS_LAMBDA(int i, int j) { gpuSimData(i, j) += 10; });
        Kokkos::fence();
    }

    bpWriter.Close();
    return 0;
}

template <class Layout>
int BPRead(const std::string fname, int rank, int size)
{
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("ReadIO");
    adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

    unsigned int step = 0;
    for (; bpReader.BeginStep() == adios2::StepStatus::OK; ++step)
    {
        auto data = io.InquireVariable<float>("bpFloats");
        if (data.Shape().size() != 2)
        {
            std::cout << "Error, the bpFloats variable on step " << step
                      << "needs to have two dimensions" << std::endl;
            break;
        }

        adios2::ArrayOrdering layout = adios2::ArrayOrdering::RowMajor;
        if (std::is_same<Layout, Kokkos::LayoutLeft>::value)
            layout = adios2::ArrayOrdering::ColumnMajor;
        auto dims = data.Shape(layout);
        size_t Nx = dims[0];
        size_t Ny = dims[1];

        Kokkos::View<float **, Layout> gpuSimData("simBuffer", Nx, Ny);
        const adios2::Dims start{0, 0};
        const adios2::Box<adios2::Dims> sel(start, {Nx, Ny});
        data.SetSelection(sel);

        bpReader.Get(data, gpuSimData);
        bpReader.EndStep();

        auto cpuData = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gpuSimData);
        // all processes are reading the entire global array, only rank 0 plots it
        if (rank == 0)
        {
            std::cout << "Rank " << rank << " receieved data:" << std::endl;
            for (size_t i = 0; i < Nx; i++)
            {
                for (size_t j = 0; j < Ny; j++)
                    std::cout << cpuData(i, j) << " ";
                std::cout << std::endl;
            }
        }
    }

    bpReader.Close();
    return 0;
}

int main(int argc, char **argv)
{
    int rank, size;
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    rank = 0;
    size = 1;
#endif
    const size_t Nx = 3, Ny = 4, nSteps = 1;

    const std::string filename = "KokkosViewWR";
    Kokkos::initialize(argc, argv);
    {
        std::cout << "Write on LayoutRight, read on LayoutRight" << std::endl;
        BPWrite<Kokkos::LayoutLeft>(filename + "_RR.bp", Nx, Ny, nSteps, rank, size);
        BPRead<Kokkos::LayoutLeft>(filename + "_RR.bp", rank, size);

        std::cout << "Write on LayoutRight, read on LayoutLeft" << std::endl;
        BPWrite<Kokkos::LayoutLeft>(filename + "_RL.bp", Nx, Ny, nSteps, rank, size);
        BPRead<Kokkos::LayoutRight>(filename + "_RL.bp", rank, size);

        std::cout << "Write on LayoutLeft, read on LayoutLeft" << std::endl;
        BPWrite<Kokkos::LayoutLeft>(filename + "_LL.bp", Nx, Ny, nSteps, rank, size);
        BPRead<Kokkos::LayoutLeft>(filename + "_LL.bp", rank, size);

        std::cout << "Write on LayoutLeft, read on LayoutRight" << std::endl;
        BPWrite<Kokkos::LayoutLeft>(filename + "_LR.bp", Nx, Ny, nSteps, rank, size);
        BPRead<Kokkos::LayoutRight>(filename + "_LR.bp", rank, size);
    }
    Kokkos::finalize();
#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
