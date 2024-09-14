/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpStepsWriteReadKokkos.cpp  Simple example of writing and reading bpFloats through ADIOS2 BP
 * engine using Kokkos on any memory space
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

template <class MemSpace, class ExecSpace>
int BPWrite(const std::string fname, const size_t Nx, const size_t Ny, const size_t nSteps,
            const std::string engine, int rank, int size)
{
    // Initialize the simulation data
    Kokkos::View<float **, MemSpace> gpuSimData("simBuffer", Nx, Ny);
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemSpace>::accessible, "");
    Kokkos::parallel_for(
        "initBuffer", Kokkos::RangePolicy<ExecSpace>(0, Nx), KOKKOS_LAMBDA(int i) {
            for (int j = 0; j < Ny; j++)
                gpuSimData(i, j) = static_cast<float>(rank + i);
        });
    Kokkos::fence();

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("WriteIO");
    io.SetEngine(engine);

    const adios2::Dims shape{static_cast<size_t>(size * Nx), Ny};
    const adios2::Dims start{static_cast<size_t>(rank * Nx), 0};
    const adios2::Dims count{Nx, Ny};
    auto data = io.DefineVariable<float>("bpFloats", shape, start, count);

    adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

    // Simulation steps
    for (unsigned int step = 0; step < nSteps; ++step)
    {
        adios2::Box<adios2::Dims> sel({static_cast<size_t>(rank * Nx), 0}, {Nx, Ny});
        data.SetSelection(sel);

        bpWriter.BeginStep();
        bpWriter.Put(data, gpuSimData.data());
        bpWriter.EndStep();

        // Update values in the simulation data
        Kokkos::parallel_for(
            "updateBuffer", Kokkos::RangePolicy<ExecSpace>(0, Nx), KOKKOS_LAMBDA(int i) {
                for (int j = 0; j < Ny; j++)
                    gpuSimData(i, j) += 10;
            });
        Kokkos::fence();
    }

    bpWriter.Close();
    ExecSpace exe_space;
    if (rank == 0)
        std::cout << "Done writing on memory space: " << exe_space.name() << std::endl;
    return 0;
}

template <class MemSpace, class ExecSpace>
int BPRead(const std::string fname, const std::string engine, int rank, int size)
{
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("ReadIO");
    io.SetEngine(engine);

    ExecSpace exe_space;
    if (rank == 0)
        std::cout << "Read on memory space: " << exe_space.name() << std::endl;

    adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

    unsigned int step = 0;
    for (; bpReader.BeginStep() == adios2::StepStatus::OK; ++step)
    {
        auto data = io.InquireVariable<float>("bpFloats");
        if (data.Shape().size() != 2)
        {
            std::cout << "Error, the bpFloats variable in the BP file "
                         " on step "
                      << step << "needs to have two dimensions" << std::endl;
            break;
        }
        adios2::MemorySpace adiosMemSpace = adios2::MemorySpace::Host;
#ifdef ADIOS2_HAVE_GPU_SUPPORT
        if (!std::is_same<MemSpace, Kokkos::HostSpace>::value)
            adiosMemSpace = adios2::MemorySpace::GPU;
#endif
        auto dims = data.Shape(adiosMemSpace);
        size_t Nx = dims[0];
        size_t Ny = dims[1];

        if (Nx > Ny)
        {
            Nx = Nx / size;
            const adios2::Dims start{static_cast<size_t>(rank * Nx), 0};
            const adios2::Box<adios2::Dims> sel(start, {Nx, Ny});
            data.SetSelection(sel);
        }
        else
        {
            Ny = Ny / size;
            const adios2::Dims start{0, static_cast<size_t>(rank * Ny)};
            const adios2::Box<adios2::Dims> sel(start, {Nx, Ny});
            data.SetSelection(sel);
        }
        Kokkos::View<float **, MemSpace> gpuSimData("simBuffer", Nx, Ny);

        bpReader.Get(data, gpuSimData);
        bpReader.EndStep();

        auto cpuData = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gpuSimData);
        if (rank == 0)
        {
            std::cout << "Rank 0: " << std::endl;
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
    const std::string engine = argv[1] ? argv[1] : "BPFile";
    if (rank == 0)
        std::cout << "Using engine " << engine << std::endl;
    const size_t Nx = 3, Ny = 4, nSteps = 1;
    const std::string fromMemSpace = "Default";
    const std::string toMemSpace = "Default";

    const std::string filename = engine + "StepsWriteReadKokkos";
    Kokkos::initialize(argc, argv);
    {
        if (fromMemSpace == "Default")
        {
            using mem_space = Kokkos::DefaultExecutionSpace::memory_space;
            if (rank == 0)
                std::cout << "Writing on memory space: DefaultMemorySpace" << std::endl;
            BPWrite<mem_space, Kokkos::DefaultExecutionSpace>(filename + ".bp", Nx, Ny, nSteps,
                                                              engine, rank, size);
        }
        else
        {
            if (rank == 0)
                std::cout << "Writing on memory space: HostSpace" << std::endl;
            BPWrite<Kokkos::HostSpace, Kokkos::Serial>(filename + ".bp", nSteps, Nx, Ny, engine,
                                                       rank, size);
        }
        if (toMemSpace == "Default")
        {
            using mem_space = Kokkos::DefaultExecutionSpace::memory_space;
            if (rank == 0)
                std::cout << "Reading on memory space: DefaultMemorySpace" << std::endl;
            BPRead<mem_space, Kokkos::DefaultExecutionSpace>(filename + ".bp", engine, rank, size);
        }
        else
        {
            if (rank == 0)
                std::cout << "Reading on memory space: HostSpace" << std::endl;
            BPRead<Kokkos::HostSpace, Kokkos::Serial>(filename + ".bp", engine, rank, size);
        }
    }
    Kokkos::finalize();
#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
