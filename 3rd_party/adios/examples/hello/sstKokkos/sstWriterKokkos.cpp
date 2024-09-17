/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * sstWriterKokkos.cpp  Simple example of writing bpFloats through ADIOS2 SST
 * engine with multiple simulations steps for every IO step using Kokkos
 */
#include <ios>
#include <iostream>
#include <vector>

#include <adios2.h>
#include <adios2/cxx11/KokkosView.h>

#include <Kokkos_Core.hpp>

template <class MemSpace, class ExecSpace>
int BPWrite(adios2::ADIOS &adios, const std::string fname, const size_t Nx, const size_t Ny,
            const size_t nSteps, const std::string engine)
{
    // Initialize the simulation data
    Kokkos::View<float **, MemSpace> gpuSimData("simBuffer", Nx, Ny);
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemSpace>::accessible, "");
    Kokkos::parallel_for(
        "initBuffer", Kokkos::RangePolicy<ExecSpace>(0, Nx), KOKKOS_LAMBDA(int i) {
            for (int j = 0; j < Ny; j++)
                gpuSimData(i, j) = static_cast<float>(i);
        });
    Kokkos::fence();

    adios2::IO io = adios.DeclareIO("WriteIO");
    io.SetEngine(engine);

    const adios2::Dims shape{Nx, Ny};
    const adios2::Dims start{0, 0};
    const adios2::Dims count{Nx, Ny};
    auto data = io.DefineVariable<float>("bpFloats", shape, start, count);

    adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

    // Simulation steps
    for (int step = 0; step < nSteps; ++step)
    {
        adios2::Box<adios2::Dims> sel({0, 0}, {Nx, Ny});
        data.SetSelection(sel);

        bpWriter.BeginStep();
        // var.SetMemorySpace(adios2::MemorySpace::GPU);
        bpWriter.Put(data, gpuSimData);
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
    std::cout << "Done writing on memory space: " << exe_space.name() << std::endl;
    return 0;
}

int main(int argc, char **argv)
{
    const std::string engine = argv[1] ? argv[1] : "SST";
    std::cout << "Using engine " << engine << std::endl;
    const size_t Nx = 600, Ny = 100, nSteps = 2;
    const std::string memorySpace = "Device";

    const std::string filename = engine + "StepsWriteReadKokkos";
    Kokkos::initialize(argc, argv);
    {
        adios2::ADIOS adios;

        std::cout << "Using engine " << engine << std::endl;
        if (memorySpace == "Device")
        {
            using mem_space = Kokkos::DefaultExecutionSpace::memory_space;
            std::cout << "Memory space: DefaultMemorySpace" << std::endl;
            BPWrite<mem_space, Kokkos::DefaultExecutionSpace>(adios, filename + "_DD.bp", Nx, Ny,
                                                              nSteps, engine);
        }
        else
        {
            std::cout << "Memory space: HostSpace" << std::endl;
            BPWrite<Kokkos::HostSpace, Kokkos::Serial>(adios, filename + "_HH.bp", Nx, Ny, nSteps,
                                                       engine);
        }
    }
    Kokkos::finalize();
    return 0;
}
