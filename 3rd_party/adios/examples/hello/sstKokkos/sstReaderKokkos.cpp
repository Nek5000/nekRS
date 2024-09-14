/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * sstReaderKokkos.cpp  Simple example of reading bpFloats through ADIOS2 SST
 * engine with multiple simulations steps for every IO step using Kokkos
 */
#include <ios>
#include <iostream>
#include <vector>

#include <adios2.h>
#include <adios2/cxx11/KokkosView.h>

#include <Kokkos_Core.hpp>

template <class MemSpace, class ExecSpace>
int BPRead(adios2::ADIOS &adios, const std::string fname, const size_t Nx, const size_t Ny,
           const size_t nSteps, const std::string engine)
{
    adios2::IO io = adios.DeclareIO("ReadIO");
    io.SetEngine(engine);

    ExecSpace exe_space;
    std::cout << "Read on memory space: " << exe_space.name() << std::endl;

    adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

    unsigned int step = 0;
    bool correctValues = true;
    Kokkos::View<float **, MemSpace> gpuSimData("simBuffer", Nx, Ny);
    for (; bpReader.BeginStep() == adios2::StepStatus::OK; ++step)
    {
        auto data = io.InquireVariable<float>("bpFloats");
        const adios2::Dims start{0, 0};
        const adios2::Dims count{Nx, Ny};
        const adios2::Box<adios2::Dims> sel(start, count);
        data.SetSelection(sel);

        // var.SetMemorySpace(adios2::MemorySpace::GPU);
        bpReader.Get(data, gpuSimData);
        bpReader.EndStep();

        auto cpuData = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gpuSimData);
        if (cpuData(0, 0) != step * 10)
        {
            std::cout << "Value mismatch at step " << step << std::endl;
            correctValues = false;
            break;
        }
    }
    if (correctValues)
        std::cout << "Read " << step << " steps successfully" << std::endl;

    bpReader.Close();
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
            BPRead<mem_space, Kokkos::DefaultExecutionSpace>(adios, filename + "_DD.bp", Nx, Ny,
                                                             nSteps, engine);
        }
        else
        {
            std::cout << "Memory space: HostSpace" << std::endl;
            BPRead<Kokkos::HostSpace, Kokkos::Serial>(adios, filename + "_HH.bp", Nx, Ny, nSteps,
                                                      engine);
        }
    }
    Kokkos::finalize();
    return 0;
}
