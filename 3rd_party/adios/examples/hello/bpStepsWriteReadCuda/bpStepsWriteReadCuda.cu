/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpStepsWriteReadCuda.cu  Simple example of writing and reading data through ADIOS2 BP engine with
 * multiple simulations steps for every IO step using CUDA
 */

#include <ios>
#include <iostream>
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <vector>

#include <adios2.h>

#include <cuda_runtime.h>

__global__ void update_array(float *vect, int val) { vect[blockIdx.x] += val; }

void writer(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
            const size_t Nx, unsigned int nSteps)
{
    // Initialize the simulation data
    float *gpuSimData;
    cudaMalloc(&gpuSimData, Nx * sizeof(float));
    cudaMemset(gpuSimData, 0, Nx);

    // Set up the ADIOS structures
    adios2::IO bpIO = adios.DeclareIO("WriteIO");
    bpIO.SetEngine(engine);

    // Declare an array for the ADIOS data of size (NumOfProcesses * Nx)
    const adios2::Dims shape{static_cast<size_t>(Nx)};
    const adios2::Dims start{static_cast<size_t>(0)};
    const adios2::Dims count{Nx};
    auto bpFloats = bpIO.DefineVariable<float>("bpFloats", shape, start, count);
    auto bpStep = bpIO.DefineVariable<unsigned int>("bpStep");

    adios2::Engine bpWriter = bpIO.Open(fname, adios2::Mode::Write);

    // Simulation steps
    for (unsigned int step = 0; step < nSteps; ++step)
    {
        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces
        const adios2::Box<adios2::Dims> sel({0}, {Nx});
        bpFloats.SetSelection(sel);

        // Start IO step every write step
        bpWriter.BeginStep();
        bpFloats.SetMemorySpace(adios2::MemorySpace::GPU);
        bpWriter.Put(bpFloats, gpuSimData);
        bpWriter.Put(bpStep, step);
        bpWriter.EndStep();

        // Update values in the simulation data
        update_array<<<Nx, 1>>>(gpuSimData, 10);
    }

    bpWriter.Close();
}

void reader(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
            const size_t Nx, unsigned int /*nSteps*/)
{
    // Create ADIOS structures
    adios2::IO bpIO = adios.DeclareIO("ReadIO");
    bpIO.SetEngine(engine);

    adios2::Engine bpReader = bpIO.Open(fname, adios2::Mode::Read);

    unsigned int inStep = 0;
    float *gpuSimData;
    cudaMalloc(&gpuSimData, Nx * sizeof(float));
    cudaMemset(gpuSimData, 0, Nx);
    for (unsigned int step = 0; bpReader.BeginStep() == adios2::StepStatus::OK; ++step)
    {
        auto bpFloats = bpIO.InquireVariable<float>("bpFloats");
        if (bpFloats)
        {
            const adios2::Dims start{0};
            const adios2::Dims count{Nx};
            const adios2::Box<adios2::Dims> sel(start, count);
            bpFloats.SetSelection(sel);

            bpFloats.SetMemorySpace(adios2::MemorySpace::GPU);
            bpReader.Get(bpFloats, gpuSimData); //, adios2::Mode::Deferred);
        }
        auto bpStep = bpIO.InquireVariable<unsigned int>("bpStep");
        if (bpStep)
        {
            bpReader.Get(bpStep, &inStep);
        }

        bpReader.EndStep();
        if (inStep != step)
        {
            std::cout << "ERROR: step mismatch\n";
            return;
        }
    }
    bpReader.Close();
}

int main(int argc, char **argv)
{
    const int device_id = 1;
    cudaSetDevice(device_id);

    const std::string engine = argv[1] ? argv[1] : "BPFile";
    std::cout << "Using engine " << engine << std::endl;

    const std::string filename = engine + "StepsWriteReadCuda.bp";
    const unsigned int nSteps = 10;
    const unsigned int Nx = 6000;
    try
    {
        /** ADIOS class factory of IO class objects */
        adios2::ADIOS adios;

        writer(adios, engine, filename, Nx, nSteps);
        reader(adios, engine, filename, Nx, nSteps);
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM\n";
        std::cout << e.what() << "\n";
    }

    return 0;
}
