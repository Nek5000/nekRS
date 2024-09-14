/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpStepsWriteReadHip.cpp  Simple example of writing and reading bpFloats through ADIOS2 BP engine
 * with multiple simulations steps for every IO step using HIP
 */
#include <algorithm>
#include <ios>
#include <iostream>
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>

#include <hip/hip_runtime.h>

__global__ void hip_initialize(float *vec) { vec[hipBlockIdx_x] = hipBlockIdx_x; }

__global__ void hip_increment(float *vec, float val) { vec[hipBlockIdx_x] += val; }

void writer(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
            const size_t Nx, unsigned int nSteps)
{
    hipError_t hipExit;
    float *gpuSimData;
    hipExit = hipMalloc((void **)&gpuSimData, Nx * sizeof(float));
    if (hipExit != hipSuccess)
    {
        std::cout << "[BPWrite] error: " << hipGetErrorString(hipExit) << std::endl;
        return;
    }
    hipLaunchKernelGGL(hip_initialize, dim3(Nx), dim3(1), 0, 0, gpuSimData);
    hipExit = hipDeviceSynchronize();
    if (hipExit != hipSuccess)
    {
        std::cout << "[BPWrite] error: " << hipGetErrorString(hipExit) << std::endl;
        return;
    }

    adios2::IO bpIO = adios.DeclareIO("WriteIO");
    bpIO.SetEngine(engine);

    const adios2::Dims shape{static_cast<size_t>(Nx)};
    const adios2::Dims start{static_cast<size_t>(0)};
    const adios2::Dims count{Nx};
    auto bpFloats = bpIO.DefineVariable<float>("bpFloats", shape, start, count);
    auto bpStep = bpIO.DefineVariable<unsigned int>("bpStep");

    adios2::Engine bpWriter = bpIO.Open(fname, adios2::Mode::Write);

    for (unsigned int step = 0; step < nSteps; ++step)
    {
        const adios2::Box<adios2::Dims> sel({0}, {Nx});
        bpFloats.SetSelection(sel);

        bpWriter.BeginStep();
        bpWriter.Put(bpFloats, gpuSimData);
        bpWriter.Put(bpStep, step);
        bpWriter.EndStep();

        hipLaunchKernelGGL(hip_increment, dim3(Nx), dim3(1), 0, 0, gpuSimData, 10);
        hipExit = hipDeviceSynchronize();
        if (hipExit != hipSuccess)
        {
            std::cout << "[BPWrite] error: " << hipGetErrorString(hipExit) << std::endl;
            return;
        }
    }

    bpWriter.Close();
}

void reader(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
            const size_t Nx, unsigned int /*nSteps*/)
{
    hipError_t hipExit;
    adios2::IO bpIO = adios.DeclareIO("ReadIO");
    bpIO.SetEngine(engine);

    adios2::Engine bpReader = bpIO.Open(fname, adios2::Mode::Read);

    unsigned int inStep = 0;
    float *gpuSimData;
    hipExit = hipMalloc((void **)&gpuSimData, Nx * sizeof(float));
    if (hipExit != hipSuccess)
    {
        std::cout << "[BPWrite] error: " << hipGetErrorString(hipExit) << std::endl;
        return;
    }
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
            bpReader.Get(bpFloats, gpuSimData);
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
    hipError_t hipExit;
    const int device_id = 0;
    hipExit = hipSetDevice(device_id);
    if (hipExit != hipSuccess)
    {
        std::cout << "[BPWrite] error: " << hipGetErrorString(hipExit) << std::endl;
        return 1;
    }

    const std::string engine = argv[1] ? argv[1] : "BPFile";
    std::cout << "Using engine " << engine << std::endl;

    const std::string filename = engine + "StepsWriteReadHip.bp";
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
