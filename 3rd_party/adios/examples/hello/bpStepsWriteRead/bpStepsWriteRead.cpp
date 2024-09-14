/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpStepsWriteRead.cpp  Simple example of writing and reading data through ADIOS2 BP engine with
 * multiple simulations steps for every IO step.
 *
 *  Created on: Feb 16, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <algorithm> // std::for_each
#include <ios>       // std::ios_base::failure
#include <iostream>  // std::cout
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>

void update_array(std::vector<float> &array, int val)
{
    std::transform(array.begin(), array.end(), array.begin(),
                   [val](float v) -> float { return v + static_cast<float>(val); });
}

void writer(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
            const size_t Nx, unsigned int nSteps, int rank, int size)
{
    std::vector<float> simData(Nx, 0);

    adios2::IO bpIO = adios.DeclareIO("WriteIO");
    bpIO.SetEngine(engine);

    const adios2::Dims shape{static_cast<size_t>(size * Nx)};
    const adios2::Dims start{static_cast<size_t>(rank * Nx)};
    const adios2::Dims count{Nx};
    auto bpFloats = bpIO.DefineVariable<float>("bpFloats", shape, start, count);

    auto bpStep = bpIO.DefineVariable<unsigned int>("bpStep");

    adios2::Engine bpWriter = bpIO.Open(fname, adios2::Mode::Write);

    for (unsigned int step = 0; step < nSteps; ++step)
    {
        const adios2::Box<adios2::Dims> sel({0}, {Nx});
        bpFloats.SetSelection(sel);

        bpWriter.BeginStep();
        bpWriter.Put(bpFloats, simData.data());
        bpWriter.Put(bpStep, step);
        bpWriter.EndStep();

        // Update values in the simulation data
        update_array(simData, 10);
    }

    bpWriter.Close();
}

void reader(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
            const size_t Nx, unsigned int /*nSteps*/, int rank, int /*size*/)
{
    adios2::IO bpIO = adios.DeclareIO("ReadIO");
    bpIO.SetEngine(engine);

    adios2::Engine bpReader = bpIO.Open(fname, adios2::Mode::Read);

    std::vector<float> simData(Nx, 0);
    unsigned int inStep = 0;
    for (unsigned int step = 0; bpReader.BeginStep() == adios2::StepStatus::OK; ++step)
    {
        auto bpFloats = bpIO.InquireVariable<float>("bpFloats");
        if (bpFloats)
        {
            const adios2::Box<adios2::Dims> sel({{Nx * rank}, {Nx}});
            bpFloats.SetSelection(sel);
            bpReader.Get(bpFloats, simData.data());
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

int main(int argc, char *argv[])
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
    std::cout << "Using engine " << engine << std::endl;

    const std::string filename = engine + "StepsWriteRead.bp";
    const unsigned int nSteps = 10;
    const unsigned int Nx = 60000;
    try
    {
        /** ADIOS class factory of IO class objects */
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        writer(adios, engine, filename, Nx, nSteps, rank, size);
        reader(adios, engine, filename, Nx, nSteps, rank, size);
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM "
                     "from rank "
                  << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
