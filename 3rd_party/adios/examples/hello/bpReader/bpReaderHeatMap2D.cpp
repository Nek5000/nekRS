/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpReaderHeatMap2D.cpp : Writes a heat map in a regular 2D mesh,
 * values grow from 0 in increments of 1
 *
 * temperature[gNx, Ny]
 * where: gNx = MPI_size_x * Nx
 *
 * 0                1       2  ...     Ny-1
 * Ny            Ny+1    Ny+2  ...   2*Ny-1
 * ...
 * ...
 * (gNx-1)*Ny   ...                  gNx*Ny-1
 *
 *
 *  Created on: Nov 1, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <algorithm>
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <mpi.h>
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>

int main(int argc, char *argv[])
{
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /** Application variable dimensions */
    constexpr std::size_t Nx = 10;
    constexpr std::size_t Ny = 10;

    const adios2::Dims count{Nx, Ny};
    const adios2::Dims start{rank * Nx, 0};
    const adios2::Dims shape{size * Nx, Ny};

    // populate local temperature values
    std::vector<unsigned int> temperatures(Nx * Ny);
    for (unsigned int i = 0; i < Nx; ++i)
    {
        const unsigned int iGlobal = static_cast<unsigned int>(start[0] + i);

        for (unsigned int j = 0; j < Ny; ++j)
        {
            const unsigned int value = static_cast<unsigned int>(iGlobal * shape[1] + j);
            temperatures[i * Ny + j] = value;
        }
    }

    try
    {
        /** ADIOS class factory of IO class objects, Debug is ON by default */
        adios2::ADIOS adios(MPI_COMM_WORLD);

        // ************************** WRITE
        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO putHeatMap = adios.DeclareIO("HeatMapWriter");

        adios2::Variable<unsigned int> outTemperature = putHeatMap.DefineVariable<unsigned int>(
            "temperature", shape, start, count, adios2::ConstantDims);

        /** Will create HeatMap.bp */
        adios2::Engine bpWriter = putHeatMap.Open("HeatMap2D.bp", adios2::Mode::Write);

        bpWriter.BeginStep();
        bpWriter.Put(outTemperature, temperatures.data());
        bpWriter.EndStep();
        bpWriter.Close();

        // ************************** READ
        if (rank == 0)
        {
            adios2::IO getHeatMap = adios.DeclareIO("HeatMapReader");
            adios2::Engine bpReader =
                getHeatMap.Open("HeatMap2D.bp", adios2::Mode::Read, MPI_COMM_SELF);

            bpReader.BeginStep();
            // this just discovers in the metadata file that the variable exists
            adios2::Variable<unsigned int> inTemperature =
                getHeatMap.InquireVariable<unsigned int>("temperature");
            // now read the variable
            if (inTemperature)
            {
                inTemperature.SetSelection({{2, 2}, {6, 1}});
                size_t elementsSize = inTemperature.SelectionSize();
                std::vector<unsigned int> inTemperatures(elementsSize);
                std::cout << "Pre-allocated " << elementsSize << " elements, "
                          << elementsSize * sizeof(unsigned int) << " bytes\n";

                bpReader.Get(inTemperature, inTemperatures.data(), adios2::Mode::Sync);

                std::cout << "Incoming temperature map:\n";

                for (size_t i = 0; i < inTemperatures.size(); ++i)
                {
                    std::cout << inTemperatures[i] << " ";
                    if ((i + 1) % inTemperature.Count().back() == 0)
                    {
                        std::cout << "\n";
                    }
                }
                std::cout << "\n";
            }
            bpReader.EndStep();

            bpReader.Close();
        }
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

    MPI_Finalize();

    return 0;
}
