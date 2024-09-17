/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CppReader.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <adios2.h>
#include <mpi.h>

#include <iostream>

int main(int argc, char *argv[])
{
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
    {
        /** ADIOS class factory of IO class objects */
        adios2::ADIOS adios(MPI_COMM_SELF);

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO bpIO = adios.DeclareIO("CppReader");

        /** Engine derived class, spawned to start IO operations */
        adios2::Engine bpReader = bpIO.Open("FWriter.bp", adios2::Mode::Read);

        bpReader.BeginStep();

        const std::map<std::string, adios2::Params> variables = bpIO.AvailableVariables();

        for (const auto &variablePair : variables)
        {
            std::cout << "Name: " << variablePair.first;

            for (const auto &parameter : variablePair.second)
            {
                std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
            }
        }

        auto bpData = bpIO.InquireVariable<float>("data2D");
        bpData.SetSelection({{1, 1}, {3, 2}});
        std::vector<float> data(bpData.SelectionSize());

        bpReader.Get(bpData, data.data());

        bpReader.EndStep();

        bpReader.Close();

        std::cout << "Selection size: " << bpData.SelectionSize() << "\n";
        for (size_t i = 0; i < bpData.Count()[0]; ++i)
        {
            for (size_t j = 0; j < bpData.Count()[1]; ++j)
            {
                std::cout << data[i * bpData.Count()[1] + j] << " ";
            }

            std::cout << "\n";
        }
        std::cout << "\n";
    }

    MPI_Finalize();

    return 0;
}
