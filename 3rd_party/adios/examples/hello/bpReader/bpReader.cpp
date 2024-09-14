/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpReader.cpp: Simple self-descriptive example of how to read a variable
 * from a BP File.
 *
 *  Created on: Feb 16, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>

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
    std::cout << "rank " << rank << " size " << size << "\n";
    try
    {
#if ADIOS2_USE_MPI
        /** ADIOS class factory of IO class objects */
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");

        /** Engine derived class, spawned to start IO operations */
        adios2::Engine bpReader = bpIO.Open("myVector_cpp.bp", adios2::Mode::Read);

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

        /** Write variable for buffering */
        adios2::Variable<float> bpFloats = bpIO.InquireVariable<float>("bpFloats");
        adios2::Variable<int> bpInts = bpIO.InquireVariable<int>("bpInts");

        const std::size_t Nx = 10;
        if (bpFloats) // means found
        {
            std::vector<float> myFloats;

            // read only the chunk corresponding to our rank
            bpFloats.SetSelection({{Nx * rank}, {Nx}});
            // myFloats.data is pre-allocated
            bpReader.Get(bpFloats, myFloats, adios2::Mode::Sync);

            if (rank == 0)
            {
                std::cout << "MyFloats: \n";
                for (const auto number : myFloats)
                {
                    std::cout << number << " ";
                }
                std::cout << "\n";
            }
        }

        if (bpInts) // means not found
        {
            std::vector<int> myInts;
            // read only the chunk corresponding to our rank
            bpInts.SetSelection({{Nx * rank}, {Nx}});

            bpReader.Get(bpInts, myInts, adios2::Mode::Sync);

            if (rank == 0)
            {
                std::cout << "myInts: \n";
                for (const auto number : myInts)
                {
                    std::cout << number << " ";
                }
                std::cout << "\n";
            }
        }
        bpReader.EndStep();

        /** Close bp file, engine becomes unreachable after this*/
        bpReader.Close();
    }
    catch (std::invalid_argument &e)
    {
        if (rank == 0)
        {
            std::cerr << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
            std::cerr << e.what() << "\n";
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::ios_base::failure &e)
    {
        if (rank == 0)
        {
            std::cerr << "IO System base failure exception, STOPPING PROGRAM "
                         "from rank "
                      << rank << "\n";
            std::cerr << e.what() << "\n";
            std::cerr << "The file myVector_cpp.bp does not exist."
                      << " Presumably this is because adios2_hello_bpWriter has not "
                         "been run."
                      << " Run ./adios2_hello_bpWriter before running this program.\n";
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cerr << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
            std::cerr << e.what() << "\n";
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
