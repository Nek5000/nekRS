/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpSZ.cpp : example passing runtime compression arguments
 *
 *  Created on: Aug 3, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <algorithm> //std::transform
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <numeric>   //std::iota
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>
#include <mpi.h>

void Usage()
{
    std::cout << "\n";
    std::cout << "USAGE:\n";
    std::cout << "./adios2_hello_bpSZ Nx sz_accuracy\n";
    std::cout << "\t Nx: size of float and double arrays to be compressed\n";
    std::cout << "\t sz_accuracy: absolute accuracy e.g. 0.1, 0.001, to skip "
                 "compression: -1\n\n";
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        Usage();
        return EXIT_SUCCESS;
    }

    int rank, size;
    int provided;

    // Add code to init MPI

    try
    {
        // Add code to get command line arguments

        // Add code to create arrays

        // Add code to create ADIOS object

        // Add code to create IO object

        // Add code to create variables

        // Add code to add SZ compressor operation

        // Add code to add attribute

        // Add code to open file

        // Add code to write variables for 3 time steps and edit them

        // Add code to close file
    }
    catch (std::invalid_argument &e)
    {
        std::cerr << "Invalid argument exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::ios_base::failure &e)
    {
        std::cerr << "IO System base failure exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }

    // Add code to finalize MPI

    return 0;
}
