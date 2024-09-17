/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpAttributeWriteRead.cpp: Simple self-descriptive example of how to write/read attributes and
 * a variable to a BP File that lives in several MPI processes.
 *
 *  Created on: Feb 16, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <mpi.h>
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <vector>

#include <adios2.h>

void writer(adios2::ADIOS &adios, int rank, int size, std::vector<float> &myFloats)
{
    // Add code to create IO object

    // Add code to create variable

    // Add code to create attributes

    // Add code to open file

    // Add code to write variables

    // Add code to close file
}

void reader(adios2::ADIOS &adios, int rank, int /*size*/)
{
    // Add code to create IO object

    // Add code to open file

    // add code to check available attributes

    // Add code to read attributes

    // Add code to read variables

    // Add code to close file
}

int main(int argc, char *argv[])
{
    int rank, size;
    int provided;

    // Add code to init MPI

    // Add code to create array
    try
    {
        // Add code to create ADIOS object

        // Call writer and reader functions
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM from rank " << rank
                  << "\n";
        std::cout << e.what() << "\n";
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Add code to finalize MPI

    return 0;
}
