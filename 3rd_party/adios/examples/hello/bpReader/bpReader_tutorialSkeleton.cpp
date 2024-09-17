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
#include <mpi.h>
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>

int main(int argc, char *argv[])
{
    int rank, size;
    int provided;

    // Add code to init MPI
    try
    {
        // Add code to create ADIOS object

        // Add code to create IO object

        // Add code to open file

        // Add code to inquire variables and optionally check all available variables

        // Add code to read variables

        // Add code to close file
    }
    catch (std::invalid_argument &e)
    {
        if (rank == 0)
        {
            std::cerr << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
            std::cerr << e.what() << "\n";
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
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
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cerr << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
            std::cerr << e.what() << "\n";
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();

    return 0;
}
