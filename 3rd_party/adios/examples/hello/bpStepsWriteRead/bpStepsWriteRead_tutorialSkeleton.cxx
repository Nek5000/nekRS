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

#include <algorithm> //std::for_each
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <mpi.h>
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
    // Add code to create the simulation data

    // Add code to create ADIOS io and set engine type

    // Add code to define sim data variable

    // Add code to define step variable

    // Add code to open file

    // Add code to write data across multiple steps, and update the simulation data

    // Add code to close file
}

void reader(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
            const size_t Nx, unsigned int /*nSteps*/, int rank, int /*size*/)
{
    // Add code to create ADIOS io and set engine type

    // Add code to open file

    // Add code to create variable for sim data and step

    // Add code to read data across multiple steps

    // Add code to close file
}

int main(int argc, char *argv[])
{
    int rank, size;
    int provided;

    // Add code to initialize MPI

    const std::string engine = argv[1] ? argv[1] : "BPFile";
    std::cout << "Using engine " << engine << std::endl;

    // Add Code to set filename, nSteps, Nx
    try
    {
        // Add code to create ADIOS object

        // Add code to call writer
        // Add code to call reader
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

    // Add code to finalize MPI

    return 0;
}
