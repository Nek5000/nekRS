/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 *  Created on: Sept 2018
 *      Author: Norbert Podhorszki
 */

#include <chrono>
#include <cstdlib>
#include <errno.h>
#include <fstream>
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <thread>
#include <vector>

#include <adios2.h>

const std::string streamname = "skeleton_stream";
const unsigned int ndx = 10;  // Local array size  per process
const unsigned int steps = 2; // number of output steps

int main(int argc, char *argv[])
{
    int rank = 0;
    int nproc = 1;
    int retval = 0;

    std::ofstream out("TestSkeletonWriterOutput.txt");
    auto coutbuf = std::cout.rdbuf(out.rdbuf()); // save and redirect

#if ADIOS2_USE_MPI
    int wrank = 0, wnproc = 1;
    MPI_Comm mpiWriterComm;
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wnproc);

    const unsigned int color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &mpiWriterComm);

    MPI_Comm_rank(mpiWriterComm, &rank);
    MPI_Comm_size(mpiWriterComm, &nproc);
#endif

    try
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(mpiWriterComm);
#else
        adios2::ADIOS adios;
#endif

        std::vector<float> myArray(ndx);
        adios2::IO io = adios.DeclareIO("writer");
        io.SetEngine("Skeleton");
        io.SetParameter("VerBose", "5");

        adios2::Variable<float> varArray =
            io.DefineVariable<float>("myArray", {(unsigned int)nproc * ndx},
                                     {(unsigned int)rank * ndx}, {ndx}, adios2::ConstantDims);

        adios2::Variable<std::string> varSyncString =
            io.DefineVariable<std::string>("mySyncString");

        adios2::Engine writer = io.Open(streamname, adios2::Mode::Write);

        for (size_t step = 0; step < steps; ++step)
        {
            int idx = 0;
            for (size_t i = 0; i < ndx; ++i)
            {
                myArray[idx] = rank + (step / 100.0f);
                ++idx;
            }
            writer.BeginStep(adios2::StepMode::Append);
            writer.Put<std::string>(varSyncString, "Hello", adios2::Mode::Sync);
            writer.Put<float>(varArray, myArray.data());
            writer.EndStep();
        }
        writer.Close();
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
        retval = 1;
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM "
                     "from rank "
                  << rank << "\n";
        std::cout << e.what() << "\n";
        retval = 2;
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
        retval = 3;
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    std::cout.rdbuf(coutbuf); // reset to standard output again

    return retval;
}
