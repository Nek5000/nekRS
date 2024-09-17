/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * dataSpacesWriter.cpp
 *
 *  Created on: Feb 06, 2019
 *      Author: Pradeep Subedi
 *      		pradeep.subedi@rutgers.edu
 */

#include <iostream>
#include <vector>

#include <adios2.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{

    int rank;
    int size;

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

    std::vector<float> myFloats = {(float)10.0 * rank + 0, (float)10.0 * rank + 1,
                                   (float)10.0 * rank + 2, (float)10.0 * rank + 3,
                                   (float)10.0 * rank + 4, (float)10.0 * rank + 5,
                                   (float)10.0 * rank + 6, (float)10.0 * rank + 7,
                                   (float)10.0 * rank + 8, (float)10.0 * rank + 9};
    const std::size_t Nx = myFloats.size();

    try
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO dataSpacesIO = adios.DeclareIO("myIO");
        dataSpacesIO.SetEngine("DATASPACES");

        // Define variable and local size
        auto bpFloats =
            dataSpacesIO.DefineVariable<float>("bpFloats", {size * Nx}, {rank * Nx}, {Nx});

        // Create engine smart pointer to Sst Engine due to polymorphism,
        // Open returns a smart pointer to Engine containing the Derived class
        adios2::Engine dataSpacesWriter = dataSpacesIO.Open("helloDataSpaces", adios2::Mode::Write);

        dataSpacesWriter.BeginStep();
        dataSpacesWriter.Put<float>(bpFloats, myFloats.data());
        dataSpacesWriter.EndStep();
        dataSpacesWriter.Close();
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM from rank " << rank
                  << "\n";
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
