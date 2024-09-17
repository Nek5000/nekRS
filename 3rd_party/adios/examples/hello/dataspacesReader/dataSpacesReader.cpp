/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * dataSpacesReader.cpp
 *
 *  Created on: Feb 06, 2019
 *      Author: Pradeep Subedi
 *      		pradeep.subedi@rutgers.edu
 */

#include <chrono>
#include <iostream>
#include <numeric>
#include <thread>
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

    std::vector<float> myFloats(10);

    try
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        adios2::IO dataSpacesIO = adios.DeclareIO("myIO");
        dataSpacesIO.SetEngine("DATASPACES");

        adios2::Engine dataSpacesReader = dataSpacesIO.Open("helloDataSpaces", adios2::Mode::Read);
        dataSpacesReader.BeginStep();
        adios2::Variable<float> bpFloats = dataSpacesIO.InquireVariable<float>("bpFloats");
        std::cout << "Incoming variable is of size " << bpFloats.Shape()[0] << "\n";
        const std::size_t total_size = bpFloats.Shape()[0];
        const std::size_t my_start = (total_size / size) * rank;
        const std::size_t my_count = (total_size / size);
        std::cout << "Reader rank " << rank << " reading " << my_count
                  << " floats starting at element " << my_start << "\n";

        const adios2::Dims start{my_start};
        const adios2::Dims count{my_count};

        const adios2::Box<adios2::Dims> sel(start, count);

        std::vector<float> myFloats;
        myFloats.resize(my_count);

        bpFloats.SetSelection(sel);
        dataSpacesReader.Get(bpFloats, myFloats.data());
        dataSpacesReader.EndStep();

        dataSpacesReader.Close();
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
