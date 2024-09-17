/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * sstWriter.cpp
 *
 *  Created on: Aug 17, 2017
 *      Author: Greg Eisenhauer
 */

#include <iostream>
#include <vector>

#include <adios2.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

template <class T>
void PrintData(const std::vector<T> &data, const int rank, const size_t step)
{
    std::cout << "Rank: " << rank << " Step: " << step << " [";
    for (size_t i = 0; i < data.size(); ++i)
    {
        std::cout << data[i] << " ";
    }
    std::cout << "]" << std::endl;
}

int main(int argc, char *argv[])
{

    int rank;
    int size;

#if ADIOS2_USE_MPI
    int provide;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provide);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    rank = 0;
    size = 1;
#endif

    std::vector<float> myFloats = {
        static_cast<float>(10.0 * rank + 0), static_cast<float>(10.0 * rank + 1),
        static_cast<float>(10.0 * rank + 2), static_cast<float>(10.0 * rank + 3),
        static_cast<float>(10.0 * rank + 4), static_cast<float>(10.0 * rank + 5),
        static_cast<float>(10.0 * rank + 6), static_cast<float>(10.0 * rank + 7),
        static_cast<float>(10.0 * rank + 8), static_cast<float>(10.0 * rank + 9)};
    const std::size_t Nx = myFloats.size();
    const float increment = Nx * size * 1.0;

    try
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO sstIO = adios.DeclareIO("myIO");
        sstIO.SetEngine("Sst");

        // Define variable and local size
        auto bpFloats = sstIO.DefineVariable<float>("bpFloats", {size * Nx}, {rank * Nx}, {Nx});

        // Create engine smart pointer to Sst Engine due to polymorphism,
        // Open returns a smart pointer to Engine containing the Derived class
        adios2::Engine sstWriter = sstIO.Open("helloSst", adios2::Mode::Write);

        for (size_t i = 0; i < 4; ++i)
        {
            PrintData(myFloats, rank, i);
            sstWriter.BeginStep();
            sstWriter.Put<float>(bpFloats, myFloats.data());
            sstWriter.EndStep();
            for (size_t k = 0; k < myFloats.size(); ++k)
            {
                myFloats[k] += increment;
            }
        }
        sstWriter.Close();
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
