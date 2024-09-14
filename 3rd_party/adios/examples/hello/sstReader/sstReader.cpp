/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * sstReader.cpp
 *
 *  Created on: Aug 17, 2017
 *      Author: Greg Eisenhauer
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

    std::vector<float> myFloats(10);

    try
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        adios2::IO sstIO = adios.DeclareIO("myIO");
        sstIO.SetEngine("Sst");

        adios2::Engine sstReader = sstIO.Open("helloSst", adios2::Mode::Read);

        while (true)
        {
            auto status = sstReader.BeginStep();
            if (status == adios2::StepStatus::EndOfStream)
            {
                break;
            }
            else if (status == adios2::StepStatus::OtherError)
            {
                std::cout << "ERROR in stream processing when calling BeginStep(). Quit"
                          << std::endl;
                break;
            }

            adios2::Variable<float> bpFloats = sstIO.InquireVariable<float>("bpFloats");
            // std::cout << "Incoming variable is of size " << bpFloats.Shape()[0] << "\n";
            const std::size_t total_size = bpFloats.Shape()[0];
            const std::size_t my_start = (total_size / size) * rank;
            std::size_t my_count = (total_size / size);
            if (rank == size - 1)
            {
                my_count += (total_size % size);
            }

            // std::cout << "Reader rank " << rank << " reading " << my_count
            //           << " floats starting at element " << my_start << "\n";

            const adios2::Dims start{my_start};
            const adios2::Dims count{my_count};

            const adios2::Box<adios2::Dims> sel(start, count);

            std::vector<float> myFloats;
            myFloats.resize(my_count);

            bpFloats.SetSelection(sel);
            sstReader.Get(bpFloats, myFloats.data());
            sstReader.EndStep();
            // myfloats is filled ONLY AFTER EndStep!!!
            // Use adios2.Mode.Sync extra parameter in Get() to get the data immediately
            PrintData(myFloats, rank, sstReader.CurrentStep());
        }

        sstReader.Close();
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
