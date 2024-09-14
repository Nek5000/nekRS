/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * dataManWriter.cpp
 *
 *  Created on: Feb 16, 2017
 *      Author: Jason Wang
 */

#include <adios2.h>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <thread>
#include <vector>

size_t Nx = 10;
size_t Ny = 10;
size_t steps = 10000;
adios2::Dims shape;
adios2::Dims start;
adios2::Dims count;

int mpiRank, mpiSize;

template <class T>
void PrintData(std::vector<T> &data, const size_t step)
{
    std::cout << "Rank: " << mpiRank << " Step: " << step << " [";
    for (const auto i : data)
    {
        std::cout << i << " ";
    }
    std::cout << "]" << std::endl;
}

template <class T>
std::vector<T> GenerateData(const size_t step)
{
    size_t datasize = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    std::vector<T> myVec(datasize);
    for (size_t i = 0; i < datasize; ++i)
    {
        myVec[i] = i + mpiRank * 10000 + step;
    }
    return myVec;
}

int main(int argc, char *argv[])
{
    // initialize MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    // initialize data dimensions
    count = {Nx, Ny};
    start = {mpiRank * Nx, 0};
    shape = {mpiSize * Nx, Ny};

    // initialize adios2
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io = adios.DeclareIO("whatever");
    io.SetEngine("DataMan");
    io.SetParameters({{"IPAddress", "127.0.0.1"},
                      {"Port", "12306"},
                      {"Timeout", "5"},
                      {"RendezvousReaderCount", "1"}});

    // open stream
    adios2::Engine engine = io.Open("HelloDataMan", adios2::Mode::Write);

    // define variable
    auto floatArrayVar = io.DefineVariable<float>("FloatArray", shape, start, count);

    // write data
    for (size_t i = 0; i < steps; ++i)
    {
        auto floatVector = GenerateData<float>(i);
        engine.BeginStep();
        engine.Put(floatArrayVar, floatVector.data(), adios2::Mode::Sync);
        PrintData(floatVector, engine.CurrentStep());
        engine.EndStep();
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    }

    engine.Close();
    MPI_Finalize();

    return 0;
}
