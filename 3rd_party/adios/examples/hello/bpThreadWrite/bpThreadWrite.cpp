
/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpThreadWrite.cpp : adios2 low-level API example to write in a threaded
 *                    application using C++11 thread and having adios2 calls
 * inside mutex regions adios2 API are not thread-safe:
 *                    1. launching MPI from a thread is not possible on many
 * supercomputers
 *                    2. I/O is highly serialized (buffering and low-level I/O
 * calls), therefore users must be aware that adios2 might introduce
 * bottlenecks. To run: Do not use MPI, just run the executable
 * ./adios2_hello_bpThreadWrite
 *
 *  Created on: Nov 14, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <adios2.h>

#include <cstddef> //std::size_t
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace
{

std::mutex mutex;

// tasks that runs on thread, each section of the vector is covered
template <class T>
void ThreadTask(const std::size_t threadID, std::vector<T> &data, const std::size_t startIndex,
                const std::size_t localSize, const std::string &variableName, adios2::IO io,
                adios2::Engine engine)
{
    (void)threadID; // unused variable
    // populate vector data, but simply adding step to index
    for (std::size_t i = 0; i < localSize; ++i)
    {
        const std::size_t index = startIndex + i;
        data[index] = static_cast<T>(index);
    }

    // I/O write region in a locked mutex
    {
        mutex.lock();

        adios2::Variable<T> variable = io.InquireVariable<T>(variableName);
        variable.SetSelection({{startIndex}, {localSize}});

        engine.Put(variable, &data[startIndex]);
        // PerformPuts must be called to collect memory per buffer
        engine.PerformPuts();

        mutex.unlock();
    }
}

} // end namespace

int main(int argc, char *argv[])
{
    try
    {
        constexpr std::size_t nx = 100;
        // data to be populated and written per thread
        std::vector<double> data(nx);

        // initialize adios2 objects serially
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("thread-write");
        // populate shape, leave start and count empty as
        // they will come from each thread SetSelection
        const std::string variableName = "data";
        io.DefineVariable<double>(variableName, adios2::Dims{nx}, adios2::Dims(), adios2::Dims());

        adios2::Engine engine = io.Open("thread-writes.bp", adios2::Mode::Write);

        // set up thread tasks
        // just grab maximum number of threads to simplify things
        const auto nthreads = static_cast<std::size_t>(std::thread::hardware_concurrency());
        std::vector<std::thread> threadTasks;
        threadTasks.reserve(nthreads);

        // launch threaded tasks (this is what OpenMP would simplify)
        // elements per thread
        const std::size_t stride = nx / nthreads;
        // elements for last thread, add remainder
        const std::size_t last = stride + nx % nthreads;

        engine.BeginStep();
        // launch threads
        for (std::size_t t = 0; t < nthreads; ++t)
        {
            const std::size_t startIndex = stride * t;
            // non-inclusive endIndex
            const std::size_t localSize = (t == nthreads - 1) ? last : stride;

            // use std::ref to pass things by reference
            // adios2 objects can be passed by value

            threadTasks.emplace_back(ThreadTask<double>, t, std::ref(data), startIndex, localSize,
                                     std::ref(variableName), io, engine);
        }

        for (auto &threadTask : threadTasks)
        {
            threadTask.join();
        }
        engine.EndStep();

        engine.Close();
    }
    catch (std::exception &e)
    {
        std::cout << "ERROR: ADIOS2 exception: " << e.what() << "\n";
    }

    return 0;
}
