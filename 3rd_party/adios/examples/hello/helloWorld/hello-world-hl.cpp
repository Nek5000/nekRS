/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * hello-world-hl.cpp : adios2 high-level API example to write and read a
 *                      std::string Variable with a greeting
 *
 *  Created on: Nov 14, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <iostream>
#include <stdexcept>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

void writer(const std::string &greeting)
{
#if ADIOS2_USE_MPI
    adios2::fstream out("hello-world-hl-cpp.bp", adios2::fstream::out, MPI_COMM_WORLD);
#else
    adios2::fstream out("hello-world-hl-cpp.bp", adios2::fstream::out);
#endif

    out.write("Greeting", greeting);
    out.close();
}

std::string reader()
{
#if ADIOS2_USE_MPI
    adios2::fstream in("hello-world-hl-cpp.bp", adios2::fstream::in, MPI_COMM_WORLD);
#else
    adios2::fstream in("hello-world-hl-cpp.bp", adios2::fstream::in);
#endif

    for (adios2::fstep iStep; adios2::getstep(in, iStep);)
    {
        const std::vector<std::string> greetings = in.read<std::string>("Greeting");
        return greetings.front();
    }
    return "";
}

int main(int argc, char *argv[])
{
#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    try
    {
        const std::string greeting = "Hello World from ADIOS2";
        writer(greeting);

        const std::string message = reader();
        std::cout << message << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "ERROR: ADIOS2 exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
