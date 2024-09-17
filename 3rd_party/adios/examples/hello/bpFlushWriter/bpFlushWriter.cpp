/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpFlushWriter.cpp: Example that tests buffer overflow forcing a flush to
 * transports when writing a large variable in independent N-to-N mode. This
 * will have performance penalties, but it's safer.
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
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /** Application variable */
    std::vector<float> myFloats(1000000); //~ 4 MB
    const std::size_t Nx = myFloats.size();

    try
    {
        /** ADIOS class factory of IO class objects */
        adios2::ADIOS adios(MPI_COMM_WORLD);

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO bpIO = adios.DeclareIO("BPFile_N2N_Flush");
        bpIO.SetEngine("BPFile");

        //        bpIO.SetParameters({{"MaxBufferSize", "9Mb"},
        //                            {"BufferGrowthFactor", "1.5"},
        //                            {"Threads", "2"}});
        //        bpIO.AddTransport("File", {{"ProfileUnits", "Microseconds"}});

        /** global array : name, { shape (total) }, { start (local) }, { count
         * (local) }, all are constant dimensions */
        adios2::Variable<float> bpFloats = bpIO.DefineVariable<float>(
            "bpFloats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

        /** Engine derived class, spawned to start IO operations */
        adios2::Engine bpWriter = bpIO.Open("myVectorFlush.bp", adios2::Mode::Write);

        bpWriter.BeginStep();
        for (unsigned int t = 0; t < 100; ++t)
        {
            /** values to time step */
            myFloats.assign(myFloats.size(), static_cast<float>(t));
            /** Write variable for buffering */
            bpWriter.Put<float>(bpFloats, myFloats.data());
        }
        bpWriter.EndStep();

        /** Create bp file, engine becomes unreachable after this*/
        bpWriter.Close();
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

    MPI_Finalize();

    return 0;
}
