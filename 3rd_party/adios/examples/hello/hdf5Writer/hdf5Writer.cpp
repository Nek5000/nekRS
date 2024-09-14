/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * hdf5Writer.cpp: Simple self-descriptive example of how to write a
 * variable to a parallel HDF5 File using MPI processes.
 *
 *  Created on: March 20, 2017
 *      Author: Junmin
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
    std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<int> myInts = {0, -1, -2, -3, -4, -5, -6, -7, -8, -9};
    double myScalar = 1.234;
    const std::size_t Nx = myFloats.size();

    try
    {
        /** ADIOS class factory of IO class objects */
        adios2::ADIOS adios(MPI_COMM_WORLD);

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO hdf5IO = adios.DeclareIO("HDFFileIO");
        hdf5IO.SetEngine("HDF5");
        hdf5IO.SetParameter("IdleH5Writer",
                            "true"); // set this if not all ranks are writting

        /** global array : name, { shape (total) }, { start (local) }, { count
         * (local) }, all are constant dimensions */
        adios2::Variable<float> h5Floats = hdf5IO.DefineVariable<float>(
            "h5Floats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

        adios2::Variable<int> h5Ints = hdf5IO.DefineVariable<int>(
            "h5Ints", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

        adios2::Variable<double> h5ScalarDouble = hdf5IO.DefineVariable<double>("h5ScalarDouble");
        /** Engine derived class, spawned to start IO operations */
        adios2::Engine hdf5Writer = hdf5IO.Open("myVector.h5", adios2::Mode::Write);
#ifdef ALL_RANKS_WRITE
        // all Ranks must call Put
        /** Write variable for buffering */
        hdf5Writer.Put<float>(h5Floats, myFloats.data());
        hdf5Writer.Put(h5Ints, myInts.data());
        hdf5Writer.Put(h5ScalarDouble, &myScalar);
#else
        // using collective Begin/EndStep() to run the
        // collective HDF5 calls. Now Ranks can skip writting if no data
        // presented
        hdf5Writer.BeginStep();
        if (rank == 0)
        {
            hdf5Writer.Put<float>(h5Floats, myFloats.data());
            hdf5Writer.Put(h5Ints, myInts.data());
            hdf5Writer.Put(h5ScalarDouble, &myScalar);
        }
        hdf5Writer.EndStep();
#endif
        std::vector<int64_t> m_globalDims = {10, 20, 30, 40};
        hdf5IO.DefineAttribute<std::string>("adios2_schema/version_major",
                                            std::to_string(ADIOS2_VERSION_MAJOR));
        hdf5IO.DefineAttribute<std::string>("adios2_schema/version_minor",
                                            std::to_string(ADIOS2_VERSION_MINOR));
        hdf5IO.DefineAttribute<std::string>("/adios2_schema/mesh/type", "explicit");
        hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension0", m_globalDims[0]);
        hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension1", m_globalDims[1]);
        hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension2", m_globalDims[2]);
        hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension3", m_globalDims[3]);
        hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension-num",
                                             m_globalDims.size());

        hdf5Writer.Close();
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
