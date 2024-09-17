/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Write local arrays from multiple processors.
 *
 * If one cannot or does not want to organize arrays present on each process
 * as one global array, still one can write them out with the same name.
 * Reading, however, needs to be handled differently: each process' array has
 * to be read separately, using Block selections. The size of each process
 * block should be discovered by the reading application by inquiring per-block
 * size information of the variable, and allocate memory for reading
 * accordingly.
 *
 * In this example we write v0, v1, v2 and v3, in 5 output steps, where
 * v0 has the same size on every process at every step
 * v1 has different size on each process but fixed over time
 * v2 has different size on each process and that is changing over time
 * v3 is like v2 but also the number of processes writing it changes over time
 *
 * bpls can show the size of each block of the variable:
 *
 * $ cd <adios build directory>
 * $ make
 * $ mpirun -n 4 ./bin/localArray
 * $ bpls -l localArray.bp
 * double   v0    5*[4]*{6} = 0 / 3.4
 * double   v1    5*[4]*{__} = 0 / 3.4
 * double   v2    5*[4]*{__} = 0 / 3.4
 * double   v3    5*[__]*{__} = 0 / 3.3
 *
 * Study the decomposition of each variable
 * $ bpls -l localArray.bp -D v0
 * and notice the progression in the changes.
 *
 * Created on: Jun 2, 2017
 *      Author: Norbert Podhorszki <pnorbert@ornl.gov
 */

#include <iostream>
#include <vector>

#include <adios2.h>

int main(int argc, char *argv[])
{
    int rank = 0;
#if ADIOS2_USE_MPI
    int nproc = 1;
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
    const int NSTEPS = 5;

    // generate different random numbers on each process,
    // but always the same sequence at each run
    srand(rank * 32767);

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    // v0 has the same size on every process at every step
    const size_t Nglobal = 6;
    std::vector<double> v0(Nglobal);

    // Application variables for output
    // random size per process, 5..10 each
    // v1 has different size on each process (but fixed over time)
    const unsigned int Nx = rand() % 6 + 5;
    // Local array, size is fixed over time on each process
    std::vector<double> v1(Nx);

    // random size per process, a different size at each step
    unsigned int Nelems;
    // Local array, size is changing over time on each process
    std::vector<double> v2;

    // Local array, size is changing over time on each process
    // Also, random number of processes will write it at each step
    std::vector<double> &v3 = v2;

    try
    {
        // Get io settings from the config file or
        // create one with default settings here
        adios2::IO io = adios.DeclareIO("Output");
        io.SetEngine("BPFile");
        io.SetParameters({{"verbose", "4"}});

        /*
         * Define local array: type, name, local size
         * Global dimension and starting offset must be an empty vector
         * Here the size of the local array is the same on every process
         */
        adios2::Variable<double> varV0 = io.DefineVariable<double>("v0", {}, {}, {Nglobal});

        /*
         * v1 is similar to v0 but on every process the local size
         * is a different value
         */
        adios2::Variable<double> varV1 = io.DefineVariable<double>("v1", {}, {}, {Nx});

        /*
         * Define local array: type, name
         * Global dimension and starting offset must be an empty vector
         * but local size CANNOT be an empty vector.
         * We can use {adios2::UnknownDim} for this purpose or any number
         * actually since we will modify it before writing
         */
        adios2::Variable<double> varV2 =
            io.DefineVariable<double>("v2", {}, {}, {adios2::UnknownDim});

        /*
         * v3 is just like v2
         */
        adios2::Variable<double> varV3 =
            io.DefineVariable<double>("v3", {}, {}, {adios2::UnknownDim});

        // Open file. "w" means we overwrite any existing file on disk,
        // but Advance() will append steps to the same file.
        adios2::Engine writer = io.Open("localArray.bp", adios2::Mode::Write);

        for (int step = 0; step < NSTEPS; step++)
        {
            writer.BeginStep();

            // v0
            for (size_t i = 0; i < Nglobal; i++)
            {
                v0[i] = rank * 1.0 + step * 0.1;
            }
            writer.Put<double>(varV0, v0.data());

            // v1
            for (size_t i = 0; i < Nx; i++)
            {
                v1[i] = rank * 1.0 + step * 0.1;
            }
            writer.Put<double>(varV1, v1.data());

            // v2

            // random size per process per step, 5..10 each
            Nelems = rand() % 6 + 5;
            v2.reserve(Nelems);
            for (size_t i = 0; i < Nelems; i++)
            {
                v2[i] = rank * 1.0 + step * 0.1;
            }

            // Set the size of the array now because we did not know
            // the size at the time of definition
            varV2.SetSelection(adios2::Box<adios2::Dims>({}, {Nelems}));
            writer.Put<double>(varV2, v2.data());

            // v3

            // random chance who writes it
            unsigned int chance = rand() % 100;
            /*if (step == 2)
            {
                chance = 0;
            }*/
            bool doWrite = (chance > 60);
            if (doWrite)
            {
                varV3.SetSelection(adios2::Box<adios2::Dims>({}, {Nelems}));
                writer.Put<double>(varV3, v3.data());
            }

            writer.EndStep();
        }

        writer.Close();
    }
    catch (std::invalid_argument &e)
    {
        if (rank == 0)
        {
            std::cout << "Invalid argument exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::ios_base::failure &e)
    {
        if (rank == 0)
        {
            std::cout << "System exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cout << "Exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
