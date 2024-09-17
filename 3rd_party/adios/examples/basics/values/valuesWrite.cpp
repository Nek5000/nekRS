/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Write single values to a file. There are four different cases:
 * 1. Global constant - same on all processes, constant over time
 * 2. Global value - same on all processes, may change over time
 * 3. Local constants - different across processes, constant over time
 * 4. Local value - different across processes, may change over time
 *
 * Constants are not handled separately from time-varying values in ADIOS.
 * Simply write them only in the first step.
 *
 * Writing a global value from multiple processes does not hurt but it is
 * useless.
 *
 * Created on: Jun 2, 2017
 *      Author: pnorbert
 */

#include <iostream>
#include <string>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
    int rank = 0, nproc = 1;
#if ADIOS2_USE_MPI
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

    // Application variables for output
    // 1. Global value, constant across processes, constant over time
    // This is 'nproc'

    // 2. Global value, constant across processes, varying value over time
    // This is 'step'

    // 3. Local value, varying across processes, constant over time
    // It will appear in reading as a 1D array of nproc elements.
    // This is 'rank'

    // 4. Local value, varying across processes, varying over time
    unsigned int Nparts; // random size per process, 5..10 each

    try
    {
        // Get io settings from the config file or
        // create one with default settings here
        adios2::IO io = adios.DeclareIO("Output");
        io.SetEngine("BPFile");
        io.SetParameters({{"verbose", "4"}});
        /*
         * Define variables
         */
        // 1. Global constant, same value across processes, constant over time
        adios2::Variable<int> varNproc = io.DefineVariable<int>("Nproc");
        (void)varNproc; // For the sake of the example we create an unused
                        // variable

        // 2. Global value, same value across processes, varying value over time
        adios2::Variable<int> varStep = io.DefineVariable<int>("Step");
        adios2::Variable<std::string> varGlobalString =
            io.DefineVariable<std::string>("GlobalString");

        // 3. Local value, varying across processes, constant over time
        adios2::Variable<int> varProcessID =
            io.DefineVariable<int>("ProcessID", {adios2::LocalValueDim});

        // 4. Local value, varying across processes, varying over time
        adios2::Variable<unsigned int> varNparts =
            io.DefineVariable<unsigned int>("Nparts", {adios2::LocalValueDim});

        // Open file. "w" means we overwrite any existing file on disk,
        // but Advance() will append steps to the same file.
        adios2::Engine writer = io.Open("values.bp", adios2::Mode::Write);

        for (int step = 0; step < NSTEPS; step++)
        {
            writer.BeginStep();

            // random size per process, 5..10 each
            Nparts = rand() % 6 + 5;

            // 1. and 2. Writing a global value from only one process
            if (rank == 0)
            {
                // 1. Writing a global constant value only once
                if (step == 0)
                {
                    adios2::Variable<int> varNproc = io.InquireVariable<int>("Nproc");
                    writer.Put<int>(varNproc, nproc);
                }
                writer.Put<int>(varStep, step);

                std::string str = "This is step " + std::to_string(step);
                // str will go out of scope before EndStep(), so we must use
                // Sync mode in Put()
                writer.Put<std::string>(varGlobalString, str, adios2::Mode::Sync);
            }

            // 3. and 4. Writing a local value on every process. Will be shown
            // at reading as a 1D array
            if (step == 0)
            {
                writer.Put<int>(varProcessID, rank);
            }
            writer.Put<unsigned int>(varNparts, Nparts);

            // Indicate we are done for this step.
            // Disk I/O will be performed during this call unless
            // time aggregation postpones all of that to some later step
            writer.EndStep();
        }

        // Called once: indicate that we are done with this output for the run
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
