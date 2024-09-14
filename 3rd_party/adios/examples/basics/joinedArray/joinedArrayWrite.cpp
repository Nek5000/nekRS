/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Write local arrays from multiple processors and make ADIOS join them
 * at reading to show a global array
 *
 * If every process has an array that is different only in one dimension
 * it can be presented as a global array by joining the arrays together.
 * E.g. if every process has a table with a different number of rows,
 * and one does not want to do a global communication to calculate the offsets
 * in the global table, one can just write the local arrays and let ADIOS
 * calculate the offsets at read time (when all sizes are known by any process).
 *
 * bpls can show the size of each block of the table:
 * bpls -D <file> <variable>
 *
 * Note: only one dimension can be joinable, every other dimension must be the
 * same on each process.
 *
 * Note: the local dimension size in the joinable dimension is allowed to change
 * over time within each processor. However, if the sum of all local sizes
 * changes over time, the result will look like a local array.
 * (Because global arrays with changing global dimension over time can only be
 * handled as local arrays in ADIOS)
 *
 *
 * Created on: Jun 2, 2017
 *      Author: pnorbert
 */

#include <iostream>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

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
    const int NSTEPS = 3;

    // generate different random numbers on each process,
    // but always the same sequence at each run
    srand(rank * 32767);

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    // Application variables for output
    // random size per process, 5..10 each
    const unsigned int Nrows = rand() % 6 + 5;
    const unsigned int Ncols = 4;
    // Local array, size is fixed over time on each process
    std::vector<double> mytable(Nrows * Ncols);

    try
    {
        // Get io settings from the config file or
        // create one with default settings here
        adios2::IO io = adios.DeclareIO("Output");
        io.SetEngine("BPFile");

        /*
         * Define joinable local array: type, name, global and local size
         * Starting offset can be an empty vector
         * Only one global dimension can be joined
         */
        adios2::Variable<double> varTable =
            io.DefineVariable<double>("table", {adios2::JoinedDim, Ncols}, {}, {Nrows, Ncols});

        // adios2::Operator op = adios.DefineOperator("blosc", "blosc");
        // varTable.AddOperation(op, {{"clevel", std::to_string(1)}});

        // Open file. "w" means we overwrite any existing file on disk,
        // but Advance() will append steps to the same file.
        adios2::Engine writer = io.Open("joinedArray.bp", adios2::Mode::Write);

        for (int step = 0; step < NSTEPS; step++)
        {
            writer.BeginStep();

            for (unsigned int row = 0; row < Nrows; row++)
            {
                for (unsigned int col = 0; col < Ncols; col++)
                {
                    mytable[row * Ncols + col] =
                        static_cast<double>(step * 10.0 + rank * 1.0 + row * 0.1 + col * 0.01);
                }
            }

            writer.Put<double>(varTable, mytable.data());

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
