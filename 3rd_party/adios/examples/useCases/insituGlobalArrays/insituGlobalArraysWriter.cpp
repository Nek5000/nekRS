/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * A Use Case for In Situ visulization frameworks (Conduit, SENSEI)
 *
 * Each processor contributes to a subset of all variables
 * Each variable is written by a subset of all processes
 * The per-writer-blocks in a variable have different sizes
 * The variable cannot be nicely defined as a global array
 *
 * We still define the variables in this writer as global arrays but
 * with k+1 dimensions for each variable, where the extra dimension is
 * used as an index to the writer blocks. The other dimensions are
 * defined as huge numbers to cover all possible dimension sizes on each
 * process.
 *
 * The reader needs to discover the content of each variable block-by-block.
 * It cannot rely on the global shape of the variable as most of it is empty.
 *
 * The global dimensions of a global array MUST NOT change over time.
 * The decomposition of the array across the processes, however, can change
 * between output steps.
 *
 * Created on: Jul 11, 2017
 *      Author: pnorbert
 */

#include <algorithm> // std::transform
#include <iostream>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
MPI_Comm writerComm;
#endif

const size_t NSTEPS = 5;
const size_t BIGDIM = 1000;

/* Variables:
   a: size 8,  written by rank 0 (5 elements) and rank 1 (3 elements)
   b: size 9,  written by rank 0 (5 elements) and rank 2 (4 elements)
   c: size 10, written by rank 0 (5 elements) and rank 3 (5 elements)
   d: size 7,  written by rank 1 (3 elements) and rank 2 (4 elements)

   Variables in the output:
   2D arrays, 4 x BIGDIM
*/

// Which process writes which variables
std::vector<std::vector<std::string>> VarTree = {{"a", "b", "c"}, {"a", "d"}, {"b", "d"}, {"c"}};

// What size of data do they write
std::vector<std::vector<size_t>> SizesTree = {{5, 5, 5}, {3, 3}, {4, 4}, {5, 5}};

std::string argEngine = "BPFile";
adios2::Params engineParams;
std::map<std::string, adios2::Params> engineTransports;

void ProcessArgs(int rank, int argc, char *argv[])
{
    if (argc > 1)
    {
        argEngine = argv[1];
    }
    std::string elc = argEngine;
    std::transform(elc.begin(), elc.end(), elc.begin(), ::tolower);
    if (elc == "sst")
    {
        engineParams["MarshalMethod"] = "BP";
    }
    else if (elc == "insitumpi")
    {
        engineParams["verbose"] = "1";
    }
    else if (elc == "dataman")
    {
        engineParams["WorkflowMode"] = "p2p";
        engineTransports["WAN"] = {
            {"Library", "ZMQ"}, {"IPAddress", "127.0.0.1"}, {"Port", "25600"}};
    }
}

int main(int argc, char *argv[])
{
    int rank = 0, nproc = 1;
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int wrank, wnproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wnproc);
    const unsigned int color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &writerComm);
    MPI_Comm_rank(writerComm, &rank);
    MPI_Comm_size(writerComm, &nproc);
#endif

    const size_t maxProc = VarTree.size();
    if (static_cast<size_t>(nproc) > maxProc)
    {
        if (!rank)
        {
            std::cout << "ERROR: Maximum number of processors for this example is " << maxProc
                      << std::endl;
        }
        exit(1);
    }

    ProcessArgs(rank, argc, argv);
    if (!rank)
    {
        std::cout << "Writer: ADIOS2 Engine set to: " << argEngine << "   Parameters:";
        for (auto &p : engineParams)
        {
            std::cout << "    " << p.first << " = " << p.second;
        }
        std::cout << "  Transports: ";
        for (auto &t : engineTransports)
        {
            std::cout << "  " << t.first << " : {";
            for (auto &p : t.second)
            {
                std::cout << " {" << p.first << ", " << p.second << "}";
            }
            std::cout << " }";
        }
        std::cout << std::endl;
    }

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(writerComm);
#else
    adios2::ADIOS adios;
#endif

    // We have a varying number of vars on each processor
    const size_t nvars = VarTree[rank].size();
    // A 1D array for each variable
    std::vector<std::vector<double>> Vars(nvars);
    for (size_t i = 0; i < nvars; i++)
    {
        Vars[i].resize(SizesTree[rank][i]);
    }

    std::vector<adios2::Variable<double>> ADIOSVars(nvars);

    try
    {
        adios2::IO io = adios.DeclareIO("Output");
        io.SetEngine(argEngine);
        io.SetParameters(engineParams);
        for (auto &t : engineTransports)
        {
            io.AddTransport(t.first, t.second);
        }

        for (size_t i = 0; i < nvars; i++)
        {
            size_t nelems = SizesTree[rank][i];
            Vars[i].resize(nelems);
            ADIOSVars[i] =
                io.DefineVariable<double>(VarTree[rank][i], {(unsigned int)nproc, BIGDIM});
        }

        adios2::Engine writer = io.Open("output.bp", adios2::Mode::Write);

        for (size_t step = 0; step < NSTEPS; step++)
        {
            writer.BeginStep();

            for (size_t i = 0; i < nvars; i++)
            {
                size_t nelems = SizesTree[rank][i];
                for (size_t j = 0; j < nelems; j++)
                {
                    Vars[i][j] = ((double)step + 1.0) / 100.0 + (double)rank;
                }

                // Make a 2D selection to describe the local dimensions of the
                // variable we write and its offsets in the global spaces
                // adios2::SelectionBoundingBox sel();
                ADIOSVars[i].SetSelection(adios2::Box<adios2::Dims>(
                    {static_cast<size_t>(rank), 0}, {1, static_cast<size_t>(nelems)}));
                writer.Put<double>(ADIOSVars[i], Vars[i].data());
            }

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
