/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * main.cpp
 *
 * Recreates heat_transfer.f90 (Fortran) ADIOS tutorial example in C++
 * This version is adapted from heatTransfer/write/main.cpp and
 * heatTransfer/read/heatRead.cpp for use with the inline engine,
 * which requires the writer and reader to be created in the same process
 * and the same ADIOS IO object.
 *
 * Created on: May 2020
 *     Author: Norbert Podhorszki
 *             Caitlin Ross
 *
 */
#include <mpi.h>

#include "adios2.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "../write/HeatTransfer.h"
#include "../write/Settings.h"
#include "InlineIO.h"

void printUsage()
{
    std::cout << "Usage: heatTransfer  config   output  N  M   nx  ny   steps "
                 "iterations\n"
              << "  config: XML config file to use\n"
              << "  output: name of output data file/stream\n"
              << "  N:      number of processes in X dimension\n"
              << "  M:      number of processes in Y dimension\n"
              << "  nx:     local array size in X dimension per processor\n"
              << "  ny:     local array size in Y dimension per processor\n"
              << "  steps:  the total number of steps to output\n"
              << "  iterations: one step consist of this many iterations\n\n";
}

void Compute(const double *Tin, std::vector<double> &Tout, std::vector<double> &dT, bool firstStep)
{
    /* Compute dT and
     * copy Tin into Tout as it will be used for calculating dT in the
     * next step
     */
    if (firstStep)
    {
        for (size_t i = 0; i < dT.size(); ++i)
        {
            dT[i] = 0;
            Tout[i] = Tin[i];
        }
    }
    else
    {
        for (size_t i = 0; i < dT.size(); ++i)
        {
            dT[i] = Tout[i] - Tin[i];
            Tout[i] = Tin[i];
        }
    }
}

std::vector<double> Tout;
std::vector<double> dT;
adios2::Variable<double> vTout;
adios2::Variable<double> vdT;
adios2::IO outIO;
adios2::ADIOS ad;
adios2::Engine writer;

void writeOutput()
{
    writer.BeginStep();
    if (vTout)
        writer.Put<double>(vTout, Tout.data());
    if (vdT)
        writer.Put<double>(vdT, dT.data());
    writer.EndStep();
}

void setupOutputIO(const Settings &s, MPI_Comm comm)
{
    ad = adios2::ADIOS(s.configfile, comm);
    outIO = ad.DeclareIO("readerOutput");
    Tout.resize(s.ndx * s.ndy);
    dT.resize(s.ndx * s.ndy);

    /* Create output variables and open output stream */
    // For inline engine, there's no exchange of data between processes,
    // so the shape of variables to be written out for validation
    // is the same as the writer
    vTout = outIO.DefineVariable<double>("T", {s.gndx, s.gndy}, {s.offsx, s.offsy}, {s.ndx, s.ndy});
    vdT = outIO.DefineVariable<double>("dT", {s.gndx, s.gndy}, {s.offsx, s.offsy}, {s.ndx, s.ndy});
    writer = outIO.Open(s.outputfile, adios2::Mode::Write, comm);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    /* When writer and reader is launched together with a single mpirun command,
       the world comm spans all applications. We have to split and create the
       local 'world' communicator mpiHeatTransferComm for the writer only.
       When writer and reader is launched separately, the mpiHeatTransferComm
       communicator will just equal the MPI_COMM_WORLD.
     */

    int wrank, wnproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wnproc);

    const unsigned int color = 1;
    MPI_Comm mpiHeatTransferComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &mpiHeatTransferComm);

    int rank, nproc;
    MPI_Comm_rank(mpiHeatTransferComm, &rank);
    MPI_Comm_size(mpiHeatTransferComm, &nproc);

    try
    {
        double timeStart = MPI_Wtime();
        Settings settings(argc, argv, rank, nproc);
        HeatTransfer ht(settings);
        InlineIO io(settings, mpiHeatTransferComm);
        setupOutputIO(settings, mpiHeatTransferComm);

        ht.init(false);
        // ht.printT("Initialized T:", mpiHeatTransferComm);
        ht.heatEdges();
        ht.exchange(mpiHeatTransferComm);
        // ht.printT("Heated T:", mpiHeatTransferComm);

        io.write(ht);
        const double *Tin = io.read(true);
        Compute(Tin, Tout, dT, true);
        writeOutput();

        for (unsigned int t = 1; t <= settings.steps; ++t)
        {
            if (rank == 0)
                std::cout << "Step " << t << ":\n";
            for (unsigned int iter = 1; iter <= settings.iterations; ++iter)
            {
                ht.iterate();
                ht.exchange(mpiHeatTransferComm);
                ht.heatEdges();
            }

            io.write(ht);
            Tin = io.read(false);
            Compute(Tin, Tout, dT, false);
            writeOutput();
        }
        MPI_Barrier(mpiHeatTransferComm);

        double timeEnd = MPI_Wtime();
        if (rank == 0)
            std::cout << "Total runtime = " << timeEnd - timeStart << "s\n";

        writer.Close();
    }
    catch (std::invalid_argument &e) // command-line argument errors
    {
        std::cout << e.what() << std::endl;
        printUsage();
    }
    catch (std::ios_base::failure &e) // I/O failure (e.g. file not found)
    {
        std::cout << "I/O base exception caught\n";
        std::cout << e.what() << std::endl;
    }
    catch (std::exception &e) // All other exceptions
    {
        std::cout << "Exception caught\n";
        std::cout << e.what() << std::endl;
    }

    MPI_Finalize();
    return 0;
}
