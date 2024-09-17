/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <mpi.h>

#include "adios2.h"

#include <cstdint>
#include <iostream>
#include <vector>

#include "PrintData.h"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if (argc < 2)
    {
        std::cout << "Not enough arguments: need an input file\n";
        return 1;
    }
    const char *inputfile = argv[1];

    /* World comm spans all applications started with the same aprun command
     on a Cray XK6. So we have to split and create the local
     'world' communicator for the reader only.
     In normal start-up, the communicator will just equal the MPI_COMM_WORLD.
     */

    int wrank, wnproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wnproc);

    const unsigned int color = 2;
    MPI_Comm mpiReaderComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &mpiReaderComm);

    int rank, nproc;
    MPI_Comm_rank(mpiReaderComm, &rank);
    MPI_Comm_size(mpiReaderComm, &nproc);

    adios2::ADIOS ad(mpiReaderComm);

    // Define method for engine creation
    // 1. Get method def from config file or define new one

    adios2::IO bpReaderIO = ad.DeclareIO("input");
    if (!bpReaderIO.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        bpReaderIO.SetEngine("BP");
        bpReaderIO.SetParameters({{"num_threads", "2"}});
        bpReaderIO.SetParameter("OpenAsFile", "true");

        // ISO-POSIX file is the default transport
        // Passing parameters to the transport
        bpReaderIO.AddTransport("File", {{"verbose", "4"}});
    }

    adios2::Engine bpReader = bpReaderIO.Open(inputfile, adios2::Mode::Read, mpiReaderComm);

    unsigned int gndx = 0;
    unsigned int gndy = 0;
    // bpReader->Read<unsigned int>("gndx", &gndx);
    // bpReader->Read<unsigned int>("gndy", &gndy);

    adios2::Variable<unsigned int> vgndx = bpReaderIO.InquireVariable<unsigned int>("gndx");
    (void)vgndx;

    // gndx = vgndx.GetData()[0];

    adios2::Variable<unsigned int> vgndy = bpReaderIO.InquireVariable<unsigned int>("gndy");
    // gndy = vgndy.GetData()[0];

    if (rank == 0)
    {
        std::cout << "gndx       = " << gndx << std::endl;
        std::cout << "gndy       = " << gndy << std::endl;
        std::cout << "# of steps = " << vgndy.Steps() << std::endl;
    }

    // 1D decomposition of the columns, which is inefficient for reading!
    adios2::Dims readsize({gndx, gndy / nproc});
    adios2::Dims offset({0LL, rank * readsize[1]});
    if (rank == nproc - 1)
    {
        // last process should read all the rest of columns
        readsize[1] = gndy - readsize[1] * (nproc - 1);
    }

    std::cout << "rank " << rank << " reads " << readsize[1] << " columns from offset " << offset[1]
              << std::endl;

    adios2::Variable<double> vT = bpReaderIO.InquireVariable<double>("T");

    double *T = new double[vT.Steps() * readsize[0] * readsize[1]];

    // Create a 2D selection for the subset
    vT.SetSelection(adios2::Box<adios2::Dims>(offset, readsize));
    vT.SetStepSelection(adios2::Box<std::size_t>(0, vT.Steps()));

    // Arrays are read by scheduling one or more of them
    // and performing the reads at once
    // bpReader->ScheduleRead<double>(*vT, T);
    // bpReader->PerformReads(adios2::ReadMode::Blocking);

    printData(T, readsize.data(), offset.data(), rank, vT.Steps());
    bpReader.Close();
    delete[] T;
    MPI_Finalize();
    return 0;
}
