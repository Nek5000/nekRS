/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_ADIOS2.cpp
 *
 *  Created on: Nov 2017
 *      Author: Norbert Podhorszki
 *
 */
#include <mpi.h>

#include "adios2.h"

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "PrintDataStep.h"
#include "ReadSettings.h"

void printUsage()
{
    std::cout << "Usage: heatRead  config  input  output N  M \n"
              << "  config:  XML config file to use\n"
              << "  input:   name of input data file/stream\n"
              << "  output:  name of output data file/stream\n"
              << "  N:       number of processes in X dimension\n"
              << "  M:       number of processes in Y dimension\n\n";
}

void Compute(const std::vector<double> &Tin, std::vector<double> &Tout, std::vector<double> &dT,
             bool firstStep)
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

int main(int argc, char *argv[])
{
    int provided;
    std::string engineName = std::string(argv[argc - 1]);

    int threadSupportLevel = MPI_THREAD_SINGLE;
    if (engineName == "SST")
    {
        threadSupportLevel = MPI_THREAD_MULTIPLE;
    }

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, threadSupportLevel, &provided);

    /* When writer and reader is launched together with a single mpirun command,
       the world comm spans all applications. We have to split and create the
       local 'world' communicator for the reader only.
       When writer and reader is launched separately, the mpiReaderComm
       communicator will just equal the MPI_COMM_WORLD.
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

    try
    {
        ReadSettings settings(argc, argv, rank, nproc);
        adios2::ADIOS ad(settings.configfile, mpiReaderComm);

        // Define method for engine creation
        // 1. Get method def from config file or define new one

        adios2::IO inIO = ad.DeclareIO("readerInput");
        if (!inIO.InConfigFile())
        {
            // if not defined by user, we can change the default settings
            // BPFile is the default engine
            inIO.SetEngine("BPFile");
            inIO.SetParameters({{"num_threads", "1"}});

            // ISO-POSIX file output is the default transport (called "File")
            // Passing parameters to the transport
            inIO.AddTransport("File", {{"verbose", "4"}});
        }

        adios2::IO outIO = ad.DeclareIO("readerOutput");

        adios2::Engine reader = inIO.Open(settings.inputfile, adios2::Mode::Read, mpiReaderComm);

        std::vector<double> Tin;
        std::vector<double> Tout;
        std::vector<double> dT;
        adios2::Variable<double> vTin;
        adios2::Variable<double> vTout;
        adios2::Variable<double> vdT;
        adios2::Engine writer;
        bool firstStep = true;
        int step = 0;

        while (true)
        {
            adios2::StepStatus status = reader.BeginStep(adios2::StepMode::Read);
            if (status != adios2::StepStatus::OK)
            {
                break;
            }

            // Variable objects disappear between steps so we need this every
            // step
            vTin = inIO.InquireVariable<double>("T");

            if (!vTin)
            {
                std::cout << "Error: NO variable T found. Unable to proceed. "
                             "Exiting. "
                          << std::endl;
                break;
            }

            if (firstStep)
            {
                // Promise that we are not going to change the variable sizes
                // nor add new variables
                reader.LockReaderSelections();

                unsigned int gndx = static_cast<unsigned int>(vTin.Shape()[0]);
                unsigned int gndy = static_cast<unsigned int>(vTin.Shape()[1]);

                if (rank == 0)
                {
                    std::cout << "gndx       = " << gndx << std::endl;
                    std::cout << "gndy       = " << gndy << std::endl;
                }

                settings.DecomposeArray(gndx, gndy);
                Tin.resize(settings.readsize[0] * settings.readsize[1]);
                Tout.resize(settings.readsize[0] * settings.readsize[1]);
                dT.resize(settings.readsize[0] * settings.readsize[1]);

                /* Create output variables and open output stream */
                vTout = outIO.DefineVariable<double>("T", {gndx, gndy}, settings.offset,
                                                     settings.readsize);
                vdT = outIO.DefineVariable<double>("dT", {gndx, gndy}, settings.offset,
                                                   settings.readsize);
                writer = outIO.Open(settings.outputfile, adios2::Mode::Write, mpiReaderComm);

                MPI_Barrier(mpiReaderComm); // sync processes just for stdout
            }

            if (!rank)
            {
                std::cout << "Processing step " << step << std::endl;
            }

            // Create a 2D selection for the subset
            vTin.SetSelection(adios2::Box<adios2::Dims>(settings.offset, settings.readsize));

            // Arrays are read by scheduling one or more of them
            // and performing the reads at once
            reader.Get<double>(vTin, Tin.data());
            /*printDataStep(Tin.data(), settings.readsize.data(),
                          settings.offset.data(), rank, step); */
            reader.EndStep();

            /* Compute dT from current T (Tin) and previous T (Tout)
             * and save Tin in Tout for output and for future computation
             */
            Compute(Tin, Tout, dT, firstStep);

            /* Output Tout and dT */
            writer.BeginStep();

            if (vTout)
                writer.Put<double>(vTout, Tout.data());
            if (vdT)
                writer.Put<double>(vdT, dT.data());
            writer.EndStep();

            step++;
            firstStep = false;
        }
        reader.Close();
        if (writer)
            writer.Close();
    }
    catch (std::invalid_argument &e) // command-line argument errors
    {
        std::cout << e.what() << std::endl;
        printUsage();
    }

    MPI_Finalize();
    return 0;
}
