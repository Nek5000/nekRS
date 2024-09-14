/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Write various outputs, that the ADIOS built-in campaign manager records
 * and one can read all data using a single reader engine on the campaign
 * output instead of the individual files
 *
 * Campaign managment handles multiple restarts of the writer, where
 * some existing outputs get appended to (at some specific step) and some
 * other new outputs are created.
 *
 * Campaign management handles BP4/BP5  outputs
 *
 * Campaign management needs to be turned on explicitly at configuration time
 *   -DADIOS2_USE_Campaign=ON
 *
 * Outputs:
 *   dataAll.bp:      written by every process
 *   // dataFirstRank.bp:    written by rank 0 only
 *   // dataLastRank.bp:     written by last rank only
 *   dataStep{N}.bp   a file series every step (written by all)
 *   dataAnother.h5:  written by rank 1..N from another ADIOS object
 *   dataNew{N}.bp :  a new file when restarting from step N
 *   dataFinal.bp:    written by rank 1 from a third ADIOS object at the end
 *
 * After running this example, the adios2_campaign_manager scripts must be used
 * to create/update/delete a campaign archive from this run
 * Then, the Campaign engine can be used in reading codes to use the campaign
 * archive to read the content.
 *
 * E.g.
 *   adios2_campaign_manager -n example_campaign-write_101 -c ~/.campaign-store create
 *   adios2_campaign_manager -n example_campaign-write_101 -c ~/.campaign-store info
 *   bpls -E campaign -P "cachepath=/tmp/campaign" \
 *        ~/.campaign-store/example_campaign-write_101.aca -la
 *
 * Created on: May 17, 2023
 *      Author: Norbert Podhorszki <pnorbert@ornl.gov>
 */

#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include <adios2.h>
#include <adios2/common/ADIOSConfig.h> // ADIOS2_USE_HDF5 macro
#include <adios2/helper/adiosString.h> // StringToSizeT()
#include <mpi.h>

/* Arguments */
size_t nSteps = 5;
size_t startStep = 0;

/* MPI variables */
int rank, nproc;
MPI_Comm comm;
int wrank, wnproc; // MPI_COMM_WORLD rank and nproc

MPI_Comm commAnother;
MPI_Comm commFirstRank;
MPI_Comm commLastRank;

constexpr size_t Nx = 6;
constexpr size_t Ny = 4;
constexpr size_t Nz = 5;

// imitating having a physical iteration counter and time value
uint64_t physicalStep = 0;
double physicalTime = 0.0;
constexpr uint64_t physicalStep_dt = 100;
constexpr double physicalTime_dt = 0.01;

void printUsage()
{
    std::cout << "Usage: campaign_writer  steps  [start_step] \n"
              << "  steps:     the total number of steps to output\n"
              << "  start_step: restart from this step\n\n";
}

unsigned int convertToUint(std::string varName, char *arg)
{
    char *end;
    unsigned int retval = std::strtoul(arg, &end, 10);
    if (end[0] || errno == ERANGE)
    {
        throw std::invalid_argument("Invalid value given for " + varName + ": " + std::string(arg));
    }
    return retval;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printUsage();
        return -1;
    }
    nSteps = adios2::helper::StringToSizeT(argv[1], "when parsing argument steps");
    if (argc > 2)
    {
        startStep = adios2::helper::StringToSizeT(argv[2], "when parsing argument start_step");
    }
    /*nSteps = convertToUint("steps", argv[1]);
    if (argc > 2)
    {
        startStep = convertToUint("start_step", argv[2]);
    }*/

    int provided;
    // MPI_THREAD_MULTIPLE is only required if you use SST with MPI DataPlane
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wnproc);
    if (wnproc < 2)
    {
        std::cout << "This MPI example needs at least 2 processes" << std::endl;
        return -2;
    }

    const unsigned int color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);
    MPI_Group group_world;
    MPI_Comm_group(comm, &group_world);

    // communicator for ranks 1..N
    const unsigned int colorAnother = (rank ? 1 : 0);
    MPI_Comm_split(comm, colorAnother, rank, &commAnother);

    // communicator for first rank
    MPI_Group group_first;
    int first_rank = 0;
    MPI_Group_incl(group_world, 1, &first_rank, &group_first);
    MPI_Comm_create(comm, group_first, &commFirstRank);

    // communicator for last rank
    MPI_Group group_last;
    int last_rank = nproc - 1;
    MPI_Group_incl(group_world, 1, &last_rank, &group_last);
    MPI_Comm_create(comm, group_last, &commLastRank);

    // some data
    std::vector<double> x(Nx);
    std::vector<double> y(Ny);
    std::vector<double> z(Nz);

    {
        adios2::ADIOS adiosAll(comm);
        adios2::ADIOS adiosAnother(commAnother);

        adios2::IO ioAll;
        adios2::IO ioFirstRank;
        adios2::IO ioLastRank;
        adios2::IO ioStep;
        adios2::IO ioNew;
        adios2::IO ioAnother;

        adios2::Engine writerAll;
        adios2::Engine writerFirstRank;
        adios2::Engine writerLastRank;
        adios2::Engine writerAnother;

        ioAll = adiosAll.DeclareIO("OutputAll");
        ioAll.SetEngine("BP5");

        if (rank == 0)
        {
            ioFirstRank = adiosAll.DeclareIO("OutputFirst");
            ioFirstRank.SetEngine("BP4");
        }
        if (rank == nproc - 1)
        {
            ioLastRank = adiosAll.DeclareIO("OutputLast");
            ioLastRank.SetEngine("BP5");
        }
        ioStep = adiosAll.DeclareIO("OutputStep");
        ioStep.SetEngine("BP4");

        ioNew = adiosAll.DeclareIO("OutputNew");
#ifdef ADIOS2_USE_HDF5
        ioNew.SetEngine("HDF5");
#else
        ioNew.SetEngine("BP5");
#endif

        if (rank > 0)
        {
            ioAnother = adiosAnother.DeclareIO("OutputAnother");
#ifdef ADIOS2_USE_HDF5
            ioAnother.SetEngine("HDF5");
#else
            ioAnother.SetEngine("BP5");
#endif
        }

        adios2::Mode mode = adios2::Mode::Write;
        if (startStep > 0)
        {
            mode = adios2::Mode::Append;
            auto v = std::to_string(startStep);
            ioAll.SetParameter("AppendAfterSteps", v);
            if (rank == 0)
            {
                ioFirstRank.SetParameter("AppendAfterSteps", v);
            }
            if (rank == nproc - 1)
            {
                ioLastRank.SetParameter("AppendAfterSteps", v);
            }
            if (rank > 0)
            {
                ioAnother.SetParameter("AppendAfterSteps", v);
            }
        }

        adios2::Variable<double> varXAll = ioAll.DefineVariable<double>(
            "x", {static_cast<size_t>(nproc), Nx}, {static_cast<size_t>(rank), 0}, {1, Nx});
        adios2::Variable<double> varYAll = ioAll.DefineVariable<double>(
            "y", {static_cast<size_t>(nproc), Ny}, {static_cast<size_t>(rank), 0}, {1, Ny});
        adios2::Variable<double> varZAll = ioAll.DefineVariable<double>(
            "z", {static_cast<size_t>(nproc), Nz}, {static_cast<size_t>(rank), 0}, {1, Nz});
        adios2::Variable<size_t> varStepAll = ioAll.DefineVariable<size_t>("AdiosStep");
        adios2::Variable<size_t> varPhysStepAll = ioAll.DefineVariable<size_t>("iteration");
        adios2::Variable<double> varPhysTimeAll = ioAll.DefineVariable<double>("time");
        ioAll.DefineAttribute<std::string>("comment", "Written by all processes");

        adios2::Variable<double> varZStep = ioStep.DefineVariable<double>(
            "z", {(unsigned int)nproc, Nz}, {static_cast<size_t>(rank), 0}, {1, Nz});
        adios2::Variable<size_t> varStepStep = ioStep.DefineVariable<size_t>("AdiosStep");
        adios2::Variable<size_t> varPhysStepStep = ioStep.DefineVariable<size_t>("iteration");
        adios2::Variable<double> varPhysTimeStep = ioStep.DefineVariable<double>("time");

        adios2::Variable<double> varYNew = ioStep.DefineVariable<double>(
            "y", {(unsigned int)nproc, Ny}, {static_cast<size_t>(rank), 0}, {1, Ny});
        adios2::Variable<size_t> varStepNew = ioNew.DefineVariable<size_t>("AdiosStep");
        adios2::Variable<size_t> varPhysStepNew = ioNew.DefineVariable<size_t>("iteration");
        adios2::Variable<double> varPhysTimeNew = ioNew.DefineVariable<double>("time");

        writerAll = ioAll.Open("dataAll.bp", mode);

        adios2::Variable<double> varXFirstRank;
        if (rank == 0)
        {
            varXFirstRank = ioFirstRank.DefineVariable<double>("x", {Nx}, {0}, {Nx});
            ioFirstRank.DefineAttribute<std::string>("comment", "Written by rank 0");
            writerFirstRank = ioFirstRank.Open("dataFirstRank.bp", mode, commFirstRank);
        }

        adios2::Variable<double> varXLastRank;
        if (rank == nproc - 1)
        {
            varXLastRank = ioLastRank.DefineVariable<double>("x", {Nx}, {0}, {Nx});
            ioLastRank.DefineAttribute<std::string>("comment",
                                                    "Written by rank " + std::to_string(nproc - 1));
            writerLastRank = ioLastRank.Open("dataLastRank.bp", mode, commLastRank);
        }

        // the other ADIOS object (valid on rank 1..N, not on rank 0)
        adios2::Variable<double> varXAnother;
        if (rank > 0)
        {
            varXAnother =
                ioAnother.DefineVariable<double>("x", {static_cast<size_t>(nproc - 1), Nx});
            ioAnother.DefineAttribute<std::string>("comment", "Written by ranks 1.." +
                                                                  std::to_string(nproc - 1));
            writerAnother = ioAnother.Open("dataAnother.bp", mode);
        }

        physicalStep = physicalStep_dt * startStep;
        physicalTime = physicalTime_dt * startStep;
        for (size_t step = startStep; step < startStep + nSteps; step++)
        {

            for (size_t i = 0; i < Nx; i++)
            {
                x[i] = step * Nx * nproc * 1.0 + rank * Nx * 1.0 + (double)i;
            }
            for (size_t i = 0; i < Ny; i++)
            {
                y[i] = step * Ny * nproc * 1.0 + rank * Ny * 1.0 + (double)i;
            }
            for (size_t i = 0; i < Nz; i++)
            {
                z[i] = step * Nz * nproc * 1.0 + rank * Nz * 1.0 + (double)i;
            }

            writerAll.BeginStep();
            writerAll.Put(varXAll, x.data());
            writerAll.Put(varYAll, y.data());
            writerAll.Put(varZAll, z.data());
            writerAll.Put(varStepAll, step);
            writerAll.Put(varPhysStepAll, (size_t)physicalStep);
            writerAll.Put(varPhysTimeAll, physicalTime);

            writerAll.EndStep();

            adios2::Engine writerStep =
                ioStep.Open("dataStep" + std::to_string(step) + ".bp", adios2::Mode::Write);
            writerStep.Put(varZStep, z.data());
            writerStep.Put(varStepStep, step);
            writerStep.Put(varPhysStepStep, (size_t)physicalStep);
            writerStep.Put(varPhysTimeStep, physicalTime);
            writerStep.Close();

            if (rank > 0)
            {
                writerAnother.BeginStep();
                varXAnother.SetSelection(adios2::Box<adios2::Dims>(
                    {static_cast<size_t>(rank - 1), 0}, {1, static_cast<size_t>(Nx)}));
                writerAnother.Put(varXAnother, x.data());
                writerAnother.EndStep();
            }

            if (rank == 0)
            {
                writerFirstRank.BeginStep();
                writerFirstRank.Put(varXFirstRank, x.data());
                writerFirstRank.EndStep();
            }

            if (rank == nproc - 1)
            {
                writerLastRank.BeginStep();
                writerLastRank.Put(varXLastRank, x.data());
                writerLastRank.EndStep();
            }

            if (startStep > 0 && step == startStep)
            {
                adios2::Engine writerNew =
                    ioNew.Open("dataNew" + std::to_string(step) + ".bp", adios2::Mode::Write);
                writerNew.Put(varYNew, y.data());
                writerNew.Put(varStepNew, step);
                writerNew.Put(varPhysStepNew, (size_t)physicalStep);
                writerNew.Put(varPhysTimeNew, physicalTime);
                writerNew.Close();
            }

            std::this_thread::sleep_for(std::chrono::duration<double>(1.0));

            physicalStep += physicalStep_dt;
            physicalTime += physicalTime_dt;
        }

        writerAll.Close();
        if (rank == 0)
        {
            writerFirstRank.Close();
        }
        if (rank == nproc - 1)
        {
            writerLastRank.Close();
        }
        if (rank > 0)
        {
            writerAnother.Close();
        }
    }

    // Final: to test creating yet another ADIOS object later
    {
        if (rank == 1)
        {
            adios2::ADIOS adiosFinal;
            adios2::IO ioFinal = adiosFinal.DeclareIO("OutputFinal");
            ioFinal.SetEngine("BP5");
            adios2::Mode mode = adios2::Mode::Write;

            adios2::Variable<double> varXFinal =
                ioFinal.DefineVariable<double>("x_on_rank_1", {Nx}, {0}, {Nx});
            adios2::Variable<size_t> varStepFinal = ioFinal.DefineVariable<size_t>("AdiosStep");
            adios2::Variable<size_t> varPhysStepFinal = ioFinal.DefineVariable<size_t>("iteration");
            adios2::Variable<double> varPhysTimeFinal = ioFinal.DefineVariable<double>("time");
            ioFinal.DefineAttribute<std::string>("comment", "Written by rank 1 at end of run");

            adios2::Engine writerFinal = ioFinal.Open("dataFinal.bp", mode);
            writerFinal.Put(varXFinal, x.data());
            writerFinal.Put(varStepFinal, startStep + nSteps - 1);
            writerFinal.Put(varPhysStepFinal, (size_t)physicalStep);
            writerFinal.Put(varPhysTimeFinal, physicalTime);
            writerFinal.Close();
        }
    }

    MPI_Finalize();

    return 0;
}
