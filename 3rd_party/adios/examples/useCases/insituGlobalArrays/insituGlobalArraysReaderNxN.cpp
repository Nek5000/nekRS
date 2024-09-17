/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * A Use Case for In Situ visulization frameworks (Conduit, SENSEI)
 *
 * Read in the variables that the Writer wrote.
 * Every process should read only what the corresponding Writer wrote
 * This is an N to N case
 *
 * Created on: Jul 11, 2017
 *      Author: pnorbert
 */

#include <algorithm> // std::transform
#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread> // sleep_for
#include <vector>

#include <adios2.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
MPI_Comm readerComm;
#endif

typedef struct
{
    std::string varName;
    std::string type;
    adios2::Dims shape;
    adios2::Dims start;
    adios2::Dims count;
    std::vector<double> data;
} VarInfo;

std::string DimsToString(adios2::Dims &dims)
{
    std::string s = "\"";
    for (size_t i = 0; i < dims.size(); i++)
    {
        if (i > 0)
        {
            s += ", ";
        }
        s += std::to_string(dims[i]);
    }
    s += "\"";
    return s;
}

/* Find the block written by rank for a given variable
 * Appends the block info to the passed Varinfo vector if found
 */
template <class T>
void ProcessVariableMetadata(int rank, const std::string &name, const std::string &type,
                             const adios2::Engine &reader, adios2::IO &io,
                             std::vector<VarInfo> &varinfos)
{
    const adios2::Variable<T> variable = io.InquireVariable<T>(name);
    std::vector<typename adios2::Variable<T>::Info> blocks =
        reader.BlocksInfo(variable, reader.CurrentStep());
    for (auto &block : blocks)
    {
        /* offset in first dimension is encoding writer's rank */
        if (block.Start.size() > 0 && block.Start[0] == static_cast<size_t>(rank))
        {
            /*std::cout << "        Rank " << rank << "     Variable '" << name
                      << "' found a block dimensions = "
                      << DimsToString(block.Count)
                      << " offset = " << DimsToString(block.Start) <<
               std::endl;*/
            varinfos.push_back({name, type, variable.Shape(), block.Start, block.Count});
        }
    }
}

std::vector<VarInfo> ProcessMetadata(int rank, const adios2::Engine &reader, adios2::IO &io,
                                     const std::map<std::string, adios2::Params> &varNameList)
{
    std::vector<VarInfo> varinfos;
    for (auto &var : varNameList)
    {
        const std::string &name(var.first);
        auto it = var.second.find("Type");
        const std::string &type = it->second;
        it = var.second.find("Shape");
        const std::string &shape = it->second;
        if (!rank)
        {
            std::cout << "    Variable '" << name << "' type " << type << " dimensions = " << shape
                      << std::endl;
        }
        if (type == "struct")
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (type == adios2::GetType<T>())                                                         \
    {                                                                                              \
        ProcessVariableMetadata<T>(rank, name, type, reader, io, varinfos);                        \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }
#if ADIOS2_USE_MPI
    // sync printouts
    MPI_Barrier(readerComm);
#endif
    return varinfos;
}

void ReadVariables(int rank, adios2::Engine &reader, adios2::IO &io, std::vector<VarInfo> &varinfos)
{
    for (auto &vi : varinfos)
    {
        /*std::cout << "       Read Variable Name: " << vi.varName
                  << " dimensions = " << DimsToString(vi.count)
                  << " offset = " << DimsToString(vi.start) << std::endl;*/
        if (vi.type == "double")
        {
            adios2::Variable<double> variable = io.InquireVariable<double>(vi.varName);
            variable.SetSelection({vi.start, vi.count});
            reader.Get(variable, vi.data);
        }
        else
        {
            std::cout << "ERROR: This example does not support reading "
                         "variables of type: "
                      << vi.type << ". Skip reading variable " << vi.varName << std::endl;
        }
    }
    reader.PerformGets();
}

void SerialPrintout(std::vector<VarInfo> &varinfos, int rank, int nproc)
{
// Serialize printout
#if ADIOS2_USE_MPI
    MPI_Barrier(readerComm);
    int token = 0;
    MPI_Status st;
    if (rank > 0)
    {
        MPI_Recv(&token, 1, MPI_INT, rank - 1, 0, readerComm, &st);
    }
#endif

    std::cout << "    Rank " << rank << " variables:" << varinfos.size() << std::endl;
    for (auto &vi : varinfos)
    {
        std::cout << "       Name: " << vi.varName << " dimensions = " << DimsToString(vi.count)
                  << " offset = " << DimsToString(vi.start) << " = [";
        for (auto d : vi.data)
        {
            std::cout << "  " << std::fixed << std::setprecision(2) << d;
        }
        std::cout << "  ]" << std::endl;
    }

#if ADIOS2_USE_MPI
    if (rank < nproc - 1)
    {
        std::chrono::milliseconds timespan(100);
        std::this_thread::sleep_for(timespan);
        MPI_Send(&token, 1, MPI_INT, rank + 1, 0, readerComm);
    }
    MPI_Barrier(readerComm);
#endif
}

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
            {"Library", "ZMQ"}, {"Timeout", "2000"}, {"IPAddress", "127.0.0.1"}, {"Port", "25600"}};
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
    const unsigned int color = 2;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &readerComm);
    MPI_Comm_rank(readerComm, &rank);
    MPI_Comm_size(readerComm, &nproc);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(readerComm);
#else
    adios2::ADIOS adios;
#endif

    ProcessArgs(rank, argc, argv);
    if (!rank)
    {
        std::cout << "Reader: ADIOS2 Engine set to: " << argEngine << "   Parameters:";
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

    try
    {
        adios2::IO io = adios.DeclareIO("Input");
        io.SetEngine(argEngine);
        io.SetParameters(engineParams);
        for (auto &t : engineTransports)
        {
            io.AddTransport(t.first, t.second);
        }
        adios2::Engine reader = io.Open("output.bp", adios2::Mode::Read);

        while (true)
        {
            adios2::StepStatus status = reader.BeginStep(adios2::StepMode::Read, 60.0f);
            if (status == adios2::StepStatus::NotReady)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                continue;
            }
            else if (status != adios2::StepStatus::OK)
            {
                break;
            }

            std::map<std::string, adios2::Params> varNameList = io.AvailableVariables();
            const size_t nTotalVars = varNameList.size();
            if (!rank)
            {
                std::cout << "File info:" << std::endl;
                std::cout << "  Current step:   " << reader.CurrentStep() << std::endl;
                std::cout << "  Total number of variables = " << nTotalVars << std::endl;
            }

            std::vector<VarInfo> varinfos = ProcessMetadata(rank, reader, io, varNameList);
            ReadVariables(rank, reader, io, varinfos);
            SerialPrintout(varinfos, rank, nproc);

            reader.EndStep();
        }

        // Called once: indicate that we are done with this output for the run
        reader.Close();
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
