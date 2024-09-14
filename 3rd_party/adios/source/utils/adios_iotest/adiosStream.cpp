/*
 * adiosStream.cpp
 *
 *  Created on: Nov 2018
 *      Author: Norbert Podhorszki
 */

#include "adiosStream.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

adiosStream::adiosStream(const std::string &streamName, adios2::IO &io, const adios2::Mode mode,
                         MPI_Comm comm, bool iotimer, size_t appid)
: Stream(streamName, mode), io(io), comm(comm)
{
    int myRank, totalRanks;
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &totalRanks);
    double timeStart, timeEnd;

    // double maxOpenTime, minOpenTime;
    hasIOTimer = iotimer;
    appID = appid;

    if (mode == adios2::Mode::Write)
    {
        timeStart = MPI_Wtime();
        engine = io.Open(streamName, adios2::Mode::Write, comm);
        timeEnd = MPI_Wtime();
    }
    else
    {
        timeStart = MPI_Wtime();
        engine = io.Open(streamName, adios2::Mode::Read, comm);
        timeEnd = MPI_Wtime();
    }
    if (hasIOTimer)
    {
        openTime = timeEnd - timeStart;
        std::vector<double> opentime_all_ranks;
        if (myRank == 0)
        {
            std::string logfilename = "app" + std::to_string(appID) + "_perf.csv";
            perfLogFP = fopen(logfilename.c_str(), "w");
            fputs("step,rank,operation,time\n", perfLogFP);
            std::cout << "performance log file open succeeded!" << std::endl;
            opentime_all_ranks.reserve(totalRanks);
        }
        MPI_Gather(&openTime, 1, MPI_DOUBLE, opentime_all_ranks.data(), 1, MPI_DOUBLE, 0, comm);
        if (myRank == 0)
        {
            for (int i = 0; i < totalRanks; i++)
            {
                std::string content = std::to_string(engine.CurrentStep()) + "," +
                                      std::to_string(i) + ",open," +
                                      std::to_string(opentime_all_ranks[i]) + "\n";
                // std::cout << content;
                fputs(content.c_str(), perfLogFP);
            }
        }
    }
    // openTime = timeEnd - timeStart;
    // MPI_Allreduce(&openTime, &maxOpenTime, 1, MPI_DOUBLE, MPI_MAX, comm);
    // MPI_Allreduce(&openTime, &minOpenTime, 1, MPI_DOUBLE, MPI_MIN, comm);
    // if (myRank == 0)
    // {
    //     std::cout << "        Max open time = " << maxOpenTime << std::endl;
    //     std::cout << "        Min open time = " << minOpenTime << std::endl;
    //     std::ofstream open_perf_log;
    //     open_perf_log.open("open_perf.txt", std::ios::app);
    //     open_perf_log << std::to_string(maxOpenTime) + ", " +
    //                          std::to_string(minOpenTime) + "\n";
    //     open_perf_log.close();
    // }
}

adiosStream::~adiosStream() {}

void adiosStream::defineADIOSArray(const std::shared_ptr<VariableInfo> ov)
{
    if (ov->type == "double")
    {
        adios2::Variable<double> v =
            io.DefineVariable<double>(ov->name, ov->shape, ov->start, ov->count, true);
        (void)v;
        // v = io->InquireVariable<double>(ov->name);
    }
    else if (ov->type == "float")
    {
        adios2::Variable<float> v =
            io.DefineVariable<float>(ov->name, ov->shape, ov->start, ov->count, true);
        (void)v;
    }
    else if (ov->type == "int")
    {
        adios2::Variable<int> v =
            io.DefineVariable<int>(ov->name, ov->shape, ov->start, ov->count, true);
        (void)v;
    }
}

void adiosStream::putADIOSArray(const std::shared_ptr<VariableInfo> ov)
{
    if (ov->type == "double")
    {
        const double *a = reinterpret_cast<const double *>(ov->data.data());
        engine.Put<double>(ov->name, a);
    }
    else if (ov->type == "float")
    {
        const float *a = reinterpret_cast<const float *>(ov->data.data());
        engine.Put<float>(ov->name, a);
    }
    else if (ov->type == "int")
    {
        const int *a = reinterpret_cast<const int *>(ov->data.data());
        engine.Put<int>(ov->name, a);
    }
}

void adiosStream::getADIOSArray(std::shared_ptr<VariableInfo> ov)
{
    // Allocate memory on first access
    if (!ov->data.size())
    {
        ov->data.resize(ov->datasize);
    }
    if (ov->type == "double")
    {
        adios2::Variable<double> v = io.InquireVariable<double>(ov->name);
        if (!v)
        {
            ov->readFromInput = false;
            return;
        }
        v.SetSelection({ov->start, ov->count});
        double *a = reinterpret_cast<double *>(ov->data.data());
        engine.Get<double>(v, a);
        ov->readFromInput = true;
    }
    else if (ov->type == "float")
    {
        adios2::Variable<float> v = io.InquireVariable<float>(ov->name);
        if (!v)
        {
            ov->readFromInput = false;
            return;
        }
        v.SetSelection({ov->start, ov->count});
        float *a = reinterpret_cast<float *>(ov->data.data());
        engine.Get<float>(v, a);
        ov->readFromInput = true;
    }
    else if (ov->type == "int")
    {
        adios2::Variable<int> v = io.InquireVariable<int>(ov->name);
        if (!v)
        {
            ov->readFromInput = false;
            return;
        }
        v.SetSelection({ov->start, ov->count});
        int *a = reinterpret_cast<int *>(ov->data.data());
        engine.Get<int>(v, a);
        ov->readFromInput = true;
    }
}

/* return true if read-in completed */
adios2::StepStatus adiosStream::readADIOS(CommandRead *cmdR, Config &cfg, const Settings &settings,
                                          size_t step)
{
    if (!settings.myRank && settings.verbose)
    {
        std::cout << "    Read " << cmdR->streamName << " with timeout value " << cmdR->timeout_sec
                  << " using the group " << cmdR->groupName;
        if (!cmdR->variables.empty())
        {
            std::cout << " with selected variables:  ";
            for (const auto &v : cmdR->variables)
            {
                std::cout << v->name << " ";
            }
        }
        std::cout << "... " << std::endl;
    }
    double timeStart, timeEnd;
    double readTime;
    // double maxReadTime, minReadTime;
    MPI_Barrier(comm);
    timeStart = MPI_Wtime();
    adios2::StepStatus status = engine.BeginStep(cmdR->stepMode, cmdR->timeout_sec);
    if (status != adios2::StepStatus::OK)
    {
        return status;
    }

    if (!settings.myRank && settings.verbose && step == 1)
    {
        const auto varmap = io.AvailableVariables();
        std::cout << "    Variables in input for reading: " << std::endl;
        for (const auto &v : varmap)
        {
            std::cout << "        " << v.first << std::endl;
        }
    }

    if (!settings.myRank && settings.verbose)
    {
        std::cout << "    Read data " << std::endl;
    }

    for (auto ov : cmdR->variables)
    {
        getADIOSArray(ov);
    }

    if (step == 1 && settings.fixedPattern)
    {
        if (!settings.myRank && settings.verbose)
        {
            std::cout << "        Lock Reader Selections for Fixed Pattern before "
                         "first EndStep"
                      << std::endl;
        }
        engine.LockWriterDefinitions();
    }

    engine.EndStep();
    timeEnd = MPI_Wtime();
    int myRank, totalRanks;
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &totalRanks);
    if (hasIOTimer)
    {
        readTime = timeEnd - timeStart;
        // double *writetime_all_ranks = NULL;
        std::vector<double> readtime_all_ranks;
        if (myRank == 0)
        {
            // writetime_all_ranks = (double *)malloc(totalRanks);
            //	if (!writetime_all_ranks)
            //	{
            //	    std::cout << "rank " << myRank << ": malloc failed!" <<
            // std::endl;
            //	}
            //	else
            //	{
            //	    std::cout << "rank " << myRank << ": malloc succeeded!" <<
            // std::endl;
            //	}
            readtime_all_ranks.reserve(totalRanks);
        }

        MPI_Gather(&readTime, 1, MPI_DOUBLE, readtime_all_ranks.data(), 1, MPI_DOUBLE, 0, comm);
        if (myRank == 0)
        {
            for (int i = 0; i < totalRanks; i++)
            {
                std::string content = std::to_string(engine.CurrentStep()) + "," +
                                      std::to_string(i) + ",read," +
                                      std::to_string(readtime_all_ranks[i]) + "\n";
                // std::cout << content;
                fputs(content.c_str(), perfLogFP);
            }
        }
    }

    // for (auto ov : cmdR->variables)
    // {
    //     if (settings.myRank == 1)
    //     {
    //         size_t varsize = std::accumulate(ov->count.begin(),
    //         ov->count.end(), 1,
    //                             std::multiplies<std::size_t>());
    //         std::cout << ov->name << ", " << varsize << std::endl;
    //         const double *dd = reinterpret_cast<double *>(ov->data.data());
    //         for (int j = 0; j < varsize; j++)
    //             std::cout << dd[j] << ", ";
    //         std::cout << std::endl;
    //     }
    // }
    return status;
}

void adiosStream::writeADIOS(CommandWrite *cmdW, Config &cfg, const Settings &settings, size_t step)
{
    if (!settings.myRank && settings.verbose)
    {
        std::cout << "    Write to output " << cmdW->streamName << " the group " << cmdW->groupName;
        if (!cmdW->variables.empty())
        {
            std::cout << " with selected variables:  ";
            for (const auto &v : cmdW->variables)
            {
                std::cout << v->name << " ";
            }
        }
        std::cout << std::endl;
    }

    size_t s = (cfg.nSteps > 0 ? cfg.nSteps : 1000);
    const double div = pow(10.0, static_cast<double>(settings.ndigits(s - 1)));
    double myValue = static_cast<double>(settings.myRank) + static_cast<double>(step - 1) / div;

    std::map<std::string, adios2::Params> definedVars = io.AvailableVariables();
    for (auto ov : cmdW->variables)
    {
        // Allocate memory on first access
        if (!ov->data.size())
        {
            ov->data.resize(ov->datasize);
        }
        // if the variable is not in the IO group it means
        // we have not defined it yet (e.g. a write-only variable or a linked
        // variable defined in another read group)
        const auto it = definedVars.find(ov->name);
        if (it == definedVars.end())
        {
            if (!settings.myRank && settings.verbose)
            {
                std::cout << "        Define array  " << ov->name << "  for output" << std::endl;
            }
            defineADIOSArray(ov);
        }

        // if we read the variable, use the read values otherwise generate data
        // now
        if (!ov->readFromInput)
        {
            if (!settings.myRank && settings.verbose)
            {
                std::cout << "        Fill array  " << ov->name << "  for output" << std::endl;
            }
            fillArray(ov, myValue);
        }
    }

    if (!settings.myRank && settings.verbose)
    {
        std::cout << "        Write data " << std::endl;
    }
    double timeStart, timeEnd;
    double writeTime;
    // double maxWriteTime, minWriteTime;
    MPI_Barrier(comm);
    timeStart = MPI_Wtime();
    engine.BeginStep();
    for (const auto &ov : cmdW->variables)
    {
        putADIOSArray(ov);
    }

    if (step == 1 && settings.fixedPattern)
    {
        if (!settings.myRank && settings.verbose)
        {
            std::cout << "        Lock Writer Definitions for Fixed Pattern before "
                         "first EndStep"
                      << std::endl;
        }
        engine.LockWriterDefinitions();
    }

    engine.EndStep();
    timeEnd = MPI_Wtime();
    int myRank, totalRanks;
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &totalRanks);
    if (hasIOTimer)
    {
        writeTime = timeEnd - timeStart;
        // double *writetime_all_ranks = NULL;
        std::vector<double> writetime_all_ranks;
        if (myRank == 0)
        {
            // writetime_all_ranks = (double *)malloc(totalRanks);
            //	if (!writetime_all_ranks)
            //	{
            //	    std::cout << "rank " << myRank << ": malloc failed!" <<
            // std::endl;
            //	}
            //	else
            //	{
            //	    std::cout << "rank " << myRank << ": malloc succeeded!" <<
            // std::endl;
            //	}
            writetime_all_ranks.reserve(totalRanks);
        }

        MPI_Gather(&writeTime, 1, MPI_DOUBLE, writetime_all_ranks.data(), 1, MPI_DOUBLE, 0, comm);
        if (myRank == 0)
        {
            for (int i = 0; i < totalRanks; i++)
            {
                std::string content = std::to_string(engine.CurrentStep()) + "," +
                                      std::to_string(i) + ",write," +
                                      std::to_string(writetime_all_ranks[i]) + "\n";
                // std::cout << content;
                fputs(content.c_str(), perfLogFP);
            }
        }
        // std::cout << "rank " << myRank << ": before
        // free(writetime_all_ranks)..." << std::endl;
        // free(writetime_all_ranks);
        // std::cout << "rank " << myRank << ": after
        // free(writetime_all_ranks)..." << std::endl; writetime_all_ranks =
        // NULL;
    }
    //    if (settings.ioTimer)
    //    {
    //        writeTime = timeEnd - timeStart;
    //        MPI_Allreduce(&writeTime, &maxWriteTime, 1, MPI_DOUBLE, MPI_MAX,
    //        comm); MPI_Allreduce(&writeTime, &minWriteTime, 1, MPI_DOUBLE,
    //        MPI_MIN, comm); if (settings.myRank == 0)
    //        {
    //            std::cout << "    App " << settings.appId
    //                      << ": Max write time = " << maxWriteTime <<
    //                      std::endl;
    //            std::cout << "    App " << settings.appId
    //                      << ": Min write time = " << minWriteTime <<
    //                      std::endl;
    //            std::ofstream wr_perf_log;
    //            wr_perf_log.open("write_perf_" +
    //            std::to_string(settings.appId) +
    //                                 ".txt",
    //                             std::ios::app);
    //            wr_perf_log << std::to_string(maxWriteTime) + ", " +
    //                               std::to_string(minWriteTime) + "\n";
    //            wr_perf_log.close();
    //        }
    //    }
}

void adiosStream::Write(CommandWrite *cmdW, Config &cfg, const Settings &settings, size_t step)
{

    writeADIOS(cmdW, cfg, settings, step);
}

adios2::StepStatus adiosStream::Read(CommandRead *cmdR, Config &cfg, const Settings &settings,
                                     size_t step)
{
    return readADIOS(cmdR, cfg, settings, step);
}

void adiosStream::Close()
{
    int myRank, totalRanks;
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &totalRanks);
    if (hasIOTimer)
    {
        if (myRank == 0)
        {
            fclose(perfLogFP);
        }
    }
    engine.Close();
}
