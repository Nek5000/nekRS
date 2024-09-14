#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "adios2.h"
#include "mpi.h"

#include "decomp.h"
#include "processConfig.h"
#include "settings.h"
#include "stream.h"

int main(int argc, char *argv[])
{
    Settings settings;

    /* Check input arguments. Quit if something is wrong. */
    if (settings.processArguments(argc, argv))
    {
        return 1;
    }

    int provided;
    int threadSupportLevel = MPI_THREAD_SINGLE;
    if (settings.multithreadedMPI)
    {
        threadSupportLevel = MPI_THREAD_MULTIPLE;
    }
    MPI_Init_thread(&argc, &argv, threadSupportLevel, &provided);

    settings.initDecomp(MPI_COMM_WORLD);

    // MPI-dependent argument checks
    if (settings.extraArgumentChecks())
    {
        return 1;
    }

    adios2::ADIOS adios;
    if (settings.adiosConfigFileName.empty())
    {
        if (!settings.myRank && settings.verbose)
        {
            std::cout << "Use ADIOS without XML configuration " << std::endl;
        }
        adios = adios2::ADIOS(settings.appComm);
    }
    else
    {
        if (!settings.myRank && settings.verbose)
        {
            std::cout << "Use ADIOS xml file " << settings.adiosConfigFileName << std::endl;
        }
        adios = adios2::ADIOS(settings.adiosConfigFileName, settings.appComm);
    }
    Config cfg;
    size_t currentConfigLineNumber = 0;

    try
    {
        cfg = processConfig(settings, &currentConfigLineNumber);
    }
    catch (std::invalid_argument &e) // config file processing errors
    {
        if (!settings.myRank)
        {
            if (!currentConfigLineNumber)
            {
                std::cout << "Config file error: " << e.what() << std::endl;
            }
            else
            {
                std::cout << "Config file error in line " << currentConfigLineNumber << ": "
                          << e.what() << std::endl;
            }
        }

        /* Quit calmly */
        MPI_Finalize();
        return 1;
    }

    double timeStart, timeEnd;
    MPI_Barrier(settings.appComm);
    timeStart = MPI_Wtime();

    uint64_t actualBusyTime_usec = 0;

    try
    {
        /* writing to one stream using two groups is not supported.
         * FIXME: we need to check for this condition and raise error
         */
        if (!settings.myRank && settings.verbose)
        {
            std::cout << "Start App " + std::to_string(settings.appId) + ": " << std::endl;
        }
        /* 1. Assign stream names with group names that appear in
           commands */
        // map of <streamName, groupName>
        std::map<std::string, std::string> groupMap;
        // a vector of streams in the order they appear
        std::vector<std::pair<std::string, Operation>> streamsInOrder;
        for (const auto &cmd : cfg.commands)
        {
            if (cmd->op == Operation::Write)
            {
                auto cmdW = dynamic_cast<CommandWrite *>(cmd.get());
                groupMap[cmdW->streamName] = cmdW->groupName;
                streamsInOrder.push_back(std::make_pair(cmdW->streamName, Operation::Write));
            }
            else if (cmd->op == Operation::Read)
            {
                auto cmdR = dynamic_cast<CommandRead *>(cmd.get());
                groupMap[cmdR->streamName] = cmdR->groupName;
                streamsInOrder.push_back(std::make_pair(cmdR->streamName, Operation::Read));
            }
        }

        std::map<std::string, std::shared_ptr<ioGroup>> ioMap;

        /* 2. Declare/define groups and open streams in the order they
         * appear */
        std::map<std::string, std::shared_ptr<Stream>> readStreamMap;
        std::map<std::string, std::shared_ptr<Stream>> writeStreamMap;

        for (const auto &st : streamsInOrder)
        {
            std::string streamName = st.first;
            std::shared_ptr<ioGroup> io;
            auto &groupName = groupMap[streamName];
            auto it = ioMap.find(groupName);
            if (it == ioMap.end())
            {
                io = createGroup(groupName, settings.iolib, adios);
                ioMap[groupName] = io;
            }
            else
            {
                io = it->second;
            }
            const bool isWrite = (st.second == Operation::Write);
            if (isWrite)
            {
                auto it = writeStreamMap.find(streamName);
                if (it == writeStreamMap.end())
                {
                    if (!settings.myRank && settings.verbose)
                    {
                        std::cout << "    Create Output Stream " << streamName << "... "
                                  << std::endl;
                    }
                    if (!settings.outputPath.empty())
                    {
                        std::string outputPath = settings.outputPath;
                        if (settings.outputPath.back() != '/')
                        {
                            outputPath += '/';
                        }

                        streamName = outputPath + streamName;
                    }
                    std::shared_ptr<Stream> writer =
                        openStream(streamName, io, adios2::Mode::Write, settings.iolib,
                                   settings.appComm, settings.ioTimer, settings.appId);
                    writeStreamMap[st.first] = writer;
                }
            }
            else /* Read */
            {
                auto it = readStreamMap.find(streamName);
                if (it == readStreamMap.end())
                {
                    std::cout << "    Open Input Stream " << streamName << "... " << std::endl;
                    if (!settings.outputPath.empty())
                    {
                        std::string outputPath = settings.outputPath;
                        if (settings.outputPath.back() != '/')
                        {
                            outputPath += '/';
                        }

                        streamName = outputPath + streamName;
                    }
                    std::shared_ptr<Stream> reader =
                        openStream(streamName, io, adios2::Mode::Read, settings.iolib,
                                   settings.appComm, settings.ioTimer, settings.appId);
                    readStreamMap[st.first] = reader;
                }
            }
        }

        /* Execute commands */
        bool exitLoop = false;
        size_t step = 1;
        while (!exitLoop)
        {
            if (!settings.myRank)
            {
                std::cout << "App " + std::to_string(settings.appId) + " Step " << step << ": "
                          << std::endl;
            }
            for (const auto &cmd : cfg.commands)
            {
                if (!cmd->conditionalStream.empty() &&
                    cfg.condMap.at(cmd->conditionalStream) != adios2::StepStatus::OK)
                {
                    if (!settings.myRank && settings.verbose)
                    {
                        std::cout << "    Skip command because of status "
                                     "of stream "
                                  << cmd->conditionalStream << std::endl;
                    }
                    continue;
                }

                switch (cmd->op)
                {
                case Operation::Sleep: {
                    std::chrono::high_resolution_clock::time_point start =
                        std::chrono::high_resolution_clock::now();
                    adios.EnterComputationBlock();
                    auto cmdS = dynamic_cast<const CommandSleep *>(cmd.get());
                    if (!settings.myRank && settings.verbose)
                    {
                        double t = static_cast<double>(cmdS->sleepTime_us) / 1000000.0;
                        std::cout << "    Sleep for " << t << "  seconds ";
                    }
                    std::this_thread::sleep_for(std::chrono::microseconds(cmdS->sleepTime_us));
                    adios.ExitComputationBlock();
                    if (!settings.myRank && settings.verbose)
                    {
                        std::chrono::high_resolution_clock::time_point end =
                            std::chrono::high_resolution_clock::now();
                        double t = static_cast<double>((end - start).count()) / 1000000000.0;
                        std::cout << " -> Slept for " << t << "  seconds " << std::endl;
                    }
                    break;
                }
                case Operation::Busy: {
                    auto cmdS = dynamic_cast<const CommandBusy *>(cmd.get());
                    std::chrono::high_resolution_clock::time_point start =
                        std::chrono::high_resolution_clock::now();
                    // auto sleeptime =
                    //     std::chrono::microseconds(cmdS->busyTime_us);
                    if (!settings.myRank && settings.verbose)
                    {
                        double t = static_cast<double>(cmdS->busyTime_us) / 1000000.0;
                        std::cout << "    Busy for " << cmdS->cycles << " cycles with " << t
                                  << " seconds work in each cycle";
                    }
                    const size_t N = 1048576;
                    double *f = (double *)calloc(N, sizeof(double));
                    double *g = (double *)malloc(N * sizeof(double));
                    for (size_t c = 0; c < cmdS->cycles; ++c)
                    {
                        auto end = std::chrono::high_resolution_clock::now() +
                                   std::chrono::microseconds(cmdS->busyTime_us);
                        if (cmdS->busyTime_us > 0.1)
                        {
                            adios.EnterComputationBlock();
                        }
                        while (std::chrono::high_resolution_clock::now() < end)
                        {
                            for (size_t i = 0; i < N; ++i)
                            {
                                f[i] = f[i] * 2.0 + 0.000001;
                            }
                        }
                        if (cmdS->busyTime_us > 0.1)
                        {
                            adios.ExitComputationBlock();
                        }
                        MPI_Allreduce(f, g, N, MPI_DOUBLE, MPI_SUM, settings.appComm);
                    }
                    std::chrono::high_resolution_clock::time_point end =
                        std::chrono::high_resolution_clock::now();
                    actualBusyTime_usec += (end - start).count() / 1000;
                    if (!settings.myRank && settings.verbose)
                    {
                        double t = static_cast<double>((end - start).count()) / 1000000000.0;
                        std::cout << " -> Was busy for " << t << "  seconds " << std::endl;
                    }
                    break;
                }
                case Operation::Write: {
                    auto cmdW = dynamic_cast<CommandWrite *>(cmd.get());
                    auto stream = writeStreamMap[cmdW->streamName];
                    // auto io = ioMap[cmdW->groupName];
                    stream->Write(cmdW, cfg, settings, step);
                    break;
                }
                case Operation::Read: {
                    auto cmdR = dynamic_cast<CommandRead *>(cmd.get());
                    auto statusIt = cfg.condMap.find(cmdR->streamName);
                    if (statusIt->second == adios2::StepStatus::OK ||
                        statusIt->second == adios2::StepStatus::NotReady)
                    {
                        auto stream = readStreamMap[cmdR->streamName];
                        // auto io = ioMap[cmdR->groupName];
                        adios2::StepStatus status = stream->Read(cmdR, cfg, settings, step);
                        statusIt->second = status;
                        switch (status)
                        {
                        case adios2::StepStatus::OK:
                            break;
                        case adios2::StepStatus::NotReady:
                            if (!settings.myRank && settings.verbose)
                            {
                                std::cout << "    Nonblocking read status: "
                                             "Not Ready "
                                          << std::endl;
                            }
                            break;
                        case adios2::StepStatus::EndOfStream:
                        case adios2::StepStatus::OtherError:
                            cfg.stepOverStreams.erase(cmdR->streamName);
                            if (!settings.myRank && settings.verbose)
                            {
                                std::cout << "    Nonblocking read status: "
                                             "Terminated "
                                          << std::endl;
                            }
                            break;
                        }
                    }
                    break;
                }
                }
                if (!settings.myRank && settings.verbose)
                {
                    std::cout << std::endl;
                }
            }
            if (!cfg.stepOverStreams.size() && step >= cfg.nSteps)
            {
                exitLoop = true;
            }
            ++step;
        }

        /* Close all streams in order of opening */
        for (const auto &st : streamsInOrder)
        {
            const std::string &streamName = st.first;
            const bool isWrite = (st.second == Operation::Write);
            if (isWrite)
            {
                auto writerIt = writeStreamMap.find(streamName);
                if (writerIt != writeStreamMap.end())
                {
                    auto writer = writeStreamMap[streamName];
                    writerIt->second->Close();
                    writeStreamMap.erase(writerIt);
                }
            }
            else /* Read */
            {
                auto readerIt = readStreamMap.find(streamName);
                if (readerIt != readStreamMap.end())
                {
                    auto reader = readStreamMap[streamName];
                    readerIt->second->Close();
                    readStreamMap.erase(readerIt);
                }
            }
        }
    }
    catch (std::exception &e) // if some unknown error occurs
    {
        if (!settings.myRank)
        {
            std::cout << "ADIOS " << e.what() << std::endl;
        }

        /* Yell and quit */
        MPI_Abort(settings.appComm, -1);
    }

    MPI_Barrier(settings.appComm);
    timeEnd = MPI_Wtime();
    if (!settings.myRank)
    {
        if (actualBusyTime_usec > 0)
        {
            std::cout << "  Total Busy time on Rank 0 was "
                      << (double)actualBusyTime_usec / 1000000.0 << " seconds " << std::endl;
        }
        std::cout << "ADIOS IOTEST App " << settings.appId << " total time " << timeEnd - timeStart
                  << " seconds " << std::endl;
    }

    MPI_Finalize();
    return 0;
}
