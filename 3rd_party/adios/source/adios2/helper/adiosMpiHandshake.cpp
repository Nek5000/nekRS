/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosMpiHandshake.cpp
 *
 *  Created on: Mar 1, 2020
 *      Author: Jason Wang
 */

#include "adiosMpiHandshake.h"
#include "adiosLog.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <unordered_set>

namespace adios2
{
namespace helper
{

void HandshakeComm(const std::string &filename, const char mode, const int timeoutSeconds,
                   MPI_Comm localComm, MPI_Group &streamGroup, MPI_Group &writerGroup,
                   MPI_Group &readerGroup, MPI_Comm &streamComm, MPI_Comm &writerComm,
                   MPI_Comm &readerComm, int verbosity)
{
    auto appRankMaps = HandshakeRank(filename, mode, timeoutSeconds, localComm, verbosity);
    MPI_Group worldGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
    MPI_Group_incl(worldGroup, static_cast<int>(appRankMaps[0].size()), appRankMaps[0].data(),
                   &streamGroup);
    MPI_Group_incl(worldGroup, static_cast<int>(appRankMaps[1].size()), appRankMaps[1].data(),
                   &writerGroup);
    MPI_Group_incl(worldGroup, static_cast<int>(appRankMaps[2].size()), appRankMaps[2].data(),
                   &readerGroup);
#ifdef _WIN32
    MPI_Comm_create(MPI_COMM_WORLD, streamGroup, &streamComm);
    MPI_Comm_create(MPI_COMM_WORLD, writerGroup, &writerComm);
    MPI_Comm_create(MPI_COMM_WORLD, readerGroup, &readerComm);
#else
    MPI_Comm_create_group(MPI_COMM_WORLD, streamGroup, 0, &streamComm);
    MPI_Comm_create_group(MPI_COMM_WORLD, writerGroup, 0, &writerComm);
    MPI_Comm_create_group(MPI_COMM_WORLD, readerGroup, 0, &readerComm);
#endif
}

const std::vector<std::vector<int>> HandshakeRank(const std::string &filename, const char mode,
                                                  const int timeoutSeconds, MPI_Comm localComm,
                                                  int verbosity)
{
    std::vector<std::vector<int>> ret(3);

    int localRank;
    int localSize;
    int worldRank;
    int worldSize;

    MPI_Comm_rank(localComm, &localRank);
    MPI_Comm_size(localComm, &localSize);

    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    std::vector<int> allLocalRanks(localSize);

    MPI_Gather(&worldRank, 1, MPI_INT, allLocalRanks.data(), 1, MPI_INT, 0, localComm);

    if (localRank == 0)
    {
        std::ofstream fs;
        fs.open(filename + "." + mode);
        for (auto rank : allLocalRanks)
        {
            fs << rank << std::endl;
        }
        fs.close();

        if (mode == 'r')
        {

            for (auto i : allLocalRanks)
            {
                ret[0].push_back(i);
                ret[2].push_back(i);
            }

            std::ofstream fsc;
            fsc.open(filename + ".r.c");
            fsc << "completed";
            fsc.close();

            auto startTime = std::chrono::system_clock::now();
            while (true)
            {
                std::ifstream fs;
                try
                {
                    auto nowTime = std::chrono::system_clock::now();
                    auto duration =
                        std::chrono::duration_cast<std::chrono::seconds>(nowTime - startTime);
                    if (duration.count() > timeoutSeconds)
                    {
                        helper::Throw<std::runtime_error>(
                            "Helper", "adiosMpiHandshake", "HandshakeRank",
                            "Mpi handshake timeout for Stream " + filename);
                    }

                    fs.open(filename + ".w.c");
                    std::string line;
                    std::getline(fs, line);
                    if (line != "completed")
                    {
                        continue;
                    }
                    fs.close();
                    remove((filename + ".w.c\0").c_str());
                    break;
                }
                catch (...)
                {
                    continue;
                }
            }

            std::ifstream fs;
            fs.open(filename + ".w");
            for (std::string line; std::getline(fs, line);)
            {
                ret[0].push_back(std::stoi(line));
                ret[1].push_back(std::stoi(line));
            }
            fs.close();
            remove((filename + ".w\0").c_str());
        }
        else if (mode == 'w')
        {
            for (auto i : allLocalRanks)
            {
                ret[0].push_back(i);
                ret[1].push_back(i);
            }

            std::ofstream fsc;
            fsc.open(filename + ".w.c");
            fsc << "completed";
            fsc.close();

            while (true)
            {
                std::ifstream fs;
                try
                {
                    fs.open(filename + ".r.c");
                    std::string line;
                    std::getline(fs, line);
                    if (line != "completed")
                    {
                        continue;
                    }
                    fs.close();
                    remove((filename + ".r.c\0").c_str());
                    break;
                }
                catch (...)
                {
                    continue;
                }
            }

            std::ifstream fs;
            fs.open(filename + ".r");
            for (std::string line; std::getline(fs, line);)
            {
                ret[0].push_back(std::stoi(line));
                ret[2].push_back(std::stoi(line));
            }
            fs.close();
            remove((filename + ".r\0").c_str());
        }
    }

    int dims[3];

    if (localRank == 0)
    {
        for (int i = 0; i < 3; ++i)
        {
            dims[i] = static_cast<int>(ret[i].size());
            std::sort(ret[i].begin(), ret[i].end());
        }
    }

    MPI_Bcast(dims, 3, MPI_INT, 0, localComm);

    if (localRank != 0)
    {
        for (int i = 0; i < 3; ++i)
        {
            ret[i].resize(dims[i]);
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        MPI_Bcast(ret[i].data(), static_cast<int>(ret[i].size()), MPI_INT, 0, localComm);
    }

    if (verbosity >= 5)
    {
        std::stringstream output;
        output << "World Rank " << worldRank << ": " << std::endl;
        int s = 0;
        for (const auto &i : ret)
        {
            output << "    " << s << ": ";
            for (const auto &j : i)
            {
                output << j << ", ";
            }
            output << std::endl;
            ++s;
        }
        std::cout << output.str();
    }

    return ret;
}

} // end namespace helper
} // end namespace adios2
