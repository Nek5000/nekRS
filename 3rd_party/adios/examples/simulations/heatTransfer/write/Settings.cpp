/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Settings.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "Settings.h"

#include <errno.h>

#include <cstdlib>

#include <stdexcept>

static unsigned int convertToUint(std::string varName, char *arg)
{
    char *end;
    long retval = std::strtol(arg, &end, 10);
    if (end[0] || errno == ERANGE)
    {
        throw std::invalid_argument("Invalid value given for " + varName + ": " + std::string(arg));
    }
    if (retval < 0)
    {
        throw std::invalid_argument("Negative value given for " + varName + ": " +
                                    std::string(arg));
    }
    return static_cast<unsigned int>(retval);
}

Settings::Settings(int argc, char *argv[], int rank, int nproc) : rank{rank}
{
    if (argc < 9)
    {
        throw std::invalid_argument("Not enough arguments");
    }
    this->nproc = (unsigned int)nproc;

    configfile = argv[1];
    outputfile = argv[2];
    npx = convertToUint("N", argv[3]);
    npy = convertToUint("M", argv[4]);
    ndx = convertToUint("nx", argv[5]);
    ndy = convertToUint("ny", argv[6]);
    steps = convertToUint("steps", argv[7]);
    iterations = convertToUint("iterations", argv[8]);

    if (npx * npy != this->nproc)
    {
        throw std::invalid_argument("N*M must equal the number of processes");
    }

    // calculate global array size and the local offsets in that global space
    gndx = npx * ndx;
    gndy = npy * ndy;
    posx = rank % npx;
    posy = rank / npx;
    offsx = posx * ndx;
    offsy = posy * ndy;

    // determine neighbors
    if (posx == 0)
        rank_up = -1;
    else
        rank_up = rank - 1;

    if (posx == npx - 1)
        rank_down = -1;
    else
        rank_down = rank + 1;

    if (posy == 0)
        rank_left = -1;
    else
        rank_left = rank - npx;

    if (posy == npy - 1)
        rank_right = -1;
    else
        rank_right = rank + npx;
}
