/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 *  Created on: Jan 2018
 *      Author: Norbert Podhorszki
 */

#include "SkeletonArgs.h"

#include <errno.h>

#include <cstdlib>
#include <iostream>

#include <stdexcept>

static void printUsage(bool isWriter)
{
    if (isWriter)
    {
        std::cout << "Usage: Usage: helloSkeletonWriter  config  N  M   nx  ny   "
                     "steps  sleeptime\n";
    }
    else
    {
        std::cout << "Usage: Usage: helloSkeletonReader  config  N  M \n  ";
    }
    std::cout << "  config:    XML config file to use\n"
              << "  N:         number of processes in X dimension\n"
              << "  M:         number of processes in Y dimension\n";
    if (isWriter)
    {
        std::cout << "  nx:        local array size in X dimension per processor\n"
                  << "  ny:        local array size in Y dimension per processor\n"
                  << "  steps:     the total number of steps to output\n"
                  << "  sleeptime: wait this many milliseconds between output "
                     "steps\n\n";
    }
}

static unsigned int convertToUint(std::string varName, char *arg)
{
    char *end;
    const auto retval = static_cast<unsigned int>(std::strtoul(arg, &end, 10));
    if (end[0] || errno == ERANGE)
    {
        throw std::invalid_argument("Invalid value given for " + varName + ": " + std::string(arg));
    }
    return retval;
}

SkeletonArgs::SkeletonArgs(bool isWriter, int argc, char *argv[], int rank, int nproc) : rank{rank}
{
    npx = npy = ndx = ndy = steps = sleeptime = 0;
    gndx = gndy = posx = posy = offsx = offsy = 0;
    int expargs = (isWriter ? 8 : 4);
    this->nproc = static_cast<unsigned int>(nproc);

    try
    {
        if (argc < expargs)
        {
            throw std::invalid_argument("Not enough arguments");
        }

        configfile = argv[1];
        npx = convertToUint("N", argv[2]);
        npy = convertToUint("M", argv[3]);

        if (isWriter)
        {
            ndx = convertToUint("nx", argv[4]);
            ndy = convertToUint("ny", argv[5]);
            steps = convertToUint("steps", argv[6]);
            sleeptime = convertToUint("sleeptime", argv[7]);

            // calculate global array size and the local offsets in that global
            // space
            gndx = npx * ndx;
            gndy = npy * ndy;
            posx = rank % npx;
            posy = rank / npx;
            offsx = posx * ndx;
            offsy = posy * ndy;
        }

        if (npx * npy != static_cast<size_t>(nproc))
        {
            throw std::invalid_argument("N*M must equal the number of processes");
        }
    }
    catch (std::invalid_argument &e)
    {
        printUsage(isWriter);
        throw e;
    }
}

void SkeletonArgs::DecomposeArray(size_t NX, size_t NY)
{
    gndx = static_cast<unsigned int>(NX);
    gndy = static_cast<unsigned int>(NY);
    posx = rank % npx;
    posy = rank / npx;

    // 2D decomposition of global array reading
    ndx = gndx / npx;
    ndy = gndy / npy;
    offsx = ndx * posx;
    offsy = ndy * posy;
    if (posx == npx - 1)
    {
        // right-most processes need to read all the rest of rows
        ndx = gndx - ndx * (npx - 1);
    }

    if (posy == npy - 1)
    {
        // bottom processes need to read all the rest of columns
        ndy = gndy - ndy * (npy - 1);
    }

    std::cout << "rank " << rank << " reads 2D slice " << ndx << " x " << ndy << " from offset ("
              << offsx << "," << offsy << ")" << std::endl;
}
