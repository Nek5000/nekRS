/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Settings.cpp
 *
 *  Created on: Dec 2017
 *      Author: Norbert Podhorszki
 */

#include "ReadSettings.h"

#include <cstdlib>
#include <errno.h>
#include <iomanip>
#include <iostream>
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

ReadSettings::ReadSettings(int argc, char *argv[], int rank, int nproc) : rank{rank}
{
    if (argc < 6)
    {
        throw std::invalid_argument("Not enough arguments");
    }
    this->nproc = (unsigned int)nproc;

    configfile = argv[1];
    inputfile = argv[2];
    outputfile = argv[3];
    npx = convertToUint("N", argv[4]);
    npy = convertToUint("M", argv[5]);

    if (npx * npy != static_cast<unsigned int>(this->nproc))
    {
        throw std::invalid_argument("N*M must equal the number of processes");
    }
    posx = rank % npx;
    posy = rank / npx;
}

void ReadSettings::DecomposeArray(int gndx, int gndy)
{
    // 2D decomposition of global array reading
    size_t ndx = gndx / npx;
    size_t ndy = gndy / npy;
    size_t offsx = ndx * posx;
    size_t offsy = ndy * posy;
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
    readsize.push_back(ndx);
    readsize.push_back(ndy);
    offset.push_back(offsx);
    offset.push_back(offsy);

    std::cout << "rank " << rank << " reads 2D slice " << ndx << " x " << ndy << " from offset ("
              << offsx << "," << offsy << ")" << std::endl;
}
