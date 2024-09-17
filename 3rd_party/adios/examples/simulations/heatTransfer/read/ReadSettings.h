/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Settings.h
 *
 *  Created on: Dec 2017
 *      Author: Norbert Podhorszki
 */

#ifndef READSETTINGS_H_
#define READSETTINGS_H_

#include <string>
#include <vector>

class ReadSettings
{

public:
    // user arguments
    std::string configfile;
    std::string inputfile;
    std::string outputfile;
    unsigned int npx; // Number of processes in X (slow) dimension
    unsigned int npy; // Number of processes in Y (fast) dimension

    int rank;
    int nproc;

    // Calculated in constructor
    unsigned int posx; // Position of this process in X dimension
    unsigned int posy; // Position of this process in Y dimension

    // Calculated in DecomposeArray
    std::vector<size_t> readsize; // Local array size in X-Y dimensions per process
    std::vector<size_t> offset;   // Offset of local array in X-Y dimensions on this process

    ReadSettings(int argc, char *argv[], int rank, int nproc);
    void DecomposeArray(int gndx, int gndy);
};

#endif /* READSETTINGS_H_ */
