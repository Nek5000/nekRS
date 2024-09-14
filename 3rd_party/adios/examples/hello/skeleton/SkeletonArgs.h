/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 *  Created on: Jan 2018
 *      Author: Norbert Podhorszki
 */

#ifndef HELLOSKELETONARGS_H_
#define HELLOSKELETONARGS_H_

#include <string>

class SkeletonArgs
{

public:
    const std::string streamname = "skeleton_stream";

    // user arguments
    std::string configfile;
    unsigned int npx;       // Number of processes in X (slow) dimension
    unsigned int npy;       // Number of processes in Y (fast) dimension
    unsigned int ndx;       // Local array size in X dimension per process
    unsigned int ndy;       // Local array size in y dimension per process
    unsigned int steps;     // Number of output steps
    unsigned int sleeptime; // Time to wait between steps (in millisec)

    // calculated values from those arguments and number of processes
    unsigned int gndx; // Global array size in slow dimension
    unsigned int gndy; // Global array size in fast dimension
    // X dim positions: rank 0, npx, 2npx... are in the same X position
    // Y dim positions: npx number of consecutive processes belong to one row
    // (npx
    // columns)
    unsigned int posx;  // Position of this process in X dimension
    unsigned int posy;  // Position of this process in Y dimension
    unsigned int offsx; // Offset of local array in X dimension on this process
    unsigned int offsy; // Offset of local array in Y dimension on this process

    int rank;           // MPI rank
    unsigned int nproc; // number of processors

    SkeletonArgs(bool isWriter, int argc, char *argv[], int rank, int nproc);

    void DecomposeArray(size_t NX, size_t NY);
};

#endif /* HELLOSKELETONARGS_H_ */
