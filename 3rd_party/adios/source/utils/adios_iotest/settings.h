/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * settings.h
 *
 *  Created on: Oct 2018
 *      Author: Norbert Podhorszki
 */

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "adios2/common/ADIOSConfig.h"

#include <fstream>
#include <string>
#include <vector>

#include <mpi.h>

enum class IOLib
{
    ADIOS
#ifdef ADIOS2_HAVE_HDF5_PARALLEL
    ,
    HDF5
#endif
};

class Settings
{

public:
    /* user arguments */
    std::string configFileName;
    std::string adiosConfigFileName;
    std::string outputPath;
    unsigned int verbose = 0;
    size_t appId = 0;
    bool isStrongScaling = true; // strong or weak scaling
    bool ioTimer = false;        // used to measure io time
    bool fixedPattern = false;   // should Lock definitions?
    bool isRatioDecomp = false;
    bool multithreadedMPI = false; // turn on MT-enabled MPI
    IOLib iolib = IOLib::ADIOS;
    //   process decomposition
    std::vector<size_t> processDecomp = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    /* public variables */
    MPI_Comm appComm = MPI_COMM_WORLD; // will change to split communicator
    size_t myRank = 0;
    size_t nProc = 1;
    std::ifstream configFile;
    size_t nDecomp = 0;

    Settings() = default;
    ~Settings() = default;
    int processArguments(int argc, char *argv[]);
    int initDecomp(MPI_Comm worldComm);
    int extraArgumentChecks();
    size_t stringToNumber(const std::string &varName, const char *arg) const;
    int parseCSDecomp(const char *arg);
    int rescaleDecomp();
    size_t ndigits(size_t n) const;

private:
    void displayHelp();
    int processArgs(int argc, char *argv[]);
};

#endif /* SETTINGS_H_ */
