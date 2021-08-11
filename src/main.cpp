/*---------------------------------------------------------------------------*\
   Copyright (c) 2019-2021, UCHICAGO ARGONNE, LLC.

   The UChicago Argonne, LLC as Operator of Argonne National
   Laboratory holds copyright in the Software. The copyright holder
   reserves all rights except those expressly granted to licensees,
   and U.S. Government license rights.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the disclaimer below.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the disclaimer (as noted below)
   in the documentation and/or other materials provided with the
   distribution.

   3. Neither the name of ANL nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
   UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
   ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Additional BSD Notice
   ---------------------
   1. This notice is required to be provided under our contract with
   the U.S. Department of Energy (DOE). This work was produced at
   Argonne National Laboratory under Contract
   No. DE-AC02-06CH11357 with the DOE.

   2. Neither the United States Government nor UCHICAGO ARGONNE,
   LLC nor any of their employees, makes any warranty,
   express or implied, or assumes any liability or responsibility for the
   accuracy, completeness, or usefulness of any information, apparatus,
   product, or process disclosed, or represents that its use would not
   infringe privately-owned rights.

   3. Also, reference herein to any specific commercial products, process,
   or services by trade name, trademark, manufacturer or otherwise does
   not necessarily constitute or imply its endorsement, recommendation,
   or favoring by the United States Government or UCHICAGO ARGONNE LLC.
   The views and opinions of authors expressed
   herein do not necessarily state or reflect those of the United States
   Government or UCHICAGO ARGONNE, LLC, and shall
   not be used for advertising or product endorsement purposes.

\*---------------------------------------------------------------------------*/

#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <cstring>
#include <getopt.h>
#include <cfenv>
#include <limits>
#include <math.h>
#include <unistd.h>
#include <vector>

#include "nekrs.hpp"

#define DEBUG

static MPI_Comm globalComm;
static MPI_Comm comm;

struct cmdOptions
{
  int buildOnly = 0;
  int ciMode = 0;
  int debug = 0;
  int sizeTarget = 0;
  std::string setupFile;
  std::string deviceID;
  std::string backend;
  bool redirectOutput;
  int neknekSessions = 1;
  std::vector<std::string> neknekSetupFiles;
  std::vector<int> neknekProcs;
  bool neknekConnected = true;
};

static cmdOptions* processCmdLineOptions(int argc, char** argv);

int main(int argc, char** argv)
{
  {
    int request = MPI_THREAD_SINGLE;
    const char* env_val = std::getenv ("NEKRS_MPI_THREAD_MULTIPLE");
    if (env_val)
      if (std::stoi(env_val)) request = MPI_THREAD_MULTIPLE;

    int provided;
    int retval =  MPI_Init_thread(&argc, &argv, request, &provided);
    if (retval != MPI_SUCCESS) {
      std::cout << "FATAL ERROR: Cannot initialize MPI!" << "\n";
      exit(EXIT_FAILURE);
    }
  }

  MPI_Comm_dup(MPI_COMM_WORLD, &globalComm);
  cmdOptions* cmdOpt = processCmdLineOptions(argc, argv);

  neknek_t *neknek = new neknek_t();
  neknek->nsessions = cmdOpt->neknekSessions;
  neknek->globalComm = globalComm;
  if (neknek->nsessions != 1) {
    neknek->connected = cmdOpt->neknekConnected;
    int grank, sessionID = -1, nextRoot = 0;
    MPI_Comm_rank(globalComm, &grank);
    for(int i = 0; i < neknek->nsessions; ++ i) {
      nextRoot += cmdOpt->neknekProcs[i];
      if (grank < nextRoot) {
        sessionID = i;
        break;
      }
    }
    MPI_Comm_split(globalComm, sessionID, 0, &comm);
    cmdOpt->setupFile = cmdOpt->neknekSetupFiles[sessionID];
    neknek->sessionID = sessionID;
  } else {
    neknek->connected = false;
    comm = globalComm;
  }
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (rank == 0 && cmdOpt->redirectOutput) {
    std::string logfile = cmdOpt->setupFile + ".log." + std::to_string(size);
    printf("redirecting stdout to %s\n",logfile.c_str());
    freopen(logfile.c_str(), "w+", stdout);
    setvbuf(stdout, NULL, _IONBF, 0);
  }

  if (cmdOpt->debug) {
    if (rank == 0) {
      std::cout << "Attach debugger, then press enter to continue\n";
    }
    printf("\tRank\tpid\n");
    for(int currRank = 0; currRank < size; ++currRank)
    {
      if(rank == currRank){
        printf("\t%d\t%d\n", rank, ::getpid());
      }
    }
    if (rank == 0) {
      std::cin.get();
    }
    MPI_Barrier(comm);
  }
  if (cmdOpt->debug) feraiseexcept(FE_ALL_EXCEPT);

  MPI_Barrier(comm);
  const double time0 = MPI_Wtime();

  std::string cacheDir;
  nekrs::setup(comm, cmdOpt->buildOnly, cmdOpt->sizeTarget,
               cmdOpt->ciMode, cacheDir, cmdOpt->setupFile,
               cmdOpt->backend, cmdOpt->deviceID,
               neknek);

  if (cmdOpt->buildOnly) {
    nekrs::finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  MPI_Barrier(comm);
  double elapsedTime = (MPI_Wtime() - time0);

  const int runTimeStatFreq = 500;
  int tStep = 0;
  double time = nekrs::startTime();
  int lastStep = nekrs::lastStep(time, tStep, elapsedTime);

  nekrs::udfExecuteStep(time, tStep, /* outputStep */ 0);

  if (rank == 0 && !lastStep) {
    if (nekrs::endTime() > nekrs::startTime())
      std::cout << "\ntimestepping to time " << nekrs::endTime() << " ...\n";
    else
      std::cout << "\ntimestepping for " << nekrs::numSteps() << " steps ...\n";
  }
  MPI_Pcontrol(1);
  while (!lastStep) {
    MPI_Barrier(comm);
    const double timeStart = MPI_Wtime();

    ++tStep;
    lastStep = nekrs::lastStep(time, tStep, elapsedTime);

    double dt;
    if (lastStep && nekrs::endTime() > 0)
      dt = nekrs::endTime() - time;
    else
      dt = nekrs::dt();

    int outputStep = nekrs::outputStep(time+dt, tStep);
    if (nekrs::writeInterval() == 0) outputStep = 0;
    if (lastStep) outputStep = 1;
    if (nekrs::writeInterval() < 0) outputStep = 0;
    nekrs::outputStep(outputStep);

    nekrs::runStep(time, dt, tStep);
    time += dt;

    if (outputStep) nekrs::outfld(time);

    if (tStep%runTimeStatFreq == 0 || lastStep) nekrs::printRuntimeStatistics();

    MPI_Barrier(comm);
    elapsedTime += (MPI_Wtime() - timeStart);
  }
  MPI_Pcontrol(0);

  if (rank == 0) {
    std::cout << "elapsedTime: " << elapsedTime << " s\n";
    std::cout << "End\n";
  }
  fflush(stdout);

  nekrs::finalize();
  MPI_Finalize();
  return EXIT_SUCCESS;
}

static cmdOptions* processCmdLineOptions(int argc, char** argv)
{
  int rank,size;
  MPI_Comm_rank(globalComm, &rank);
  MPI_Comm_size(globalComm, &size);

  cmdOptions* cmdOpt = new cmdOptions();

  int err = 0;

  if (rank == 0) {
    while(1) {
      static struct option long_options[] =
      {
        {"setup", required_argument, 0, 's'},
        {"cimode", required_argument, 0, 'c'},
	{"build-only", required_argument, 0, 'b'},
        {"debug", no_argument, 0, 'd'},
        {"backend", required_argument, 0, 't'},
        {"device-id", required_argument, 0, 'i'},
        {"neknek", required_argument, 0, 'n'},
        {"neknek-procs", required_argument, 0, 'p'},
        {"neknek-unconnected", no_argument, 0, 'u'},
        {0, 0, 0, 0}
      };
      int option_index = 0;
      int c = getopt_long (argc, argv, "s:", long_options, &option_index);

      if (c == -1)
        break;

      switch(c) {
      case 's':
        cmdOpt->setupFile.assign(optarg);
        break;
      case 'b':
        cmdOpt->buildOnly = 1;
	cmdOpt->sizeTarget = atoi(optarg);
        break;
      case 'c':
        cmdOpt->ciMode = atoi(optarg);
        if (cmdOpt->ciMode < 1) {
          std::cout << "ERROR: ci test id has to be >0!\n";
          err = 1;
        }
        break;
      case 'd':
        cmdOpt->debug = 1;
        break;
      case 'i':
        cmdOpt->deviceID.assign(optarg);
        break;
      case 't':
        cmdOpt->backend.assign(optarg);
        break;
      case 'n': {
        cmdOpt->neknekSessions = atoi(optarg);

        cmdOpt->neknekSetupFiles.resize(cmdOpt->neknekSessions, cmdOpt->setupFile);
        std::string caseName;
        std::istringstream stream (cmdOpt->setupFile);
        int i = 0;
        for(int i = 0; i < cmdOpt->neknekSessions; ++i) {
          std::getline(stream, caseName, ',');
          cmdOpt->neknekSetupFiles[i] = caseName;
        }

        cmdOpt->redirectOutput = true;
        break;
      }
      case 'p': {
        std::string proc;
        std::string input;
        input.assign(optarg);
        std::istringstream stream (input);
        int i = 0;
        while(std::getline(stream, proc, ',')) {
          cmdOpt->neknekProcs.push_back(atoi(proc.c_str()));
        }
        break;
      }
      case 'u':
        cmdOpt->neknekConnected = false;
        break;
      default:
        err = 1;
      }
    }

    if (cmdOpt->neknekSessions > 1
        && cmdOpt->neknekSessions != cmdOpt->neknekProcs.size()) {
      err = 1;
    }
  }

  char buf[FILENAME_MAX];
  strcpy(buf, cmdOpt->setupFile.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, globalComm);
  cmdOpt->setupFile.assign(buf);
  strcpy(buf, cmdOpt->deviceID.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, globalComm);
  cmdOpt->deviceID.assign(buf);
  strcpy(buf, cmdOpt->backend.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, globalComm);
  cmdOpt->backend.assign(buf);
  MPI_Bcast(&cmdOpt->buildOnly, sizeof(cmdOpt->buildOnly), MPI_BYTE, 0, globalComm);
  MPI_Bcast(&cmdOpt->sizeTarget, sizeof(cmdOpt->sizeTarget), MPI_BYTE, 0, globalComm);
  MPI_Bcast(&cmdOpt->ciMode, sizeof(cmdOpt->ciMode), MPI_BYTE, 0, globalComm);
  MPI_Bcast(&cmdOpt->debug, sizeof(cmdOpt->debug), MPI_BYTE, 0, globalComm);
  MPI_Bcast(&cmdOpt->redirectOutput, sizeof(cmdOpt->redirectOutput), MPI_BYTE, 0, globalComm);

  MPI_Bcast(&cmdOpt->neknekSessions, sizeof(cmdOpt->neknekSessions), MPI_BYTE, 0, globalComm);
  if(cmdOpt->neknekSessions > 1) {
    if (rank != 0) {
      cmdOpt->neknekProcs.resize(cmdOpt->neknekSessions);
      cmdOpt->neknekSetupFiles.resize(cmdOpt->neknekSessions);
    }
    MPI_Bcast(cmdOpt->neknekProcs.data(), cmdOpt->neknekSessions, MPI_INT, 0, globalComm);
    for (int i = 0; i < cmdOpt->neknekSessions; ++i) {
      strcpy(buf, cmdOpt->neknekSetupFiles[i].c_str());
      MPI_Bcast(buf, sizeof(buf), MPI_CHAR, 0, globalComm);
      cmdOpt->neknekSetupFiles[i].assign(buf);
    }
  }

  if(cmdOpt->setupFile.empty()){
    err++;
  } else {
    std::string casepath, casename;
    size_t last_slash = cmdOpt->setupFile.rfind('/') + 1;
    casepath = cmdOpt->setupFile.substr(0,last_slash);
    casename = cmdOpt->setupFile.substr(last_slash, cmdOpt->setupFile.length() - last_slash);
    if(casepath.length() > 0) chdir(casepath.c_str());
    cmdOpt->setupFile.assign(casename);
  }

  MPI_Bcast(&err, sizeof(err), MPI_BYTE, 0, globalComm);
  if (err) {
    if (rank == 0)
      std::cout << "usage: ./nekrs --setup <case name> "
                << "[ --build-only <#procs> ] [ --cimode <id> ] [ --debug ] "
                << "[ --backend <CPU|CUDA|HIP|OPENCL> ] [ --device-id <id|LOCAL-RANK> ] "
                << "[ --neknek <# sessions> --neknek-procs <#procs list> [ --neknek-unconnected ] ]"
                << "\n";
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  return cmdOpt;
}
