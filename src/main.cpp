
/*---------------------------------------------------------------------------*\
   Copyright (c) 2019-2022, UCHICAGO ARGONNE, LLC.

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
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <cstring>
#include <getopt.h>
#include <cfenv>
#include <limits>
#include <math.h>
#include <unistd.h>
#include <libgen.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fcntl.h>
#include <chrono>
#include <csignal>
#include "stacktrace.hpp"
#include <filesystem>

#include "nekrs.hpp"

#ifdef ENABLE_SENSEI
#include "sensei/Bridge.h"
#endif

#define DEBUG

namespace {

int worldRank;

struct cmdOptions
{
  int buildOnly = 0;
  int ciMode = 0;
  int debug = 0;
  int sizeTarget = 0;
  std::string multiSessionFile;
  std::string setupFile;
  std::string deviceID;
  std::string backend;
  int nSessions = 1;
  int sessionID = 0;
  bool redirectOutput;
};

struct session_t {
  int size;
  std::string setupFile;
};

std::vector<std::string> serializeString(const std::string sin, char dlim)
{
  std::vector<std::string> slist;
  std::string s(sin);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  std::stringstream ss;
  ss.str(s);
  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, dlim);
    if(!substr.empty()) slist.push_back(substr);
  }
  return slist;
}

cmdOptions* processCmdLineOptions(int argc, char** argv, const MPI_Comm &comm)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  cmdOptions* cmdOpt = new cmdOptions();

  int err = 0;
  int printHelp = 0;
  std::string helpCat;

  if (rank == 0) {
    while(1) {
      static struct option long_options[] =
      {
        {"setup", required_argument, 0, 's'},
        {"cimode", required_argument, 0, 'c'},
        {"build-only", optional_argument, 0, 'b'},
        {"debug", no_argument, 0, 'd'},
        {"backend", required_argument, 0, 't'},
        {"device-id", required_argument, 0, 'i'},
        {"help", optional_argument, 0, 'h'},
        {0, 0, 0, 0}
      };
      int option_index = 0;
      int c = getopt_long (argc, argv, "", long_options, &option_index);

      if (c == -1)
        break;

      switch(c) {
      case 's':
        cmdOpt->setupFile.assign(optarg);
        if (cmdOpt->setupFile.find(".par") != std::string::npos)
          cmdOpt->setupFile.erase(cmdOpt->setupFile.find(".par"), std::string::npos);
        if (cmdOpt->setupFile.substr(cmdOpt->setupFile.find_last_of(".") + 1) == "sess") {
          cmdOpt->multiSessionFile = cmdOpt->setupFile;
          cmdOpt->setupFile.clear();
        }
        break;
      case 'b':
        cmdOpt->buildOnly = 1;
        cmdOpt->sizeTarget = size;
        if(!optarg && argv[optind] != NULL && argv[optind][0] != '-')
          cmdOpt->sizeTarget = std::stoi(argv[optind++]);
        break;
      case 'c':
        cmdOpt->ciMode = atoi(optarg);
        if (cmdOpt->ciMode < 0) {
          std::cout << "ERROR: ci test id has to be >= 0!\n";
          printHelp = 1;
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
      case 'h':
        printHelp++;
        if(!optarg && argv[optind] != NULL && argv[optind][0] != '-')
          helpCat.assign(argv[optind++]);
        break;
      default:
        err = 1;
      }
    }

    if(cmdOpt->setupFile.empty() && cmdOpt->multiSessionFile.empty()){
      int cnt = 0;
      for (auto &p : std::filesystem::directory_iterator{"."})
      {
        if (p.path().extension() == ".par") {
          cmdOpt->setupFile.assign(p.path().stem().string());
          cnt++; 
        }
      }
      if(cnt > 1) {
        std::cout << "Multiple .par files found!\n"; 
        err++; 
      }
    }
  }

  for(auto opt: {&cmdOpt->multiSessionFile, &cmdOpt->setupFile, &cmdOpt->deviceID, &cmdOpt->backend})
  {
    int bufSize = opt->size() + 1;
    MPI_Bcast(&bufSize, 1, MPI_INT, 0, comm);
    auto buf = (char*) std::calloc(bufSize, sizeof(char));
    if(rank == 0) std::strcpy(buf, opt->c_str());
    MPI_Bcast(buf, bufSize, MPI_BYTE, 0, comm);
    opt->assign(buf);
    free(buf);
  }

  MPI_Bcast(&cmdOpt->buildOnly, sizeof(cmdOpt->buildOnly), MPI_BYTE, 0, comm);
  MPI_Bcast(&cmdOpt->sizeTarget, sizeof(cmdOpt->sizeTarget), MPI_BYTE, 0, comm);
  MPI_Bcast(&cmdOpt->ciMode, sizeof(cmdOpt->ciMode), MPI_BYTE, 0, comm);
  MPI_Bcast(&cmdOpt->debug, sizeof(cmdOpt->debug), MPI_BYTE, 0, comm);

  if(cmdOpt->setupFile.empty() && cmdOpt->multiSessionFile.empty())
    printHelp++;

  MPI_Bcast(&printHelp, sizeof(printHelp), MPI_BYTE, 0, comm);
  MPI_Bcast(&err, sizeof(err), MPI_BYTE, 0, comm);
  if (err | printHelp) {
    if (rank == 0) {
      if (helpCat == "par") {
        std::string installDir;
        installDir.assign(getenv("NEKRS_HOME"));
        std::ifstream f(installDir + "/doc/parHelp.txt");
        if (f.is_open()) std::cout << f.rdbuf();
        f.close();
      } else {
        std::cout << "usage: ./nekrs [--help <par>] "
                  << "--setup <par|sess file> "
                  << "[ --build-only <#procs> ] [ --cimode <id> ] [ --debug ] "
                  << "[ --backend <CPU|CUDA|HIP|DPCPP|OPENCL> ] [ --device-id <id|LOCAL-RANK> ]"
                  << "\n";
      }
    }
    MPI_Finalize();
    exit((err) ? EXIT_FAILURE : EXIT_SUCCESS);
  }

  return cmdOpt;
}

MPI_Comm setupSession(cmdOptions *cmdOpt, const MPI_Comm &comm)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  MPI_Comm newComm = comm;

  if(cmdOpt->multiSessionFile.size()) {
    std::string multiSessionFileContent;

    if(rank == 0) {
      std::ifstream f(cmdOpt->multiSessionFile);
      if (!f) {
        std::cout << "FATAL ERROR: Cannot find sess file "
                  << cmdOpt->multiSessionFile << "!\n";
        fflush(stdout);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      std::ostringstream ss;
      if (f.is_open()) ss << f.rdbuf();
      f.close();
      multiSessionFileContent = ss.str();
    }
    int bufSize = multiSessionFileContent.size() + 1;
    MPI_Bcast(&bufSize, sizeof(bufSize), MPI_BYTE, 0, comm);
    char* buf = (char*) malloc(bufSize * sizeof(char));
    strcpy(buf, multiSessionFileContent.c_str());
    MPI_Bcast(buf, bufSize * sizeof(char), MPI_BYTE, 0, comm);
    multiSessionFileContent = std::string(buf);
    free(buf);

    auto list = serializeString(multiSessionFileContent, ';');
    auto sessionList = new session_t[list.size()];

    int nSessions = 0;
    int rankSum = 0;
    for(std::string s : list) {
      auto items = serializeString(s,':');
      if(items.size() != 2) {
        if(rank == 0) std::cout << "FATAL ERROR: invalid sess file entry!\n";
        fflush(stdout);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      sessionList[nSessions].setupFile = items[0];
      sessionList[nSessions].size = std::stoi(items[1]);
      rankSum += sessionList[nSessions].size;
      nSessions++;
    }

    int err = 0;
    if(rankSum != size) err = 1;
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, comm);
    if(err) {
      if(rank == 0) std::cout << "FATAL ERROR: size of sub-communicators does not match parent!\n";
      fflush(stdout);
      MPI_Abort(comm, EXIT_FAILURE);
    }

    int color = MPI_UNDEFINED;
    int rankOffsetSession = 0;
    for(int i = 0; i < nSessions; i++) {
      if(rank - rankOffsetSession < sessionList[i].size) {
        color = i;
        break;
      }
      rankOffsetSession += sessionList[i].size;
    }

    int rankGlobal, sizeGlobal;
    MPI_Comm_rank(comm, &rankGlobal);
    MPI_Comm_size(comm, &sizeGlobal);

    MPI_Comm_split(comm, color, rankGlobal, &newComm);

    MPI_Comm_rank(newComm, &rank);
    MPI_Comm_size(newComm, &size);

    cmdOpt->setupFile = sessionList[color].setupFile;
    cmdOpt->sizeTarget = size;
    cmdOpt->sessionID = color;
    cmdOpt->nSessions = nSessions;

    delete[] sessionList;

    if(cmdOpt->debug) {
      std::cout << "globalRank:" << rankGlobal << " localRank: " << rank << " pid: " << ::getpid()
                << " commSize: " << size << " setupFile:" << cmdOpt->setupFile << "\n";
    }
    fflush(stdout);
    MPI_Barrier(comm);

    if(rank == 0) {
      const std::string outputFile = cmdOpt->setupFile + ".log";
      std::cout << "redirecting output to " << outputFile << " ...\n";
      const int fd = open(outputFile.c_str(), O_WRONLY|O_CREAT|O_TRUNC, S_IWUSR|S_IRUSR);
      dup2(fd, fileno(stderr));
      dup2(fd, fileno(stdout));
    }
  }
  return newComm;
}

void signalHandler(int signum) 
{
   { // needs to be refactored as this is not async-signal-safe
     std::cerr << "generating stacktrace ...\n";

     const std::string fileName = "stacktrace." + std::to_string(worldRank);

     FILE *fp;
     fp = fopen (fileName.c_str(), "w");
     print_stacktrace(fp);
     fclose(fp);
   }
}

} // namespace


int main(int argc, char** argv)
{

  const auto timeStart = std::chrono::high_resolution_clock::now();
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
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  }

  {
    const char* env_val = std::getenv("NEKRS_SIGNUM_BACKTRACE");
    if (env_val)
      std::signal(std::atoi(env_val), signalHandler);  
  }

  MPI_Barrier(MPI_COMM_WORLD);
  const double time0 = MPI_Wtime(); 

  MPI_Comm commGlobal;
  MPI_Comm_dup(MPI_COMM_WORLD, &commGlobal);

  {
    if(!getenv("NEKRS_HOME")) {
      std::cout << "FATAL ERROR: Cannot find env variable NEKRS_HOME!" << "\n";
      fflush(stdout);
      MPI_Abort(commGlobal, EXIT_FAILURE);
    }

    std::string bin(getenv("NEKRS_HOME"));
    bin += "/bin/nekrs";
    const char* ptr = realpath(bin.c_str(), NULL);
    if(!ptr) {
      std::cout << "FATAL ERROR: Cannot find " << bin << "!\n";
      fflush(stdout);
      MPI_Abort(commGlobal, EXIT_FAILURE);
    }
  }

  cmdOptions* cmdOpt = processCmdLineOptions(argc, argv, commGlobal);

  MPI_Comm comm = setupSession(cmdOpt, commGlobal);

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (rank == 0) {
     time_t now = time(0);
     tm *gmtm = gmtime(&now);
     char *dt = asctime(gmtm);
     std::cout << "UTC time: " << dt << std::endl;
     fflush(stdout);
  }

  if (cmdOpt->debug) {
    MPI_Barrier(comm);
    for(int currRank = 0; currRank < size; ++currRank)
      if(rank == currRank) printf("rank %d: pid<%d>\n", rank, ::getpid());
    fflush(stdout);
    MPI_Barrier(comm);
    if (rank == 0) std::cout << "Attach debugger, then press enter to continue\n";
    if (rank == 0) std::cin.get();
    MPI_Barrier(comm);
  }

  try{

  if (cmdOpt->debug) feraiseexcept(FE_ALL_EXCEPT);

  { // change working dir
    const size_t last_slash = cmdOpt->setupFile.rfind('/') + 1;
    const std::string casepath = cmdOpt->setupFile.substr(0,last_slash);
    chdir(casepath.c_str());
    const std::string casename = cmdOpt->setupFile.substr(last_slash, cmdOpt->setupFile.length() - last_slash);
    if(casepath.length() > 0) chdir(casepath.c_str());
    cmdOpt->setupFile.assign(casename);
  }

  nekrs::setup(commGlobal, comm, 
    	       cmdOpt->buildOnly, cmdOpt->sizeTarget,
               cmdOpt->ciMode, cmdOpt->setupFile,
               cmdOpt->backend, cmdOpt->deviceID,
               cmdOpt->nSessions, cmdOpt->sessionID,
               cmdOpt->debug);

  if (cmdOpt->buildOnly) {
    nekrs::finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  int tStep = 0;
  double time = nekrs::startTime();

  double elapsedTime = 0;
  {
    MPI_Barrier(comm);
    const auto timeStop = std::chrono::high_resolution_clock::now();
    elapsedTime += std::chrono::duration<double, std::milli>(timeStop - timeStart).count() / 1e3;
    MPI_Allreduce(MPI_IN_PLACE, &elapsedTime, 1, MPI_DOUBLE, MPI_MAX, comm);
    nekrs::updateTimer("setup", elapsedTime);
    if (rank == 0)
      std::cout << "initialization took " << elapsedTime << " s" << std::endl;
  }

  int isLastStep = 1;
  if (nekrs::endTime() > nekrs::startTime() || nekrs::numSteps() > tStep) {
    isLastStep = 0;
  }
  nekrs::lastStep(isLastStep);

  nekrs::udfExecuteStep(time, tStep, /* outputStep */ 0);
  nekrs::resetTimer("udfExecuteStep");

  double elapsedStepSum = 0;

  double tSolveStepMin = std::numeric_limits<double>::max();
  double tSolveStepMax = std::numeric_limits<double>::min();

  if (rank == 0) {
    if (isLastStep)
      std::cout << "endTime or numSteps reached already -> skip timestepping\n"; 
    else if (nekrs::endTime() > nekrs::startTime())
      std::cout << "\ntimestepping to time " << nekrs::endTime() << " ...\n";
    else
      std::cout << "\ntimestepping for " << nekrs::numSteps() << " steps ...\n";
  }

  fflush(stdout);
  MPI_Pcontrol(1);
  while (!isLastStep) {
    MPI_Barrier(comm);
    const double timeStartStep = MPI_Wtime();

    ++tStep;
    isLastStep = nekrs::lastStep(time, tStep, elapsedTime);

    double dt;
    if (isLastStep && nekrs::endTime() > 0)
      dt = nekrs::endTime() - time;
    else
      dt = nekrs::dt(tStep);

    int outputStep = nekrs::outputStep(time + dt, tStep);
    if (nekrs::writeInterval() == 0) outputStep = 0;
    if (isLastStep) outputStep = 1;
    if (nekrs::writeInterval() < 0) outputStep = 0;
    nekrs::outputStep(outputStep);

    if (tStep <= 1000) nekrs::verboseInfo(true); 

    nekrs::initStep(time, dt, tStep);
    
    int corrector = 1;
    bool converged = false;
    do {
      converged = nekrs::runStep(corrector++);
    } while (!converged);
 
    time = nekrs::finishStep();

    if(nekrs::updateFileCheckFreq()) {
      if(tStep % nekrs::updateFileCheckFreq()) 
        nekrs::processUpdFile();
    }

    if (nekrs::printInfoFreq()) {
      if (tStep % nekrs::printInfoFreq() == 0)
        nekrs::printInfo(time, tStep, false, true);
    }

    if (outputStep) nekrs::outfld(time, tStep);

#ifdef ENABLE_SENSEI
    nekrs::runSensei(time, dt, tStep);
#endif

    MPI_Barrier(comm);
    const double elapsedStep = MPI_Wtime() - timeStartStep;
    tSolveStepMin = std::min(elapsedStep, tSolveStepMin);
    tSolveStepMax = std::max(elapsedStep, tSolveStepMax);
    nekrs::updateTimer("minSolveStep", tSolveStepMin);
    nekrs::updateTimer("maxSolveStep", tSolveStepMax);

    elapsedStepSum += elapsedStep;
    elapsedTime += elapsedStep;
    nekrs::updateTimer("elapsedStep", elapsedStep);
    nekrs::updateTimer("elapsedStepSum", elapsedStepSum);
    nekrs::updateTimer("elapsed", elapsedTime);

    if (nekrs::printInfoFreq()) {
      if (tStep % nekrs::printInfoFreq() == 0)
        nekrs::printInfo(time, tStep, true, false);
    }

    if(nekrs::runTimeStatFreq()) {
      if (tStep % nekrs::runTimeStatFreq() == 0 || isLastStep)
        nekrs::printRuntimeStatistics(tStep);
    }

    if (tStep % 100 == 0) fflush(stdout);
  }
  MPI_Pcontrol(0);

  delete cmdOpt;

  const int exitValue = nekrs::finalize();

#ifdef ENABLE_SENSEI
  sensei_bridge_finalize();
#endif

  MPI_Barrier(commGlobal);
  MPI_Finalize();

  if(exitValue)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;

  } // try
  catch(const std::exception& ex)
  {
    std::cerr << ex.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
}
