/*---------------------------------------------------------------------------*\
Copyright (c) 2019-2020, UCHICAGO ARGONNE, LLC. 

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
#include <cstring>
#include <getopt.h>
#include <cfenv>
#include "nekrs.hpp"

#define DEBUG

static MPI_Comm comm;

struct cmdOptions {
  int buildOnly = 0;
  int ciMode = 0;
  int sizeTarget = 0;
  int debug = 0;
  std::string setupFile;
  std::string deviceID;
  std::string backend;
};

static cmdOptions *processCmdLineOptions(int argc, char **argv); 


int main(int argc, char **argv)
{
  int retval;
  retval =  MPI_Init(&argc, &argv);
  if (retval != MPI_SUCCESS) {
    std::cout << "FATAL ERROR: Cannot initialize MPI!" << "\n";
    exit(1);
  }

  int rank, size;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  cmdOptions *cmdOpt = processCmdLineOptions(argc, argv);

  if (cmdOpt->debug) { 
    if (rank == 0) {
      std::cout << "Attach debugger, then press enter to continue\n";
      std::cin.get(); 
    }
    MPI_Barrier(comm);
  } 

  if(cmdOpt->debug) feraiseexcept(FE_ALL_EXCEPT);  

  std::string cacheDir; 
  nekrs::setup(comm, cmdOpt->buildOnly, cmdOpt->sizeTarget,
               cmdOpt->ciMode, cacheDir, cmdOpt->setupFile, 
               cmdOpt->backend, cmdOpt->deviceID);

  if (cmdOpt->buildOnly) {
    MPI_Finalize(); 
    return EXIT_SUCCESS;
  }

  const int outputStep = nekrs::outputStep();
  const int NtimeSteps = nekrs::NtimeSteps(); 
  const double startTime = nekrs::startTime();
  const double finalTime = nekrs::finalTime();

  if (rank == 0) std::cout << "\nstarting time loop" << "\n";
    
  double time = startTime;
  int tStep = 1;
  MPI_Pcontrol(1);
  while ((finalTime-time)/finalTime > 1e-6*nekrs::dt()) {

    nekrs::runStep(time, nekrs::dt(), tStep);
    time += nekrs::dt();

    int isOutputStep = 0;
    if (outputStep > 0) {
      if (tStep%outputStep == 0 || tStep == NtimeSteps) isOutputStep = 1;
    }

    nekrs::udfExecuteStep(time, tStep, isOutputStep);
    if (isOutputStep) {
      nekrs::copyToNek(time, tStep); 
      nekrs::nekOutfld();
    }

    ++tStep;
  }
  MPI_Pcontrol(0);

  nekrs::printRuntimeStatistics();

  if(rank == 0) std::cout << "\nEnd." << "\n";

  MPI_Finalize();
  return EXIT_SUCCESS;
}

static cmdOptions *processCmdLineOptions(int argc, char **argv) 
{
  int rank;
  MPI_Comm_rank(comm, &rank);
 
  cmdOptions *cmdOpt = new cmdOptions();
 
  int err = 0;

  if (rank == 0){
    while(1){
      static struct option long_options[] =
      {
          {"setup", required_argument, 0, 's'},
          {"cimode", required_argument, 0, 'c'},
          {"build-only", required_argument, 0, 'b'},
          {"debug", no_argument, 0, 'd'},
          {"backend", required_argument, 0, 't'},
          {"device-id", required_argument, 0, 'i'},
        {0, 0, 0, 0}
      };
      int option_index = 0;
      int c = getopt_long (argc, argv, "s:", long_options, &option_index);
    
      if (c == -1)
        break;
    
      switch(c){ 
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
          default:  
              err = 1;
      }
    }  
  } 

  char buf[FILENAME_MAX];
  strcpy(buf, cmdOpt->setupFile.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, comm);
  cmdOpt->setupFile.assign(buf);
  strcpy(buf, cmdOpt->deviceID.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, comm);
  cmdOpt->deviceID.assign(buf);
  strcpy(buf, cmdOpt->backend.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, comm);
  cmdOpt->backend.assign(buf);
  MPI_Bcast(&cmdOpt->buildOnly, sizeof(cmdOpt->buildOnly), MPI_BYTE, 0, comm);
  MPI_Bcast(&cmdOpt->sizeTarget, sizeof(cmdOpt->sizeTarget), MPI_BYTE, 0, comm);
  MPI_Bcast(&cmdOpt->ciMode, sizeof(cmdOpt->ciMode), MPI_BYTE, 0, comm);
  MPI_Bcast(&cmdOpt->debug, sizeof(cmdOpt->debug), MPI_BYTE, 0, comm);

  if(cmdOpt->setupFile.empty()) err++;

  MPI_Bcast(&err, sizeof(err), MPI_BYTE, 0, comm);
  if (err) {
    if (rank == 0)
      std::cout << "usage: ./nekrs --setup <case name> "
                << "[ --build-only <#procs> ] [ --cimode <id> ] [ --debug ] "
                << "[ --backend <CPU|CUDA|HIP|OPENCL> ] [ --device-id <id|LOCAL-RANK> ]"
                << "\n";
    MPI_Finalize(); 
    exit(1);
  }



 return cmdOpt;
}
