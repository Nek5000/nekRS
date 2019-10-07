/*---------------------------------------------------------------------------*\
Copyright (c) 2019, UCHICAGO ARGONNE, LLC. 

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
#include <cstdio>
#include <string>
#include <cstring>
#include <getopt.h>
#include "nrs.hpp"

static MPI_Comm comm;

struct cmdOptions {
  int buildOnly = 0;
  int ciMode = 0;
  int sizeTarget = 0;
  string setupFile;
};

static cmdOptions *processCmdLineOptions(int argc, char **argv); 


int main(int argc, char **argv)
{
  int retval;
  retval =  MPI_Init(&argc, &argv);
  if (retval != MPI_SUCCESS) {
    cout << "FATAL ERROR: Cannot initialize MPI!" << endl;
    exit(1);
  }

  int rank, size;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
 
  cmdOptions *cmdOpt = processCmdLineOptions(argc, argv);

  string cacheDir; 
  setupAide options;
  options = nrs::setup(comm, cmdOpt->buildOnly, cmdOpt->sizeTarget,
                       cmdOpt->ciMode, cacheDir, cmdOpt->setupFile);

  if (cmdOpt->buildOnly) {
    MPI_Finalize(); 
    return 0;
  }

  int outputStep, isOutputStep; 
  int NtimeSteps; 
  double startTime;
  double dt;
  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", outputStep);
  options.getArgs("NUMBER TIMESTEPS", NtimeSteps); 
  options.getArgs("START TIME", startTime);
  options.getArgs("DT", dt);

  if (rank == 0)
    cout << "\nstarting time loop" << endl;

  double time = startTime;
  MPI_Pcontrol(1);
  for(int tstep=1; tstep <= NtimeSteps; ++tstep){

    isOutputStep = 0;
    if (outputStep > 0) {
      if ((tstep%outputStep)==0 ||  tstep == NtimeSteps) {
        isOutputStep = 1;
      }
    }

    nrs::runStep(time, tstep);

    if (isOutputStep) nrs::copyToNek(time+dt, tstep);
    nrs::udfExecuteStep(time+dt, tstep, isOutputStep);
    if (isOutputStep) nrs::nekOutfld();

    time += dt;

  }
  MPI_Pcontrol(0);

  if(rank == 0)
    cout << "\nEnd." << endl;
  MPI_Finalize();
  return 0;
}

static cmdOptions *processCmdLineOptions(int argc, char **argv) 
{
  int rank;
  MPI_Comm_rank(comm, &rank);
 
  cmdOptions *cmdOpt = new cmdOptions();
 
  int err = 0;
  int foundSetup = 0;

  if (rank == 0){
    while(1){
      static struct option long_options[] =
      {
          {"setup", required_argument, 0, 's'},
          {"cimode", required_argument, 0, 'c'},
          {"build-only", required_argument, 0, 'b'},
          {0, 0, 0, 0}
      };
      int option_index = 0;
      int c = getopt_long (argc, argv, "s:d:", long_options, &option_index);
    
      if (c == -1)
        break;
    
      switch(c){ 
          case 's':
              cmdOpt->setupFile.assign(optarg);  
              cmdOpt->setupFile = cmdOpt->setupFile + ".par";
              if (const char *ptr = realpath(cmdOpt->setupFile.c_str(), NULL)) 
                foundSetup = 1;
              else 
                cout << "ERROR: Cannot find " << cmdOpt->setupFile << "!\n";
              break;  
          case 'b':  
              cmdOpt->buildOnly = 1;
              cmdOpt->sizeTarget = atoi(optarg);  
              break; 
           case 'c':  
              cmdOpt->ciMode = atoi(optarg);
              if (cmdOpt->ciMode < 1) {
                cout << "ERROR: ci test id has to be >0!\n";
                err = 1;
              }
              break;  
          default:  
              err = 1;
      }
    }  
    if (!foundSetup) err = 1;
  } 

  MPI_Bcast(&err, sizeof(err), MPI_BYTE, 0, comm);
  if (err) {
    if (rank == 0)
      cout << "usage: ./nekrs --setup <case name> [ --build-only <#procs> ] [ --cimode <id> ]"
           << endl;
    MPI_Finalize(); 
    exit(1);
  }

  char buf[FILENAME_MAX];
  strcpy(buf, cmdOpt->setupFile.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, comm);
  cmdOpt->setupFile.assign(buf);
  MPI_Bcast(&cmdOpt->sizeTarget, sizeof(cmdOpt->sizeTarget), MPI_BYTE, 0, comm);
  MPI_Bcast(&cmdOpt->ciMode, sizeof(cmdOpt->ciMode), MPI_BYTE, 0, comm);

 return cmdOpt;
}
