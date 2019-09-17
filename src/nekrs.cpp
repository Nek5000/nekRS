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

#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

#include "nekrs.hpp"
#include "meshNekSetupHex3D.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "header.h"

#define NEKLDIMT 2

int rank, size;
MPI_Comm comm;
int buildOnly = 0;  // hack for libParanumal::meshPhysicalNodes() 
int ciMode = 0;

void dryRun(libParanumal::setupAide &options, string udfFile, 
            string casename,int N, int npTarget)
{
  if (rank == 0) 
    cout << "performing dry-run for "
         << npTarget
         << " MPI ranks ...\n" << endl;

  // jit compile udf
  if (!udfFile.empty()) {
    if(rank == 0) udfBuild(udfFile.c_str()); 
    MPI_Barrier(comm);
    *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
    *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
  } 

  if(udf.setup0) udf.setup0(comm, options);

  // jit compile nek
  if (rank == 0) buildNekInterface(casename.c_str(), NEKLDIMT, N, npTarget);
  MPI_Barrier(comm);

  // jit compile libP kernels
  mesh_t *mesh = new mesh_t[1];
  mesh->comm = comm;
  mesh->rank = rank;
  mesh->size = size;
  meshBoxSetupHex3D(N, mesh);
  ins_t *ins = setup(mesh, options);

  // jit compile udf kernels
  if (udf.loadKernels) {
    if (rank == 0) cout << "building udf kernels ...";
    udf.loadKernels(ins);
    if (rank == 0) cout << " done" << endl;
  }

  if (mesh->rank == 0) cout << "\nBuild successful." << endl;
}

int processCmdLineOptions(int argc, char **argv, string &setupFile, 
                          int &build, int &sizeTarget)
{
  int foundSetup = 0;

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
            setupFile.assign(optarg);  
            setupFile = setupFile + ".par";
            if (const char *ptr = realpath(setupFile.c_str(), NULL)) 
              foundSetup = 1;
            else 
              cout << "ERROR: Cannot find " << setupFile << "!\n";
            break;  
        case 'b':  
            build = 1;
            sizeTarget = atoi(optarg);  
            break; 
         case 'c':  
            ciMode = atoi(optarg);
            if (ciMode < 1) {
              cout << "ERROR: ci test id has to be >0!\n";
              return 1;
            }
            break;  
        default:  
            return 1;
    }  
  }  
  if (!foundSetup) return 1; 

  return 0;
}

void setDataFile(libParanumal::setupAide &options)
{
  std::string oklFile;
  options.getArgs("UDF OKL FILE",oklFile);

  char buf[FILENAME_MAX];

  char *ptr = realpath(oklFile.c_str(), NULL);
  if(!ptr) {
    if (rank ==0) cout << "ERROR: Cannot find " << oklFile << "!\n";
    EXIT(-1);
  }
  free(ptr);

  std::string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
  string casename;
  options.getArgs("CASENAME", casename);
  string dataFile = cache_dir + "/" + casename + ".okl";

  if (rank == 0) {

    std::ifstream in;
    in.open(oklFile);
    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();

    std::ofstream out;
    out.open(dataFile, std::ios::trunc);

    out << buffer.str();

    out << "// automatically added \n" 
        << "void insFlowField3D(bcData *bc){}\n"
        << "void insPressureNeumannConditions3D(bcData *bc){}\n"; 

    std::size_t found;
    found = buffer.str().find("void insVelocityDirichletConditions");
    if (found == std::string::npos)
      out << "void insVelocityDirichletConditions3D(bcData *bc){}\n"; 

    found = buffer.str().find("void insVelocityNeumannConditions");
    if (found == std::string::npos)
      out << "void insVelocityNeumannConditions3D(bcData *bc){}\n"; 

    found = buffer.str().find("void insPressureDirichletConditions");
    if (found ==std::string::npos)
      out << "void insPressureDirichletConditions3D(bcData *bc){}\n"; 

    out.close();
  }

  options.setArgs("DATA FILE", dataFile);
}

int main(int argc, char **argv)
{
  int retval;
  int provided, required = MPI_THREAD_MULTIPLE;
  //retval =  MPI_Init_thread(&argc, &argv, required, &provided);
  retval =  MPI_Init(&argc, &argv);

  if (retval != MPI_SUCCESS) {
    cout << "FATAL ERROR: Cannot initialize MPI!" << endl;
    exit(-1);
  }
/*
  if (provided != required) {
    cout << "FATAL ERROR: MPI does not support MPI_THREAD_MULTIPLE!\n";
    exit(-1);
  }
*/

  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // parse cmdline
  char buf[FILENAME_MAX];
  string setupFile;
  int sizeTarget;
  if (rank == 0)
    retval = processCmdLineOptions(argc, argv, setupFile, buildOnly, sizeTarget);

  MPI_Bcast(&retval, sizeof(retval), MPI_BYTE, 0, comm);
  if (retval) {
    if (rank == 0)
      cout << "usage: ./nekrs --setup <case name> [ --build-only <#procs> ] [ --cimode <id> ]"
           << endl;
    EXIT(-1);
  }

  strcpy(buf, setupFile.c_str());
  MPI_Bcast(buf, sizeof(buf), MPI_BYTE, 0, comm);
  setupFile.assign(buf);
  MPI_Bcast(&sizeTarget, sizeof(sizeTarget), MPI_BYTE, 0, comm);
  MPI_Bcast(&ciMode, sizeof(ciMode), MPI_BYTE, 0, comm);

  // set cache dir
  getcwd(buf, sizeof(buf));
  string cwd;
  cwd.assign(buf);
  sprintf(buf,"%s/.cache", cwd.c_str());
  setenv("NEKRS_CACHE_DIR", buf, 1);
  string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
  if (rank == 0) mkdir(cache_dir.c_str(), S_IRWXU);
  MPI_Barrier(comm);

  // print header
  if (rank == 0) {
    printHeader();
    
    cout << "MPI ranks: " << size << endl << endl;

    if (const char *env_p = getenv("OCCA_CACHE_DIR"))
      cout << "using OCCA_CACHE_DIR: " << env_p << endl << endl;
    else
      cout << "OCCA_CACHE_DIR undefined fallback to default" << endl << endl;
  }

#ifdef DEBUG
  if (rank == 0) {
    char str[10];
    cout << "Press any key to continue" << endl;
    gets(str);
  }
  MPI_Barrier(comm);
#endif

  // read settings
  libParanumal::setupAide options = parRead(setupFile, comm);

  int N;
  string udfFile, casename;
  int readRestartFile;
  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("UDF FILE", udfFile);
  options.getArgs("CASENAME", casename);
  options.getArgs("RESTART FROM FILE", readRestartFile);

  setDataFile(options);

  if (buildOnly){
    dryRun(options, udfFile, casename, N, sizeTarget);
    EXIT(0);
  }

  MPI_Barrier(comm); double t0 = MPI_Wtime();

  // jit compile udf
  if (!udfFile.empty()){
    if(rank == 0) udfBuild(udfFile.c_str()); 
    MPI_Barrier(comm);
    udfLoad();
  } 

  if(udf.setup0) udf.setup0(comm, options);

  // jit compile nek
  if(rank == 0) buildNekInterface(casename.c_str(), NEKLDIMT, N, size);
  MPI_Barrier(comm);

  if(rank == 0 && ciMode) 
    cout << "performing continous integration tests\n" << endl;

  // init nek
  nek_setup(comm, options);
  nek_setic();
  nek_userchk();

  // enfore num steps
  dfloat startTime;
  options.getArgs("START TIME", startTime);
  if(startTime > 0.0) {
    int numSteps;
    if(options.getArgs("NUMBER TIMESTEPS", numSteps)) {
      dfloat endTime;
      options.getArgs("FINAL TIME", endTime);
      endTime += startTime;
      options.setArgs("FINAL TIME", to_string_f(endTime));
    }
  }

  // setup libP
  mesh_t *mesh = new mesh_t[1];
  mesh->comm = comm;
  mesh->rank = rank;
  mesh->size = size;
  meshNekSetupHex3D(N, mesh);
  ins_t *ins = setup(mesh, options);

  // jit compile udf kernels
  if (udf.loadKernels) {
    if (rank == 0) cout << "building udf kernels ...";
    udf.loadKernels(ins);
    if (rank == 0) cout << " done" << endl;
  }

  // set initial condition
  for (int n = 0; n < mesh->Np*mesh->Nelements; ++n) {
    ins->U[0*ins->fieldOffset + n] = 0;
    ins->U[1*ins->fieldOffset + n] = 0;
    ins->U[2*ins->fieldOffset + n] = 0;
    ins->P[n] = 0;
  }
  if(readRestartFile) {
    dlong Nlocal = mesh->Nelements*mesh->Np;
    if (*(nekData.ifgetu)) {
      dfloat *vx = ins->U + 0*ins->fieldOffset;
      dfloat *vy = ins->U + 1*ins->fieldOffset;
      dfloat *vz = ins->U + 2*ins->fieldOffset;
      memcpy(vx, nekData.vx, sizeof(dfloat)*Nlocal);
      memcpy(vy, nekData.vy, sizeof(dfloat)*Nlocal);
      memcpy(vz, nekData.vz, sizeof(dfloat)*Nlocal);
    }
    if (*(nekData.ifgetp)) memcpy(ins->P, nekData.pr, sizeof(dfloat)*Nlocal);
  }
  if(udf.setup) udf.setup(ins);
  ins->o_U.copyFrom(ins->U);
  ins->o_P.copyFrom(ins->P);

  if (udf.executeStep) udf.executeStep(ins, ins->startTime, 0);
  nek_copyFrom(ins, ins->startTime, 0);

  if(rank == 0) {
    cout << "\nsettings:\n" << endl;
    cout << ins->vOptions << endl;
    cout << "\ninitialization took " << MPI_Wtime() - t0 << " seconds" << endl; 
    cout << "starting time loop" << endl;
  }
  fflush(stdout);

  // run time stepper
  MPI_Barrier(mesh->comm); t0 = MPI_Wtime();
  MPI_Pcontrol(1);

  if(ins->options.compareArgs("TIME INTEGRATOR", "TOMBO")) runTombo(ins); 

  MPI_Pcontrol(0);
  MPI_Barrier(mesh->comm); double tElapsed = MPI_Wtime() - t0;


  // finalize
  if(rank == 0) { 
    cout << "\nreached final time " << ins->finalTime << " in " 
         << tElapsed << " seconds" << endl
         << "\nEnd." << endl;
  }
  MPI_Finalize();

  return 0;
}
