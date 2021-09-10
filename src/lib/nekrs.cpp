#include <stdlib.h>
#include "nrs.hpp"
#include "meshSetup.hpp"
#include "setup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "parReader.hpp"
#include "configReader.hpp"
#include "timeStepper.hpp"
#include "platform.hpp"
#include "nrssys.hpp"
#include "linAlg.hpp"
#include "cfl.hpp"
#include "amgx.h"

// extern variable from nrssys.hpp
platform_t* platform;

static int rank, size;
static MPI_Comm comm;
static nrs_t* nrs;
static setupAide options;
static dfloat lastOutputTime = 0;
static int enforceLastStep = 0;
static int enforceOutputStep = 0;

static void setOccaVars();

void printHeader()
{
  std::cout << R"(                 __    ____  _____)" << std::endl
            << R"(   ____   ___   / /__ / __ \/ ___/)" << std::endl
            << R"(  / __ \ / _ \ / //_// /_/ /\__ \ )" << std::endl
            << R"( / / / //  __// ,<  / _, _/___/ / )" << std::endl
            << R"(/_/ /_/ \___//_/|_|/_/ |_|/____/  )"
            << "v" << NEKRS_VERSION << "." << NEKRS_SUBVERSION 
            << " (" << GITCOMMITHASH << ")" << std::endl
            << std::endl
            << "COPYRIGHT (c) 2019-2021 UCHICAGO ARGONNE, LLC" << std::endl
            << std::endl;
}

namespace nekrs
{
double startTime(void)
{
  double val = 0;
  platform->options.getArgs("START TIME", val);
  return val;
}

void setup(MPI_Comm comm_in, int buildOnly, int commSizeTarget,
           int ciMode, std::string cacheDir, std::string _setupFile,
           std::string _backend, std::string _deviceID)
{
  MPI_Comm_dup(comm_in, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (rank == 0) {
    printHeader();
    std::cout << "MPI tasks: " << size << std::endl << std::endl;
  }

  srand48((long int) rank);

  configRead(comm);

  oogs::gpu_mpi(std::stoi(getenv("NEKRS_GPU_MPI")));

  if (rank == 0) std::cout << "reading par file ...\n"; 
  auto par = new inipp::Ini();	  
  std::string setupFile = _setupFile + ".par";
  options = parRead((void*) par, setupFile, comm);

  {
    char buf[FILENAME_MAX];
    char * ret = getcwd(buf, sizeof(buf));
    if(!ret) ABORT(EXIT_FAILURE);;
    std::string cwd;
    cwd.assign(buf);
 
    std::string dir(cacheDir);
    if (cacheDir.empty()) dir = cwd + "/.cache";
    if(getenv("NEKRS_CACHE_DIR")) dir.assign(getenv("NEKRS_CACHE_DIR"));
    setenv("NEKRS_CACHE_DIR", dir.c_str(), 1);
  }

  setOccaVars();

  if (rank == 0) {
    std::string install_dir;
    install_dir.assign(getenv("NEKRS_HOME"));
    std::cout << std::endl;
    std::cout << "using NEKRS_HOME: " << install_dir << std::endl;

    std:: string cache_dir;
    cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
    std::cout << "using NEKRS_CACHE_DIR: " << cache_dir << std::endl;

    std::cout << "using OCCA_CACHE_DIR: " << occa::env::OCCA_CACHE_DIR << std::endl << std::endl;
  }

  options.setArgs("BUILD ONLY", "FALSE");
  if(buildOnly) {
    options.setArgs("BUILD ONLY", "TRUE");
    options.setArgs("NP TARGET", std::to_string(commSizeTarget));
    if(rank == 0){
      std::cout << "jit-compiling for >="
           << commSizeTarget 
           << " MPI tasks ...\n" << std::endl;
    }
    fflush(stdout);	
  }

  if (options.getArgs("THREAD MODEL").length() == 0) 
    options.setArgs("THREAD MODEL", getenv("NEKRS_OCCA_MODE_DEFAULT"));
  if(!_backend.empty()) options.setArgs("THREAD MODEL", _backend);
  if(!_deviceID.empty()) options.setArgs("DEVICE NUMBER", _deviceID);

  // setup device
  platform_t* _platform = platform_t::getInstance(options, comm);
  platform = _platform;
  platform->par = par;

  platform->timer.tic("setup", 1);

  int buildRank = rank;
  if(getenv("NEKRS_BUILD_NODE_LOCAL"))
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);    

  if(buildRank == 0) {
    std::string cache_dir;
    cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
    mkdir(cache_dir.c_str(), S_IRWXU);
    std::string udf_dir;
    udf_dir.assign(getenv("NEKRS_UDF_DIR"));
    mkdir(udf_dir.c_str(), S_IRWXU);
  }

  oudfInit(options);

  // jit compile udf
  std::string udfFile;
  options.getArgs("UDF FILE", udfFile);
  if (!udfFile.empty()) {
    int err = 0;
    if(buildRank == 0) err = udfBuild(udfFile.c_str(), options);

    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, comm);
    if(err) ABORT(EXIT_FAILURE);

    if(buildOnly) {
      *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
      *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
    } else {
      udfLoad();
    }
  }

  options.setArgs("CI-MODE", std::to_string(ciMode));
  if(rank == 0 && ciMode)
    std::cout << "enabling continous integration mode " << ciMode << "\n";

  nek::bootstrap();

  if(udf.setup0) udf.setup0(comm, options);

  compileKernels();

  if(buildOnly) {
    MPI_Barrier(platform->comm.mpiComm);
    if(buildRank == 0) {
      std::string cache_dir;
      cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
      std::string file = cache_dir + "/build-only.timestamp";
      remove(file.c_str());
      std::ofstream ofs;
      ofs.open(file, std::ofstream::out);
      ofs.close();
      std::cout << "\nBuild successful." << std::endl;
    }
    return;
  }

  platform->linAlg = linAlg_t::getInstance();

  nrs = new nrs_t();

  nrsSetup(comm, options, nrs);

  nrs->o_U.copyFrom(nrs->U);
  nrs->o_P.copyFrom(nrs->P);
  nrs->o_prop.copyFrom(nrs->prop);
  if(nrs->Nscalar) {
    nrs->cds->o_S.copyFrom(nrs->cds->S);
    nrs->cds->o_prop.copyFrom(nrs->cds->prop);
  }

  evaluateProperties(nrs, startTime());
  nrs->o_prop.copyTo(nrs->prop);
  if(nrs->Nscalar) nrs->cds->o_prop.copyTo(nrs->cds->prop);

  nek::ocopyToNek(startTime(), 0);

  platform->timer.toc("setup");
  const double setupTime = platform->timer.query("setup", "DEVICE:MAX");
  if(rank == 0) {
    std::cout << "\nsettings:\n" << std::endl << options << std::endl;
    std::cout << "occa memory usage: " << platform->device.memoryAllocated()/1e9 << " GB" << std::endl;
    std::cout << "initialization took " << setupTime << " s" << std::endl;
  }
  fflush(stdout);

  platform->timer.set("setup", setupTime);
}

void runStep(double time, double dt, int tstep)
{
  timeStepper::step(nrs, time, dt, tstep);
}

void copyFromNek(double time, int tstep)
{
  nek::ocopyToNek(time, tstep);
}

void udfExecuteStep(double time, int tstep, int isOutputStep)
{
  platform->timer.tic("udfExecuteStep", 1);
  if (isOutputStep) {
    nek::ifoutfld(1);
    nrs->isOutputStep = 1;
  }

  if (udf.executeStep) udf.executeStep(nrs, time, tstep);

  nek::ifoutfld(0);
  nrs->isOutputStep = 0;
  platform->timer.toc("udfExecuteStep");
}

void nekUserchk(void)
{
  nek::userchk();
}

double dt(int tstep)
{
  if(platform->options.compareArgs("VARIABLE DT", "TRUE")){

    if(tstep == 1)
    {
      // use user-specified initial dt on startup
      double initialDt = 0.0;
      platform->options.getArgs("DT", initialDt);
      if(initialDt > 0.0){
        nrs->dt[0] = initialDt;
        return nrs->dt[0];
      }
    }
    timeStepper::computeTimeStepFromCFL(nrs, tstep);
  }
  
  {
    double maxDt = std::numeric_limits<double>::max();
    platform->options.getArgs("MAX DT", maxDt);
    nrs->dt[0] = (nrs->dt[0] < maxDt) ? nrs->dt[0] : maxDt;
  }

  return nrs->dt[0];
}

double writeInterval(void)
{
  double val = -1;
  platform->options.getArgs("SOLUTION OUTPUT INTERVAL", val);
  return val;
}

int writeControlRunTime(void)
{
  return platform->options.compareArgs("SOLUTION OUTPUT CONTROL", "RUNTIME");
}

int outputStep(double time, int tStep)
{
  int outputStep = 0;
  if (writeControlRunTime()) {
    double val;
    platform->options.getArgs("START TIME", val);
    if(lastOutputTime == 0 && val > 0) lastOutputTime = val;
    outputStep = ((time - lastOutputTime) + 1e-10) > nekrs::writeInterval();
  } else {
    if (writeInterval() > 0) outputStep = (tStep%(int)writeInterval() == 0);
  }

  if (enforceOutputStep) {
    enforceOutputStep = 0;
    return 1;
  }
  return outputStep;
}

void outputStep(int val)
{
  nrs->isOutputStep = val;
}

void outfld(double time)
{
  writeFld(nrs, time);
  lastOutputTime = time;
}

void outfld(double time, std::string suffix)
{
  writeFld(nrs, time, suffix);
  lastOutputTime = time;
}

double endTime(void)
{
  double endTime = -1;
  platform->options.getArgs("END TIME", endTime);
  return endTime;
}

int numSteps(void)
{
  int numSteps = -1;
  platform->options.getArgs("NUMBER TIMESTEPS", numSteps);
  return numSteps;
}

int lastStep(double time, int tstep, double elapsedTime)
{
  if(!platform->options.getArgs("STOP AT ELAPSED TIME").empty()) {
    double maxElaspedTime;
    platform->options.getArgs("STOP AT ELAPSED TIME", maxElaspedTime);
    if(elapsedTime > 60.0*maxElaspedTime) nrs->lastStep = 1; 
  } else if (endTime() > 0) { 
     const double eps = 1e-12;
     nrs->lastStep = fabs((time+nrs->dt[0]) - endTime()) < eps || (time+nrs->dt[0]) > endTime();
  } else {
    nrs->lastStep = tstep == numSteps();
  }
 
  if(enforceLastStep) return 1;
  return nrs->lastStep;
}

void* nekPtr(const char* id)
{
  return nek::ptr(id);
}

void* nrsPtr(void)
{
  return nrs;
}

void finalize(void)
{
  AMGXfree();
}

void printRuntimeStatistics(int step)
{
  platform_t* platform = platform_t::getInstance(options, comm);
  platform->timer.printRunStat(step);
}

void processUpdFile()
{
  char* rbuf = nullptr;
  long fsize = 0;

  if (rank == 0) {
    const std::string cmdFile = "nekrs.upd";
    const char* ptr = realpath(cmdFile.c_str(), NULL);
    if (ptr) {
      if(rank == 0) std::cout << "processing " << cmdFile << " ...\n";
      FILE* f = fopen(cmdFile.c_str(), "rb");
      fseek(f, 0, SEEK_END);
      fsize = ftell(f);
      fseek(f, 0, SEEK_SET);
      rbuf = new char[fsize];
      fread(rbuf, 1, fsize, f);
      fclose(f);
      remove(cmdFile.c_str());
    }
  }
  MPI_Bcast(&fsize, sizeof(fsize), MPI_BYTE, 0, comm);

  if (fsize) {
    if(rank != 0) rbuf = new char[fsize];
    MPI_Bcast(rbuf, fsize, MPI_CHAR, 0, comm);
    std::stringstream is;
    is.write(rbuf, fsize);
    inipp::Ini ini;
    ini.parse(is, false);

    std::string end;
    ini.extract("", "end", end);
    if (end == "true") {
      enforceLastStep = 1; 
      platform->options.setArgs("END TIME", "-1");
    }

    std::string checkpoint;
    ini.extract("", "checkpoint", checkpoint);
    if (checkpoint == "true") enforceOutputStep = 1; 

    std::string endTime;
    ini.extract("general", "endtime", endTime);
    if (!endTime.empty()) {
      if (rank == 0) std::cout << "  set endTime = " << endTime << "\n";
      platform->options.setArgs("END TIME", endTime);
    }

    std::string numSteps;
    ini.extract("general", "numsteps", numSteps);
    if (!numSteps.empty()) {
      if (rank == 0) std::cout << "  set numSteps = " << numSteps << "\n";
      platform->options.setArgs("NUMBER TIMESTEPS", numSteps);
    }

    std::string writeInterval;
    ini.extract("general", "writeinterval", writeInterval);
    if(!writeInterval.empty()) {
      if(rank == 0) std::cout << "  set writeInterval = " << writeInterval << "\n";
      platform->options.setArgs("SOLUTION OUTPUT INTERVAL", writeInterval);
    }

    delete[] rbuf;
  }
}

} // namespace

static void setOccaVars()
{
  std::string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));

  if (!getenv("OCCA_CACHE_DIR"))
    occa::env::OCCA_CACHE_DIR = cache_dir + "/occa/";

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_HOME"));

  if (!getenv("OCCA_DIR"))
    occa::env::OCCA_DIR = install_dir + "/";

  occa::env::OCCA_INSTALL_DIR = occa::env::OCCA_DIR;
}
