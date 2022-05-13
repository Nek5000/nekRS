#include <stdlib.h>
#include "nrs.hpp"
#include "meshSetup.hpp"
#include "setup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "printHeader.hpp"
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
static MPI_Comm commg, comm;
static nrs_t* nrs;
static setupAide options;
static dfloat lastOutputTime = 0;
static int firstOutfld = 1;
static int enforceLastStep = 0;
static int enforceOutputStep = 0;

static void setOccaVars();

bool useNodeLocalCache(){
  int buildNodeLocal = 0;
  if (getenv("NEKRS_CACHE_LOCAL"))
    buildNodeLocal = std::stoi(getenv("NEKRS_CACHE_LOCAL"));
  return (buildNodeLocal > 0);
}

namespace nekrs
{
double startTime(void)
{
  double val = 0;
  platform->options.getArgs("START TIME", val);
  return val;
}

void setup(MPI_Comm commg_in, MPI_Comm comm_in, 
    	   int buildOnly, int commSizeTarget,
           int ciMode, std::string _setupFile,
           std::string _backend, std::string _deviceID,
           int debug)
{
  MPI_Comm_dup(commg_in, &commg);
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
    std::string dir = std::string(buf) + "/.cache";
    if(getenv("NEKRS_CACHE_DIR")) dir.assign(getenv("NEKRS_CACHE_DIR"));
    setenv("NEKRS_CACHE_DIR", dir.c_str(), 1);
  }

  setOccaVars();

  if (rank == 0) {
    std::string installDir;
    installDir.assign(getenv("NEKRS_HOME"));
    std::cout << std::endl;
    std::cout << "using NEKRS_HOME: " << installDir << std::endl;

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

  if(options.getArgs("THREAD MODEL").length() == 0) 
    options.setArgs("THREAD MODEL", getenv("NEKRS_OCCA_MODE_DEFAULT"));
  if(!_backend.empty()) options.setArgs("THREAD MODEL", _backend);
  if(!_deviceID.empty()) options.setArgs("DEVICE NUMBER", _deviceID);

  // setup device
  platform_t* _platform = platform_t::getInstance(options, commg, comm);
  platform = _platform;
  platform->par = par;

  if(debug) platform->options.setArgs("VERBOSE","TRUE");

  int buildRank = rank;
  const bool buildNodeLocal = useNodeLocalCache();
  if(buildNodeLocal)
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);    

  if(buildRank == 0) {
    std::string cache_dir;
    cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
    mkdir(cache_dir.c_str(), S_IRWXU);
    std::string udf_cache_dir = cache_dir + "/udf";
    mkdir(udf_cache_dir.c_str(), S_IRWXU);
  }

  oudfInit(options);

  // jit compile udf
  std::string udfFile;
  options.getArgs("UDF FILE", udfFile);
  if (!udfFile.empty()) {
    udfBuild(udfFile.c_str(), options);
    udfLoad();
  }

  options.setArgs("CI-MODE", std::to_string(ciMode));
  if(rank == 0 && ciMode)
    std::cout << "enabling continous integration mode " << ciMode << "\n";

  nek::bootstrap();

  if(udf.setup0) udf.setup0(comm, options);

  {
    int overlap = 1;
    if(options.compareArgs("GS OVERLAP", "FALSE")) overlap = 0;
    oogs::overlap(overlap);
  }

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
      if(rank == 0) std::cout << "\nBuild successful." << std::endl;
    }
    return;
  }

  platform->linAlg = linAlg_t::getInstance();

  nrs = new nrs_t();

  {
    int result = 0;
    MPI_Comm_compare(commg, comm, &result);

    nrs->multiSession = (result == MPI_UNEQUAL);
  }

  nrsSetup(comm, options, nrs);

  const double setupTime = platform->timer.query("setup", "DEVICE:MAX");
  if(rank == 0) {
    std::cout << "\nsettings:\n" << std::endl << options << std::endl;
    std::cout << "occa memory usage: " << platform->device.occaDevice().memoryAllocated() / 1e9 << " GB"
              << std::endl;
  }
  fflush(stdout);

  platform->flopCounter->clear();
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
    if(tstep == 1) {
      double initialDt = 0.0;
      platform->options.getArgs("DT", initialDt);
      if(initialDt > 0.0){
        nrs->dt[0] = initialDt;
        return nrs->dt[0];
      }
    }
    const double dtOld = nrs->dt[0];
    timeStepper::adjustDt(nrs, tstep);

    // limit dt to control introduced error
    if(tstep > 1) nrs->dt[0] = std::min(nrs->dt[0], 1.25*dtOld);
    double maxDt = 0;
    platform->options.getArgs("MAX DT", maxDt);
    if(maxDt > 0) nrs->dt[0] = std::min(nrs->dt[0], maxDt);
  }
  
  if(nrs->dt[0] < 1e-10 || std::isnan(nrs->dt[0]) || std::isinf(nrs->dt[0])) {
    if(platform->comm.mpiRank == 0) printf("Invalid time step size %.2e\n", nrs->dt[0]);
    ABORT(EXIT_FAILURE);
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

void outfld(double time, std::string suffix)
{
  std::string oldValue;
  platform->options.getArgs("CHECKPOINT OUTPUT MESH", oldValue);

  if(firstOutfld)
    platform->options.setArgs("CHECKPOINT OUTPUT MESH", "TRUE");

  if(platform->options.compareArgs("MOVING MESH", "TRUE"))
    platform->options.setArgs("CHECKPOINT OUTPUT MESH", "TRUE");

  writeFld(nrs, time, suffix);
  lastOutputTime = time;
  firstOutfld = 0;

  platform->options.setArgs("CHECKPOINT OUTPUT MESH", oldValue);
}

void outfld(double time)
{
  outfld(time, "");
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
  if(options.compareArgs("BUILD ONLY", "FALSE")) {
    AMGXfree();
    nek::end();
  }
}

int runTimeStatFreq()
{
  int freq = 500;
  platform->options.getArgs("RUNTIME STATISTICS FREQUENCY", freq);
  return freq;
}

void printRuntimeStatistics(int step) { platform->timer.printRunStat(step); }

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

void printInfo(double time, int tstep) { timeStepper::printInfo(nrs, time, tstep); }

void verboseInfo(bool enabled)
{
  platform->options.setArgs("VERBOSE SOLVER INFO", "FALSE");
  if(enabled) platform->options.setArgs("VERBOSE SOLVER INFO", "TRUE");
}

void updateTimer(const std::string &key, double time) { platform->timer.set(key, time); }

} // namespace

static void setOccaVars()
{
  std::string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));

  if (!getenv("OCCA_CACHE_DIR")) {
    const std::string path= cache_dir + "/occa/";
    occa::env::OCCA_CACHE_DIR = path;
    setenv("OCCA_CACHE_DIR", path.c_str(), 1);
  }

  std::string installDir;
  installDir.assign(getenv("NEKRS_HOME"));

  if (!getenv("OCCA_DIR"))
    occa::env::OCCA_DIR = installDir + "/";

  occa::env::OCCA_INSTALL_DIR = occa::env::OCCA_DIR;
}
