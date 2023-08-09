#include <stdlib.h>
#include <filesystem>
#include <functional>
#include "nrs.hpp"
#include "setup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "printHeader.hpp"
#include "udf.hpp"
#include "bcMap.hpp"
#include "parReader.hpp"
#include "re2Reader.hpp"
#include "configReader.hpp"
#include "re2Reader.hpp"
#include "timeStepper.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "cfl.hpp"
#include "AMGX.hpp"
#include "hypreWrapper.hpp"
#include "hypreWrapperDevice.hpp"

namespace fs = std::filesystem;

#ifdef ENABLE_SENSEI
#include "sensei/Bridge.h"
#endif

// extern variable from nrssys.hpp
platform_t *platform;

static nrs_t *nrs;
static setupAide options;

static int rank, size;
static MPI_Comm commg, comm;

static dfloat lastOutputTime = 0;
static int firstOutfld = 1;
static int enforceLastStep = 0;
static int enforceOutputStep = 0;
static bool initialized = false;

namespace nekrs {

#ifdef ENABLE_SENSEI
dfloat* vel_x;
dfloat* vel_y;
dfloat* vel_z;
dfloat* temp;
dfloat* pr;
#endif

void reset()
{
  lastOutputTime = 0;
  firstOutfld = 1;
  enforceLastStep = 0;
  enforceOutputStep = 0;
}

double startTime(void)
{
  double val = 0;
  platform->options.getArgs("START TIME", val);
  return val;
}

void setup(MPI_Comm commg_in,
           MPI_Comm comm_in,
           int buildOnly,
           int commSizeTarget,
           int ciMode,
           std::string _setupFile,
           std::string _backend,
           std::string _deviceID,
           int nSessions,
           int sessionID,
           int debug)
{
  nrsCheck(initialized, comm_in, EXIT_FAILURE, "%s\n", "Calling setup twice is erroneous!");

  commg = commg_in;
  comm = comm_in;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (rank == 0) {
    printHeader();
    std::cout << "MPI tasks: " << size << std::endl << std::endl;
  }

  configRead(comm);

  if(nSessions > 1) {
    options.setArgs("NEKNEK NUMBER OF SESSIONS", std::to_string(nSessions));
    options.setArgs("NEKNEK SESSION ID", std::to_string(sessionID));
  }

  options.setArgs("BUILD ONLY", "FALSE");
  if (buildOnly) {
    options.setArgs("BUILD ONLY", "TRUE");
    options.setArgs("NP TARGET", std::to_string(commSizeTarget));
    if (rank == 0) {
      std::cout << "jit-compiling for >=" << commSizeTarget << " MPI tasks ...\n" << std::endl;
    }
    fflush(stdout);
  }

  auto par = new inipp::Ini();
  if (rank == 0)
    std::cout << "reading par file ...\n";
  parRead(par, _setupFile + ".par", comm, options);

  // precedence: cmd arg, par, env-var
  if (options.getArgs("THREAD MODEL").length() == 0)
    options.setArgs("THREAD MODEL", getenv("NEKRS_OCCA_MODE_DEFAULT"));
  if (!_backend.empty())
    options.setArgs("THREAD MODEL", _backend);
  if (!_deviceID.empty())
    options.setArgs("DEVICE NUMBER", _deviceID);

  // setup platform (requires THREAD MODEL)
  platform_t *_platform = platform_t::getInstance(options, commg, comm);
  platform = _platform;
  platform->par = par;

  if (debug)
    platform->options.setArgs("VERBOSE", "TRUE");

  int buildRank = rank;
  if (platform->cacheLocal)
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);

  if (rank == 0) {
    std::cout << "using NEKRS_HOME: " << getenv("NEKRS_HOME") << std::endl;
    std::cout << "using NEKRS_CACHE_DIR: " << getenv("NEKRS_CACHE_DIR") << std::endl;
    std::cout << "using OCCA_CACHE_DIR: " << occa::env::OCCA_CACHE_DIR << std::endl << std::endl;
  }

  options.setArgs("CI-MODE", std::to_string(ciMode));
  if (rank == 0 && ciMode)
    std::cout << "enabling continous integration mode " << ciMode << "\n";

  {
    int nelgt, nelgv;
    re2::nelg(options.getArgs("MESH FILE"), nelgt, nelgv, comm);
    nrsCheck(size > nelgv, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "MPI tasks > number of elements!");
  }

  bcMap::setup();

  nek::bootstrap();

  // jit compile udf
  std::string udfFile;
  options.getArgs("UDF FILE", udfFile);
  if (!udfFile.empty()) {
    udfBuild(udfFile, options);
    udfLoad();
  }

  // here we might access some nek variables
  if (udf.setup0)
    udf.setup0(comm, options);

  if (rank == 0) {
    if (!buildOnly && options.compareArgs("STDOUT PAR", "TRUE"))
      parEcho();
    if (!buildOnly && options.compareArgs("STDOUT UDF", "TRUE"))
      udfEcho();
  }

  compileKernels();

  oogs::overlap(options.compareArgs("ENABLE GS COMM OVERLAP", "FALSE") ? 0 : 1);

  if (buildOnly) {
    MPI_Barrier(platform->comm.mpiComm);
    if (buildRank == 0) {
      std::string cache_dir;
      cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
      std::string file = cache_dir + "/build-only.timestamp";
      remove(file.c_str());
      std::ofstream ofs;
      ofs.open(file, std::ofstream::out);
      ofs.close();
      if (rank == 0)
        std::cout << "\nBuild successful." << std::endl;
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
  if (neknekCoupled()) {
    new neknek_t(nrs, nSessions, sessionID);
  }

  const double setupTime = platform->timer.query("setup", "DEVICE:MAX");
  if (rank == 0) {
    std::cout << "\nsettings:\n" << std::endl << options << std::endl;
    std::cout << "occa memory usage: " << platform->device.occaDevice().memoryAllocated() / 1e9 << " GB"
              << std::endl;
  }
  fflush(stdout);

#ifdef ENABLE_SENSEI
  std::string casename;
  options.getArgs("CASENAME", casename);
  
  dfloat* vx = nrs->U + 0 * nrs->fieldOffset;
  dfloat* vy = nrs->U + 1 * nrs->fieldOffset;
  dfloat* vz = nrs->U + 2 * nrs->fieldOffset;


  mesh_t *mesh = nrs->meshV;

  double mi = 0.0f;
  double ma = 100.0f;
  double x_min = -0.5;
  double x_max = 0.5f;
  double z_min = 0.0f;
  double z_max = 20.0f;
  
  dfloat global_min_mesh_x = 0, global_min_mesh_y = 0, global_min_mesh_z = 0, global_min_vx = 0, global_min_vy = 0, global_min_vz = 0, global_min_p = 0, global_min_t = 0;
  dfloat global_max_mesh_x = 0, global_max_mesh_y = 0, global_max_mesh_z = 0, global_max_vx = 0, global_max_vy = 0, global_max_vz = 0, global_max_p = 0, global_max_t = 0;
  
  dfloat local_min_mesh_x = *std::min_element(mesh->x, mesh->x + nrs->fieldOffset);
  dfloat local_min_mesh_y = *std::min_element(mesh->y, mesh->y + nrs->fieldOffset);
  dfloat local_min_mesh_z = *std::min_element(mesh->z, mesh->z + nrs->fieldOffset);

  dfloat local_min_vx = *std::min_element(vx, vx + nrs->fieldOffset);
  dfloat local_min_vy = *std::min_element(vy, vy + nrs->fieldOffset);
  dfloat local_min_vz = *std::min_element(vz, vz + nrs->fieldOffset);

  dfloat local_min_p = *std::min_element(nrs->P, nrs->P + nrs->fieldOffset);
  dfloat local_min_t = nrs->cds ? *std::min_element(nrs->cds->S, nrs->cds->S + nrs->fieldOffset) : 0;

  dfloat local_max_mesh_x = *std::max_element(mesh->x, mesh->x + nrs->fieldOffset);
  dfloat local_max_mesh_y = *std::max_element(mesh->y, mesh->y + nrs->fieldOffset);
  dfloat local_max_mesh_z = *std::max_element(mesh->z, mesh->z + nrs->fieldOffset);

  dfloat local_max_vx = *std::max_element(vx, vx + nrs->fieldOffset);
  dfloat local_max_vy = *std::max_element(vy, vy + nrs->fieldOffset);
  dfloat local_max_vz = *std::max_element(vz, vz + nrs->fieldOffset);

  dfloat local_max_p = *std::max_element(nrs->P, nrs->P + nrs->fieldOffset);
  dfloat local_max_t = nrs->cds ? *std::max_element(nrs->cds->S, nrs->cds->S + nrs->fieldOffset) : 0;

  MPI_Allreduce(&local_min_vx, &global_min_vx, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
  MPI_Allreduce(&local_min_vy, &global_min_vy, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
  MPI_Allreduce(&local_min_vz, &global_min_vz, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
 
  MPI_Allreduce(&local_min_mesh_x, &global_min_mesh_x, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
  MPI_Allreduce(&local_min_mesh_y, &global_min_mesh_y, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
  MPI_Allreduce(&local_min_mesh_z, &global_min_mesh_z, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);

  MPI_Allreduce(&local_min_p, &global_min_p, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
  MPI_Allreduce(&local_min_t, &global_min_t, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);

  MPI_Allreduce(&local_max_vx, &global_max_vx, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
  MPI_Allreduce(&local_max_vy, &global_max_vy, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
  MPI_Allreduce(&local_max_vz, &global_max_vz, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
 
  MPI_Allreduce(&local_max_mesh_x, &global_max_mesh_x, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
  MPI_Allreduce(&local_max_mesh_y, &global_max_mesh_y, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
  MPI_Allreduce(&local_max_mesh_z, &global_max_mesh_z, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);

  MPI_Allreduce(&local_max_p, &global_min_p, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
  MPI_Allreduce(&local_max_t, &global_min_t, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);

	dlong Nlocal = mesh->Nelements * mesh->Np;
	vel_x = (dfloat*)malloc(Nlocal * sizeof(dfloat));
	vel_y = (dfloat*)malloc(Nlocal * sizeof(dfloat));
	vel_z = (dfloat*)malloc(Nlocal * sizeof(dfloat));

	temp = (dfloat*)malloc(nrs->Nscalar * Nlocal * sizeof(dfloat));
	pr = (dfloat*)malloc(Nlocal * sizeof(dfloat));

  sensei_bridge_initialize(&comm, (char*)casename.c_str(), &mesh->Nelements,
      &mesh->Nq, &mesh->Nq, &mesh->Nq, mesh->x, mesh->y, mesh->z,
      &mesh->Nq, &mesh->Nq, &mesh->Nq, mesh->x, mesh->y, mesh->z,
      vel_x, vel_y, vel_z,
      nullptr, nullptr, nullptr,
      pr, temp, nullptr, &size,
      &global_min_mesh_x, &global_max_mesh_x, &global_min_mesh_y, &global_max_mesh_y, &global_min_mesh_z, &global_max_mesh_z,
      &global_min_vx, &global_max_vx, &global_min_vy, &global_max_vy, &global_min_vz, &global_max_vz,
      &mi, &ma, &mi, &ma, &mi, &ma, //double* vort_x_min, double* vort_x_max, double* vort_y_min, double* vort_y_max, double* vort_z_min, double* vort_z_max,
      &global_min_p, &global_max_p, &global_min_t, &global_max_t, &mi, &ma); //double* pr_min, double* pr_max, double* temp_min, double* temp_max, double* jac_min, double* jac_max);

#endif

  platform->flopCounter->clear();

#if 1
  if (platform->cacheBcast) {
    MPI_Barrier(platform->comm.mpiComm);

    int rankLocal;
    MPI_Comm_rank(platform->comm.mpiCommLocal, &rankLocal);
    if (rankLocal == 0) {
      for (auto &entry : std::filesystem::directory_iterator(platform->tmpDir))
        fs::remove_all(entry.path());
    }
  }
#endif

  initialized = true;
}

void copyFromNek(double time, int tstep) { nek::ocopyToNek(time, tstep); }

void udfExecuteStep(double time, int tstep, int isOutputStep)
{
  platform->timer.tic("udfExecuteStep", 1);
  if (isOutputStep) {
    nek::ifoutfld(1);
    nrs->isOutputStep = 1;
  }

  if (udf.executeStep)
    udf.executeStep(nrs, time, tstep);

  nek::ifoutfld(0);
  nrs->isOutputStep = 0;
  platform->timer.toc("udfExecuteStep");
}

void nekUserchk(void) { nek::userchk(); }

double dt(int tstep)
{
  if (platform->options.compareArgs("VARIABLE DT", "TRUE")) {
    if (tstep == 1) {
      double initialDt = 0.0;
      platform->options.getArgs("DT", initialDt);
      if (initialDt > 0.0) {
        nrs->dt[0] = initialDt;
        return nrs->dt[0];
      }
    }
    const double dtOld = nrs->dt[0];
    timeStepper::adjustDt(nrs, tstep);

    double maxDt = 0;
    platform->options.getArgs("MAX DT", maxDt);
    if (maxDt > 0)
      nrs->dt[0] = std::min(nrs->dt[0], maxDt);
  }

  nrsCheck(nrs->dt[0] < 1e-10 || std::isnan(nrs->dt[0]) || std::isinf(nrs->dt[0]),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "Invalid time step size %.2e\n",
           nrs->dt[0]);

  // during a neknek simulation, sync dt across all ranks
  if (nrs->neknek) {
    MPI_Allreduce(MPI_IN_PLACE, &nrs->dt[0], 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiCommParent);
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
  return platform->options.compareArgs("SOLUTION OUTPUT CONTROL", "SIMULATIONTIME");
}

int outputStep(double time, int tStep)
{
  int outputStep = 0;
  if (writeControlRunTime()) {
    double val;
    platform->options.getArgs("START TIME", val);
    if (lastOutputTime == 0 && val > 0)
      lastOutputTime = val;
    outputStep = ((time - lastOutputTime) + 1e-10) > nekrs::writeInterval();
  }
  else {
    if (writeInterval() > 0)
      outputStep = (tStep % (int)writeInterval() == 0);
  }

  if (enforceOutputStep) {
    enforceOutputStep = 0;
    return 1;
  }
  return outputStep;
}

void outputStep(int val) { nrs->isOutputStep = val; }

void outfld(double time, int step, std::string suffix)
{
  std::string oldValue;
  platform->options.getArgs("CHECKPOINT OUTPUT MESH", oldValue);

  if (firstOutfld)
    platform->options.setArgs("CHECKPOINT OUTPUT MESH", "TRUE");

  if (platform->options.compareArgs("MOVING MESH", "TRUE"))
    platform->options.setArgs("CHECKPOINT OUTPUT MESH", "TRUE");

  writeFld(nrs, time, step, suffix);
  lastOutputTime = time;
  firstOutfld = 0;

  platform->options.setArgs("CHECKPOINT OUTPUT MESH", oldValue);
}

void outfld(double time, int step) { outfld(time, step, ""); }

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

void lastStep(int val) { nrs->lastStep = val; }

int lastStep(double time, int tstep, double elapsedTime)
{
  if (!platform->options.getArgs("STOP AT ELAPSED TIME").empty()) {
    double maxElaspedTime;
    platform->options.getArgs("STOP AT ELAPSED TIME", maxElaspedTime);
    if (elapsedTime > 60.0 * maxElaspedTime)
      nrs->lastStep = 1;
  }
  else if (endTime() >= 0) {
    const double eps = 1e-12;
    nrs->lastStep = fabs((time + nrs->dt[0]) - endTime()) < eps || (time + nrs->dt[0]) > endTime();
  }
  else {
    nrs->lastStep = (tstep == numSteps());
  }

  if (enforceLastStep)
    return 1;
  return nrs->lastStep;
}

void *nekPtr(const char *id) { return nek::ptr(id); }

void *nrsPtr(void) { return nrs; }

int finalize(void) { 
#ifdef ENABLE_SENSEI
	free(vel_x);
	free(vel_y);
	free(vel_z);
	free(temp);
	free(pr);
#endif
  return nrsFinalize(nrs); 
}

int runTimeStatFreq()
{
  int freq = 500;
  platform->options.getArgs("RUNTIME STATISTICS FREQUENCY", freq);
  return freq;
}

int printInfoFreq()
{
  int freq = 1;
  platform->options.getArgs("PRINT INFO FREQUENCY", freq);
  return freq;
}

int updateFileCheckFreq()
{
  int freq = 20;
  platform->options.getArgs("UPDATE FILE CHECK FREQUENCY", freq);
  return freq;
}

void printRuntimeStatistics(int step) { platform->timer.printRunStat(step); }

void processUpdFile()
{
  char *rbuf = nullptr;
  long long int fsize = 0;
  const std::string updFile = "nekrs.upd";

  if (rank == 0) {
    if (fs::exists(updFile)) {
      FILE *f = fopen(updFile.c_str(), "r");
      fseek(f, 0, SEEK_END);
      fsize = ftell(f);
      fseek(f, 0, SEEK_SET);
      rbuf = new char[fsize];
      fread(rbuf, 1, fsize, f);
      fclose(f);
      remove(updFile.c_str());
    }
  }

  MPI_Bcast(&fsize, 1, MPI_LONG_LONG_INT, 0, comm);

  if (fsize) {
    exit(1);
    if (rank == 0)
      std::cout << "processing " << updFile << " ...\n";

    if (rank != 0)
      rbuf = new char[fsize];
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
    if (checkpoint == "true")
      enforceOutputStep = 1;

    std::string endTime;
    ini.extract("general", "endtime", endTime);
    if (!endTime.empty()) {
      if (rank == 0)
        std::cout << "  set endTime = " << endTime << "\n";
      platform->options.setArgs("END TIME", endTime);
    }

    std::string numSteps;
    ini.extract("general", "numsteps", numSteps);
    if (!numSteps.empty()) {
      if (rank == 0)
        std::cout << "  set numSteps = " << numSteps << "\n";
      platform->options.setArgs("NUMBER TIMESTEPS", numSteps);
    }

    std::string writeInterval;
    ini.extract("general", "writeinterval", writeInterval);
    if (!writeInterval.empty()) {
      if (rank == 0)
        std::cout << "  set writeInterval = " << writeInterval << "\n";
      platform->options.setArgs("SOLUTION OUTPUT INTERVAL", writeInterval);
    }

    delete[] rbuf;
  }
}

void printInfo(double time, int tstep, bool printStepInfo, bool printVerboseInfo)
{
  timeStepper::printInfo(nrs, time, tstep, printStepInfo, printVerboseInfo);
}

void verboseInfo(bool enabled)
{
  platform->options.setArgs("VERBOSE SOLVER INFO", "FALSE");
  if (enabled)
    platform->options.setArgs("VERBOSE SOLVER INFO", "TRUE");
}

void updateTimer(const std::string &key, double time) { platform->timer.set(key, time); }

void resetTimer(const std::string &key) { platform->timer.reset(key); }

int exitValue() { return platform->exitValue; }

void initStep(double time, double dt, int tstep) { timeStepper::initStep(nrs, time, dt, tstep); }

bool runStep(std::function<bool(int)> convergenceCheck, int corrector)
{
  return timeStepper::runStep(nrs, convergenceCheck, corrector);
}

bool runStep(int corrector)
{

  auto _nrs = &nrs;
  auto _udf = &udf;

  std::function<bool(int)> convergenceCheck = [](int corrector) -> bool {
    if (udf.timeStepConverged)
      return udf.timeStepConverged(nrs, corrector);
    else
      return true;
  };

  return timeStepper::runStep(nrs, convergenceCheck, corrector);
}

double finishStep()
{
  timeStepper::finishStep(nrs);
  return nrs->timePrevious + nrs->dt[0];
}

bool stepConverged() { return nrs->timeStepConverged; }

#ifdef ENABLE_SENSEI
void runSensei(double time, double dt, int tstep) {
	sensei::DataAdaptor *daOut = nullptr;

  // Copy to CPU
  // Code from nekInterface/nekInterfaceAdapter.cpp
	mesh_t* mesh = nrs->meshV;
	
	dlong Nlocal = mesh->Nelements * mesh->Np;
	occa::memory o_u, o_p, o_s;
	void* o_uu = &nrs->o_U;
	void* o_pp = &nrs->o_P;
	void* o_ss = &nrs->cds->o_S;
	if(o_uu) o_u = *((occa::memory *) o_uu);
	if(o_pp) o_p = *((occa::memory *) o_pp);
	if(o_ss && nrs->Nscalar) o_s = *((occa::memory *) o_ss);
	
	if(o_u.ptr()) {
    occa::memory o_vx = o_u + 0 * nrs->fieldOffset * sizeof(dfloat);
    occa::memory o_vy = o_u + 1 * nrs->fieldOffset * sizeof(dfloat);
    occa::memory o_vz = o_u + 2 * nrs->fieldOffset * sizeof(dfloat);
    o_vx.copyTo(vel_x, Nlocal * sizeof(dfloat));
		o_vy.copyTo(vel_y, Nlocal * sizeof(dfloat));
		o_vz.copyTo(vel_z, Nlocal * sizeof(dfloat));
  }
  if(nrs->o_P.ptr()) {
    o_p.copyTo(pr, Nlocal * sizeof(dfloat));
  }
  if(nrs->cds) {
    if(nrs->cds->o_S.ptr()) {
      const dlong nekFieldOffset = nekData.lelt * mesh->Np;
      int NSfields = nrs->Nscalar;
      for(int is = 0; is < NSfields; is++) {
        mesh_t* mesh = nrs->meshV;
        if(nrs->cds)
          (is) ? mesh = nrs->meshV: mesh = nrs->cds->mesh[0];
        const dlong Nlocal = mesh->Nelements * mesh->Np;
        dfloat* Ti = temp + is * nekFieldOffset;
        occa::memory o_Si = o_s + is * nrs->fieldOffset * sizeof(dfloat);
        o_Si.copyTo(Ti, Nlocal * sizeof(dfloat));
      }
    }
  }
	sensei_bridge_update(&tstep, &time, &daOut);
}
#endif
} // namespace nekrs

int nrsFinalize(nrs_t *nrs)
{
  auto exitValue = nekrs::exitValue();
  if (platform->options.compareArgs("BUILD ONLY", "FALSE")) {
    if (nrs->uSolver)
      delete nrs->uSolver;
    if (nrs->vSolver)
      delete nrs->vSolver;
    if (nrs->wSolver)
      delete nrs->wSolver;
    if (nrs->uvwSolver)
      delete nrs->uvwSolver;
    if (nrs->pSolver)
      delete nrs->pSolver;
    for (int is; is < nrs->Nscalar; is++) {
      if (nrs->cds->solver[is])
        delete nrs->cds->solver[is];
    }
    if (nrs->cvode)
      delete nrs->cvode;
    if (nrs->meshSolver)
      delete nrs->meshSolver;

    hypreWrapper::finalize();
    hypreWrapperDevice::finalize();
    AMGXfinalize();
    nek::finalize();
  }

  if (platform->comm.mpiRank == 0)
    std::cout << "finished with exit code " << exitValue << std::endl;

  return exitValue;
}


