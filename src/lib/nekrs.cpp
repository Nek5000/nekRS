#include "nrs.hpp"
#include <stdlib.h>
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
#include "amgx.h"

// extern variable from nrssys.hpp
platform_t* platform;

static int rank, size;
static MPI_Comm comm;
static nrs_t* nrs;
static setupAide options;
static dfloat lastOutputTime = 0;

static void setOccaVars(string dir);
static void setOUDF(setupAide &options);
static void dryRun(setupAide &options, int npTarget);

void printHeader()
{
  cout << R"(                 __    ____  _____)" << endl
       << R"(   ____   ___   / /__ / __ \/ ___/)" << endl
       << R"(  / __ \ / _ \ / //_// /_/ /\__ \ )" << endl
       << R"( / / / //  __// ,<  / _, _/___/ / )" << endl
       << R"(/_/ /_/ \___//_/|_|/_/ |_|/____/  )"
       << "v" << NEKRS_VERSION << "." << NEKRS_SUBVERSION
       << " (" << GITCOMMITHASH << ")" << endl
       << endl
       << "COPYRIGHT (c) 2019-2021 UCHICAGO ARGONNE, LLC" << endl
       << endl;
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
           int ciMode, string cacheDir, string _setupFile,
           string _backend, string _deviceID,
           neknek_t *neknek)
{
  if(buildOnly) {
    int rank, size;
    MPI_Comm_rank(comm_in, &rank);
    MPI_Comm_size(comm_in, &size);
    int color = MPI_UNDEFINED;
    if (rank == 0) color = 1;
    MPI_Comm_split(comm_in, color, 0, &comm);
    if (rank != 0) return;
  } else {
    MPI_Comm_dup(comm_in, &comm);
  }

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  srand48((long int) rank);

  configRead(comm);
  oogs::gpu_mpi(std::stoi(getenv("NEKRS_GPU_MPI")));
  setOccaVars(cacheDir);

  if (rank == 0) {
    printHeader();
    cout << "MPI tasks: " << size << endl << endl;
    string install_dir;
    install_dir.assign(getenv("NEKRS_HOME"));
    cout << "using NEKRS_HOME: " << install_dir << endl;
    cout << "using OCCA_CACHE_DIR: " << occa::env::OCCA_CACHE_DIR << endl << endl;
  }

  nrs = new nrs_t();
  nrs->neknek = neknek;

  nrs->par = new inipp::Ini<char>();
  string setupFile = _setupFile + ".par";
  options = parRead((void*) nrs->par, setupFile, comm);

  options.setArgs("BUILD ONLY", "FALSE");
  if(buildOnly) options.setArgs("BUILD ONLY", "TRUE");

  if (options.getArgs("THREAD MODEL").length() == 0)
    options.setArgs("THREAD MODEL", getenv("NEKRS_OCCA_MODE_DEFAULT"));
  if(!_backend.empty()) options.setArgs("THREAD MODEL", _backend);
  if(!_deviceID.empty()) options.setArgs("DEVICE NUMBER", _deviceID);

  setOUDF(options);

  // setup device
  platform_t* _platform = platform_t::getInstance(options, comm);
  platform = _platform;

  if (buildOnly) {
    dryRun(options, commSizeTarget);
    return;
  }

  platform->timer.tic("setup", 1);

  platform->linAlg = linAlg_t::getInstance();

  // jit compile udf
  string udfFile;
  options.getArgs("UDF FILE", udfFile);
  string casename;
  options.getArgs("CASENAME", casename);
  if (!udfFile.empty()) {
    int err = 0;
    if(rank == 0) err = udfBuild(casename.c_str(), udfFile.c_str(), 0);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, comm);
    if(err) ABORT(EXIT_FAILURE);;
    udfLoad(casename.c_str());
  }

  options.setArgs("CI-MODE", std::to_string(ciMode));
  if(rank == 0 && ciMode)
    cout << "enabling continous integration mode " << ciMode << "\n";

  if(udf.setup0) udf.setup0(comm, options);

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

  neknekSetup(nrs);

  nek::ocopyToNek(startTime(), 0);

  platform->timer.toc("setup");
  const double setupTime = platform->timer.query("setup", "DEVICE:MAX");
  if(rank == 0) {
    cout << "\nsettings:\n" << endl << options << endl;
    cout << "occa memory usage: " << platform->device.memoryAllocated()/1e9 << " GB" << endl;
    cout << "initialization took " << setupTime << " s" << endl;
  }
  fflush(stdout);

  platform->timer.reset();
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

double dt(void)
{
  // TODO: adjust dt for target CFL
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
  return outputStep;
}

void outputStep(int val)
{
  nrs->isOutputStep = val;
}

void outfld(double time)
{
  writeFld(nrs, time, 0);
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

void printRuntimeStatistics()
{
  platform_t* platform = platform_t::getInstance(options, comm);
  platform->timer.printRunStat();
}
} // namespace

static void dryRun(setupAide &options, int npTarget)
{
  cout << "performing dry-run to jit-compile for >="
       << npTarget
       << " MPI tasks ...\n" << endl;
  fflush(stdout);

  options.setArgs("NP TARGET", std::to_string(npTarget));
  options.setArgs("BUILD ONLY", "TRUE");

  platform->linAlg = linAlg_t::getInstance();

  // jit compile udf
  string udfFile;
  options.getArgs("UDF FILE", udfFile);
  string casename;
  options.getArgs("CASENAME", casename);
  if (!udfFile.empty()) {
    int err = 0;
    if(rank == 0) err = udfBuild(casename.c_str(), udfFile.c_str(), 1);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, comm);
    if(err) ABORT(EXIT_FAILURE);
    MPI_Barrier(comm);
    *(void**)(&udf.loadKernels) = udfLoadFunction(casename.c_str(), "UDF_LoadKernels",0);
    *(void**)(&udf.setup0) = udfLoadFunction(casename.c_str(), "UDF_Setup0",0);
  }

  if(udf.setup0) udf.setup0(comm, options);

  // init solver
  platform_t* platform = platform_t::getInstance();
  nrsSetup(comm, options, nrs);
  neknekSetup(nrs);

  cout << "\nBuild successful." << endl;
}

static void setOUDF(setupAide &options)
{
  std::string oklFile;
  options.getArgs("UDF OKL FILE",oklFile);

  // char buf[FILENAME_MAX];

  char* ptr = realpath(oklFile.c_str(), NULL);
  if(!ptr) {
    if (rank == 0) cout << "ERROR: Cannot find " << oklFile << "!\n";
    ABORT(EXIT_FAILURE);;
  }
  free(ptr);

  std::string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
  string casename;
  options.getArgs("CASENAME", casename);
  const string dataFileDir = cache_dir + "/udf/";
  const string dataFile = dataFileDir + "udf-" + casename + ".okl";

  if (rank == 0) {
    mkdir(dataFileDir.c_str(), S_IRWXU);

    std::ifstream in;
    in.open(oklFile);
    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();

    std::ofstream out;
    out.open(dataFile, std::ios::trunc);

    out << buffer.str();

    std::size_t found;
    found = buffer.str().find("void nrsVelocityDirichletConditions");
    if (found == std::string::npos) found = buffer.str().find("void insVelocityDirichletConditions");
    if (found == std::string::npos) found = buffer.str().find("void velocityDirichletConditions");
    if (found == std::string::npos)
      out << "void velocityDirichletConditions(bcData *bc){}\n";

    found = buffer.str().find("void nrsVelocityNeumannConditions");
    if (found == std::string::npos) found = buffer.str().find("void insVelocityNeumannConditions");
    if (found == std::string::npos) found = buffer.str().find("void velocityNeumannConditions");
    if (found == std::string::npos)
      out << "void velocityNeumannConditions(bcData *bc){}\n";

    found = buffer.str().find("void nrsPressureDirichletConditions");
    if (found == std::string::npos) found = buffer.str().find("void insPressureDirichletConditions");
    if (found == std::string::npos) found = buffer.str().find("void pressureDirichletConditions");
    if (found == std::string::npos)
      out << "void pressureDirichletConditions(bcData *bc){}\n";

    found = buffer.str().find("void cdsNeumannConditions");
    if (found == std::string::npos) found = buffer.str().find("void scalarNeumannConditions");
    if (found == std::string::npos)
      out << "void scalarNeumannConditions(bcData *bc){}\n";

    found = buffer.str().find("void cdsDirichletConditions");
    if (found == std::string::npos) found = buffer.str().find("void scalarDirichletConditions");
    if (found == std::string::npos)
      out << "void scalarDirichletConditions(bcData *bc){}\n";

    out <<
      "@kernel void __dummy__(int N) {"
      "  for (int i = 0; i < N; ++i; @tile(16, @outer, @inner)) {}"
      "}";

    out.close();
  }

  options.setArgs("DATA FILE", dataFile);
}

static void setOccaVars(string dir)
{
  char buf[FILENAME_MAX];
  char * ret = getcwd(buf, sizeof(buf));
  if(!ret) ABORT(EXIT_FAILURE);;
  string cwd;
  cwd.assign(buf);

  if (dir.empty())
    sprintf(buf,"%s/.cache", cwd.c_str());
  else
    sprintf(buf,"%s/%s", cwd.c_str(), dir.c_str());

  setenv("NEKRS_CACHE_DIR", buf, 1);
  string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
  if (rank == 0) mkdir(cache_dir.c_str(), S_IRWXU);
  MPI_Barrier(comm);

  if (!getenv("OCCA_CACHE_DIR"))
    occa::env::OCCA_CACHE_DIR = cache_dir + "/occa/";

  string install_dir;
  install_dir.assign(getenv("NEKRS_HOME"));

  if (!getenv("OCCA_DIR"))
    occa::env::OCCA_DIR = install_dir + "/";

  occa::env::OCCA_INSTALL_DIR = occa::env::OCCA_DIR;
}
