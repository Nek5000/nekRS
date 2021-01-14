#include "nrs.hpp"
#include "meshSetup.hpp"
#include "setup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "parReader.hpp"
#include "configReader.hpp"
#include "runTime.hpp"

static int rank, size;
static MPI_Comm comm;
static occa::device device;
static nrs_t* nrs;
static setupAide options;
static int ioStep;

static void setOccaVars(string dir);
static void setOUDF(setupAide &options);
static void dryRun(setupAide &options, int npTarget);

namespace nekrs
{
const double startTime(void)
{
  double val = 0;
  nrs->options.getArgs("START TIME", val);
  return val;
}

void setup(MPI_Comm comm_in, int buildOnly, int sizeTarget,
           int ciMode, string cacheDir, string _setupFile,
           string _backend, string _deviceID)
{
  MPI_Comm_dup(comm_in, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  configRead(comm);

  string setupFile = _setupFile + ".par";
  setOccaVars(cacheDir);

  if (rank == 0) {
#include "printHeader.inc"
    cout << "MPI tasks: " << size << endl << endl;
    cout << "using OCCA_DIR: " << occa::env::OCCA_DIR << endl;
    cout << "using OCCA_CACHE_DIR: " << occa::env::OCCA_CACHE_DIR << endl << endl;
  }

  nrs = new nrs_t();

  nrs->par = new inipp::Ini<char>();	   
  options = parRead((void*) nrs->par, setupFile, comm);

  if(buildOnly) 
    options.setArgs("BUILD ONLY", "TRUE");
  else
    options.setArgs("BUILD ONLY", "FALSE");
  if(!_backend.empty()) options.setArgs("THREAD MODEL", _backend);
  if(!_deviceID.empty()) options.setArgs("DEVICE NUMBER", _deviceID);

  setOUDF(options);

  // configure device
  device = occaDeviceConfig(options, comm);
  timer::init(comm, device, 0);

  if (buildOnly) {
    dryRun(options, sizeTarget);
    return;
  }

  timer::tic("setup", 1);

  // jit compile udf
  string udfFile;
  options.getArgs("UDF FILE", udfFile);
  if (!udfFile.empty()) {
    if(rank == 0) udfBuild(udfFile.c_str());
    MPI_Barrier(comm);
    udfLoad();
  }

  options.setArgs("CI-MODE", std::to_string(ciMode));
  if(rank == 0 && ciMode)
    cout << "enabling continous integration mode " << ciMode << "\n";

  if(udf.setup0) udf.setup0(comm, options);

  nrsSetup(comm, device, options, buildOnly, nrs);

  nrs->o_U.copyFrom(nrs->U);
  nrs->o_P.copyFrom(nrs->P);
  nrs->o_prop.copyFrom(nrs->prop);
  if(nrs->Nscalar) {
    nrs->cds->o_S.copyFrom(nrs->cds->S);
    nrs->cds->o_prop.copyFrom(nrs->cds->prop);
  }

  if(udf.properties) {
    occa::memory o_S = nrs->o_wrk0;
    occa::memory o_SProp = nrs->o_wrk0;
    if(nrs->Nscalar) {
      o_S = nrs->cds->o_S;
      o_SProp = nrs->cds->o_prop;
    }
    udf.properties(nrs, startTime(), nrs->o_U, o_S,
                   nrs->o_prop, o_SProp);
    nrs->o_prop.copyTo(nrs->prop);
    if(nrs->Nscalar) nrs->cds->o_prop.copyTo(nrs->cds->prop);
  }

  if(udf.executeStep) udf.executeStep(nrs, startTime(), 0);
  nek_ocopyFrom(startTime(), 0);

  timer::toc("setup");
  const double setupTime = timer::query("setup", "DEVICE:MAX");
  if(rank == 0) {
    cout << "\nsettings:\n" << endl << options << endl;
    size_t dMB = nrs->mesh->device.memoryAllocated() / 1e6;
    cout << "device memory allocation: " << dMB << " MB" << endl;
    cout << "initialization took " <<  setupTime << " s" << endl;
  }
  fflush(stdout);

  timer::reset();
  timer::set("setup", setupTime);
}

void runStep(double time, double dt, int tstep)
{
  runStep(nrs, time, dt, tstep);
}

void copyToNek(double time, int tstep)
{
  nek_ocopyFrom(time, tstep);
}

void udfExecuteStep(double time, int tstep, int isOutputStep)
{
  timer::tic("udfExecuteStep", 1);
  if (isOutputStep) {
    nek_ifoutfld(1);
    nrs->isOutputStep = 1;
  }

  if (udf.executeStep) udf.executeStep(nrs, time, tstep);

  nek_ifoutfld(0);
  nrs->isOutputStep = 0;
  timer::toc("udfExecuteStep");
}

void nekUserchk(void)
{
  nek_userchk();
}

const double dt(void)
{
  // TODO: adjust dt for target CFL
  return nrs->dt[0];
}

const double writeInterval(void)
{
  double val = -1;
  nrs->options.getArgs("SOLUTION OUTPUT INTERVAL", val);
  return val;
}

const int writeControlRunTime(void) 
{
  return nrs->options.compareArgs("SOLUTION OUTPUT CONTROL", "RUNTIME");
}

void outfld(double time, double outputTime)
{
  writeFld(nrs, time, 0);
}

const double endTime(void)
{
  double endTime = -1;
  nrs->options.getArgs("END TIME", endTime);
  return endTime;
}

const int numSteps(void)
{
  int numSteps = -1;
  nrs->options.getArgs("NUMBER TIMESTEPS", numSteps);
  return numSteps;
}

const int lastStep(double time, int tstep, double elapsedTime)
{
  if(!nrs->options.getArgs("STOP AT ELAPSED TIME").empty()) {
    double maxElaspedTime;
    nrs->options.getArgs("STOP AT ELAPSED TIME", maxElaspedTime);
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
  return nek_ptr(id);
}

void* nrsPtr(void)
{
  return nrs;
}



void printRuntimeStatistics()
{
  timer::printRunStat();
}
} // namespace

static void dryRun(setupAide &options, int npTarget)
{
  if (rank == 0)
    cout << "performing dry-run for "
         << npTarget
         << " MPI ranks ...\n" << endl;

  options.setArgs("NP TARGET", std::to_string(npTarget));

  // jit compile udf
  string udfFile;
  options.getArgs("UDF FILE", udfFile);
  if (!udfFile.empty()) {
    if(rank == 0) udfBuild(udfFile.c_str());
    MPI_Barrier(comm);
    *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",0);
    *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
  }

  if(udf.setup0) udf.setup0(comm, options);

  // init solver
  nrsSetup(comm, device, options, 1, nrs);

  if (rank == 0) cout << "\nBuild successful." << endl;
}

static void setOUDF(setupAide &options)
{
  std::string oklFile;
  options.getArgs("UDF OKL FILE",oklFile);

  char buf[FILENAME_MAX];

  char* ptr = realpath(oklFile.c_str(), NULL);
  if(!ptr) {
    if (rank == 0) cout << "ERROR: Cannot find " << oklFile << "!\n";
    EXIT(1);
  }
  free(ptr);

  std::string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
  string casename;
  options.getArgs("CASENAME", casename);
  const string dataFileDir = cache_dir + "/udf/";
  const string dataFile = dataFileDir + "udf.okl";

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
  getcwd(buf, sizeof(buf));
  string cwd;
  cwd.assign(buf);

  string install_dir;
  install_dir.assign(getenv("NEKRS_HOME"));

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

  if (!getenv("OCCA_DIR"))
    occa::env::OCCA_DIR = install_dir + "/";

  occa::env::OCCA_INSTALL_DIR = occa::env::OCCA_DIR;
}
