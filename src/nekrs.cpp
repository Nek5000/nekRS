#include "nrs.hpp"
#include "meshSetup.hpp"
#include "insSetup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "parReader.hpp"
#include "configReader.hpp"
#include "runTime.hpp"

static int rank, size;
static MPI_Comm comm;
static occa::device device;
static ins_t* ins;
static libParanumal::setupAide options;
static int ioStep;

int nrsBuildOnly = 0;  // hack for meshPhysicalNodes()

static void setCache(string dir);
static void setOUDF(libParanumal::setupAide &options);
static void dryRun(libParanumal::setupAide &options, int npTarget);

namespace nekrs
{
void setup(MPI_Comm comm_in, int buildOnly, int sizeTarget,
           int ciMode, string cacheDir, string _setupFile,
           string _backend, string _deviceID)
{
  MPI_Comm_dup(comm_in, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  configRead(comm);

  nrsBuildOnly = buildOnly;
  string setupFile = _setupFile + ".par";
  setCache(cacheDir);

  if (rank == 0) {
#include "printHeader.inc"
    cout << "MPI tasks: " << size << endl << endl;
    cout << "using OCCA_CACHE_DIR: " << occa::env::OCCA_CACHE_DIR << endl << endl;
  }

  options = parRead(setupFile, comm);

  if(!_backend.empty()) options.setArgs("THREAD MODEL", _backend);
  if(!_deviceID.empty()) options.setArgs("DEVICE NUMBER", _deviceID);

  setOUDF(options);

  // configure device
  device = occaDeviceConfig(options, comm);
  timer::init(comm, device, 0);

  if (nrsBuildOnly) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    dryRun(options, sizeTarget);
    return;
  }

  MPI_Barrier(comm);
  double t0 = MPI_Wtime();

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

  int nscal;
  options.getArgs("NUMBER OF SCALARS", nscal);

  // jit compile nek
  int N;
  string casename;
  options.getArgs("CASENAME", casename);
  options.getArgs("POLYNOMIAL DEGREE", N);
  if(rank == 0) buildNekInterface(casename.c_str(), mymax(1,nscal), N, size);
  MPI_Barrier(comm);

  // init nek
  nek_setup(comm, options, &ins);
  nek_setic();
  nek_userchk();

  // init solver
  ins = insSetup(comm, device, options, nrsBuildOnly);

  // set initial condition
  int readRestartFile;
  options.getArgs("RESTART FROM FILE", readRestartFile);
  if(readRestartFile) nek_copyRestart();
  if(udf.setup) udf.setup(ins);
  if(options.compareArgs("VARIABLEPROPERTIES", "TRUE")) {
    if(!udf.properties) {
      if (rank ==
          0) cout << "ERROR: variableProperties requires assigned udf.properties pointer" << "!\n";
      EXIT(1);
    }
  }
  ins->o_U.copyFrom(ins->U);
  ins->o_P.copyFrom(ins->P);
  ins->o_prop.copyFrom(ins->prop);
  if(ins->Nscalar) {
    ins->cds->o_S.copyFrom(ins->cds->S);
    ins->cds->o_prop.copyFrom(ins->cds->prop);
  }

  if(udf.properties) {
    udf.properties(ins, ins->startTime, ins->o_U, ins->cds->o_S,
                   ins->o_prop, ins->cds->o_prop);
    ins->o_prop.copyTo(ins->prop);
    if(ins->Nscalar) ins->cds->o_prop.copyTo(ins->cds->prop);
  }

  if(udf.executeStep) udf.executeStep(ins, ins->startTime, 0);
  nek_ocopyFrom(ins->startTime, 0);

  if(rank == 0) {
    cout << "\nsettings:\n" << endl << options << endl;
    size_t dMB = ins->mesh->device.memoryAllocated() / 1e6;
    cout << "device memory allocation: " << dMB << " MB" << endl;
    cout << "initialization took " << MPI_Wtime() - t0 << " seconds" << endl;
  }
  fflush(stdout);

  timer::reset();
}

void runStep(double time, double dt, int tstep)
{
  runStep(ins, time, dt, tstep);
}

void copyToNek(double time, int tstep)
{
  nek_ocopyFrom(time, tstep);
}

void udfExecuteStep(double time, int tstep, int isOutputStep)
{
  if (isOutputStep) {
    nek_ifoutfld(1);
    ins->isOutputStep = 1;
  }

  if (udf.executeStep) udf.executeStep(ins, time, tstep);

  nek_ifoutfld(0);
  ins->isOutputStep = 0;
}

void nekUserchk(void)
{
  nek_userchk();
}

void nekOutfld(void)
{
  nek_outfld();
}

const double dt(void)
{
  return ins->dt;
}

const int outputStep(void)
{
  return ins->outputStep;
}

const int NtimeSteps(void)
{
  return ins->NtimeSteps;
}

const double startTime(void)
{
  return ins->startTime;
}

const double finalTime(void)
{
  return ins->finalTime;
}

void* nekPtr(const char* id)
{
  return nek_ptr(id);
}

void printRuntimeStatistics()
{
  timer::printRunStat();
}
} // namespace

static void dryRun(libParanumal::setupAide &options, int npTarget)
{
  if (rank == 0)
    cout << "performing dry-run for "
         << npTarget
         << " MPI ranks ...\n" << endl;

  string udfFile;
  options.getArgs("UDF FILE", udfFile);

  // jit compile udf
  if (!udfFile.empty()) {
    if(rank == 0) udfBuild(udfFile.c_str());
    MPI_Barrier(comm);
    *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
    *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
  }

  if(udf.setup0) udf.setup0(comm, options);

  int N;
  string casename;
  options.getArgs("CASENAME", casename);
  options.getArgs("POLYNOMIAL DEGREE", N);

  // jit compile nek
  int nscal;
  options.getArgs("NUMBER OF SCALARS", nscal);
  if (rank == 0) buildNekInterface(casename.c_str(), nscal, N, npTarget);
  MPI_Barrier(comm);

  // init solver
  ins = insSetup(comm, device, options, nrsBuildOnly);

  if (rank == 0) cout << "\nBuild successful." << endl;
}

static void setOUDF(libParanumal::setupAide &options)
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
    if (found == std::string::npos)
      out << "void insPressureDirichletConditions3D(bcData *bc){}\n";

    found = buffer.str().find("void cdsNeumannConditions");
    if (found == std::string::npos)
      out << "void cdsNeumannConditions3D(bcData *bc){}\n";

    found = buffer.str().find("void cdsDirichletConditions");
    if (found == std::string::npos)
      out << "void cdsDirichletConditions3D(bcData *bc){}\n";

    out <<
      "@kernel void __dummy__(int N) {"
      "  for (int i = 0; i < N; ++i; @tile(16, @outer, @inner)) {}"
      "}";

    out.close();
  }

  options.setArgs("DATA FILE", dataFile);
}

static void setCache(string dir)
{
  char buf[FILENAME_MAX];
  getcwd(buf, sizeof(buf));
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
}
