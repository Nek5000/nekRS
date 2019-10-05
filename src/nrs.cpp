#include "nekrs.hpp"
#include "meshSetup.hpp"
#include "insSetup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "parReader.hpp"
#include "runTime.hpp"

static int rank, size;
static MPI_Comm comm;

static double timeLast = -1;
static ins_t *ins;

static libParanumal::setupAide options;

static int ioStep;

int nrsBuildOnly = 0;  // hack for meshPhysicalNodes() 

static void setCache(string dir);
static void setOUDF(libParanumal::setupAide &options);
static void dryRun(libParanumal::setupAide &options, int npTarget); 

namespace nrs {

setupAide setup(MPI_Comm comm_in, int buildOnly, int sizeTarget, 
                int ciMode, string cacheDir, string setupFile)
{
  MPI_Comm_dup(comm_in, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
 
  nrsBuildOnly = buildOnly; 
  setCache(cacheDir);

  if (rank == 0) {
    #include "printHeader.inc"
    cout << "MPI ranks: " << size << endl << endl;
    cout << "using OCCA_CACHE_DIR: " << occa::env::OCCA_CACHE_DIR << endl << endl;
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
  int N;
  int nscal;
  string udfFile, casename;
  int readRestartFile;

  options = parRead(setupFile, comm);

  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("UDF FILE", udfFile);
  options.getArgs("CASENAME", casename);
  options.getArgs("RESTART FROM FILE", readRestartFile);
  options.getArgs("NUMBER OF SCALARS", nscal);

  setOUDF(options);

  if (nrsBuildOnly){
    dryRun(options, sizeTarget);
    return options;
  }

  MPI_Barrier(comm); double t0 = MPI_Wtime();

  // jit compile udf
  if (!udfFile.empty()){
    if(rank == 0) udfBuild(udfFile.c_str()); 
    MPI_Barrier(comm);
    udfLoad();
  } 

  options.setArgs("CI-MODE", std::to_string(ciMode));
  if(rank == 0 && ciMode) 
    cout << "enabling continous integration mode\n" << endl;

  if(udf.setup0) udf.setup0(comm, options);

  // jit compile nek
  if(rank == 0) buildNekInterface(casename.c_str(), nscal+3, N, size);
  MPI_Barrier(comm);

  // init nek
  nek_setup(comm, options);
  nek_setic();
  nek_userchk();

  // init solver
  mesh_t *mesh = meshSetup(comm, options, nrsBuildOnly);
  ins = insSetup(mesh, options);

  // set initial condition
  if(readRestartFile) nek_copyRestart(ins);
  if(udf.setup) udf.setup(ins);
  ins->o_U.copyFrom(ins->U);
  ins->o_P.copyFrom(ins->P);
  if(ins->Nscalar) ins->cds->o_S.copyFrom(ins->cds->S);    

  if (udf.executeStep) udf.executeStep(ins, ins->startTime, 0);
  nek_ocopyFrom(ins, ins->startTime, 0);

  if(rank == 0) {
    cout << "\nsettings:\n" << endl;
    cout << ins->vOptions << endl;
    cout << "initialization took " << MPI_Wtime() - t0 << " seconds" << endl; 
    size_t dMB = ins->mesh->device.memoryAllocated() / 1e6;
    cout << "device memory allocation: " << dMB << " MB" << endl;
  }
  fflush(stdout);

  return options;
}

void runStep(double time, int tstep)
{
  runStep(ins, time, tstep);
}

void copyToNek(double time, int tstep)
{
  nek_ocopyFrom(ins, time, tstep);
  timeLast = time;
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

void nekOutfld(void)
{
  nek_outfld();
}

} // namespace

static void dryRun(libParanumal::setupAide &options, int npTarget) 
{
  if (rank == 0) 
    cout << "performing dry-run for "
         << npTarget
         << " MPI ranks ...\n" << endl;

  int N;
  string udfFile, casename;
  options.getArgs("UDF FILE", udfFile);
  options.getArgs("CASENAME", casename);
  options.getArgs("POLYNOMIAL DEGREE", N);

  // jit compile udf
  if (!udfFile.empty()) {
    if(rank == 0) udfBuild(udfFile.c_str()); 
    MPI_Barrier(comm);
    *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
    *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
  } 

  if(udf.setup0) udf.setup0(comm, options);

  // jit compile nek
  int nscal;
  options.getArgs("NUMBER OF SCALARS", nscal);
  if (rank == 0) buildNekInterface(casename.c_str(), nscal+3, N, npTarget);
  MPI_Barrier(comm);

  mesh_t *mesh = meshSetup(comm, options, nrsBuildOnly);

  ins_t *ins = insSetup(mesh, options);

  if (mesh->rank == 0) cout << "\nBuild successful." << endl;
}

static void setOUDF(libParanumal::setupAide &options)
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

    found = buffer.str().find("void cdsNeumannConditions");
    if (found == std::string::npos)
      out << "void cdsNeumannConditions3D(bcData *bc){}\n"; 

    found = buffer.str().find("void cdsDirichletConditions");
    if (found ==std::string::npos)
      out << "void cdsDirichletConditions3D(bcData *bc){}\n"; 

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

