#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <regex>

#include "udf.hpp"
#include "io.hpp"
#include "platform.hpp"

UDF udf = {NULL, NULL, NULL, NULL};

static int velocityDirichletConditions = 0;
static int velocityNeumannConditions = 0;
static int pressureDirichletConditions = 0;
static int scalarDirichletConditions = 0;
static int scalarNeumannConditions = 0;

uint32_t fchecksum(std::ifstream& file) 
{
    uint32_t checksum = 0;
    unsigned shift = 0;
    for (uint32_t ch = file.get(); file; ch = file.get()) {
        checksum += (ch << shift);
        shift += 8;
        if (shift == 32) {
            shift = 0;
        }
    }
    return checksum;
}

void oudfFindDirichlet(std::string &field)
{
  if (field.find("velocity") != std::string::npos && !velocityDirichletConditions) {
    if (platform->comm.mpiRank == 0) std::cout << "Cannot find oudf function: velocityDirichletConditions!\n";
    ABORT(EXIT_FAILURE);
  }
  if (field.find("scalar") != std::string::npos && !scalarDirichletConditions) {
    if (platform->comm.mpiRank == 0) std::cout << "Cannot find oudf function: scalarDirichletConditions!\n";
    ABORT(EXIT_FAILURE);
  }
  if(field == "pressure" && !pressureDirichletConditions) {
    if (platform->comm.mpiRank == 0) std::cout << "Cannot find oudf function: pressureDirichletConditions!\n";
    // ABORT(EXIT_FAILURE); this bc is optional 
  }
}

void oudfFindNeumann(std::string &field)
{
  if (field.find("velocity") != std::string::npos && !velocityNeumannConditions) {
    if (platform->comm.mpiRank == 0) std::cout << "Cannot find oudf function: velocityNeumannConditions!\n";
    ABORT(EXIT_FAILURE);
  }
  if (field.find("scalar") != std::string::npos && !scalarNeumannConditions) {
    if (platform->comm.mpiRank == 0) std::cout << "Cannot find oudf function: scalarNeumannConditions!\n";
    ABORT(EXIT_FAILURE);
  }
}

void oudfInit(setupAide &options)
{
  std::string oklFile;
  options.getArgs("UDF OKL FILE",oklFile);

  char* ptr = realpath(oklFile.c_str(), NULL);
  if(!ptr) {
    if (platform->comm.mpiRank == 0) std::cout << "ERROR: Cannot find " << oklFile << "!\n";
    ABORT(EXIT_FAILURE);;
  }
  free(ptr);

  std::string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
  std::string casename;
  options.getArgs("CASENAME", casename);
  std::string udf_dir = cache_dir + "/udf";
  const std::string dataFile = udf_dir + "/udf.okl";

  int buildRank = platform->comm.mpiRank;
  MPI_Comm comm = platform->comm.mpiComm;
  int buildNodeLocal = 0;
  if (getenv("NEKRS_BUILD_NODE_LOCAL"))
    buildNodeLocal = std::stoi(getenv("NEKRS_BUILD_NODE_LOCAL"));
  if(buildNodeLocal) {
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);
    comm = platform->comm.mpiCommLocal;
  }

  if (buildRank == 0) {
    std::ifstream in;
    in.open(oklFile);
    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();

    std::ofstream out;
    out.open(dataFile, std::ios::trunc);

    out << buffer.str();

    bool found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+velocityDirichletConditions)"));
    velocityDirichletConditions = found;
    if(!found)
      out << "void velocityDirichletConditions(bcData *bc){}\n";

    found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+velocityNeumannConditions)"));
    velocityNeumannConditions = found;
    if(!found)
      out << "void velocityNeumannConditions(bcData *bc){}\n";

    found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+pressureDirichletConditions)"));
    pressureDirichletConditions = found;
    if(!found)
      out << "void pressureDirichletConditions(bcData *bc){}\n";

    found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+scalarNeumannConditions)"));
    scalarNeumannConditions = found;
    if(!found)
      out << "void scalarNeumannConditions(bcData *bc){}\n";

    found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+scalarDirichletConditions)"));
    scalarDirichletConditions = found;
    if(!found)
      out << "void scalarDirichletConditions(bcData *bc){}\n";

    out <<
      "@kernel void __dummy__(int N) {"
      "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {}"
      "}";

    out.close();
  }

  MPI_Bcast(&velocityDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&velocityNeumannConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&pressureDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&scalarNeumannConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&scalarDirichletConditions, 1, MPI_INT, 0, comm);

  options.setArgs("DATA FILE", dataFile);
}


int udfBuild(const char* udfFile, setupAide& options)
{
  double tStart = MPI_Wtime();

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  std::string udf_dir = install_dir + "/udf";

  std::string cache_dir;
  cache_dir.assign(getenv("NEKRS_CACHE_DIR"));

  std::string udfFileCache = cache_dir + "/udf/udf.cpp";
  std::string udfLib = cache_dir + "/udf/libUDF.so";

  char buf[FILENAME_MAX];
  char * ret = getcwd(buf, sizeof(buf));
  if(!ret) ABORT(EXIT_FAILURE);
  std::string case_dir;
  case_dir.assign(buf);

  const int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;

  if(!fileExists(udfFile)) {
    printf("\nERROR: Cannot find %s!\n", udfFile);
    return EXIT_FAILURE;
  }

  char cmd[BUFSIZ];
  printf("building udf ... \n"); fflush(stdout);
  if(isFileNewer(udfFile, udfFileCache.c_str()) || !fileExists(udfLib.c_str())) {
    char udfFileResolved[BUFSIZ];
    realpath(udfFile, udfFileResolved);
    sprintf(cmd,
            "cd %s/udf && cp -f %s udf.cpp && cp -f %s/CMakeLists.txt . && \
             rm -f *.so && cmake -Wno-dev -DCASE_DIR=\"%s\" -DCMAKE_CXX_COMPILER=\"$NEKRS_CXX\" \
	         -DCMAKE_CXX_FLAGS=\"$NEKRS_CXXFLAGS\" .",
             cache_dir.c_str(),
             udfFileResolved,
             udf_dir.c_str(),
             case_dir.c_str());
    if(verbose) printf("%s\n", cmd);
    if(system(cmd)) return EXIT_FAILURE; 
  }
  sprintf(cmd, "cd %s/udf && make", cache_dir.c_str());
  if(system(cmd)) return EXIT_FAILURE; 
  fileSync(udfLib.c_str());
  printf("done (%gs)\n", MPI_Wtime() - tStart);
  fflush(stdout);

  return 0;
}

void* udfLoadFunction(const char* fname, int errchk)
{
  char udfLib[BUFSIZ];

  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  sprintf(udfLib, "%s/udf/libUDF.so", cache_dir);

  void* h, * fptr;
  h = dlopen(udfLib, RTLD_LAZY | RTLD_GLOBAL);
  if (!h) goto errOpen;

  fptr = dlsym(h,fname);
  if (!fptr) {
    if(platform->comm.mpiRank == 0) printf("Cannot find udf function: %s!\n", fname);
    if(errchk) goto err;
  }

  return fptr;

errOpen:
  fprintf(stderr, "Error in %s(): %s\n", __func__, dlerror());

err:
  ABORT(EXIT_FAILURE);
}

void udfLoad(void)
{
  *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
  *(void**)(&udf.setup) = udfLoadFunction("UDF_Setup",1);
  *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
  *(void**)(&udf.executeStep) = udfLoadFunction("UDF_ExecuteStep",0);
}

occa::kernel udfBuildKernel(occa::properties kernelInfo, const char* function)
{
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string bcDataFile = install_dir + "/include/core/bcData.h";
  kernelInfo["includes"] += bcDataFile.c_str();

  // provide some common kernel args
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const int Nq = N+1;
  const int Np = Nq * Nq * Nq;
  kernelInfo["defines/p_Nq"] = Nq;
  kernelInfo["defines/p_Np"] = Np;

  std::string oudf;
  platform->options.getArgs("DATA FILE", oudf);

  return platform->device.buildKernel(oudf.c_str(), function, kernelInfo);
}
