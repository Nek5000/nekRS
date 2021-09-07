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
  const std::string dataFileDir = cache_dir + "/udf/";
  const std::string dataFile = dataFileDir + "udf.okl";

  if (platform->comm.mpiRank == 0) {
    mkdir(dataFileDir.c_str(), S_IRWXU);

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
      "  for (int i = 0; i < N; ++i; @tile(16, @outer, @inner)) {}"
      "}";

    out.close();
  }

  MPI_Bcast(&velocityDirichletConditions, 1, MPI_INT, 0, platform->comm.mpiComm);
  MPI_Bcast(&velocityNeumannConditions, 1, MPI_INT, 0, platform->comm.mpiComm);
  MPI_Bcast(&pressureDirichletConditions, 1, MPI_INT, 0, platform->comm.mpiComm);
  MPI_Bcast(&scalarNeumannConditions, 1, MPI_INT, 0, platform->comm.mpiComm);
  MPI_Bcast(&scalarDirichletConditions, 1, MPI_INT, 0, platform->comm.mpiComm);

  options.setArgs("DATA FILE", dataFile);
}


int udfBuild(const char* udfFile, setupAide& options)
{
  double tStart = MPI_Wtime();
  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  const char* udf_dir = getenv("NEKRS_UDF_DIR");

  const int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;

  if(!fileExists(udfFile)) {
    printf("\nERROR: Cannot find %s!\n", udfFile);
    return EXIT_FAILURE;
  }

  char udfFileCache[BUFSIZ];
  sprintf(udfFileCache,"%s/udf/udf.cpp",cache_dir);
  char udfLib[BUFSIZ];
  sprintf(udfLib, "%s/udf/libUDF.so", cache_dir);

  char cmd[BUFSIZ];
  printf("building udf ... \n"); fflush(stdout);
  if(isFileNewer(udfFile, udfFileCache) || !fileExists(udfLib)) {
    char udfFileResolved[BUFSIZ];
    realpath(udfFile, udfFileResolved);
    sprintf(cmd,
            "mkdir -p %s/udf && cd %s/udf && cp -f %s udf.cpp && cp %s/CMakeLists.txt . && \
             rm -rf *.so && cmake -Wno-dev -DCMAKE_CXX_COMPILER=\"$NEKRS_CXX\" \
	         -DCMAKE_CXX_FLAGS=\"$NEKRS_CXXFLAGS\" -DUDF_DIR=\"%s\" .",
             cache_dir,
             cache_dir,
             udfFileResolved,
             udf_dir,
             udf_dir);
    fileSync(udfFileCache);
    fileSync(udfLib);
    if(verbose) printf("%s\n", cmd);
    if(system(cmd)) return EXIT_FAILURE; 
  }
  sprintf(cmd, "cd %s/udf && make", cache_dir);
  if(system(cmd)) return EXIT_FAILURE; 
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
  if (!h) goto err;

  fptr = dlsym(h,fname);
  if (!fptr) {
    if(platform->comm.mpiRank == 0) printf("Cannot find udf function: %s!\n", fname);
    if(errchk) goto err;
  }

  return fptr;

err:
  fprintf(stderr, "Error in %s(): %s\n", __func__, dlerror());
  ABORT(EXIT_FAILURE);
}

void udfLoad(void)
{
  *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
  *(void**)(&udf.setup) = udfLoadFunction("UDF_Setup",1);
  *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",0);
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
