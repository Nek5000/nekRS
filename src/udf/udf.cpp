#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <regex>

#include "udf.hpp"
#include "ioUtils.hpp"
#include "platform.hpp"
#include "bcMap.hpp"

UDF udf = {NULL, NULL, NULL, NULL};

static int velocityDirichletConditions = 0;
static int meshVelocityDirichletConditions = 0;
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
    if (platform->comm.mpiRank == 0) 
     std::cout << "WARNING: Cannot find oudf function: pressureDirichletConditions!\n";
    // ABORT(EXIT_FAILURE); this bc is optional 
  }
  if(field.find("mesh") != std::string::npos && !meshVelocityDirichletConditions && !bcMap::useDerivedMeshBoundaryConditions()) {
    if (platform->comm.mpiRank == 0) std::cout << "Cannot find oudf function: meshVelocityDirichletConditions!\n";
    ABORT(EXIT_FAILURE);
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
  const bool buildNodeLocal = useNodeLocalCache();
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

    found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+meshVelocityDirichletConditions)"));
    meshVelocityDirichletConditions = found;
    if(!found)
      out << "void meshVelocityDirichletConditions(bcData *bc){\n"
      "  velocityDirichletConditions(bc);\n"
      "}\n";

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

    out.close();
  }

  MPI_Bcast(&velocityDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&meshVelocityDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&velocityNeumannConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&pressureDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&scalarNeumannConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&scalarDirichletConditions, 1, MPI_INT, 0, comm);

  options.setArgs("DATA FILE", dataFile);
}


void udfBuild(const char* udfFile, setupAide& options)
{
  int buildRank = platform->comm.mpiRank;
  const bool buildNodeLocal = useNodeLocalCache();
  if(buildNodeLocal)
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);    
  
  int err = [&](){
    if(buildRank == 0){
      double tStart = MPI_Wtime();

      std::string installDir;
      installDir.assign(getenv("NEKRS_INSTALL_DIR"));
      std::string udf_dir = installDir + "/udf";

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
        if(platform->comm.mpiRank == 0) printf("\nERROR: Cannot find %s!\n", udfFile);
        return EXIT_FAILURE;
      }

      char cmd[10*BUFSIZ];
      if(platform->comm.mpiRank == 0) printf("building udf ... \n"); fflush(stdout);
      std::string pipeToNull = (platform->comm.mpiRank == 0) ?
        std::string("") :
        std::string("> /dev/null 2>&1");
      if(isFileNewer(udfFile, udfFileCache.c_str()) || !fileExists(udfLib.c_str())) {
        char udfFileResolved[BUFSIZ];
        realpath(udfFile, udfFileResolved);
        sprintf(cmd,
                "cp -f %s %s/udf/udf.cpp && cp -f %s/CMakeLists.txt %s/udf && rm -f %s/udf/*.so",
                udfFileResolved,
                cache_dir.c_str(),
                udf_dir.c_str(),
                cache_dir.c_str(),
                cache_dir.c_str());
        if(verbose && platform->comm.mpiRank == 0) printf("%s\n", cmd);
        if(system(cmd)) return EXIT_FAILURE; 

        std::string cmakeFlags("-Wno-dev");
        if(verbose) cmakeFlags += " --trace-expand";
        std::string cmakeBuildDir = cache_dir + "/udf"; 
        sprintf(cmd, "cmake %s -S %s -B %s -DCASE_DIR=\"%s\" -DCMAKE_CXX_COMPILER=\"$NEKRS_CXX\" "
	            "-DCMAKE_CXX_FLAGS=\"$NEKRS_CXXFLAGS\" %s",
                 cmakeFlags.c_str(),
                 cmakeBuildDir.c_str(),
                 cmakeBuildDir.c_str(),
                 case_dir.c_str(),
                 pipeToNull.c_str());
        const int retVal = system(cmd);
        if(verbose && platform->comm.mpiRank == 0) {
          printf("%s (retVal: %d)\n", cmd, retVal);
        }
        if(retVal) return EXIT_FAILURE; 
      }

      {
        sprintf(cmd, "cd %s/udf && make %s", cache_dir.c_str(), pipeToNull.c_str());
        const int retVal = system(cmd);
        if(verbose && platform->comm.mpiRank == 0) {
          printf("%s (retVal: %d)\n", cmd, retVal);
        }
        if(retVal) return EXIT_FAILURE; 
      }

      fileSync(udfLib.c_str());
      if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStart);
      fflush(stdout);
    }

    return 0;
  }();

  MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, platform->comm.mpiComm);
  if(err) ABORT(EXIT_FAILURE);
}

void* udfLoadFunction(const char* fname, int errchk)
{
  char udfLib[BUFSIZ];

  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  sprintf(udfLib, "%s/udf/libUDF.so", cache_dir);

  void* h, * fptr;
  h = dlopen(udfLib, RTLD_NOW | RTLD_GLOBAL);
  if (!h) goto errOpen;

  fptr = dlsym(h,fname);
  if (!fptr) {
    if(platform->comm.mpiRank == 0 && errchk) printf("Cannot find udf function: %s!\n", fname);
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
  *(void**)(&udf.setup) = udfLoadFunction("UDF_Setup",0);
  *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
  *(void**)(&udf.executeStep) = udfLoadFunction("UDF_ExecuteStep",0);
}

occa::kernel oudfBuildKernel(occa::properties kernelInfo, const char *function)
{
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string bcDataFile = installDir + "/include/bdry/bcData.h";
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
