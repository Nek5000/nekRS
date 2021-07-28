#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>

#include "udf.hpp"
#include "io.hpp"
#include "platform.hpp"

UDF udf = {NULL, NULL, NULL, NULL};

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

int udfBuild(const char* casename, const char* udfFile, int buildOnly)
{
  double tStart = MPI_Wtime();
  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  const char* udf_dir = getenv("NEKRS_UDF_DIR");

  if(!fileExists(udfFile)) {
    printf("\nERROR: Cannot find %s!\n", udfFile);
    return EXIT_FAILURE;
  }

  char udfFileCache[BUFSIZ];
  sprintf(udfFileCache,"%s/udf/udf-%s.cpp",cache_dir, casename);
  char udfLib[BUFSIZ];
  sprintf(udfLib, "%s/udf/libUDF-%s.so", cache_dir, casename);

  char cmd[BUFSIZ];
  printf("building udf ... \n"); fflush(stdout);
  if(isFileNewer(udfFile, udfFileCache) || !fileExists(udfLib)) {
    char udfFileResolved[BUFSIZ];
    realpath(udfFile, udfFileResolved);
    sprintf(cmd,
            "mkdir -p %s/udf && cd %s/udf && cp -f %s udf-%s.cpp && cp %s/CMakeLists.txt . && \
             rm -rf libUDF-%s.so && cmake -Wno-dev -DCMAKE_CXX_COMPILER=\"$NEKRS_CXX\" \
	     -DCMAKE_CXX_FLAGS=\"$NEKRS_CXXFLAGS\" -DUDF_DIR=\"%s\" -DCASENAME=\"%s\" .",
             cache_dir,
             cache_dir,
             udfFileResolved,
             casename,
             udf_dir,
             casename,
             udf_dir,
             casename);
    if(system(cmd)) return EXIT_FAILURE; 
  }
  sprintf(cmd, "cd %s/udf && make", cache_dir);
  if(system(cmd)) return EXIT_FAILURE; 
  printf("done (%gs)\n", MPI_Wtime() - tStart);
  fflush(stdout);

  return 0;
}

void* udfLoadFunction(const char* casename, const char* fname, int errchk)
{
  char udfLib[BUFSIZ];

  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  sprintf(udfLib, "%s/udf/libUDF-%s.so", cache_dir, casename);

  void* h, * fptr;
  h = dlopen(udfLib, RTLD_LAZY | RTLD_GLOBAL);
  if (!h) goto err;

  fptr = dlsym(h,fname);
  if (!fptr && errchk) goto err;

  return fptr;

err:
  fprintf(stderr, "Error in %s(): %s\n", __func__, dlerror());
  ABORT(EXIT_FAILURE);
}

void udfLoad(const char* casename)
{
  *(void**)(&udf.setup0) = udfLoadFunction(casename, "UDF_Setup0",0);
  *(void**)(&udf.setup) = udfLoadFunction(casename, "UDF_Setup",1);
  *(void**)(&udf.loadKernels) = udfLoadFunction(casename, "UDF_LoadKernels",0);
  *(void**)(&udf.executeStep) = udfLoadFunction(casename, "UDF_ExecuteStep",0);
}

occa::kernel udfBuildKernel(nrs_t* nrs, const char* function)
{
  int rank;
  mesh_t* mesh = nrs->meshV;
  
  MPI_Comm_rank(platform->comm.mpiComm, &rank);

  string install_dir;
  occa::properties kernelInfo = *nrs->kernelInfo;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const string bcDataFile = install_dir + "/include/core/bcData.h";
  kernelInfo["includes"] += bcDataFile.c_str();

  string oudf;
  platform->options.getArgs("DATA FILE", oudf);

  return platform->device.buildKernel(oudf.c_str(), function, kernelInfo);
}
