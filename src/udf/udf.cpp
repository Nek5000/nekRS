#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>

#include "udf.hpp"
#include "io.hpp"

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

int udfBuild(const char* udfFile)
{
  char cmd[BUFSIZ];
  char abs_path[BUFSIZ];

  double tStart = MPI_Wtime();
  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  const char* udf_dir = getenv("NEKRS_UDF_DIR");

  if(! realpath(udfFile, abs_path)) {
    printf("\nERROR: Cannot find %s!\n", udfFile);
    return EXIT_FAILURE;
  }

  char udfFileCache[BUFSIZ];
  sprintf(udfFileCache,"%s/udf/udf.cpp",cache_dir);

  if(isFileNewer(udfFile, udfFileCache)) {
    printf("building udf ... "); fflush(stdout);
    sprintf(cmd,
            "mkdir -p %s/udf && cd %s/udf && cp %s/CMakeLists.txt . && \
             CXX=\"${NEKRS_CXX}\" CXXFLAGS=\"${NEKRS_CXXFLAGS}\" \
             cmake -DUDF_DIR=\"%s\" -DFILENAME=\"%s\" . >build.log 2>&1 && \
             make >>build.log 2>&1",
            cache_dir,
            cache_dir,
            udf_dir,
            udf_dir,
            abs_path);

    if(system(cmd)) { 
      printf("\nAn ERROR occured, see %s/udf/build.log for details!\n", cache_dir);
      return EXIT_FAILURE;
    }
    printf("done (%gs)\n", MPI_Wtime() - tStart);
    fflush(stdout);
  }

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
  if (!fptr && errchk) goto err;

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

occa::kernel udfBuildKernel(nrs_t* nrs, const char* function)
{
  int rank;
  mesh_t* mesh = nrs->mesh;
  MPI_Comm_rank(mesh->comm, &rank);

  string install_dir;
  occa::properties kernelInfo = *nrs->kernelInfo;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const string bcDataFile = install_dir + "/include/core/bcData.h";
  kernelInfo["includes"] += bcDataFile.c_str();

  string oudf;
  nrs->options.getArgs("DATA FILE", oudf);

  occa::kernel k;
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && rank == 0) || (r == 1 && rank > 0))
      k = mesh->device.buildKernel(oudf.c_str(), function, kernelInfo);
    MPI_Barrier(mesh->comm);
  }
  return k;
}
