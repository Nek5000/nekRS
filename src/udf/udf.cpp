#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>

#include "udf.hpp"

UDF udf = {NULL, NULL, NULL, NULL};

void udfBuild(const char *udfFile)
{
  char cmd[BUFSIZ];
  char abs_path[BUFSIZ];
  char *ptr;
  int retval;

  const char *cache_dir = getenv("NEKRS_CACHE_DIR");
  const char *udf_dir = getenv("NEKRS_UDF_DIR");

  printf("building udf ... "); fflush(stdout);
  sprintf(cmd, "mkdir -p %s/udf && cd %s/udf && rm -rf CMake* Makefile cmake_install.cmake", 
          cache_dir, cache_dir);
  system(cmd);

  ptr = realpath(udfFile, abs_path);
  if(!ptr) {
    printf("\nERROR: Cannot find %s!\n", udfFile);
    ABORT(EXIT_FAILURE);
  }
 
  sprintf(cmd,"cd %s/udf && cp %s/CMakeLists.txt . && CXX=\"${NEKRS_CXX}\" CXXFLAGS=\"${NEKRS_CXXFLAGS}\" \
              cmake -DUDF_DIR=\"%s\" -DFILENAME=\"%s\" . >build.log 2>&1 && \
              make >>build.log 2>&1",
              cache_dir, udf_dir, udf_dir, abs_path);
  retval = system(cmd);
  if(retval) goto err;

  printf("done\n");
  fflush(stdout);
  return;

err:
  printf("\nAn ERROR occured, see %s/udf/build.log for details!\n", cache_dir);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void *udfLoadFunction(const char *fname, int errchk)
{
  char udfLib[BUFSIZ];

  const char *cache_dir = getenv("NEKRS_CACHE_DIR");
  sprintf(udfLib, "%s/udf/libUDF.so", cache_dir);

  void *h, *fptr;
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
  *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
  *(void**)(&udf.executeStep) = udfLoadFunction("UDF_ExecuteStep",1);
}

occa::kernel udfBuildKernel(ins_t *ins, const char *function)
{
  int rank;
  mesh_t *mesh = ins->mesh;
  MPI_Comm_rank(mesh->comm, &rank);

  string install_dir;
  occa::properties kernelInfo = *ins->kernelInfo;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const string bcDataFile = install_dir + "/include/insBcData.h";
  kernelInfo["includes"] += bcDataFile.c_str();
 
  string oudf;
  ins->options.getArgs("DATA FILE", oudf);

  occa::kernel k;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
      k = mesh->device.buildKernel(oudf.c_str(), function, kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
  return k;
}
