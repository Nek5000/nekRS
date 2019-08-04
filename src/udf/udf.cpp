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
  sprintf(cmd, "mkdir -p %s/udf", cache_dir);
  system(cmd);

  ptr = realpath(udfFile, abs_path);
  if(!ptr) {
    printf("\nERROR: Cannot find %s!\n", udfFile);
    exit(EXIT_FAILURE);
  }
 
  sprintf(cmd,"cp %s/CMakeLists.txt %s/udf", udf_dir, cache_dir);
  system(cmd);
  sprintf(cmd,"cd %s/udf && CXX=\"${NEKRS_CXX}\" CXXFLAGS=\"${NEKRS_CXXFLAGS}\" \
      cmake -DUDF_DIR=\"%s\" -DFILENAME=\"%s\" . >build.log 2>&1", cache_dir, udf_dir, abs_path);
  retval = system(cmd);
  if(retval) goto err;

  sprintf(cmd,"cd %s/udf && make >>build.log 2>&1", cache_dir);
  retval = system(cmd);
  if(retval) goto err;

  printf("done\n");
  fflush(stdout);
  sync();
  return;

err:
  printf("\nAn ERROR occured, see %s/udf/build.log for details!\n", cache_dir);
  exit(EXIT_FAILURE);
}

void *udfLoadFunction(const char *fname, int errchk)
{
  char udfLib[BUFSIZ];

  const char *cache_dir = getenv("NEKRS_CACHE_DIR");
  sprintf(udfLib, "%s/udf/libUDF.so", cache_dir);

  void *h, *fptr;
  h = dlopen(udfLib, RTLD_NOW | RTLD_GLOBAL);
  if (!h) goto err;

  fptr = dlsym(h,fname);
  if (!fptr && errchk) goto err;

  return fptr;

err:
  fprintf(stderr, "Error in %s(): %s\n", __func__, dlerror());
  exit(EXIT_FAILURE);
}

void udfLoad(void)
{
  *(void**)(&udf.setup0) = udfLoadFunction("UDF_Setup0",0);
  *(void**)(&udf.setup) = udfLoadFunction("UDF_Setup",1);
  *(void**)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels",1);
  *(void**)(&udf.executeStep) = udfLoadFunction("UDF_ExecuteStep",1);
}

occa::kernel udfBuildKernel(mesh_t *mesh, const char *function, const occa::properties kernelInfo)
{
  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  const char *fname = "__udf.okl";

  // add dummy to make occa happy
  if (rank == 0){
    std::ofstream ss;
    ss.open(fname, std::ofstream::out);
    ss << "@kernel void __dummy(const int entries){ \n"
       << "  for (int group = 0; group < entries; group += 128; @outer) {\n"
       << "    for (int item = group; item < entries; ++item; @inner) {\n"
       << "      const int n = item;\n\n"
       << "    }\n"
       << "  }\n"
       << "}\n"
       << "\n";
    ss.close();
  }

  occa::kernel k;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
       k = mesh->device.buildKernel(fname, function, kernelInfo);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (rank == 0) remove(fname);
  return k;
}
