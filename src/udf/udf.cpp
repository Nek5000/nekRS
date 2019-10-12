#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>

#include "udf.hpp"
#include "environment.hpp"

UDF udf = {NULL, NULL, NULL, NULL};

void udfBuild(const std::string udfFile)
{
  int retval;
  std::stringstream cmd;

  printf("building udf ... "); fflush(stdout);
  os::makeDir(os::joinPath(env::cacheDir(),"udf"));

  if(!os::exist(udfFile) || !os::readable(udfFile)) {
    std::cout << "ERROR: Cannot find or read: " << udfFile << std::endl;
    exit(EXIT_FAILURE);
  }

  cmd << "cp " << os::joinPath(env::udfDir(),"/CMakeLists.txt") << " "
      << os::joinPath(env::cacheDir(),"udf");
  retval=system(cmd.str().c_str());
  if(retval) goto err;

  cmd.str(std::string());
  cmd << "cp " << udfFile << " " <<os::joinPath(env::cacheDir(),"udf");
  retval=system(cmd.str().c_str());
  if(retval) goto err;

  cmd.str(std::string());
  cmd << "cd "              << env::cacheDir()    <<"/udf && "
      << "CXX="             << env::cxxCompiler() << " "
      << "CXXFLAGS=\""      << env::cxxFlags()    << "\" "
      << "cmake -DUDF_DIR=" << env::udfDir()      << " " << "-DFILENAME=" << udfFile << " "
      << "-DNEKRS_INSTALL_DIR=" << env::installDir() << " "
      << "-DNEKRS_LIBP_DEFINES=" << env::libPDefines() << " . "
      << ". >build.log 2>&1";
  retval = system(cmd.str().c_str());
  if(retval) goto err;

  cmd.str(std::string());
  cmd << "cd " << env::cacheDir() << "/udf && "
      << "make >>build.log 2>&1";
  retval = system(cmd.str().c_str());
  if(retval) goto err;

  printf("done\n");
  fflush(stdout);
  sync();
  return;

err:
  std::cout << "\nAn ERROR occured, see "<<env::cacheDir() << "/udf/build.log for details!\n";
  exit(EXIT_FAILURE);
}

void *udfLoadFunction(const char *fname, int errchk)
{
  char udfLib[BUFSIZ];

  std::string cacheDir=env::cacheDir();
  sprintf(udfLib, "%s/udf/libUDF.so", cacheDir.c_str());

  void *h, *fptr;
  h = dlopen(udfLib, RTLD_LAZY | RTLD_GLOBAL);
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

occa::kernel udfBuildKernel(ins_t *ins, const char *function)
{
  int rank;
  mesh_t *mesh = ins->mesh;
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
  occa::properties& kernelInfo = *ins->kernelInfo;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
       k = mesh->device.buildKernel(fname, function, kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }

  if (rank == 0) {
    remove(fname);
  }
  return k;
}
