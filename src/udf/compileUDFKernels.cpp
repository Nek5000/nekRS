#include <compileKernels.hpp>
#include "udf.hpp"

occa::properties compileUDFKernels()
{
  const bool buildNodeLocal = platform->cacheLocal;

  std::string installDir;
  installDir.assign(getenv("NEKRS_HOME"));
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  occa::properties kernelInfo = platform->kernelInfo + meshKernelProperties(N);
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0)
    std::cout << "loading udf kernels ... " << std::endl;

  occa::properties kernelInfoBC = kernelInfo;
  if (udf.loadKernels) {
    // kernelInfoBC will include any relevant user-defined kernel props
    udf.loadKernels(kernelInfoBC);
  }
  if (udf.autoloadKernels) {
    // kernelInfoBC will include any relevant user-defined kernel props
    udf.autoloadKernels(kernelInfoBC);
  }
  udf.autoloadPlugins(kernelInfoBC);
  
  const std::string bcDataFile = installDir + "/include/bdry/bcData.h";
  kernelInfoBC["includes"] += bcDataFile.c_str();

  std::string oklFileCache;
  platform->options.getArgs("OKL FILE CACHE", oklFileCache);
  kernelInfoBC["includes"] += realpath(oklFileCache.c_str(), NULL);

  kernelInfoBC += meshKernelProperties(N);

  // just to bail out early in case included source doesn't compile 
  {
   const std::string dummyKernelName = "compileUDFKernelsTest";
   const std::string dummyKernelStr = std::string("@kernel void compileUDFKernelsTest(int N) {"
                                                  "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {}" "}");

  if (platform->comm.mpiRank == 0)
     platform->device.occaDevice().buildKernelFromString(dummyKernelStr, dummyKernelName, kernelInfoBC);
 
  }

  MPI_Barrier(platform->comm.mpiComm);
  const double loadTime = MPI_Wtime() - tStart;
  if (platform->comm.mpiRank == 0)
    printf("done (%gs)\n", loadTime);
  fflush(stdout);

  return kernelInfoBC;
}
