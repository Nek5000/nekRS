#include <compileKernels.hpp>
#include "udf.hpp"

occa::properties compileUDFKernels()
{
  const bool buildNodeLocal = useNodeLocalCache();

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
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
    printf("loading udf kernels ... ");
  fflush(stdout);

  occa::properties kernelInfoBC = kernelInfo;
  if (udf.loadKernels) {
    // kernelInfoBC will include any relevant user-defined kernel props
    udf.loadKernels(kernelInfoBC);
  }
  const std::string bcDataFile = installDir + "/include/bdry/bcData.h";
  kernelInfoBC["includes"] += bcDataFile.c_str();
  std::string boundaryHeaderFileName;
  platform->options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfoBC["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

  kernelInfoBC += meshKernelProperties(N);

  MPI_Barrier(platform->comm.mpiComm);
  const double loadTime = MPI_Wtime() - tStart;
  if (platform->comm.mpiRank == 0)
    printf("done (%gs)\n", loadTime);
  fflush(stdout);

  return kernelInfoBC;
}
