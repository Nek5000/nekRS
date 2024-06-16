#include "nekrsSys.hpp"
#include "compileKernels.hpp"
#include "udf.hpp"

occa::properties registerUDFKernels()
{
  const auto registerOnly = platform->options.compareArgs("REGISTER ONLY", "TRUE") ? true : false;

  std::string installDir;
  installDir.assign(getenv("NEKRS_HOME"));
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  deviceKernelProperties kernelInfo(platform->kernelInfo + meshKernelProperties(N));

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0 && !registerOnly) {
    std::cout << "loading udf kernels ... " << std::endl;
  }

  const std::string bcDataFile = installDir + "/include/nrs/bdry/bcData.h";
  kernelInfo.include() += bcDataFile.c_str();
  kernelInfo.okl_include_paths() += std::string(fs::current_path()).c_str();

  if (udf.loadKernels) {
    udf.loadKernels(kernelInfo);
  }

  // envoke after udf.loadKernels as the user may have modified kernelInfo
  if (udf.autoloadKernels) {
    udf.autoloadKernels(kernelInfo);
  }

  // internal BC kernels call device functions defined in this file
  std::string oklFileCache;
  platform->options.getArgs("OKL FILE CACHE", oklFileCache);
  kernelInfo.include() += realpath(oklFileCache.c_str(), NULL);

  udf.autoloadPlugins(kernelInfo);

  MPI_Barrier(platform->comm.mpiComm);
  const double loadTime = MPI_Wtime() - tStart;
  if (platform->comm.mpiRank == 0 && !registerOnly) {
    printf("done (%gs)\n", loadTime);
  }
  fflush(stdout);

  return static_cast<occa::properties>(kernelInfo);
}

occa::kernel oudfBuildKernel(occa::properties kernelInfo, const std::string& kernelName)
{
  const auto registerOnly = platform->options.compareArgs("REGISTER ONLY", "TRUE") ? true : false;
  const auto fileName = platform->options.getArgs("OKL FILE CACHE");
  static std::vector<std::string> reqNames;
 
  auto knl = [&] ()
  {
    if (registerOnly) {
      // hash reqName as the user might call oudfBuildKernel using different props
      const auto reqName = "udf::" + std::string(kernelInfo.hash().getString());
      reqNames.push_back(reqName); 
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      if (platform->verbose && platform->comm.mpiRank == 0) { 
        std::cout << kernelName << std::endl; 
      }
      static int idx = 0;  
      const auto reqName = reqNames[idx++];
      return platform->kernelRequests.load(reqName, kernelName);
    }
  }();

  return knl;
}
