#include "compileKernels.hpp"
#include "nrssys.hpp"
#include "platform.hpp"
#include "mesh.h" // for meshKernelProperties
#include "postProcessing.hpp"
void registerPostProcessingKernels()
{
  // gatherPlanarValues and scatterPlanarValues kernels require use of atomics
  if (!platform->device.deviceAtomic)
    return;

  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const int Nq = N + 1;
  const int Np = Nq * Nq * Nq;

  auto kernelInfo = platform->kernelInfo + meshKernelProperties(N);

  kernelInfo["includes"].asArray();

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  std::string oklpath = installDir + "/okl/";
  std::string kernelName, fileName;

  kernelInfo["includes"] += oklpath + "postProcessing/planarAveraging.h";

  for(const std::string dir : {"XY", "XZ", "YZ"}){
    kernelName = "gatherPlanarValues" + dir;
    fileName = oklpath + "postProcessing/" + kernelName + ".okl";
    platform->kernels.add(kernelName, fileName, kernelInfo);
  
    kernelName = "scatterPlanarValues" + dir;
    fileName = oklpath + "postProcessing/" + kernelName + ".okl";
    platform->kernels.add(kernelName, fileName, kernelInfo);
  }
}