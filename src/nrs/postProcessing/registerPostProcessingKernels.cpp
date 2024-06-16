#include "compileKernels.hpp"
#include "nekrsSys.hpp"
#include "platform.hpp"
#include "mesh.h" // for meshKernelProperties

void registerPostProcessingKernels()
{
  // gatherPlanarValues and scatterPlanarValues kernels require use of atomics
  if (!platform->device.deviceAtomic) {
    return;
  }

  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const int Nq = N + 1;
  const int Np = Nq * Nq * Nq;

  auto kernelInfo = platform->kernelInfo + meshKernelProperties(N);
  const std::string section = "nrs-";

  kernelInfo["includes"].asArray();

  const std::string oklpath = getenv("NEKRS_KERNEL_DIR");
  std::string kernelName, fileName;

  kernelName = "aeroForces";
  fileName = oklpath + "/nrs/postProcessing/" + kernelName + ".okl";
  platform->kernelRequests.add(section + kernelName, fileName, kernelInfo);

  kernelName = "Qcriterion";
  fileName = oklpath + "/nrs/postProcessing/" + kernelName + ".okl";
  platform->kernelRequests.add(section + kernelName, fileName, kernelInfo);
}
