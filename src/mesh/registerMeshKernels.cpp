#include "nrs.hpp"
#include <compileKernels.hpp>
#include "mesh.h"

void registerMeshKernels(occa::properties kernelInfoBC) {
  int N, cubN;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);
  const int Nq = N + 1;
  const int cubNq = cubN + 1;
  const int Np = Nq * Nq * Nq;
  const int cubNp = cubNq * cubNq * cubNq;

  int nAB;
  platform->options.getArgs("MESH INTEGRATION ORDER", nAB);

  auto kernelInfo = platform->kernelInfo + meshKernelProperties(N);
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  std::string oklpath = installDir + "/okl/";
  std::string kernelName;

  const std::string meshPrefix = "mesh-";
  {
    std::string fileName = oklpath + "mesh/velocityBCHex3D.okl";
    kernelName = "velocityDirichletBCHex3D";
    platform->kernels.add(meshPrefix + kernelName, fileName, kernelName, kernelInfo);
    occa::properties meshKernelInfo = kernelInfo;
    meshKernelInfo["defines/p_cubNq"] = cubNq;
    meshKernelInfo["defines/p_cubNp"] = cubNp;

    fileName = oklpath + "mesh/geometricFactorsHex3D.okl";
    kernelName = "geometricFactorsHex3D";
    platform->kernels.add(
        meshPrefix + kernelName, fileName, kernelName, meshKernelInfo);
    fileName = oklpath + "mesh/surfaceGeometricFactorsHex3D.okl";
    kernelName = "surfaceGeometricFactorsHex3D";
    platform->kernels.add(
        meshPrefix + kernelName, fileName, kernelName, meshKernelInfo);

    meshKernelInfo = kernelInfo;
    meshKernelInfo["defines/p_nAB"] = nAB;
    fileName = oklpath + "core/nStagesSum.okl";
    kernelName = "nStagesSumVector";
    platform->kernels.add(
        meshPrefix + kernelName, fileName, kernelName, meshKernelInfo);
  }
}
