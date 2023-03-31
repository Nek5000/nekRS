#include "nrs.hpp"
#include <compileKernels.hpp>
#include "mesh.h"

void registerMeshKernels(occa::properties kernelInfoBC)
{
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
  std::string oklpath = getenv("NEKRS_KERNEL_DIR");
  std::string kernelName;

  const std::string meshPrefix = "mesh-";
  std::string fileName;
  {

    kernelName = "velocityDirichletBCHex3D";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernels.add(meshPrefix + kernelName, fileName, kernelInfoBC);

    {
      int N;
      platform->options.getArgs("POLYNOMIAL DEGREE", N);
      const int Nq = N + 1;
      nrsCheck(BLOCKSIZE < Nq * Nq, platform->comm.mpiComm, EXIT_FAILURE,
               "surfaceIntegral kernel requires BLOCKSIZE >= Nq * Nq\nBLOCKSIZE = %d, Nq*Nq = %d\n",
               BLOCKSIZE, Nq * Nq);
    }

    kernelName = "surfaceIntegral";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernels.add(meshPrefix + kernelName, fileName, kernelInfo);

    occa::properties meshKernelInfo = kernelInfo;
    meshKernelInfo["defines/p_cubNq"] = cubNq;
    meshKernelInfo["defines/p_cubNp"] = cubNp;

    kernelName = "geometricFactorsHex3D";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernels.add(meshPrefix + kernelName, fileName, meshKernelInfo);
    kernelName = "surfaceGeometricFactorsHex3D";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernels.add(meshPrefix + kernelName, fileName, meshKernelInfo);

    kernelName = "cubatureGeometricFactorsHex3D";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernels.add(meshPrefix + kernelName, fileName, meshKernelInfo);

    meshKernelInfo = kernelInfo;
    meshKernelInfo["defines/p_nAB"] = nAB;
    kernelName = "nStagesSumVector";
    fileName = oklpath + "/core/" + kernelName + ".okl";
    platform->kernels.add(meshPrefix + kernelName, fileName, meshKernelInfo);
  }
}
