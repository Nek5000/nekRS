#include "platform.hpp"
#include "compileKernels.hpp"
#include "mesh.h"

void registerMeshKernels(occa::properties kernelInfoBC)
{
  int p, pCub = 0;
  platform->options.getArgs("POLYNOMIAL DEGREE", p);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", pCub);

  std::vector<int> Nlist = {p};

  for (auto &N : Nlist) {
    const int Nq = N + 1;
    const int Np = Nq * Nq * Nq;

    const std::string meshPrefix = "mesh-";
    const std::string orderSuffix = "_" + std::to_string(N);

    auto kernelInfo = platform->kernelInfo + meshKernelProperties(N);
    std::string oklpath = getenv("NEKRS_KERNEL_DIR");

    std::string fileName;
    std::string kernelName;

    occa::properties meshKernelInfo = kernelInfo;
    
    kernelName = "geometricFactorsHex3D";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernelRequests.add(meshPrefix + kernelName + orderSuffix, fileName, meshKernelInfo);

    kernelName = "surfaceGeometricFactorsHex3D";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernelRequests.add(meshPrefix + kernelName + orderSuffix, fileName, meshKernelInfo);

    const int cubNq = (N == p) ? pCub + 1 : 1;
    const int cubNp = cubNq * cubNq * cubNq;

    auto meshCubKernelInfo = meshKernelInfo;
    meshCubKernelInfo["defines/p_cubNq"] = cubNq;
    meshCubKernelInfo["defines/p_cubNp"] = cubNp;

    kernelName = "cubatureGeometricFactorsHex3D";
    fileName = oklpath + "/mesh/" + kernelName + ".okl";
    platform->kernelRequests.add(meshPrefix + kernelName + orderSuffix, fileName, meshCubKernelInfo);

    if (N == p) {
      kernelName = "velocityDirichletBCHex3D";
      fileName = oklpath + "/mesh/" + kernelName + ".okl";
      platform->kernelRequests.add(meshPrefix + kernelName, fileName, kernelInfoBC);

      auto addIntpKernels = [&] (int Nf, int Nc, std::string kernelName)
      {
        if (Nf < Nc) return;

        auto props = kernelInfo;
        props["defines/p_NqFine"] = Nf + 1;
        props["defines/p_NqCoarse"] = Nc + 1;
        props["defines/pfloat"] = dfloatString;

        props["defines/p_NpFine"] = (Nf + 1) * (Nf + 1) * (Nf + 1); 
        props["defines/p_NpCoarse"] = (Nc + 1) * (Nc + 1) * (Nc + 1);;

        const std::string ext = platform->serial ? ".c" : ".okl";
        const std::string orderSuffix =
          std::string("_Nf_") + std::to_string(Nf) + std::string("_Nc_") + std::to_string(Nc);

        fileName = oklpath + "/mesh/" + kernelName + ext;
        platform->kernelRequests.add(meshPrefix + kernelName + orderSuffix, fileName, props);
      };

      // N to M
      for (int M = 1; M < mesh_t::maxNqIntp; M++) {
        //if (M == N) continue;

        {
          auto transpose = false;
          bool condition = transpose ? (N > M) : (N <= M);
          const auto Nf = condition ? M : N;
          const auto Nc = condition ? N : M;
          kernelName = condition ? "prolongateHex3D" : "coarsenHex3D";
          addIntpKernels(Nf, Nc, kernelName);
        }
 
        {
          auto transpose = true;
          bool condition = transpose ? (M > N) : (M <= N);
          const auto Nf = condition ? N : M;
          const auto Nc = condition ? M : N;
          kernelName = condition ? "prolongateHex3D" : "coarsenHex3D";
          addIntpKernels(Nf, Nc, kernelName);
        }
      }

      // M to N
      for (int M = 1; M < mesh_t::maxNqIntp; M++) {
        if (M == N) continue;

        {
          auto transpose = false;
          bool condition = transpose ? (M > N) : (M <= N);
          const auto Nf = condition ? N : M;
          const auto Nc = condition ? M : N;
          kernelName = condition ? "prolongateHex3D" : "coarsenHex3D";
          addIntpKernels(Nf, Nc, kernelName);
        }
 
        {
          auto transpose = true;
          bool condition = transpose ? (M > N) : (M <= N);
          const auto Nf = condition ? N : M;
          const auto Nc = condition ? M : N;
          kernelName = condition ? "prolongateHex3D" : "coarsenHex3D";
          addIntpKernels(Nf, Nc, kernelName);
        }
      }

      meshKernelInfo = kernelInfo;
      int nAB = 3;
      platform->options.getArgs("MESH INTEGRATION ORDER", nAB);
      meshKernelInfo["defines/p_nAB"] = nAB;
      kernelName = "nStagesSumVector";
      fileName = oklpath + "/core/" + kernelName + ".okl";
      platform->kernelRequests.add(meshPrefix + kernelName, fileName, meshKernelInfo);

      auto prop = kernelInfo;
      prop["defines/p_ndot"] = 0;
      kernelName = "surfaceAreaNormalMultiplyIntegrateHex3D";
      fileName = oklpath + "/mesh/" + kernelName + ".okl";
      platform->kernelRequests.add(meshPrefix + kernelName, fileName, prop);
 
      prop["defines/p_ndot"] = 1;
      kernelName = "surfaceAreaNormalMultiplyIntegrateHex3D-ndot";
      platform->kernelRequests.add(meshPrefix + kernelName, fileName, prop);

      kernelName = "surfaceAreaMultiplyHex3D";
      fileName = oklpath + "/mesh/" + kernelName + ".okl";
      platform->kernelRequests.add(meshPrefix + kernelName, fileName, kernelInfo);

      kernelName = "setBIDHex3D";
      fileName = oklpath + "/mesh/" + kernelName + ".okl";
      platform->kernelRequests.add(meshPrefix + kernelName, fileName, kernelInfo);
 
      kernelName = "distanceHex3D";
      fileName = oklpath + "/mesh/" + kernelName + ".okl";
      platform->kernelRequests.add(meshPrefix + kernelName, fileName, kernelInfo);

      auto hlongSumKernelInfo = kernelInfo;
      hlongSumKernelInfo["defines/dfloat"] = hlongString;
      kernelName = "sum";
      fileName = oklpath + "/core/linAlg/" + kernelName + ".okl";
      platform->kernelRequests.add("hlong-" + meshPrefix + kernelName, fileName, hlongSumKernelInfo);

      for (const std::string dir : {"XY", "XZ", "YZ"}) {
        auto props = kernelInfo;
        props["includes"].asArray();
        props["includes"] += oklpath + "/mesh/planarAveraging.h";
 
        kernelName = "gatherPlanarValues" + dir;
        fileName = oklpath + "/mesh/" + kernelName + ".okl";
        platform->kernelRequests.add(kernelName, fileName, props);
 
        kernelName = "scatterPlanarValues" + dir;
        fileName = oklpath + "/mesh/" + kernelName + ".okl";
        platform->kernelRequests.add(kernelName, fileName, props);
      }
    }
  }
}
