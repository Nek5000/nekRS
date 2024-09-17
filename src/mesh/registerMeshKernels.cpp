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

      // N to M
      for (int Nc = 1; Nc < N; Nc++) { 
        auto props = kernelInfo;
        props["defines/p_NqFine"] = N + 1;
        props["defines/p_NqCoarse"] = Nc + 1;
        props["defines/pfloat"] = dfloatString;

        props["defines/p_NpFine"] = (N + 1) * (N + 1) * (N + 1); 
        props["defines/p_NpCoarse"] = (Nc + 1) * (Nc + 1) * (Nc + 1);;

        const std::string ext = platform->serial ? ".c" : ".okl";
        const std::string orderSuffix =
          std::string("_Nf_") + std::to_string(N) + std::string("_Nc_") + std::to_string(Nc);

        kernelName = "coarsenHex3D";
        fileName = oklpath + "/mesh/" + kernelName + ext;
        platform->kernelRequests.add(meshPrefix + kernelName + orderSuffix, fileName, props);
      }

      for (int Nf = N; Nf < mesh_t::maxNqIntp; Nf++) { 
        auto props = kernelInfo;
        props["defines/p_NqFine"] = Nf + 1;
        props["defines/p_NqCoarse"] = N + 1;
        props["defines/pfloat"] = dfloatString;

        props["defines/p_NpFine"] = (Nf + 1) * (Nf + 1) * (Nf + 1); 
        props["defines/p_NpCoarse"] = (N + 1) * (N + 1) * (N + 1);

        const std::string ext = platform->serial ? ".c" : ".okl";
        const std::string orderSuffix =
          std::string("_Nf_") + std::to_string(Nf) + std::string("_Nc_") + std::to_string(N);

        kernelName = "prolongateHex3D";
        fileName = oklpath + "/mesh/" + kernelName + ext;
        platform->kernelRequests.add(meshPrefix + kernelName + orderSuffix, fileName, props);
      }

      auto addIntpKernels = [&](bool transpose = false)
      {
        for (int M = 1; M < mesh_t::maxNqIntp; M++) {
          if (M == N) continue;
 
          bool condition = transpose ? (M > N) : (M <= N);

          const auto Nf = condition ? N : M;
          const auto Nc = condition ? M : N;
          kernelName = condition ? "prolongateHex3D" : "coarsenHex3D";

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
        }
      };

      // N to M
      addIntpKernels();
      // M to N
      addIntpKernels(true);

      // M to N
      for (int M = 1; M < mesh_t::maxNqIntp; M++) {
        if (M == N) continue;

        int Nf;
        int Nc;

        if (M > N) {
          Nf = M;
          Nc = N;
          kernelName = "coarsenHex3D";
        } else {
          Nf = N;
          Nc = M;
          kernelName = "prolongateHex3D";
        }
 
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
