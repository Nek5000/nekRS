#include "nrs.hpp"
#include "mesh.h"
#include <compileKernels.hpp>

#include "re2Reader.hpp"
#include "benchmarkAdvsub.hpp"

void registerNrsKernels(occa::properties kernelInfoBC)
{
  const bool serial = platform->serial;
  const std::string extension = serial ? ".c" : ".okl";
  const device_t &device = platform->device;
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  // build kernels
  std::string fileName, kernelName;
  const std::string suffix = "Hex3D";
  const std::string oklpath = installDir + "/kernels/";
  int N, cubN;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);
  const int Nq = N + 1;
  const int cubNq = cubN + 1;
  const int Np = Nq * Nq * Nq;
  const int cubNp = cubNq * cubNq * cubNq;
  constexpr int Nfaces{6};

  occa::properties kernelInfo = platform->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  constexpr int NVfields{3};
  kernelInfo["defines/p_NVfields"] = NVfields;

  int Nsubsteps = 0;
  int nBDF = 0;
  int nEXT = 0;
  platform->options.getArgs("SUBCYCLING STEPS", Nsubsteps);

  if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO1")) {
    nBDF = 1;
  }
  else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO2")) {
    nBDF = 2;
  }
  else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO3")) {
    nBDF = 3;
  }
  nEXT = 3;
  if (Nsubsteps)
    nEXT = nBDF;

  {
    kernelName = "nStagesSum3";
    fileName = oklpath + "core/" + kernelName + ".okl";
    const std::string section = "nrs-";
    platform->kernels.add(section + kernelName, fileName, platform->kernelInfo);

    kernelName = "computeFieldDotNormal";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, platform->kernelInfo);

    occa::properties centroidProp = kernelInfo;
    centroidProp["defines/p_Nfp"] = Nq * Nq;
    centroidProp["defines/p_Nfaces"] = Nfaces;
    {
      int N;
      platform->options.getArgs("POLYNOMIAL DEGREE", N);
      const int Nq = N + 1;
      if (BLOCKSIZE < Nq * Nq) {
        if (platform->comm.mpiRank == 0)
          printf("ERROR: computeFaceCentroid kernel requires BLOCKSIZE >= Nq * Nq."
                 "BLOCKSIZE = %d, Nq*Nq = %d\n",
                 BLOCKSIZE,
                 Nq * Nq);
        ABORT(EXIT_FAILURE);
      }
    }
    kernelName = "computeFaceCentroid";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, centroidProp);

    occa::properties meshProps = kernelInfo;
    meshProps += meshKernelProperties(N);

    {
      occa::properties prop = meshProps;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;

      kernelName = "strongAdvectionVolume" + suffix;
      fileName = oklpath + "nrs/" + kernelName + ".okl";
      platform->kernels.add(section + kernelName, fileName, prop);
      kernelName = "strongAdvectionCubatureVolume" + suffix;
      fileName = oklpath + "nrs/" + kernelName + ".okl";
      platform->kernels.add(section + kernelName, fileName, prop);
    }

    kernelName = "curl" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    kernelName = "gradientVolume" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    kernelName = "wGradientVolume" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    {
      occa::properties prop = kernelInfo;
      const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/p_nEXT"] = nEXT;
      prop["defines/p_nBDF"] = nBDF;
      prop["defines/p_MovingMesh"] = movingMesh;
      if (Nsubsteps)
        prop["defines/p_SUBCYCLING"] = 1;
      else
        prop["defines/p_SUBCYCLING"] = 0;

      kernelName = "sumMakef";
      fileName = oklpath + "nrs/" + kernelName + ".okl";
      platform->kernels.add(section + kernelName, fileName, prop);
    }

    kernelName = "wDivergenceVolume" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfoBC);
    kernelName = "divergenceVolume" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfoBC);

    kernelName = "divergenceSurface" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfoBC);

    kernelName = "advectMeshVelocityHex3D";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    kernelName = "pressureRhs" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    kernelName = "pressureStress" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    kernelName = "pressureDirichletBC" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfoBC);

    kernelName = "velocityRhs" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    auto zeroNormalProps = kernelInfoBC;
    zeroNormalProps["defines/p_ZERO_NORMAL"] = ZERO_NORMAL;
    zeroNormalProps["defines/p_NO_OP"] = NO_OP;
    kernelName = "averageNormalBcType";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, zeroNormalProps);

    kernelName = "fixZeroNormalMask";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, zeroNormalProps);

    kernelName = "applyZeroNormalMask";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, zeroNormalProps);

    kernelName = "initializeZeroNormalMask";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, zeroNormalProps);

    kernelName = "velocityDirichletBC" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfoBC);

    kernelName = "velocityNeumannBC" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfoBC);

    occa::properties prop = meshProps;
    const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
    prop["defines/p_relative"] = movingMesh && Nsubsteps;
    prop["defines/p_cubNq"] = cubNq;
    prop["defines/p_cubNp"] = cubNp;
    fileName = oklpath + "nrs/Urst" + suffix + ".okl";

    kernelName = "UrstCubature" + suffix;
    fileName = oklpath + "nrs/" + kernelName + extension;
    platform->kernels.add(section + kernelName, fileName, prop);

    kernelName = "Urst" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, prop);

    {
      occa::properties prop = meshProps;
      const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/p_MovingMesh"] = movingMesh;
      prop["defines/p_nEXT"] = nEXT;
      prop["defines/p_nBDF"] = nBDF;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;

      occa::properties subCycleStrongCubatureProps = prop;

      int nelgt, nelgv;
      const std::string meshFile = platform->options.getArgs("MESH FILE");
      re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);
      const int NelemBenchmark = nelgv / platform->comm.mpiCommSize;

      bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");
      const int verbosity = verbose ? 2 : 1;

      int Nsubsteps;
      platform->options.getArgs("SUBCYCLING STEPS", Nsubsteps);

      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE") && Nsubsteps) {
        auto subCycleKernel = benchmarkAdvsub(3,
                                              NelemBenchmark,
                                              Nq,
                                              cubNq,
                                              nEXT,
                                              true,
                                              false,
                                              verbosity,
                                              nrs_t::targetTimeBenchmark,
                                              false);

        kernelName = "subCycleStrongCubatureVolume" + suffix;
        platform->kernels.add(section + kernelName, subCycleKernel);
      }

      kernelName = "subCycleStrongVolume" + suffix;
      fileName = oklpath + "nrs/" + kernelName + ".okl";
      platform->kernels.add(section + kernelName, fileName, prop);

      kernelName = "subCycleRKUpdate";
      fileName = oklpath + "nrs/" + kernelName + ".okl";
      platform->kernels.add(section + kernelName, fileName, prop);
      kernelName = "subCycleRK";
      fileName = oklpath + "nrs/" + kernelName + ".okl";
      platform->kernels.add(section + kernelName, fileName, prop);

      kernelName = "subCycleInitU0";
      fileName = oklpath + "nrs/" + kernelName + ".okl";
      platform->kernels.add(section + kernelName, fileName, prop);
    }

    kernelName = "extrapolate";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    kernelName = "maskCopy";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfo);

    kernelName = "maskCopy2";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfo);

    kernelName = "mask";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfo);

    kernelName = "filterRT" + suffix;
    fileName = oklpath + "nrs/regularization/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    {
      int N;
      platform->options.getArgs("POLYNOMIAL DEGREE", N);
      const int Nq = N + 1;
      if (BLOCKSIZE < Nq * Nq) {
        if (platform->comm.mpiRank == 0)
          printf("ERROR: cfl kernel requires BLOCKSIZE >= Nq * Nq."
                 "BLOCKSIZE = %d, Nq*Nq = %d\n",
                 BLOCKSIZE,
                 Nq * Nq);
        ABORT(EXIT_FAILURE);
      }
    }

    occa::properties cflProps = meshProps;
    cflProps["defines/p_MovingMesh"] = movingMesh;
    kernelName = "cfl" + suffix;
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, cflProps);

    kernelName = "pressureAddQtl";
    fileName = oklpath + "nrs/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, meshProps);

    kernelName = "setEllipticCoeff";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfo);
    kernelName = "setEllipticCoeffPressure";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(section + kernelName, fileName, kernelInfo);
  }
}