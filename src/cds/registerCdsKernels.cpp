#include <nrs.hpp>
#include <compileKernels.hpp>
#include "re2Reader.hpp"
#include "benchmarkAdvsub.hpp"

void registerCdsKernels(occa::properties kernelInfoBC) {
  const bool serial = platform->serial;
  const std::string extension = serial ? ".c" : ".okl";
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  occa::properties kernelInfo = platform->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  int N, cubN;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);
  const int Nq = N + 1;
  const int cubNq = cubN + 1;
  const int Np = Nq * Nq * Nq;
  const int cubNp = cubNq * cubNq * cubNq;
  constexpr int Nfaces{6};

  constexpr int NVfields{3};
  kernelInfo["defines/p_NVfields"] = NVfields;

  int Nsubsteps = 0;
  int nBDF = 0;
  int nEXT = 0;
  platform->options.getArgs("SUBCYCLING STEPS", Nsubsteps);

  if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO1")) {
    nBDF = 1;
  } else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO2")) {
    nBDF = 2;
  } else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO3")) {
    nBDF = 3;
  }
  nEXT = 3;
  if (Nsubsteps)
    nEXT = nBDF;


  std::string fileName, kernelName;
  const std::string suffix = "Hex3D";
  const std::string oklpath = installDir + "/okl/";
  const std::string section = "cds-";
  occa::properties meshProps = kernelInfo;
  meshProps += meshKernelProperties(N);
  {
    kernelName = "relativeMassHighestMode";
    fileName = oklpath + "cds/regularization/" + kernelName + ".okl";
    platform->kernels.add(kernelName, fileName, meshProps);

    kernelName = "computeMaxVisc";
    fileName = oklpath + "cds/regularization/" + kernelName + ".okl";
    platform->kernels.add(kernelName, fileName, meshProps);

    kernelName = "interpolateP1";
    fileName = oklpath + "cds/regularization/" + kernelName + ".okl";
    platform->kernels.add(kernelName, fileName, meshProps);

    {
      occa::properties prop = meshProps;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;

      kernelName = "strongAdvectionVolume" + suffix;
      fileName = oklpath + "cds/" + kernelName + ".okl";
      platform->kernels.add(
          section + kernelName, fileName, prop);

      kernelName = "strongAdvectionCubatureVolume" + suffix;
      fileName = oklpath + "cds/" + kernelName + ".okl";
      platform->kernels.add(
          section + kernelName, fileName, prop);
    }

    kernelName = "advectMeshVelocityHex3D";
    fileName = oklpath + "cds/" + kernelName + ".okl";
    platform->kernels.add(
        section + kernelName, fileName, meshProps);

    kernelName = "maskCopy";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(
        section + kernelName, fileName, meshProps);

    {
      occa::properties prop = kernelInfo;
      const int movingMesh =
          platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/p_MovingMesh"] = movingMesh;
      prop["defines/p_nEXT"] = nEXT;
      prop["defines/p_nBDF"] = nBDF;
      if (Nsubsteps)
        prop["defines/p_SUBCYCLING"] = 1;
      else
        prop["defines/p_SUBCYCLING"] = 0;

      kernelName = "sumMakef";
      fileName = oklpath + "cds/" + kernelName + ".okl";
      platform->kernels.add(
          section + kernelName, fileName, prop);
    }

    kernelName = "helmholtzBC" + suffix;
    fileName = oklpath + "cds/" + kernelName + ".okl";
    platform->kernels.add(
        section + kernelName, fileName, kernelInfoBC);
    kernelName = "dirichletBC";
    fileName = oklpath + "cds/" + kernelName + ".okl";
    platform->kernels.add(
        section + kernelName, fileName, kernelInfoBC);

    kernelName = "setEllipticCoeff";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(
        section + kernelName, fileName, kernelInfo);

    kernelName = "filterRT" + suffix;
    fileName = oklpath + "cds/regularization/" + kernelName + ".okl";
    platform->kernels.add(
        section + kernelName, fileName, meshProps);

    kernelName = "nStagesSum3";
    fileName = oklpath + "core/" + kernelName + ".okl";
    platform->kernels.add(
        section + kernelName, fileName, platform->kernelInfo);

    {
      occa::properties prop = meshProps;
      const int movingMesh =
          platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/p_MovingMesh"] = movingMesh;
      prop["defines/p_nEXT"] = nEXT;
      prop["defines/p_nBDF"] = nBDF;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;

      occa::properties subCycleStrongCubatureProps = prop;

      int nelgt, nelgv;
      const std::string meshFile = platform->options.getArgs("MESH FILE");
      re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);
      const int NelemBenchmark = nelgv/platform->comm.mpiCommSize;
      bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");
      const int verbosity = verbose ? 2 : 1;

      int Nsubsteps;
      platform->options.getArgs("SUBCYCLING STEPS", Nsubsteps);

      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE") && Nsubsteps) {
        auto subCycleKernel =
            benchmarkAdvsub(1, NelemBenchmark, Nq, cubNq, nEXT, true, true, verbosity, 0.5, false);

        kernelName = "subCycleStrongCubatureVolume" + suffix;
        platform->kernels.add(section + kernelName, subCycleKernel);
      }

      kernelName = "subCycleStrongVolume" + suffix;
      fileName = oklpath + "cds/" + kernelName + ".okl";
      platform->kernels.add(
          section + kernelName, fileName, prop);

      kernelName = "subCycleRKUpdate";
      fileName = oklpath + "cds/" + kernelName + ".okl";
      platform->kernels.add(
          section + kernelName, fileName, prop);
      kernelName = "subCycleRK";
      fileName = oklpath + "cds/" + kernelName + ".okl";
      platform->kernels.add(
          section + kernelName, fileName, prop);

      kernelName = "subCycleInitU0";
      fileName = oklpath + "cds/" + kernelName + ".okl";
      platform->kernels.add(
          section + kernelName, fileName, prop);
    }
  }
}
