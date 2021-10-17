#include "nrs.hpp"
#include "mesh.h"
#include <compileKernels.hpp>

void registerNrsKernels(occa::properties kernelInfoBC)
{
  const device_t &device = platform->device;
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  // build kernels
  std::string fileName, kernelName;
  const std::string suffix = "Hex3D";
  const std::string oklpath = installDir + "/okl/";
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
  } else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO2")) {
    nBDF = 2;
  } else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO3")) {
    nBDF = 3;
  }
  nEXT = 3;
  if (Nsubsteps)
    nEXT = nBDF;

  {
    fileName = oklpath + "core/nStagesSum.okl";
    kernelName = "nStagesSum3";
    const std::string section = "nrs-";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, platform->kernelInfo);

    fileName = oklpath + "nrs/computeFieldDotNormal.okl";
    kernelName = "computeFieldDotNormal";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, platform->kernelInfo);

    occa::properties centroidProp = kernelInfo;
    centroidProp["defines/p_Nfp"] = Nq * Nq;
    centroidProp["defines/p_Nfaces"] = Nfaces;
    fileName = oklpath + "nrs/computeFaceCentroid.okl";
    kernelName = "computeFaceCentroid";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, centroidProp);

    occa::properties meshProps = kernelInfo;
    meshProps += meshKernelProperties(N);

    {
      occa::properties prop = meshProps;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;
      fileName = oklpath + "nrs/advection" + suffix + ".okl";
      kernelName = "strongAdvectionVolume" + suffix;
      platform->kernels.add(
          section + kernelName, fileName, kernelName, prop);
      kernelName = "strongAdvectionCubatureVolume" + suffix;
      platform->kernels.add(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "nrs/curl" + suffix + ".okl";
    kernelName = "curl" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/gradient" + suffix + ".okl";
    kernelName = "gradientVolume" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    kernelName = "nrswGradientVolume" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    {
      occa::properties prop = kernelInfo;
      const int movingMesh =
          platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/p_nEXT"] = nEXT;
      prop["defines/p_nBDF"] = nBDF;
      prop["defines/p_MovingMesh"] = movingMesh;
      if (Nsubsteps)
        prop["defines/p_SUBCYCLING"] = 1;
      else
        prop["defines/p_SUBCYCLING"] = 0;

      fileName = oklpath + "nrs/sumMakef.okl";
      kernelName = "sumMakef";
      platform->kernels.add(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "nrs/divergence" + suffix + ".okl";
    kernelName = "nrswDivergenceVolume" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfoBC);
    kernelName = "divergenceVolume" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    kernelName = "divergenceSurfaceTOMBO" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    fileName = oklpath + "nrs/advectMeshVelocityHex3D.okl";
    kernelName = "advectMeshVelocityHex3D";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/pressureRhs" + suffix + ".okl";
    kernelName = "pressureRhsTOMBO" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/pressureStress" + suffix + ".okl";
    kernelName = "pressureStress" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/pressureBC" + suffix + ".okl";
    kernelName = "pressureDirichletBC" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    fileName = oklpath + "nrs/velocityRhs" + suffix + ".okl";
    kernelName = "velocityRhsTOMBO" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/velocityBC" + suffix + ".okl";
    kernelName = "velocityDirichletBC" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    kernelName = "velocityNeumannBC" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    occa::properties prop = meshProps;
    const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
    prop["defines/p_relative"] = movingMesh && Nsubsteps;
    prop["defines/p_cubNq"] = cubNq;
    prop["defines/p_cubNp"] = cubNp;
    fileName = oklpath + "nrs/Urst" + suffix + ".okl";

    const bool serial = useSerial();
    if(serial) fileName = oklpath + "nrs/Urst" + suffix + ".c";
    kernelName = "UrstCubature" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, prop);

    fileName = oklpath + "nrs/Urst" + suffix + ".okl";
    kernelName = "Urst" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, prop);

    {
      occa::properties prop = meshProps;
      const int movingMesh =
          platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/p_MovingMesh"] = movingMesh;
      prop["defines/p_nEXT"] = nEXT;
      prop["defines/p_nBDF"] = nBDF;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;

      fileName = oklpath + "nrs/subCycle" + suffix + ".okl";
      occa::properties subCycleStrongCubatureProps = prop;

      if(serial)
        fileName = oklpath + "nrs/subCycle" + suffix + ".c";

      kernelName = "subCycleStrongCubatureVolume" + suffix;
      platform->kernels.add(section + kernelName,
          fileName,
          kernelName,
          subCycleStrongCubatureProps);
      fileName = oklpath + "nrs/subCycle" + suffix + ".okl";

      kernelName = "subCycleStrongVolume" + suffix;
      platform->kernels.add(
          section + kernelName, fileName, kernelName, prop);

      fileName = oklpath + "nrs/subCycleRKUpdate" + ".okl";
      kernelName = "subCycleERKUpdate";
      platform->kernels.add(
          section + kernelName, fileName, kernelName, prop);
      kernelName = "subCycleRK";
      platform->kernels.add(
          section + kernelName, fileName, kernelName, prop);

      kernelName = "subCycleInitU0";
      platform->kernels.add(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "nrs/extrapolate" + ".okl";
    kernelName = "multiExtrapolate";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "core/mask" + ".okl";
    kernelName = "maskCopy";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfo);
    kernelName = "mask";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfo);

    fileName = oklpath + "nrs/regularization/filterRT" + suffix + ".okl";
    kernelName = "filterRT" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    occa::properties cflProps = meshProps;
    cflProps["defines/p_MovingMesh"] = movingMesh;
    fileName = oklpath + "nrs/cfl" + suffix + ".okl";
    kernelName = "cfl" + suffix;
    platform->kernels.add(
        section + kernelName, fileName, kernelName, cflProps);

    fileName = oklpath + "nrs/pressureAddQtl" + ".okl";
    kernelName = "pressureAddQtl";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "core/setEllipticCoeff.okl";
    kernelName = "setEllipticCoeff";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfo);
    kernelName = "setEllipticCoeffPressure";
    platform->kernels.add(
        section + kernelName, fileName, kernelName, kernelInfo);
  }
}