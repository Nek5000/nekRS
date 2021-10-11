#include "elliptic.h"
#include "mesh.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "udf.hpp"
#include <map>
#include <vector>

namespace {
static occa::properties kernelInfoBC;
std::string createOptionsPrefix(const std::string &section) {
  std::string prefix = section + std::string(" ");
  if (section.find("temperature") != std::string::npos) {
    prefix = std::string("scalar00 ");
  }
  std::transform(
      prefix.begin(), prefix.end(), prefix.begin(), [](unsigned char c) {
        return std::toupper(c);
      });
  return prefix;
}

void registerGMRESKernels(const std::string &section, int Nfields) {
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  std::string fileName;
  const bool serial = (platform->device.mode() == "Serial" ||
                       platform->device.mode() == "OpenMP");

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";
  const std::string sectionIdentifier = std::to_string(Nfields) + "-";

  occa::properties gmresKernelInfo = platform->kernelInfo;
  gmresKernelInfo["defines/p_Nfields"] = Nfields;

  std::string kernelName = "gramSchmidtOrthogonalization";
  fileName = oklpath + "ellipticGramSchmidtOrthogonalization" + fileNameExtension;
  platform->kernels.add_kernel(
      sectionIdentifier + kernelName, fileName, kernelName, gmresKernelInfo);

  kernelName = "updatePGMRESSolution";
  fileName = oklpath + "ellipticUpdatePGMRES" + fileNameExtension;
  platform->kernels.add_kernel(
      sectionIdentifier + kernelName, fileName, kernelName, gmresKernelInfo);

  kernelName = "fusedResidualAndNorm";
  fileName = oklpath + "ellipticFusedResidualAndNorm" + fileNameExtension;
  platform->kernels.add_kernel(
      sectionIdentifier + kernelName, fileName, kernelName, gmresKernelInfo);
}

void registerNrsKernels() {
  const device_t &device = platform->device;
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  // build kernels
  std::string fileName, kernelName;
  const std::string suffix = "Hex3D";
  const std::string oklpath = install_dir + "/okl/";
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
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, platform->kernelInfo);

    fileName = oklpath + "nrs/computeFieldDotNormal.okl";
    kernelName = "computeFieldDotNormal";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, platform->kernelInfo);

    occa::properties centroidProp = kernelInfo;
    centroidProp["defines/p_Nfp"] = Nq * Nq;
    centroidProp["defines/p_Nfaces"] = Nfaces;
    fileName = oklpath + "nrs/computeFaceCentroid.okl";
    kernelName = "computeFaceCentroid";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, centroidProp);

    occa::properties meshProps = kernelInfo;
    meshProps += meshKernelProperties(N);

    {
      occa::properties prop = meshProps;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;
      fileName = oklpath + "nrs/advection" + suffix + ".okl";
      kernelName = "strongAdvectionVolume" + suffix;
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
      kernelName = "strongAdvectionCubatureVolume" + suffix;
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "nrs/curl" + suffix + ".okl";
    kernelName = "curl" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/gradient" + suffix + ".okl";
    kernelName = "gradientVolume" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    kernelName = "nrswGradientVolume" + suffix;
    platform->kernels.add_kernel(
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
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "nrs/divergence" + suffix + ".okl";
    kernelName = "nrswDivergenceVolume" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);
    kernelName = "divergenceVolume" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    kernelName = "divergenceSurfaceTOMBO" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    fileName = oklpath + "nrs/advectMeshVelocityHex3D.okl";
    kernelName = "advectMeshVelocityHex3D";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/pressureRhs" + suffix + ".okl";
    kernelName = "pressureRhsTOMBO" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/pressureStress" + suffix + ".okl";
    kernelName = "pressureStress" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/pressureBC" + suffix + ".okl";
    kernelName = "pressureDirichletBC" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    fileName = oklpath + "nrs/velocityRhs" + suffix + ".okl";
    kernelName = "velocityRhsTOMBO" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "nrs/velocityBC" + suffix + ".okl";
    kernelName = "velocityDirichletBC" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    kernelName = "velocityNeumannBC" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    occa::properties prop = meshProps;
    const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
    prop["defines/p_relative"] = movingMesh && Nsubsteps;
    prop["defines/p_cubNq"] = cubNq;
    prop["defines/p_cubNp"] = cubNp;
    fileName = oklpath + "nrs/Urst" + suffix + ".okl";

    const bool serial = (platform->device.mode() == "Serial" ||
      platform->device.mode() == "OpenMP");
    if(serial) fileName = oklpath + "nrs/Urst" + suffix + ".c";
    kernelName = "UrstCubature" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, prop);

    fileName = oklpath + "nrs/Urst" + suffix + ".okl";
    kernelName = "Urst" + suffix;
    platform->kernels.add_kernel(
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
      if (platform->device.mode() == "Serial" ||
          platform->device.mode() == "OpenMP") {
        fileName = oklpath + "nrs/subCycle" + suffix + ".c";
      }
      kernelName = "subCycleStrongCubatureVolume" + suffix;
      platform->kernels.add_kernel(section + kernelName,
          fileName,
          kernelName,
          subCycleStrongCubatureProps);
      fileName = oklpath + "nrs/subCycle" + suffix + ".okl";

      kernelName = "subCycleStrongVolume" + suffix;
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);

      fileName = oklpath + "nrs/subCycleRKUpdate" + ".okl";
      kernelName = "subCycleERKUpdate";
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
      kernelName = "subCycleRK";
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);

      kernelName = "subCycleInitU0";
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "nrs/extrapolate" + ".okl";
    kernelName = "multiExtrapolate";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "core/mask" + ".okl";
    kernelName = "maskCopy";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfo);
    kernelName = "mask";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfo);

    fileName = oklpath + "nrs/regularization/filterRT" + suffix + ".okl";
    kernelName = "filterRT" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    occa::properties cflProps = meshProps;
    cflProps["defines/p_MovingMesh"] = movingMesh;
    fileName = oklpath + "nrs/cfl" + suffix + ".okl";
    kernelName = "cfl" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, cflProps);

    fileName = oklpath + "nrs/pressureAddQtl" + ".okl";
    kernelName = "pressureAddQtl";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "core/setEllipticCoeff.okl";
    kernelName = "setEllipticCoeff";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfo);
    kernelName = "setEllipticCoeffPressure";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfo);
  }
}
void registerCdsKernels() {
  const device_t &device = platform->device;
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
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
  const std::string oklpath = install_dir + "/okl/";
  const std::string section = "cds-";
  occa::properties meshProps = kernelInfo;
  meshProps += meshKernelProperties(N);
  {
    {
      occa::properties prop = meshProps;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;
      fileName = oklpath + "cds/advection" + suffix + ".okl";

      kernelName = "strongAdvectionVolume" + suffix;
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);

      kernelName = "strongAdvectionCubatureVolume" + suffix;
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "cds/advectMeshVelocityHex3D.okl";
    kernelName = "advectMeshVelocityHex3D";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "core/mask.okl";
    kernelName = "maskCopy";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

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

      fileName = oklpath + "cds/sumMakef.okl";
      kernelName = "sumMakef";
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
    }

    fileName = oklpath + "cds/helmholtzBC" + suffix + ".okl";
    kernelName = "helmholtzBC" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);
    kernelName = "dirichletBC";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfoBC);

    fileName = oklpath + "core/setEllipticCoeff.okl";
    kernelName = "setEllipticCoeff";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, kernelInfo);

    fileName = oklpath + "cds/regularization/filterRT" + suffix + ".okl";
    kernelName = "filterRT" + suffix;
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, meshProps);

    fileName = oklpath + "core/nStagesSum.okl";
    kernelName = "nStagesSum3";
    platform->kernels.add_kernel(
        section + kernelName, fileName, kernelName, platform->kernelInfo);

    {
      occa::properties prop = meshProps;
      const int movingMesh =
          platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/p_MovingMesh"] = movingMesh;
      prop["defines/p_nEXT"] = nEXT;
      prop["defines/p_nBDF"] = nBDF;
      prop["defines/p_cubNq"] = cubNq;
      prop["defines/p_cubNp"] = cubNp;

      fileName = oklpath + "cds/subCycle" + suffix + ".okl";
      occa::properties subCycleStrongCubatureProps = prop;
      if (platform->device.mode() == "Serial" ||
          platform->device.mode() == "OpenMP") {
        fileName = oklpath + "cds/subCycle" + suffix + ".c";
      }
      kernelName = "subCycleStrongCubatureVolume" + suffix;
      platform->kernels.add_kernel(section + kernelName,
          fileName,
          kernelName,
          subCycleStrongCubatureProps);
      fileName = oklpath + "cds/subCycle" + suffix + ".okl";
      kernelName = "subCycleStrongVolume" + suffix;
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);

      fileName = oklpath + "cds/subCycleRKUpdate.okl";
      kernelName = "subCycleERKUpdate";
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
      kernelName = "subCycleRK";
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);

      kernelName = "subCycleInitU0";
      platform->kernels.add_kernel(
          section + kernelName, fileName, kernelName, prop);
    }
  }
}
void registerCommonMGPreconditionerKernels(int N, occa::properties kernelInfo) {
  const std::string prefix = "Hex3D";
  std::string fileName, kernelName;

  kernelInfo["defines/pfloat"] = pfloatString;

  kernelInfo["defines/p_Nfields"] = 1;

  occa::properties pfloatKernelInfo = kernelInfo;
  pfloatKernelInfo["defines/dfloat"] = pfloatString;
  pfloatKernelInfo["defines/pfloat"] = pfloatString;

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));

  const std::string orderSuffix = std::string("_") + std::to_string(N);

  {
    const std::string oklpath = install_dir + "/okl/core/";
    std::string fileName;

    fileName = oklpath + "mask.okl";
    kernelName = "mask";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    fileName = oklpath + "mask.okl";
    platform->kernels.add_kernel(kernelName + orderSuffix + "pfloat",
        fileName,
        kernelName,
        pfloatKernelInfo,
        orderSuffix + "pfloat");
    fileName = install_dir + "/okl/elliptic/ellipticLinAlg.okl";
    kernelName = "fusedCopyDfloatToPfloat";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    kernelName = "copyDfloatToPfloat";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    kernelName = "copyPfloatToDfloat";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    kernelName = "scaledAdd";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    kernelName = "dotMultiply";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    fileName = install_dir + "/okl/elliptic/chebyshev.okl";
    kernelName = "updateSmoothedSolutionVec";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    kernelName = "updateChebyshevSolutionVec";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    kernelName = "updateIntermediateSolutionVec";
    platform->kernels.add_kernel(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
  }
}

void registerSchwarzKernels(const std::string &section, int N) {
  const std::string optionsPrefix = createOptionsPrefix(section);
  const int Nq = N + 1;
  const int Nq_e = Nq + 2;
  const int Np = Nq * Nq * Nq;
  const int Np_e = Nq_e * Nq_e * Nq_e;

  bool overlap = false;
  const bool serial = (platform->device.mode() == "Serial" ||
                       platform->device.mode() == "OpenMP");
  if (Nq >= 5 && !serial)
    overlap = true;

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  std::string fileName, kernelName;

  {
    occa::properties properties = platform->kernelInfo;
    properties["defines/p_Nq"] = Nq;
    properties["defines/p_Nq_e"] = Nq_e;
    properties["defines/p_restrict"] = 0;
    const std::string suffix =
        std::string("_") + std::to_string(Nq_e - 1) + std::string("pfloat");
    properties["defines/p_overlap"] = (int)overlap;
    if (platform->options.compareArgs(
            optionsPrefix + "MULTIGRID SMOOTHER", "RAS"))
      properties["defines/p_restrict"] = 1;

    fileName = oklpath + "ellipticSchwarzSolverHex3D.okl";
    if (serial) {
      fileName = oklpath + "ellipticSchwarzSolverHex3D.c";
    }
    platform->kernels.add_kernel(
        "preFDM" + suffix, fileName, "preFDM", properties, suffix);
    platform->kernels.add_kernel(
        "fusedFDM" + suffix, fileName, "fusedFDM", properties, suffix);
    platform->kernels.add_kernel(
        "postFDM" + suffix, fileName, "postFDM", properties, suffix);
  }
}
void registerFineLevelKernels(const std::string &section, int N) {
  auto gen_suffix = [N](const char *floatString) {
    const std::string precision = std::string(floatString);
    if (precision.find(pfloatString) != std::string::npos) {
      return std::string("_") + std::to_string(N) + std::string("pfloat");
    } else {
      return std::string("_") + std::to_string(N);
    }
  };

  auto kernelInfo = platform->kernelInfo + meshKernelProperties(N);
  registerCommonMGPreconditionerKernels(N, kernelInfo);

  const std::string suffix = "Hex3D";
  constexpr int Nfields{1};

  kernelInfo["defines/p_Nfields"] = Nfields;

  std::string fileName, kernelName;

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  const bool serial = platform->device.mode() == "Serial" ||
                      platform->device.mode() == "OpenMP";
  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  {
    occa::properties AxKernelInfo = kernelInfo;
    AxKernelInfo["defines/p_poisson"] = 1;

    kernelName = "ellipticAx" + suffix;
    fileName = oklpath + kernelName + fileNameExtension;
    {
      const std::string kernelSuffix = gen_suffix(dfloatString);
      platform->kernels.add_kernel(kernelName + kernelSuffix,
          fileName,
          kernelName,
          AxKernelInfo,
          kernelSuffix);
    }

    if (!strstr(pfloatString, dfloatString)) {
      AxKernelInfo["defines/dfloat"] = pfloatString;
      const std::string kernelSuffix = gen_suffix(pfloatString);
      platform->kernels.add_kernel(kernelName + kernelSuffix,
          fileName,
          kernelName,
          AxKernelInfo,
          kernelSuffix);
      AxKernelInfo["defines/dfloat"] = dfloatString;
    }

    if (platform->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
      kernelName = "ellipticPartialAxTrilinear" + suffix;
    else
      kernelName = "ellipticPartialAx" + suffix;

    fileName = oklpath + kernelName + fileNameExtension;

    if (!serial) {
      {
        const std::string kernelSuffix = gen_suffix(dfloatString);
        platform->kernels.add_kernel(kernelName + kernelSuffix,
            fileName,
            kernelName,
            AxKernelInfo,
            kernelSuffix);
      }
      if (!strstr(pfloatString, dfloatString)) {
        AxKernelInfo["defines/dfloat"] = pfloatString;
        const std::string kernelSuffix = gen_suffix(pfloatString);
        platform->kernels.add_kernel(kernelName + kernelSuffix,
            fileName,
            kernelName,
            AxKernelInfo,
            kernelSuffix);
        AxKernelInfo["defines/dfloat"] = dfloatString;
      }
    }
  }

  registerSchwarzKernels(section, N);
}
void registerSEMFEMKernels(const std::string &section, int N);

void registerLevelKernels(const std::string &section, int Nf, int N) {
  const int Nc = N;
  auto gen_suffix = [N](const char *floatString) {
    const std::string precision = std::string(floatString);
    if (precision.find(pfloatString) != std::string::npos) {
      return std::string("_") + std::to_string(N) + std::string("pfloat");
    } else {
      return std::string("_") + std::to_string(N);
    }
  };

  occa::properties kernelInfo = platform->kernelInfo + meshKernelProperties(N);

  const std::string suffix = "Hex3D";

  std::string fileName, kernelName;

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  registerCommonMGPreconditionerKernels(N, kernelInfo);

  const bool serial = platform->device.mode() == "Serial" ||
                      platform->device.mode() == "OpenMP";

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  constexpr int elementType = HEXAHEDRA;

  {
    occa::properties AxKernelInfo = kernelInfo;
    AxKernelInfo["defines/p_poisson"] = 1;

    kernelName = "ellipticAx" + suffix;
    fileName = oklpath + kernelName + fileNameExtension;

    {
      const std::string kernelSuffix = gen_suffix(dfloatString);
      platform->kernels.add_kernel(kernelName + kernelSuffix,
          fileName,
          kernelName,
          AxKernelInfo,
          kernelSuffix);
    }
    if (!strstr(pfloatString, dfloatString)) {
      AxKernelInfo["defines/dfloat"] = pfloatString;
      {
        const std::string kernelSuffix = gen_suffix(pfloatString);
        platform->kernels.add_kernel(kernelName + kernelSuffix,
            fileName,
            kernelName,
            AxKernelInfo,
            kernelSuffix);
      }
      AxKernelInfo["defines/dfloat"] = dfloatString;
    }

    if (platform->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
      kernelName = "ellipticPartialAxTrilinear" + suffix;
    else
      kernelName = "ellipticPartialAx" + suffix;

    fileName = oklpath + kernelName + fileNameExtension;

    if (!serial) {
      {
        const std::string kernelSuffix = gen_suffix(dfloatString);
        platform->kernels.add_kernel(kernelName + kernelSuffix,
            fileName,
            kernelName,
            AxKernelInfo,
            kernelSuffix);
      }
      if (!strstr(pfloatString, dfloatString)) {
        AxKernelInfo["defines/dfloat"] = pfloatString;
        const std::string kernelSuffix = gen_suffix(pfloatString);
        platform->kernels.add_kernel(kernelName + kernelSuffix,
            fileName,
            kernelName,
            AxKernelInfo,
            kernelSuffix);
        AxKernelInfo["defines/dfloat"] = dfloatString;
      }
    }
  }

  {
    fileName = oklpath + "ellipticBlockJacobiPrecon.okl";
    kernelName = "ellipticBlockJacobiPrecon";
    // sizes for the coarsen and prolongation kernels. degree NFine to degree N
    int NqFine = (Nf + 1);
    int NqCoarse = (Nc + 1);
    occa::properties coarsenProlongateKernelInfo = kernelInfo;
    coarsenProlongateKernelInfo["defines/p_NqFine"] = Nf + 1;
    coarsenProlongateKernelInfo["defines/p_NqCoarse"] = Nc + 1;

    const int NpFine = (Nf + 1) * (Nf + 1) * (Nf + 1);
    const int NpCoarse = (Nc + 1) * (Nc + 1) * (Nc + 1);
    coarsenProlongateKernelInfo["defines/p_NpFine"] = NpFine;
    coarsenProlongateKernelInfo["defines/p_NpCoarse"] = NpCoarse;

    const std::string orderSuffix = std::string("_") + std::to_string(Nf);

    if (serial) {
      fileName = oklpath + "ellipticPreconCoarsen" + suffix + ".c";
      kernelName = "ellipticPreconCoarsen" + suffix;
      platform->kernels.add_kernel(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
      fileName = oklpath + "ellipticPreconProlongate" + suffix + ".c";
      kernelName = "ellipticPreconProlongate" + suffix;
      platform->kernels.add_kernel(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
    } else {
      fileName = oklpath + "ellipticPreconCoarsen" + suffix + ".okl";
      kernelName = "ellipticPreconCoarsen" + suffix;
      platform->kernels.add_kernel(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
      fileName = oklpath + "ellipticPreconProlongate" + suffix + ".okl";
      kernelName = "ellipticPreconProlongate" + suffix;
      platform->kernels.add_kernel(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
    }
  }
  { registerSchwarzKernels(section, N); }
}
void registerMultiGridKernels(const std::string &section) {
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const std::string optionsPrefix = createOptionsPrefix(section);

  registerFineLevelKernels(section, N);

  auto levels = determineMGLevels(section);

  for (unsigned levelIndex = 1U; levelIndex < levels.size(); ++levelIndex) {
    const int levelFine = levels[levelIndex - 1];
    const int levelCoarse = levels[levelIndex];
    registerLevelKernels(section, levelFine, levelCoarse);
  }
  const int coarseLevel = levels.back();
  if (platform->options.compareArgs(
          optionsPrefix + "MULTIGRID COARSE SOLVE", "TRUE")) {
    if (platform->options.compareArgs(
            optionsPrefix + "MULTIGRID COARSE SEMFEM", "TRUE")) {
      registerSEMFEMKernels(section, coarseLevel);
    } else {
      {
        std::string install_dir;
        install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
        const std::string oklpath = install_dir + "/okl/";
        std::string fileName = oklpath + "parAlmond/convertFP64ToFP32.okl";
        std::string kernelName = "convertFP64ToFP32";
        platform->kernels.add_kernel(
            kernelName, fileName, kernelName, platform->kernelInfo);

        fileName = oklpath + "parAlmond/convertFP32ToFP64.okl";
        kernelName = "convertFP32ToFP64";
        platform->kernels.add_kernel(
            kernelName, fileName, kernelName, platform->kernelInfo);
        fileName = oklpath + "parAlmond/vectorDotStar.okl";
        kernelName = "vectorDotStar2";
        platform->kernels.add_kernel(
            kernelName, fileName, kernelName, platform->kernelInfo);
      }
    }
  }
}
void registerSEMFEMKernels(const std::string &section, int N) {
  const int Nq = N + 1;
  const int Np = Nq * Nq * Nq;
  const std::string optionsPrefix = createOptionsPrefix(section);
  const int useFP32 = platform->options.compareArgs(
      optionsPrefix + "SEMFEM SOLVER PRECISION", "FP32");
  occa::properties SEMFEMKernelProps = platform->kernelInfo;
  if (useFP32) {
    SEMFEMKernelProps["defines/pfloat"] = "float";
  } else {
    SEMFEMKernelProps["defines/pfloat"] = "double";
  }
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  std::string fileName = oklpath + "ellipticGather.okl";
  platform->kernels.add_kernel("gather", fileName, "gather", SEMFEMKernelProps);
  fileName = oklpath + "ellipticScatter.okl";
  platform->kernels.add_kernel(
      "scatter", fileName, "scatter", SEMFEMKernelProps);
  occa::properties stiffnessKernelInfo = platform->kernelInfo;
  fileName = oklpath + "ellipticSEMFEMStiffness.okl";
  stiffnessKernelInfo["defines/p_Nq"] = Nq;
  stiffnessKernelInfo["defines/p_Np"] = Np;
  stiffnessKernelInfo["defines/p_rows_sorted"] = 1;
  stiffnessKernelInfo["defines/p_cols_sorted"] = 0;

  const bool constructOnHost =
      platform->device.mode() == std::string("OpenCL") ||
      platform->device.mode() == std::string("HIP") ||
      platform->device.mode() == std::string("Serial");

  if (!constructOnHost) {
    platform->kernels.add_kernel("computeStiffnessMatrix",
        fileName,
        "computeStiffnessMatrix",
        stiffnessKernelInfo);
  }
}
void registerJacobiKernels(const std::string &section) {
  const std::string optionsPrefix = createOptionsPrefix(section);
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/";
  std::string fileName = oklpath + "elliptic/ellipticJacobi.okl";
  std::string kernelName = "axmyzManyPfloat";
  platform->kernels.add_kernel(
    kernelName, fileName, kernelName, platform->kernelInfo);

  kernelName = "adyManyPfloat";
  platform->kernels.add_kernel(
    kernelName, fileName, kernelName, platform->kernelInfo);
}
void registerEllipticPreconditionerKernels(const std::string &section) {
  const std::string optionsPrefix = createOptionsPrefix(section);
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);

  if(platform->options.compareArgs(optionsPrefix + "PRECONDITIONER", "MULTIGRID")) {
    registerMultiGridKernels(section);
  } else if(platform->options.compareArgs(optionsPrefix + "PRECONDITIONER", "SEMFEM")) {
    registerSEMFEMKernels(section, N);
  } else if(platform->options.compareArgs(optionsPrefix + "PRECONDITIONER", "JACOBI")) {
    registerJacobiKernels(section);
  } else if(platform->options.compareArgs(optionsPrefix + "PRECONDITIONER", "NONE")) {
    // nothing 
  } else {
    printf("ERROR: Unknown preconditioner!\n");
    ABORT(EXIT_FAILURE);
  }
}
void registerEllipticKernels(const std::string &section) {
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const std::string optionsPrefix = createOptionsPrefix(section);

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  occa::properties kernelInfo = platform->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();
  kernelInfo += meshKernelProperties(N);

  const bool blockSolver = [&section]() {
    if (section.find("velocity") == std::string::npos)
      return false;
    if (platform->options.compareArgs("STRESSFORMULATION", "TRUE"))
      return true;
    if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"))
      return true;
    return false;
  }();
  const int Nfields = (blockSolver) ? 3 : 1;
  const bool stressForm = [&section]() {
    if (section.find("velocity") == std::string::npos)
      return false;
    if (platform->options.compareArgs("STRESSFORMULATION", "TRUE"))
      return true;
    return false;
  }();

  const bool serial = platform->device.mode() == "Serial" ||
                      platform->device.mode() == "OpenMP";

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  const std::string sectionIdentifier = std::to_string(Nfields) + "-";

  if (platform->options.compareArgs(
          optionsPrefix + "KRYLOV SOLVER", "PGMRES")) {
    registerGMRESKernels(section, Nfields);
  }

  // solution projection kernels
  {
    const std::string oklpath = install_dir + "/okl/elliptic/";
    std::string fileName, kernelName;

    {
      occa::properties properties = platform->kernelInfo;
      properties["defines/p_Nfields"] = Nfields;

      fileName = oklpath + "ellipticResidualProjection.okl";
      kernelName = "multiScaledAddwOffset";
      platform->kernels.add_kernel(
          sectionIdentifier + kernelName, fileName, kernelName, properties);
      kernelName = "accumulate";
      platform->kernels.add_kernel(
          sectionIdentifier + kernelName, fileName, kernelName, properties);
    }
  }

  {
    const std::string oklpath = install_dir + "/okl/core/";
    std::string fileName;

    fileName = oklpath + "mask.okl";
    platform->kernels.add_kernel("mask", fileName, "mask", kernelInfo);
  }

  kernelInfo["defines/p_Nfields"] = Nfields;

  occa::properties dfloatKernelInfo = kernelInfo;
  occa::properties floatKernelInfo = kernelInfo;
  floatKernelInfo["defines/pfloat"] = pfloatString;
  floatKernelInfo["defines/dfloat"] = pfloatString;

  constexpr bool var_coeff = true;
  constexpr int elementType{HEXAHEDRA};

  const std::string suffix = "Hex3D";

  occa::properties AxKernelInfo = dfloatKernelInfo;
  {
    const std::string oklpath = install_dir + "/okl/elliptic/";
    std::string fileName;
    std::string kernelName;

    fileName = oklpath + "ellipticBuildDiagonal" + suffix + ".okl";
    kernelName = "ellipticBlockBuildDiagonal" + suffix;
    dfloatKernelInfo["defines/dfloat"] = dfloatString;
    dfloatKernelInfo["defines/pfloat"] = pfloatString;
    platform->kernels.add_kernel(
        sectionIdentifier + kernelName, fileName, kernelName, dfloatKernelInfo);

    // Ax
    dfloatKernelInfo["defines/pfloat"] = dfloatString;

    std::string kernelNamePrefix = "elliptic";
    if (blockSolver)
      kernelNamePrefix += (stressForm) ? "Stress" : "Block";

    kernelName = "Ax";
    if (var_coeff) kernelName += "Var";
    if (platform->options.compareArgs("ELEMENT MAP", "TRILINEAR")) kernelName += "Trilinear";
    kernelName += suffix; 
    if (blockSolver && !stressForm) kernelName += "_N" + std::to_string(Nfields);

    {
      std::string _kernelName = kernelNamePrefix + kernelName;
      fileName = oklpath + _kernelName + fileNameExtension; 
      platform->kernels.add_kernel(
        _kernelName, fileName, _kernelName, AxKernelInfo);
    }

    if (!serial) {
      std::string _kernelName = kernelNamePrefix + "Partial" + kernelName;
      fileName = oklpath + _kernelName + fileNameExtension; 
      platform->kernels.add_kernel(
        _kernelName, fileName, _kernelName, AxKernelInfo);
    }

    // PCG update
    if (serial) {
      fileName = oklpath + "ellipticSerialUpdatePCG.c";
    } else {
      fileName = oklpath + "ellipticUpdatePCG.okl";
    }
    platform->kernels.add_kernel(sectionIdentifier + "ellipticBlockUpdatePCG",
        fileName,
        "ellipticBlockUpdatePCG",
        dfloatKernelInfo);
  }

  // projection
  {}
}
void registerMeshKernels() {
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
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  std::string oklpath = install_dir + "/okl/";
  std::string kernelName;

  const std::string meshPrefix = "mesh-";
  {
    std::string fileName = oklpath + "mesh/velocityBCHex3D.okl";
    kernelName = "velocityDirichletBCHex3D";
    platform->kernels.add_kernel(meshPrefix + kernelName, fileName, kernelName, kernelInfo);
    occa::properties meshKernelInfo = kernelInfo;
    meshKernelInfo["defines/p_cubNq"] = cubNq;
    meshKernelInfo["defines/p_cubNp"] = cubNp;

    fileName = oklpath + "mesh/geometricFactorsHex3D.okl";
    kernelName = "geometricFactorsHex3D";
    platform->kernels.add_kernel(
        meshPrefix + kernelName, fileName, kernelName, meshKernelInfo);
    fileName = oklpath + "mesh/surfaceGeometricFactorsHex3D.okl";
    kernelName = "surfaceGeometricFactorsHex3D";
    platform->kernels.add_kernel(
        meshPrefix + kernelName, fileName, kernelName, meshKernelInfo);

    meshKernelInfo = kernelInfo;
    meshKernelInfo["defines/p_nAB"] = nAB;
    fileName = oklpath + "core/nStagesSum.okl";
    kernelName = "nStagesSumVector";
    platform->kernels.add_kernel(
        meshPrefix + kernelName, fileName, kernelName, meshKernelInfo);
  }
}

void registerLinAlgKernels() {
  occa::properties kernelInfo = platform->kernelInfo;

  std::string oklDir;
  oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
  oklDir += "/okl/linAlg/";
  std::string fileName;
  const bool serial = (platform->device.mode() == "Serial" ||
                       platform->device.mode() == "OpenMP");

  platform->kernels.add_kernel(
      "fill", oklDir + "linAlgFill.okl", "fill", kernelInfo);
  platform->kernels.add_kernel(
      "vabs", oklDir + "linAlgAbs.okl", "vabs", kernelInfo);
  platform->kernels.add_kernel(
      "add", oklDir + "linAlgAdd.okl", "add", kernelInfo);
  platform->kernels.add_kernel(
      "scale", oklDir + "linAlgScale.okl", "scale", kernelInfo);
  platform->kernels.add_kernel(
      "scaleMany", oklDir + "linAlgScale.okl", "scaleMany", kernelInfo);
  fileName = std::string("linAlgAXPBY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel("axpby", oklDir + fileName, "axpby", kernelInfo);
  fileName = std::string("linAlgAXPBY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "axpbyMany", oklDir + fileName, "axpbyMany", kernelInfo);
  platform->kernels.add_kernel(
      "axpbyz", oklDir + "linAlgAXPBY.okl", "axpbyz", kernelInfo);
  platform->kernels.add_kernel(
      "axpbyzMany", oklDir + "linAlgAXPBY.okl", "axpbyzMany", kernelInfo);
  fileName = std::string("linAlgAXMY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel("axmy", oklDir + fileName, "axmy", kernelInfo);
  fileName = std::string("linAlgAXMY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "axmyMany", oklDir + fileName, "axmyMany", kernelInfo);
  fileName = std::string("linAlgAXMY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "axmyVector", oklDir + fileName, "axmyVector", kernelInfo);
  platform->kernels.add_kernel(
      "axmyz", oklDir + "linAlgAXMY.okl", "axmyz", kernelInfo);
  platform->kernels.add_kernel(
      "axmyzMany", oklDir + "linAlgAXMY.okl", "axmyzMany", kernelInfo);
  platform->kernels.add_kernel(
      "ady", oklDir + "linAlgAXDY.okl", "ady", kernelInfo);
  platform->kernels.add_kernel(
      "adyMany", oklDir + "linAlgAXDY.okl", "adyMany", kernelInfo);
  platform->kernels.add_kernel(
      "axdy", oklDir + "linAlgAXDY.okl", "axdy", kernelInfo);
  platform->kernels.add_kernel(
      "aydx", oklDir + "linAlgAXDY.okl", "aydx", kernelInfo);
  platform->kernels.add_kernel(
      "aydxMany", oklDir + "linAlgAXDY.okl", "aydxMany", kernelInfo);
  platform->kernels.add_kernel(
      "axdyz", oklDir + "linAlgAXDY.okl", "axdyz", kernelInfo);
  platform->kernels.add_kernel(
      "sum", oklDir + "linAlgSum.okl", "sum", kernelInfo);
  platform->kernels.add_kernel(
      "sumMany", oklDir + "linAlgSum.okl", "sumMany", kernelInfo);
  platform->kernels.add_kernel(
      "min", oklDir + "linAlgMin.okl", "min", kernelInfo);
  platform->kernels.add_kernel(
      "max", oklDir + "linAlgMax.okl", "max", kernelInfo);
  fileName = std::string("linAlgNorm2") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "norm2", oklDir + fileName, "norm2", kernelInfo);
  platform->kernels.add_kernel(
      "norm2Many", oklDir + fileName, "norm2Many", kernelInfo);
  fileName = std::string("linAlgNorm1") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "norm1", oklDir + fileName, "norm1", kernelInfo);
  platform->kernels.add_kernel(
      "norm1Many", oklDir + fileName, "norm1Many", kernelInfo);
  fileName = std::string("linAlgWeightedNorm1") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "weightedNorm1", oklDir + fileName, "weightedNorm1", kernelInfo);
  platform->kernels.add_kernel(
      "weightedNorm1Many", oklDir + fileName, "weightedNorm1Many", kernelInfo);
  fileName = std::string("linAlgWeightedNorm2") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "weightedNorm2", oklDir + fileName, "weightedNorm2", kernelInfo);
  fileName = std::string("linAlgWeightedNorm2") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "weightedNorm2Many", oklDir + fileName, "weightedNorm2Many", kernelInfo);
  platform->kernels.add_kernel(
      "innerProd", oklDir + "linAlgInnerProd.okl", "innerProd", kernelInfo);
  fileName = std::string("linAlgWeightedInnerProd") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel(
      "weightedInnerProd", oklDir + fileName, "weightedInnerProd", kernelInfo);
  fileName = std::string("linAlgWeightedInnerProd") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add_kernel("weightedInnerProdMany",
      oklDir + fileName,
      "weightedInnerProdMany",
      kernelInfo);
  platform->kernels.add_kernel("weightedInnerProdMulti",
      oklDir + "linAlgWeightedInnerProd.okl",
      "weightedInnerProdMulti",
      kernelInfo);
}
void compileUDFKernels()
{
  int buildNodeLocal = 0;
  if (getenv("NEKRS_CACHE_LOCAL"))
    buildNodeLocal = std::stoi(getenv("NEKRS_CACHE_LOCAL"));

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  occa::properties kernelInfo = platform->kernelInfo + meshKernelProperties(N);
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  auto rank = buildNodeLocal ? platform->comm.localRank : platform->comm.mpiRank;
  auto communicator = buildNodeLocal ? platform->comm.mpiCommLocal : platform->comm.mpiComm;

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0)
    printf("loading udf kernels ... ");
  fflush(stdout);

  for(int pass = 0; pass < 2; ++pass)
  {
    bool executePass = (pass == 0) && (rank == 0);
    executePass |= (pass == 1) && (rank != 0);
    if(executePass){
      kernelInfoBC = kernelInfo;
      if (udf.loadKernels) {
        // side-effect: kernelInfoBC will include any relevant user-defined kernel props
        udf.loadKernels(kernelInfoBC);
      }
      const std::string bcDataFile = install_dir + "/include/core/bcData.h";
      kernelInfoBC["includes"] += bcDataFile.c_str();
      std::string boundaryHeaderFileName;
      platform->options.getArgs("DATA FILE", boundaryHeaderFileName);
      kernelInfoBC["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

      kernelInfoBC += meshKernelProperties(N);
    }
    MPI_Barrier(communicator);
  }

  MPI_Barrier(platform->comm.mpiComm);
  const double loadTime = MPI_Wtime() - tStart;
  if (platform->comm.mpiRank == 0)
    printf("done (%gs)\n\n", loadTime);
  fflush(stdout);
}
void compileDummyKernel()
{
  int buildNodeLocal = 0;
  if (getenv("NEKRS_CACHE_LOCAL"))
    buildNodeLocal = std::stoi(getenv("NEKRS_CACHE_LOCAL"));
  auto rank = buildNodeLocal ? platform->comm.localRank : platform->comm.mpiRank;
  const std::string dummyKernelName = "myDummyKernelName";
  const std::string dummyKernelStr = std::string(
      "@kernel void myDummyKernelName(int N) {"
      "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {}"
      "}"
  );

  if(rank == 0){
    platform->device.buildKernelFromString(
      dummyKernelStr,
      dummyKernelName,
      platform->kernelInfo
    );
  }

}
} // namespace

void compileKernels() {

  compileDummyKernel(); // trigger occa's compilerVendorTest

  compileUDFKernels();

  registerLinAlgKernels();

  registerMeshKernels();

  registerNrsKernels();

  {
    int Nscalars;
    platform->options.getArgs("NUMBER OF SCALARS", Nscalars);
    if (Nscalars) {
      registerCdsKernels();
    }
  }

  {
    const std::vector<std::string> sections = {
        "pressure",
        "velocity",
    };
    for (auto &&section : sections) {
      registerEllipticKernels(section);
      registerEllipticPreconditionerKernels(section);
    }
  }

  {
    MPI_Barrier(platform->comm.mpiComm);
    const double tStart = MPI_Wtime();
    if (platform->comm.mpiRank == 0)
      printf("loading kernels ... ");
    fflush(stdout);

    {
      int buildNodeLocal = 0;
      if (getenv("NEKRS_CACHE_LOCAL"))
        buildNodeLocal = std::stoi(getenv("NEKRS_CACHE_LOCAL"));
      const bool buildOnly = platform->options.compareArgs("BUILD ONLY", "TRUE");
      auto communicator = buildNodeLocal ? platform->comm.mpiCommLocal : platform->comm.mpiComm;
      oogs::compile(
          platform->device, platform->device.mode(), communicator, buildOnly);
    }

    platform->kernels.compile();

    MPI_Barrier(platform->comm.mpiComm);
    const double loadTime = MPI_Wtime() - tStart;


    fflush(stdout);
    if (platform->comm.mpiRank == 0) {
      std::ofstream ofs;
      ofs.open(occa::env::OCCA_CACHE_DIR + "cache/compile.timestamp", 
	       std::ofstream::out | std::ofstream::trunc);
      ofs.close();
    }
 
    platform->timer.set("loadKernels", loadTime);
    if (platform->comm.mpiRank == 0)
      printf("done (%gs)\n\n", loadTime);
    fflush(stdout);
  }
}
