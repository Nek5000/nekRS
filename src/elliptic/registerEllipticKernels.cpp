#include <compileKernels.hpp>
#include "elliptic.h"

namespace{

void registerGMRESKernels(const std::string &section, int Nfields) {
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/elliptic/";
  std::string fileName;
  const bool serial = platform->serial;

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";
  const std::string sectionIdentifier = std::to_string(Nfields) + "-";

  occa::properties gmresKernelInfo = platform->kernelInfo;
  gmresKernelInfo["defines/p_Nfields"] = Nfields;

  std::string kernelName = "gramSchmidtOrthogonalization";
  fileName = oklpath + kernelName + fileNameExtension;
  platform->kernels.add(
      sectionIdentifier + kernelName, fileName, gmresKernelInfo);

  kernelName = "updatePGMRESSolution";
  fileName = oklpath + kernelName + fileNameExtension;
  platform->kernels.add(
      sectionIdentifier + kernelName, fileName, gmresKernelInfo);

  kernelName = "fusedResidualAndNorm";
  fileName = oklpath + kernelName + fileNameExtension;
  platform->kernels.add(
      sectionIdentifier + kernelName, fileName, gmresKernelInfo);
}

}

void registerEllipticKernels(std::string section) {
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const std::string optionsPrefix = createOptionsPrefix(section);

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
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

  const bool serial = platform->serial;

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  const std::string sectionIdentifier = std::to_string(Nfields) + "-";

  if (platform->options.compareArgs(
          optionsPrefix + "KRYLOV SOLVER", "PGMRES")) {
    registerGMRESKernels(section, Nfields);
  }

  // solution projection kernels
  {
    const std::string oklpath = installDir + "/okl/elliptic/";
    std::string fileName, kernelName;

    {
      const std::string extension = ".okl";
      occa::properties properties = platform->kernelInfo;
      properties["defines/p_Nfields"] = Nfields;

      kernelName = "multiScaledAddwOffset";
      fileName = oklpath + kernelName + extension;
      platform->kernels.add(
          sectionIdentifier + kernelName, fileName, properties);
      kernelName = "accumulate";
      fileName = oklpath + kernelName + extension;
      platform->kernels.add(
          sectionIdentifier + kernelName, fileName, properties);
    }
  }

  {
    const std::string oklpath = installDir + "/okl/core/";
    std::string fileName;

    fileName = oklpath + "mask.okl";
    platform->kernels.add("mask", fileName, kernelInfo);
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
    const std::string oklpath = installDir + "/okl/elliptic/";
    std::string fileName;
    std::string kernelName;

    kernelName = "ellipticBlockBuildDiagonal" + suffix;
    fileName = oklpath + kernelName + ".okl";
    dfloatKernelInfo["defines/dfloat"] = dfloatString;
    dfloatKernelInfo["defines/pfloat"] = pfloatString;
    platform->kernels.add(
        sectionIdentifier + kernelName, fileName, dfloatKernelInfo);

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
      platform->kernels.add(
        _kernelName, fileName, AxKernelInfo);
    }

    if (!serial) {
      std::string _kernelName = kernelNamePrefix + "Partial" + kernelName;
      fileName = oklpath + _kernelName + fileNameExtension; 
      platform->kernels.add(
        _kernelName, fileName, AxKernelInfo);
    }

    // PCG update
    fileName = oklpath + "ellipticBlockUpdatePCG" + fileNameExtension;
    platform->kernels.add(sectionIdentifier + "ellipticBlockUpdatePCG",
        fileName,
        kernelInfo);
  }
}