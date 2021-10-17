#include <compileKernels.hpp>
#include "nrs.hpp"
#include "elliptic.h"

namespace {

void registerJacobiKernels(const std::string &section) {
  const std::string optionsPrefix = createOptionsPrefix(section);
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/";
  std::string fileName = oklpath + "elliptic/ellipticJacobi.okl";
  std::string kernelName = "axmyzManyPfloat";
  platform->kernels.add(
    kernelName, fileName, kernelName, platform->kernelInfo);

  kernelName = "adyManyPfloat";
  platform->kernels.add(
    kernelName, fileName, kernelName, platform->kernelInfo);
}

void registerCommonMGPreconditionerKernels(int N, occa::properties kernelInfo) {
  const std::string prefix = "Hex3D";
  std::string fileName, kernelName;

  kernelInfo["defines/pfloat"] = pfloatString;

  kernelInfo["defines/p_Nfields"] = 1;

  occa::properties pfloatKernelInfo = kernelInfo;
  pfloatKernelInfo["defines/dfloat"] = pfloatString;
  pfloatKernelInfo["defines/pfloat"] = pfloatString;

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));

  const std::string orderSuffix = std::string("_") + std::to_string(N);

  {
    const std::string oklpath = installDir + "/okl/core/";
    std::string fileName;

    fileName = oklpath + "mask.okl";
    kernelName = "mask";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    fileName = oklpath + "mask.okl";
    platform->kernels.add(kernelName + orderSuffix + "pfloat",
        fileName,
        kernelName,
        pfloatKernelInfo,
        orderSuffix + "pfloat");
    fileName = installDir + "/okl/elliptic/ellipticLinAlg.okl";
    kernelName = "fusedCopyDfloatToPfloat";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    kernelName = "copyDfloatToPfloat";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    kernelName = "copyPfloatToDfloat";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    kernelName = "scaledAdd";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    kernelName = "dotMultiply";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    fileName = installDir + "/okl/elliptic/chebyshev.okl";
    kernelName = "updateSmoothedSolutionVec";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);
    kernelName = "updateChebyshevSolutionVec";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelName,
        kernelInfo,
        orderSuffix);

    kernelName = "updateIntermediateSolutionVec";
    platform->kernels.add(kernelName + orderSuffix,
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
  const bool serial = useSerial();
  if (Nq >= 5 && !serial)
    overlap = true;

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/elliptic/";
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
    platform->kernels.add(
        "preFDM" + suffix, fileName, "preFDM", properties, suffix);
    platform->kernels.add(
        "fusedFDM" + suffix, fileName, "fusedFDM", properties, suffix);
    platform->kernels.add(
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

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/elliptic/";
  const bool serial = useSerial();
  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  {
    occa::properties AxKernelInfo = kernelInfo;
    AxKernelInfo["defines/p_poisson"] = 1;

    kernelName = "ellipticAx" + suffix;
    fileName = oklpath + kernelName + fileNameExtension;
    {
      const std::string kernelSuffix = gen_suffix(dfloatString);
      platform->kernels.add(kernelName + kernelSuffix,
          fileName,
          kernelName,
          AxKernelInfo,
          kernelSuffix);
    }

    if (!strstr(pfloatString, dfloatString)) {
      AxKernelInfo["defines/dfloat"] = pfloatString;
      const std::string kernelSuffix = gen_suffix(pfloatString);
      platform->kernels.add(kernelName + kernelSuffix,
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
        platform->kernels.add(kernelName + kernelSuffix,
            fileName,
            kernelName,
            AxKernelInfo,
            kernelSuffix);
      }
      if (!strstr(pfloatString, dfloatString)) {
        AxKernelInfo["defines/dfloat"] = pfloatString;
        const std::string kernelSuffix = gen_suffix(pfloatString);
        platform->kernels.add(kernelName + kernelSuffix,
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

void registerMultigridLevelKernels(const std::string &section, int Nf, int N) {
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

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/elliptic/";
  registerCommonMGPreconditionerKernels(N, kernelInfo);

  const bool serial = useSerial();

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  constexpr int elementType = HEXAHEDRA;

  {
    occa::properties AxKernelInfo = kernelInfo;
    AxKernelInfo["defines/p_poisson"] = 1;

    kernelName = "ellipticAx" + suffix;
    fileName = oklpath + kernelName + fileNameExtension;

    {
      const std::string kernelSuffix = gen_suffix(dfloatString);
      platform->kernels.add(kernelName + kernelSuffix,
          fileName,
          kernelName,
          AxKernelInfo,
          kernelSuffix);
    }
    if (!strstr(pfloatString, dfloatString)) {
      AxKernelInfo["defines/dfloat"] = pfloatString;
      {
        const std::string kernelSuffix = gen_suffix(pfloatString);
        platform->kernels.add(kernelName + kernelSuffix,
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
        platform->kernels.add(kernelName + kernelSuffix,
            fileName,
            kernelName,
            AxKernelInfo,
            kernelSuffix);
      }
      if (!strstr(pfloatString, dfloatString)) {
        AxKernelInfo["defines/dfloat"] = pfloatString;
        const std::string kernelSuffix = gen_suffix(pfloatString);
        platform->kernels.add(kernelName + kernelSuffix,
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
      platform->kernels.add(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
      fileName = oklpath + "ellipticPreconProlongate" + suffix + ".c";
      kernelName = "ellipticPreconProlongate" + suffix;
      platform->kernels.add(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
    } else {
      fileName = oklpath + "ellipticPreconCoarsen" + suffix + ".okl";
      kernelName = "ellipticPreconCoarsen" + suffix;
      platform->kernels.add(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
      fileName = oklpath + "ellipticPreconProlongate" + suffix + ".okl";
      kernelName = "ellipticPreconProlongate" + suffix;
      platform->kernels.add(kernelName + orderSuffix,
          fileName,
          kernelName,
          coarsenProlongateKernelInfo,
          orderSuffix);
    }
  }
  registerSchwarzKernels(section, N);
}
void registerMultiGridKernels(const std::string &section) {
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const std::string optionsPrefix = createOptionsPrefix(section);

  registerFineLevelKernels(section, N);

  std::vector<int> levels = determineMGLevels(section);

  for (unsigned levelIndex = 1U; levelIndex < levels.size(); ++levelIndex) {
    const int levelFine = levels[levelIndex - 1];
    const int levelCoarse = levels[levelIndex];
    registerMultigridLevelKernels(section, levelFine, levelCoarse);
  }
  const int coarseLevel = levels.back();
  if (platform->options.compareArgs(
          optionsPrefix + "MULTIGRID COARSE SOLVE", "TRUE")) {
    if (platform->options.compareArgs(
            optionsPrefix + "MULTIGRID COARSE SEMFEM", "TRUE")) {
      registerSEMFEMKernels(section, coarseLevel);
    } else {
      {
        std::string installDir;
        installDir.assign(getenv("NEKRS_INSTALL_DIR"));
        const std::string oklpath = installDir + "/okl/";
        std::string fileName = oklpath + "parAlmond/convertFP64ToFP32.okl";
        std::string kernelName = "convertFP64ToFP32";
        platform->kernels.add(
            kernelName, fileName, kernelName, platform->kernelInfo);

        fileName = oklpath + "parAlmond/convertFP32ToFP64.okl";
        kernelName = "convertFP32ToFP64";
        platform->kernels.add(
            kernelName, fileName, kernelName, platform->kernelInfo);
        fileName = oklpath + "parAlmond/vectorDotStar.okl";
        kernelName = "vectorDotStar2";
        platform->kernels.add(
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
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/elliptic/";
  std::string fileName = oklpath + "ellipticGather.okl";
  platform->kernels.add("gather", fileName, "gather", SEMFEMKernelProps);
  fileName = oklpath + "ellipticScatter.okl";
  platform->kernels.add(
      "scatter", fileName, "scatter", SEMFEMKernelProps);
  occa::properties stiffnessKernelInfo = platform->kernelInfo;
  fileName = oklpath + "ellipticSEMFEMStiffness.okl";
  stiffnessKernelInfo["defines/p_Nq"] = Nq;
  stiffnessKernelInfo["defines/p_Np"] = Np;
  stiffnessKernelInfo["defines/p_rows_sorted"] = 1;
  stiffnessKernelInfo["defines/p_cols_sorted"] = 0;

  const bool constructOnHost = !supportsAtomicReductions();

  if (!constructOnHost) {
    platform->kernels.add("computeStiffnessMatrix",
        fileName,
        "computeStiffnessMatrix",
        stiffnessKernelInfo);
  }
}

}

void registerEllipticPreconditionerKernels(std::string section) {
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
