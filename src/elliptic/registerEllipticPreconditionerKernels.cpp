#include <compileKernels.hpp>
#include "nrs.hpp"
#include "elliptic.h"
#include "benchmarkFDM.hpp"
#include "benchmarkAx.hpp"

#include "re2Reader.hpp"

namespace {

void registerAxKernels(const std::string& section, int N, int poissonEquation)
{
  auto gen_suffix = [N](const char *floatString) {
    const std::string precision = std::string(floatString);
    if (precision.find(pfloatString) != std::string::npos) {
      return std::string("_") + std::to_string(N) + std::string("pfloat");
    } else {
      return std::string("_") + std::to_string(N);
    }
  };
  constexpr int Nfields{1};

  auto kernelInfo = platform->kernelInfo + meshKernelProperties(N);
  kernelInfo["defines/p_Nfields"] = Nfields;

  std::string fileName, kernelName;

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/elliptic/";
  const bool serial = platform->serial;
  const std::string fileNameExtension = (serial) ? ".c" : ".okl";
  const std::string poissonPrefix = poissonEquation ? "poisson-" : "";
  {
    int nelgt, nelgv;
    const std::string meshFile = platform->options.getArgs("MESH FILE");
    re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);
    const int NelemBenchmark = nelgv/platform->comm.mpiCommSize;

    occa::properties AxKernelInfo = kernelInfo;
    const auto Nq = N+1;
    for(auto&& coeffField : {true,false}){
      for(auto&& floatString : {std::string(dfloatString), std::string(pfloatString)}){

        dlong wordSize = 8;
        if(floatString.find("float") != std::string::npos){
          wordSize = 4;
        }

        bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");
        const int verbosity = verbose ? 2 : 1;
        const std::string kernelSuffix = gen_suffix(floatString.c_str());
        auto axKernel = benchmarkAx(NelemBenchmark,
                                    Nq,
                                    Nq - 1,
                                    !coeffField,
                                    poissonEquation,
                                    false,
                                    wordSize,
                                    Nfields,
                                    false, // no stress formulation in preconditioner
                                    verbosity,
                                    0.2,
                                    false,
                                    kernelSuffix);

        const std::string suffix = coeffField ? "CoeffHex3D" : "Hex3D";

        if (platform->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
          kernelName = "ellipticPartialAxTrilinear" + suffix;
        else
          kernelName = "ellipticPartialAx" + suffix;

        fileName = oklpath + kernelName + fileNameExtension;

        platform->kernels.add(poissonPrefix + kernelName + kernelSuffix,
          axKernel);
      }
    }
  }
}

void registerJacobiKernels(const std::string &section, int poissonEquation) {
  const bool serial = platform->serial;
  const std::string extension = serial ? ".c" : ".okl";
  const std::string optionsPrefix = createOptionsPrefix(section);
  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/";
  occa::properties pfloatProps = platform->kernelInfo;
  pfloatProps["defines/dfloat"] = pfloatString;

  // This kernel is needed as it used for mixed-precision Jacobi preconditioning 
  std::string kernelName = "axmyzManyPfloat";
  std::string fileName = oklpath + "elliptic/" + kernelName + extension;
  platform->kernels.add(
    kernelName, fileName, platform->kernelInfo);

  kernelName = "adyManyPfloat";
  fileName = oklpath + "linAlg/adyMany.okl";
  platform->kernels.add(
    kernelName, fileName, pfloatProps);
}

void registerCommonMGPreconditionerKernels(int N, occa::properties kernelInfo, int poissonEquation) {
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

  const bool serial = platform->serial;
  const std::string extension = serial ? ".c" : ".okl";

  {
    const std::string oklpath = installDir + "/okl/core/";
    std::string fileName;

    fileName = oklpath + "mask.okl";
    kernelName = "mask";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);

    fileName = oklpath + "mask.okl";
    platform->kernels.add(kernelName + orderSuffix + "pfloat",
        fileName,
        pfloatKernelInfo,
        orderSuffix + "pfloat");
    kernelName = "fusedCopyDfloatToPfloat";
    fileName = installDir + "/okl/elliptic/" + kernelName + ".okl";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);

    kernelName = "copyDfloatToPfloat";
    fileName = installDir + "/okl/core/" + kernelName + extension;
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);

    kernelName = "copyPfloatToDfloat";
    fileName = installDir + "/okl/core/" + kernelName + extension;
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);

    kernelName = "scaledAdd";
    fileName = installDir + "/okl/elliptic/" + kernelName + ".okl";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);
    kernelName = "dotMultiply";
    fileName = installDir + "/okl/elliptic/" + kernelName + ".okl";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);

    kernelName = "updateSmoothedSolutionVec";
    fileName = installDir + "/okl/elliptic/" + kernelName + ".okl";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);
    kernelName = "updateChebyshevSolutionVec";
    fileName = installDir + "/okl/elliptic/" + kernelName + ".okl";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);

    kernelName = "updateIntermediateSolutionVec";
    fileName = installDir + "/okl/elliptic/" + kernelName + ".okl";
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        kernelInfo,
        orderSuffix);

    occa::properties buildDiagInfo = kernelInfo;
    if(poissonEquation) buildDiagInfo["defines/p_poisson"] = 1;
    const std::string poissonPrefix = poissonEquation ? "poisson-" : "";
    kernelName = "ellipticBlockBuildDiagonalHex3D";
    fileName = installDir + "/okl/elliptic/" + kernelName + ".okl";
    platform->kernels.add(poissonPrefix + kernelName + orderSuffix, fileName, buildDiagInfo, orderSuffix);
  }
}

void registerSchwarzKernels(const std::string &section, int N) {
  const std::string optionsPrefix = createOptionsPrefix(section);
  const int Nq = N + 1;
  const int Nq_e = Nq + 2;
  const int Np = Nq * Nq * Nq;
  const int Np_e = Nq_e * Nq_e * Nq_e;

  bool overlap = false;
  const bool serial = platform->serial;
  if (Nq >= (elliptic_t::minNFDMOverlap + 1) && !serial)
    overlap = true;

  std::string installDir;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/elliptic/";
  std::string fileName, kernelName;
  const std::string extension = serial ? ".c" : ".okl";

  {
    occa::properties properties = platform->kernelInfo;
    properties["defines/p_Nq"] = Nq;
    properties["defines/p_Nq_e"] = Nq_e;
    properties["defines/p_restrict"] = 0;
    bool useRAS = platform->options.compareArgs(optionsPrefix + "MULTIGRID SMOOTHER", "RAS");
    const std::string suffix =
        std::string("_") + std::to_string(Nq_e - 1) + std::string("pfloat");
    properties["defines/p_overlap"] = (int)overlap;
    if(useRAS){
      properties["defines/p_restrict"] = 1;
    }

    fileName = oklpath + "preFDM" + extension;
    platform->kernels.add(
        "preFDM" + suffix, fileName, properties, suffix);

    int nelgt, nelgv;
    const std::string meshFile = platform->options.getArgs("MESH FILE");
    re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);
    const int NelemBenchmark = nelgv/platform->comm.mpiCommSize;

    bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");
    const int verbosity = verbose ? 2 : 1;
    auto fdmKernel = benchmarkFDM(NelemBenchmark,
                                  Nq_e,
                                  sizeof(pfloat),
                                  useRAS,
                                  static_cast<int>(overlap),
                                  verbosity,
                                  0.2,
                                  false,
                                  suffix);
    platform->kernels.add("fusedFDM" + suffix, fdmKernel);

    fileName = oklpath + "postFDM" + extension;
    platform->kernels.add(
        "postFDM" + suffix, fileName, properties, suffix);
  }
}
void registerFineLevelKernels(const std::string &section, int N, int poissonEquation) {
  auto gen_suffix = [N](const char *floatString) {
    const std::string precision = std::string(floatString);
    if (precision.find(pfloatString) != std::string::npos) {
      return std::string("_") + std::to_string(N) + std::string("pfloat");
    } else {
      return std::string("_") + std::to_string(N);
    }
  };

  auto kernelInfo = platform->kernelInfo + meshKernelProperties(N);
  registerCommonMGPreconditionerKernels(N, kernelInfo, poissonEquation);

  registerAxKernels(section, N, poissonEquation);
  registerSchwarzKernels(section, N);
}
void registerSEMFEMKernels(const std::string &section, int N, int poissonEquation);

void registerMultigridLevelKernels(const std::string &section, int Nf, int N, int poissonEquation) {
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
  registerCommonMGPreconditionerKernels(N, kernelInfo, poissonEquation);

  const bool serial = platform->serial;

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  constexpr int elementType = HEXAHEDRA;

  {
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

    const std::string orderSuffix =
        std::string("_Nf_") + std::to_string(Nf) + std::string("_Nc_") + std::to_string(Nc);
    const std::string fileExtension = serial ? ".c" : ".okl";

    fileName = oklpath + "ellipticPreconCoarsen" + suffix + fileNameExtension;
    kernelName = "ellipticPreconCoarsen" + suffix;
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        coarsenProlongateKernelInfo,
        orderSuffix);
    fileName = oklpath + "ellipticPreconProlongate" + suffix + fileNameExtension;
    kernelName = "ellipticPreconProlongate" + suffix;
    platform->kernels.add(kernelName + orderSuffix,
        fileName,
        coarsenProlongateKernelInfo,
        orderSuffix);
  }

  const std::string optionsPrefix = createOptionsPrefix(section);
  if (N == 1 && platform->options.compareArgs(optionsPrefix + "MULTIGRID COARSE SOLVE", "TRUE")) {
    return;
  }

  registerAxKernels(section, N, poissonEquation);
  registerSchwarzKernels(section, N);
}
void registerMultiGridKernels(const std::string &section, int poissonEquation) {
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const std::string optionsPrefix = createOptionsPrefix(section);

  registerFineLevelKernels(section, N, poissonEquation);

  std::vector<int> levels = determineMGLevels(section);

  if (levels.empty())
    return;

  for (unsigned levelIndex = 1U; levelIndex < levels.size(); ++levelIndex) {
    const int levelFine = levels[levelIndex - 1];
    const int levelCoarse = levels[levelIndex];
    registerMultigridLevelKernels(section, levelFine, levelCoarse, poissonEquation);
  }
  const int coarseLevel = levels.back();
  if (platform->options.compareArgs(
          optionsPrefix + "MULTIGRID COARSE SOLVE", "TRUE")) {
    if (platform->options.compareArgs(
            optionsPrefix + "MULTIGRID COARSE SEMFEM", "TRUE")) {
      registerSEMFEMKernels(section, coarseLevel, poissonEquation);
    } else {
      {
        std::string installDir;
        installDir.assign(getenv("NEKRS_INSTALL_DIR"));
        const std::string oklpath = installDir + "/okl/";
        std::string fileName = oklpath + "parAlmond/convertFP64ToFP32.okl";
        std::string kernelName = "convertFP64ToFP32";
        platform->kernels.add(
            kernelName, fileName, platform->kernelInfo);

        fileName = oklpath + "parAlmond/convertFP32ToFP64.okl";
        kernelName = "convertFP32ToFP64";
        platform->kernels.add(
            kernelName, fileName, platform->kernelInfo);
        fileName = oklpath + "parAlmond/vectorDotStar2.okl";
        kernelName = "vectorDotStar2";
        platform->kernels.add(
            kernelName, fileName, platform->kernelInfo);
      }
    }
  }
}
void registerSEMFEMKernels(const std::string &section, int N, int poissonEquation) {
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
  std::string fileName = oklpath + "gather.okl";
  platform->kernels.add("gather", fileName, SEMFEMKernelProps);
  fileName = oklpath + "scatter.okl";
  platform->kernels.add(
      "scatter", fileName, SEMFEMKernelProps);
  occa::properties stiffnessKernelInfo = platform->kernelInfo;
  fileName = oklpath + "computeStiffnessMatrix.okl";
  stiffnessKernelInfo["defines/p_Nq"] = Nq;
  stiffnessKernelInfo["defines/p_Np"] = Np;
  stiffnessKernelInfo["defines/p_rows_sorted"] = 1;
  stiffnessKernelInfo["defines/p_cols_sorted"] = 0;

  const bool constructOnHost = !platform->device.deviceAtomic;

  if (!constructOnHost) {
    platform->kernels.add("computeStiffnessMatrix",
        fileName,
        stiffnessKernelInfo);
  }
}

}

void registerEllipticPreconditionerKernels(std::string section, int poissonEquation)
{
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const std::string optionsPrefix = createOptionsPrefix(section);

  if(platform->options.compareArgs(optionsPrefix + "PRECONDITIONER", "MULTIGRID")){
    registerMultiGridKernels(section, poissonEquation);
  }
  if(platform->options.compareArgs(optionsPrefix + "PRECONDITIONER", "SEMFEM")){
    registerSEMFEMKernels(section, N, poissonEquation);
  }
  if(platform->options.compareArgs(optionsPrefix + "PRECONDITIONER", "JACOBI")){
    registerJacobiKernels(section, poissonEquation);
  }
}
