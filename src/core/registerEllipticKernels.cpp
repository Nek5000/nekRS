#include <compileKernels.hpp>
#include "re2Reader.hpp"
#include "benchmarkAx.hpp"

namespace {

void registerGMRESKernels(const std::string &section, int Nfields)
{
  const std::string oklpath = getenv("NEKRS_KERNEL_DIR") + std::string("/elliptic/");
  std::string fileName;
  const bool serial = platform->serial;

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";
  const std::string sectionIdentifier = std::to_string(Nfields) + "-";

  occa::properties gmresKernelInfo = platform->kernelInfo;
  gmresKernelInfo["defines/p_Nfields"] = Nfields;

  std::string kernelName = "gramSchmidtOrthogonalization";
  fileName = oklpath + kernelName + fileNameExtension;
  platform->kernels.add(sectionIdentifier + kernelName, fileName, gmresKernelInfo);

  kernelName = "updatePGMRESSolution";
  fileName = oklpath + kernelName + fileNameExtension;
  platform->kernels.add(sectionIdentifier + kernelName, fileName, gmresKernelInfo);

  kernelName = "fusedResidualAndNorm";
  fileName = oklpath + kernelName + fileNameExtension;
  platform->kernels.add(sectionIdentifier + kernelName, fileName, gmresKernelInfo);
}

} // namespace

void registerEllipticKernels(std::string section, int poissonEquation)
{
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const std::string optionsPrefix = createOptionsPrefix(section);

  occa::properties kernelInfo = platform->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();
  kernelInfo += meshKernelProperties(N);

  const bool blockSolver = [&section]() {
    if (section == "velocity" && 
        platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"))
      return true;
    if (section == "velocity" && 
        platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE"))
      return true;
    if (section == "mesh" && 
        platform->options.compareArgs("MESH BLOCK SOLVER", "TRUE"))
      return true;

    return false;
  }();

  const int Nfields = (blockSolver) ? 3 : 1;

  const bool stressForm = [&section]() {
    if (section == "velocity" && platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE"))
      return true;
    if (section == "mesh" && platform->options.compareArgs("MESH STRESSFORMULATION", "TRUE"))
      return true;
    return false;
  }();

  const bool serial = platform->serial;

  const std::string fileNameExtension = (serial) ? ".c" : ".okl";

  const std::string sectionIdentifier = std::to_string(Nfields) + "-";

  if (platform->options.compareArgs(optionsPrefix + "SOLVER", "PGMRES")) {
    registerGMRESKernels(section, Nfields);
  }

  {
    const std::string oklpath = getenv("NEKRS_KERNEL_DIR") + std::string("/elliptic/");
    std::string fileName, kernelName;

    {
      const std::string extension = ".okl";
      occa::properties properties = platform->kernelInfo;
      properties["defines/p_Nfields"] = Nfields;

      kernelName = "multiScaledAddwOffset";
      fileName = oklpath + kernelName + extension;
      platform->kernels.add(sectionIdentifier + kernelName, fileName, properties);
      kernelName = "accumulate";
      fileName = oklpath + kernelName + extension;
      platform->kernels.add(sectionIdentifier + kernelName, fileName, properties);

      kernelName = "fusedCopyDfloatToPfloat";
      fileName = oklpath + kernelName + extension;
      platform->kernels.add(kernelName, fileName, properties);
    }
  }

  {
    const std::string oklpath = getenv("NEKRS_KERNEL_DIR") + std::string("/core/");
    std::string fileName;

    fileName = oklpath + "mask.okl";
    platform->kernels.add("mask", fileName, kernelInfo);

    occa::properties pfloatKernelInfo = kernelInfo;
    pfloatKernelInfo["defines/dfloat"] = pfloatString;
    platform->kernels.add("maskPfloat", fileName, pfloatKernelInfo);
  }

  kernelInfo["defines/p_Nfields"] = Nfields;

  occa::properties dfloatKernelInfo = kernelInfo;
  occa::properties floatKernelInfo = kernelInfo;
  floatKernelInfo["defines/pfloat"] = pfloatString;
  floatKernelInfo["defines/dfloat"] = pfloatString;

  const std::string suffix = "Hex3D";

  occa::properties AxKernelInfo = dfloatKernelInfo;
  const std::string oklpath = getenv("NEKRS_KERNEL_DIR") + std::string("/elliptic/");
  std::string fileName;
  std::string kernelName;

  kernelName = "ellipticBlockBuildDiagonal" + suffix;
  fileName = oklpath + kernelName + ".okl";
  dfloatKernelInfo["defines/dfloat"] = dfloatString;
  dfloatKernelInfo["defines/pfloat"] = pfloatString;
  platform->kernels.add(sectionIdentifier + kernelName, fileName, dfloatKernelInfo);

  if (poissonEquation) {
    AxKernelInfo["defines/p_poisson"] = 1;
  }

  int nelgt, nelgv;
  const std::string meshFile = platform->options.getArgs("MESH FILE");
  re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);
  const int NelemBenchmark = nelgv / platform->comm.mpiCommSize;
  bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const int verbosity = verbose ? 2 : 1;

  if (section == "pressure" && platform->options.compareArgs("LOWMACH", "TRUE")) {
    platform->options.setArgs(optionsPrefix + "ELLIPTIC COEFF FIELD", "TRUE");
    platform->options.setArgs(optionsPrefix + "ELLIPTIC PRECO COEFF FIELD", "TRUE");
  }

  for (auto &&coeffField : {true, false}) {
    if (platform->options.compareArgs(optionsPrefix + "ELLIPTIC COEFF FIELD", "TRUE") != coeffField)
      continue;

    std::string kernelNamePrefix = "elliptic";
    if (blockSolver)
      kernelNamePrefix += (stressForm) ? "Stress" : "Block";

    kernelName = "Ax";
    kernelName += "Coeff";
    if (platform->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
      kernelName += "Trilinear";
    kernelName += suffix;

    const std::string _kernelName = kernelNamePrefix + "Partial" + kernelName;
    const std::string prefix = (poissonEquation) ? "poisson-" : "";
    fileName = oklpath + _kernelName + fileNameExtension;

    auto axKernel = benchmarkAx(NelemBenchmark,
                                N + 1,
                                N,
                                !coeffField,
                                poissonEquation,
                                false,
                                sizeof(dfloat),
                                Nfields,
                                stressForm,
                                verbosity,
                                elliptic_t::targetTimeBenchmark,
                                platform->options.compareArgs("KERNEL AUTOTUNING", "FALSE") ? false : true,
                                "");
    platform->kernels.add(prefix + _kernelName, axKernel);
  }

  kernelName = "ellipticBlockBuildDiagonal" + suffix;
  fileName = oklpath + kernelName + ".okl";
  dfloatKernelInfo["defines/dfloat"] = dfloatString;
  dfloatKernelInfo["defines/pfloat"] = pfloatString;
  platform->kernels.add(sectionIdentifier + kernelName, fileName, dfloatKernelInfo);
  dfloatKernelInfo["defines/pfloat"] = dfloatString;

  // PCG update
  fileName = oklpath + "ellipticBlockUpdatePCG" + fileNameExtension;
  platform->kernels.add(sectionIdentifier + "ellipticBlockUpdatePCG", fileName, kernelInfo);
}
