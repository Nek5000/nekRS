#include "compileKernels.hpp"
#include "nrs.hpp"

void registerCvodeKernels(occa::properties kernelInfoBC)
{
  int Nscalars = 0;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalars);

  int NscalarsCvode = 0;

  if (Nscalars) {
    for (int is = 0; is < Nscalars; is++) {
      const auto sid = scalarDigitStr(is);

      if (platform->options.compareArgs("SCALAR" + sid + " SOLVER", "CVODE")) {
        NscalarsCvode++;
      }
    }
  }

  if (!NscalarsCvode) {
    return;
  }

  std::string kernelName;
  std::string fileName;
  
  const std::string oklpath = getenv("NEKRS_KERNEL_DIR") + std::string("/cvode/");
  const std::string prefix = "cvode_t::";

  kernelName = "errorWeight";
  fileName = oklpath + kernelName + ".okl";
  platform->kernels.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "mapToMaskedPoint";
  fileName = oklpath + kernelName + ".okl";
  platform->kernels.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "extrapolateDirichlet";
  fileName = oklpath + kernelName + ".okl";
  platform->kernels.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "cvToNrs";
  fileName = oklpath + kernelName + ".okl";
  platform->kernels.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "nrsToCv";
  fileName = oklpath + kernelName + ".okl";
  platform->kernels.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "nrsToCv";
  fileName = oklpath + kernelName + ".okl";
  platform->kernels.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "cvToNrs";
  fileName = oklpath + kernelName + ".okl";
  platform->kernels.add(prefix + kernelName, fileName, platform->kernelInfo);

  {
    auto prop = platform->kernelInfo;
    kernelName = "fusedAddRhoDiv";
    fileName = oklpath + kernelName + ".okl";
 
    prop["defines/p_addPointSource"] = 0; 
    platform->kernels.add(prefix + "rhoDiv", fileName, prop);
 
    prop["defines/p_addPointSource"] = 1; 
    platform->kernels.add(prefix + kernelName, fileName, prop);
  }

  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);

  auto weakLaplacianKernelInfo = platform->kernelInfo;
  weakLaplacianKernelInfo["includes"].asArray();
  weakLaplacianKernelInfo += meshKernelProperties(N);
  std::string derivDataFile = std::string(getenv("NEKRS_KERNEL_DIR")) + "/mesh/constantGLLDifferentiationMatrices.h";

  weakLaplacianKernelInfo["includes"] += derivDataFile.c_str();
  kernelName = "weakLaplacianHex3D";
  fileName = oklpath + kernelName + ".okl";

  platform->kernels.add(prefix + kernelName, fileName, weakLaplacianKernelInfo);
}
