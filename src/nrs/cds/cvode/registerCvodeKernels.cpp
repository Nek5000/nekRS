#include "platform.hpp"
#include "mesh.h"
#include "compileKernels.hpp"

void registerCvodeKernels()
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

  const std::string suffix = "Hex3D";
  const std::string prefix = "cvode_t::";
  const std::string extension = platform->serial ? ".c" : ".okl";

  const std::string oklpath = getenv("NEKRS_KERNEL_DIR");

  kernelName = "extrapolate";
  fileName = oklpath + "/core/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "errorWeight";
  fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "mapToMaskedPoint";
  fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "extrapolateDirichlet";
  fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "axpby";
  fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "axmyz";
  fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "linearCombination";
  fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  kernelName = "innerProdMulti";
  fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
  platform->kernelRequests.add(prefix + kernelName, fileName, platform->kernelInfo);

  {
    auto prop = platform->kernelInfo;
    kernelName = "fusedAddRhoDiv";
    fileName = oklpath + "/nrs/cds/cvode/" + kernelName + ".okl";
 
    prop["defines/p_addPointSource"] = 0; 
    platform->kernelRequests.add(prefix + "rhoDiv", fileName, prop);
 
    prop["defines/p_addPointSource"] = 1; 
    platform->kernelRequests.add(prefix + kernelName, fileName, prop);
  }

  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  int cubN;
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);

  auto prop = platform->kernelInfo;
  prop["includes"].asArray();
  prop += meshKernelProperties(N);

  std::string derivDataFile = std::string(getenv("NEKRS_KERNEL_DIR")) + "/mesh/constantGLLDifferentiationMatrices.h";

  prop["includes"] += derivDataFile.c_str();
  prop["defines/p_weightInputAdd"] = 1;
  kernelName = "weakLaplacianHex3D";
  fileName = oklpath + "/core/" + kernelName + ".okl";

  platform->kernelRequests.add(prefix + kernelName, fileName, prop);
}
