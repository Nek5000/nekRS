#include <nrs.hpp>
#include <compileKernels.hpp>
#include <tuple>

void registerLinAlgKernels()
{
  occa::properties kernelInfo = platform->kernelInfo;

  std::string oklDir;
  oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
  oklDir += "/kernels/linAlg/";
  const bool serial = platform->serial;

  const std::string extension = serial ? ".c" : ".okl";
  const std::vector<std::pair<std::string, bool>> allKernels{
      {"fill", false},
      {"pfill", false},
      {"vabs", false},
      {"add", false},
      {"scale", false},
      {"scaleMany", false},
      {"axpby", true},
      {"paxpby", true},
      {"axpbyMany", true},
      {"paxpbyMany", true},
      {"axpbyz", false},
      {"axpbyzMany", false},
      {"axmy", true},
      {"paxmy", false},
      {"axmyMany", true},
      {"axmyVector", true},
      {"axmyz", false},
      {"paxmyz", true},
      {"axmyzMany", false},
      {"paxmyzMany", false},
      {"ady", false},
      {"adyMany", false},
      {"padyMany", false},
      {"axdy", false},
      {"aydx", false},
      {"aydxMany", false},
      {"axdyz", false},
      {"sum", false},
      {"sumMany", false},
      {"min", false},
      {"max", false},
      {"amax", false},
      {"amaxMany", false},
      {"norm2", true},
      {"norm2Many", true},
      {"norm1", true},
      {"norm1Many", true},
      {"weightedNorm1", true},
      {"weightedNorm1Many", true},
      {"weightedNorm2", true},
      {"weightedNorm2Many", true},
      {"innerProd", true},
      {"weightedInnerProd", true},
      {"weightedInnerProdMany", true},
      {"weightedInnerProdMulti", false},
      {"weightedInnerProdMultiDevice", false},
      {"crossProduct", false},
      {"unitVector", false},
  };

  std::string kernelName;
  bool nativeSerialImplementation;
  for (auto &&nameAndSerialImpl : allKernels) {
    std::tie(kernelName, nativeSerialImplementation) = nameAndSerialImpl;
    const std::string extension = (serial && nativeSerialImplementation) ? ".c" : ".okl";
    const bool pfloatKernel = (kernelName.front() == 'p') ? true : false;
    if (pfloatKernel) {
      occa::properties props = kernelInfo;
      props["defines/dfloat"] = pfloatString;
      std::string fileName = kernelName;
      fileName.erase(0, 1);
      platform->kernels.add(kernelName, oklDir + fileName + extension, props);
    }
    else {
      platform->kernels.add(kernelName, oklDir + kernelName + extension, kernelInfo);
    }
  }
}
