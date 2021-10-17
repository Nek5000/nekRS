#include <nrs.hpp>
#include <compileKernels.hpp>

void registerLinAlgKernels()
{
  occa::properties kernelInfo = platform->kernelInfo;

  std::string oklDir;
  oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
  oklDir += "/okl/linAlg/";
  std::string fileName;
  const bool serial = useSerial();

  platform->kernels.add(
      "fill", oklDir + "linAlgFill.okl", "fill", kernelInfo);
  platform->kernels.add(
      "vabs", oklDir + "linAlgAbs.okl", "vabs", kernelInfo);
  platform->kernels.add(
      "add", oklDir + "linAlgAdd.okl", "add", kernelInfo);
  platform->kernels.add(
      "scale", oklDir + "linAlgScale.okl", "scale", kernelInfo);
  platform->kernels.add(
      "scaleMany", oklDir + "linAlgScale.okl", "scaleMany", kernelInfo);
  fileName = std::string("linAlgAXPBY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add("axpby", oklDir + fileName, "axpby", kernelInfo);
  fileName = std::string("linAlgAXPBY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "axpbyMany", oklDir + fileName, "axpbyMany", kernelInfo);
  platform->kernels.add(
      "axpbyz", oklDir + "linAlgAXPBY.okl", "axpbyz", kernelInfo);
  platform->kernels.add(
      "axpbyzMany", oklDir + "linAlgAXPBY.okl", "axpbyzMany", kernelInfo);
  fileName = std::string("linAlgAXMY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add("axmy", oklDir + fileName, "axmy", kernelInfo);
  fileName = std::string("linAlgAXMY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "axmyMany", oklDir + fileName, "axmyMany", kernelInfo);
  fileName = std::string("linAlgAXMY") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "axmyVector", oklDir + fileName, "axmyVector", kernelInfo);
  platform->kernels.add(
      "axmyz", oklDir + "linAlgAXMY.okl", "axmyz", kernelInfo);
  platform->kernels.add(
      "axmyzMany", oklDir + "linAlgAXMY.okl", "axmyzMany", kernelInfo);
  platform->kernels.add(
      "ady", oklDir + "linAlgAXDY.okl", "ady", kernelInfo);
  platform->kernels.add(
      "adyMany", oklDir + "linAlgAXDY.okl", "adyMany", kernelInfo);
  platform->kernels.add(
      "axdy", oklDir + "linAlgAXDY.okl", "axdy", kernelInfo);
  platform->kernels.add(
      "aydx", oklDir + "linAlgAXDY.okl", "aydx", kernelInfo);
  platform->kernels.add(
      "aydxMany", oklDir + "linAlgAXDY.okl", "aydxMany", kernelInfo);
  platform->kernels.add(
      "axdyz", oklDir + "linAlgAXDY.okl", "axdyz", kernelInfo);
  platform->kernels.add(
      "sum", oklDir + "linAlgSum.okl", "sum", kernelInfo);
  platform->kernels.add(
      "sumMany", oklDir + "linAlgSum.okl", "sumMany", kernelInfo);
  platform->kernels.add(
      "min", oklDir + "linAlgMin.okl", "min", kernelInfo);
  platform->kernels.add(
      "max", oklDir + "linAlgMax.okl", "max", kernelInfo);
  fileName = std::string("linAlgNorm2") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "norm2", oklDir + fileName, "norm2", kernelInfo);
  platform->kernels.add(
      "norm2Many", oklDir + fileName, "norm2Many", kernelInfo);
  fileName = std::string("linAlgNorm1") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "norm1", oklDir + fileName, "norm1", kernelInfo);
  platform->kernels.add(
      "norm1Many", oklDir + fileName, "norm1Many", kernelInfo);
  fileName = std::string("linAlgWeightedNorm1") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "weightedNorm1", oklDir + fileName, "weightedNorm1", kernelInfo);
  platform->kernels.add(
      "weightedNorm1Many", oklDir + fileName, "weightedNorm1Many", kernelInfo);
  fileName = std::string("linAlgWeightedNorm2") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "weightedNorm2", oklDir + fileName, "weightedNorm2", kernelInfo);
  fileName = std::string("linAlgWeightedNorm2") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "weightedNorm2Many", oklDir + fileName, "weightedNorm2Many", kernelInfo);
  platform->kernels.add(
      "innerProd", oklDir + "linAlgInnerProd.okl", "innerProd", kernelInfo);
  fileName = std::string("linAlgWeightedInnerProd") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add(
      "weightedInnerProd", oklDir + fileName, "weightedInnerProd", kernelInfo);
  fileName = std::string("linAlgWeightedInnerProd") +
             (serial ? std::string(".c") : std::string(".okl"));
  platform->kernels.add("weightedInnerProdMany",
      oklDir + fileName,
      "weightedInnerProdMany",
      kernelInfo);
  platform->kernels.add("weightedInnerProdMulti",
      oklDir + "linAlgWeightedInnerProd.okl",
      "weightedInnerProdMulti",
      kernelInfo);
}