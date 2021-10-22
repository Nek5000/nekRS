#include <compileKernels.hpp>
#include "elliptic.h"
#include "mesh.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "udf.hpp"
#include <vector>

namespace {

void compileDummyKernel()
{
  const bool buildNodeLocal = useNodeLocalCache();
  auto rank = buildNodeLocal ? platform->comm.localRank : platform->comm.mpiRank;
  const std::string dummyKernelName = "myDummyKernelName";
  const std::string dummyKernelStr = std::string(
      "@kernel void myDummyKernelName(int N) {"
      "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {}"
      "}"
  );

  if(rank == 0){
    platform->device.occaDevice().buildKernelFromString(
      dummyKernelStr,
      dummyKernelName,
      platform->kernelInfo
    );
  }

}

} // namespace

std::string createOptionsPrefix(std::string section) {
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

void compileKernels() {

  compileDummyKernel(); // trigger occa's compilerVendorTest

  const occa::properties kernelInfoBC = compileUDFKernels();

  registerLinAlgKernels();

  registerMeshKernels(kernelInfoBC);

  registerNrsKernels(kernelInfoBC);

  int Nscalars;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalars);
  if (Nscalars) {
    registerCdsKernels(kernelInfoBC);
  }

  // Scalar section is omitted
  // as pressure section kernels are the same.
  const std::vector<std::string> sections = {
      "pressure",
      "velocity",
  };
  for (auto &&section : sections) {
    registerEllipticKernels(section);
    registerEllipticPreconditionerKernels(section);
  }

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0)
    printf("loading kernels ... ");
  fflush(stdout);

    {
      const bool buildNodeLocal = useNodeLocalCache();
      const bool buildOnly = platform->options.compareArgs("BUILD ONLY", "TRUE");
      auto communicator = buildNodeLocal ? platform->comm.mpiCommLocal : platform->comm.mpiComm;
      oogs::compile(
          platform->device.occaDevice(), platform->device.mode(), communicator, buildOnly);
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
