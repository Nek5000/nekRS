#include <compileKernels.hpp>
#include "elliptic.h"
#include "mesh.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "udf.hpp"
#include <vector>
#include <tuple>

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

  const occa::properties kernelInfoBC = compileUDFKernels();

  registerLinAlgKernels();

  registerMeshKernels(kernelInfoBC);

  registerNrsKernels(kernelInfoBC);

  int Nscalars;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalars);
  if (Nscalars) {
    registerCdsKernels(kernelInfoBC);
    for(int is = 0; is < Nscalars; is++){
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is;
      std::string sid = ss.str();
      const std::string section = "scalar" + sid;
      const int poisson = 0;
      registerEllipticKernels(section, poisson);
      registerEllipticPreconditionerKernels(section, poisson);
    }
  }

  // Scalar section is omitted
  // as pressure section kernels are the same.
  const std::vector<std::pair<std::string,int>> sections = {
      {"pressure", 1},
      {"velocity", 0}
  };

  std::string section;
  int poissonEquation;
  for (auto&& entry : sections) {
    std::tie(section, poissonEquation) = entry;
    registerEllipticKernels(section, poissonEquation);
    registerEllipticPreconditionerKernels(section, poissonEquation);
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
