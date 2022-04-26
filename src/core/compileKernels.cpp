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

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0)
    printf("loading kernels (this may take awhile) ...\n");
  fflush(stdout);

  const occa::properties kernelInfoBC = compileUDFKernels();

  registerLinAlgKernels();

  registerMeshKernels(kernelInfoBC);

  registerNrsKernels(kernelInfoBC);

  int Nscalars;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalars);
  const int scalarWidth = getDigitsRepresentation(NSCALAR_MAX - 1);

  if (Nscalars) {
    registerCdsKernels(kernelInfoBC);
    for(int is = 0; is < Nscalars; is++){
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(scalarWidth) << is;
      std::string sid = ss.str();
      const std::string section = "scalar" + sid;
      const int poisson = 0;

      if(!platform->options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")){
        registerEllipticKernels(section, poisson);
        registerEllipticPreconditionerKernels(section, poisson);
      }
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


  {
    const bool buildNodeLocal = useNodeLocalCache();
    const bool buildOnly = platform->options.compareArgs("BUILD ONLY", "TRUE");
    auto communicator = buildNodeLocal ? platform->comm.mpiCommLocal : platform->comm.mpiComm;
    oogs::compile(
        platform->device.occaDevice(), platform->device.mode(), communicator, buildOnly);
  }

  platform->kernels.compile();

  // load platform related kernels
  std::string kernelName;
  kernelName = "copyDfloatToPfloat";
  platform->copyDfloatToPfloatKernel = platform->kernels.get(kernelName);

  kernelName = "copyPfloatToDfloat";
  platform->copyPfloatToDfloatKernel = platform->kernels.get(kernelName);

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
