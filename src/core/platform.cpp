#include <cstdlib>
#include <strings.h>

#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "flopCounter.hpp"
#include "fileUtils.hpp"

namespace {

static void compileDummyKernel(platform_t &plat)
{
  const bool buildNodeLocal = plat.cacheLocal;
  auto rank = buildNodeLocal ? plat.comm.localRank : plat.comm.mpiRank;
  const std::string dummyKernelName = "myDummyKernelName";
  const std::string dummyKernelStr = std::string("@kernel void myDummyKernelName(int N) {"
                                                 "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {}"
                                                 "}");

  if (rank == 0) {
    plat.device.occaDevice().buildKernelFromString(dummyKernelStr, dummyKernelName, plat.kernelInfo);
  }
}

} // namespace

platform_t *platform_t::singleton = nullptr;

platform_t::platform_t(setupAide &_options, MPI_Comm _commg, MPI_Comm _comm)
    : options(_options), warpSize(32), comm(_commg, _comm), device(options, comm),
      timer(_comm, device.occaDevice(), 0, 0), kernels(*this)
{
  int rank;
  MPI_Comm_rank(_comm, &rank);

  exitValue = 0;

  // only relevant for SERIAL backend
  setenv("OCCA_MEM_BYTE_ALIGN", std::to_string(ALIGN_SIZE).c_str(), 1);

  cacheLocal = 0;
  if(getenv("NEKRS_CACHE_LOCAL"))
    cacheLocal = std::stoi(getenv("NEKRS_CACHE_LOCAL"));

  cacheBcast = 0;
  if(getenv("NEKRS_CACHE_BCAST"))
    cacheBcast = std::stoi(getenv("NEKRS_CACHE_BCAST"));

  // build-only mode has to use cacheBcast as well otherwise include paths 
  // change triggering a re-build 
#if 0
  if (options.compareArgs("BUILD ONLY", "TRUE"))
    cacheBcast = 0;
#endif

  nrsCheck(cacheLocal && cacheBcast,
           _comm, EXIT_FAILURE, 
           "%s\n", "NEKRS_CACHE_LOCAL=1 and NEKRS_CACHE_BCAST=1 is incompatible!");

  oogs::gpu_mpi(std::stoi(getenv("NEKRS_GPU_MPI")));

  if(getenv("OOGS_SYNC_RECV"))
    oogs::sync_recv(std::stoi(getenv("OOGS_SYNC_RECV")));

  verbose = options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  timer.enableSync();
  if (options.compareArgs("ENABLE TIMER SYNC", "FALSE"))
    timer.disableSync();

  flopCounter = std::make_unique<flopCounter_t>();

  tmpDir = "/";

  // bcast install dir 
  if(cacheBcast || cacheLocal) {
    if(getenv("NEKRS_LOCAL_TMP_DIR")) { 
      tmpDir = getenv("NEKRS_LOCAL_TMP_DIR");
    } else {
      nrsAbort(_comm, EXIT_FAILURE, "%s\n", "NEKRS_LOCAL_TMP_DIR undefined!");
    }

    int rankLocal;
    MPI_Comm_rank(comm.mpiCommLocal, &rankLocal);

    if(rankLocal == 0)
      nrsCheck(!fs::exists(tmpDir), MPI_COMM_SELF, EXIT_FAILURE,
               "Cannot find NEKRS_LOCAL_TMP_DIR %s\n", tmpDir.c_str());
  
    auto nSessions = 1; 
    options.getArgs("NEKNEK NUMBER OF SESSIONS", nSessions);
    if(nSessions > 1) {
      auto sessionID = 0; 
      options.getArgs("NEKNEK SESSION ID", sessionID);

      tmpDir = fs::path(tmpDir) / fs::path(std::string("nrs_") + std::to_string(sessionID)); 

      if(rankLocal == 0) {
        fs::create_directory(tmpDir);
        nrsCheck(!fs::exists(tmpDir), MPI_COMM_SELF, EXIT_FAILURE,
                 "Cannot create %s\n", tmpDir.c_str());
      }
    }

    const auto NEKRS_HOME_NEW = fs::path(tmpDir) / "nekrs";
    const auto srcPath = fs::path(getenv("NEKRS_HOME"));
    for (auto &entry : {fs::path("udf"),
                        fs::path("nek5000"),
                        fs::path("nekInterface"),
                        fs::path("include"),
                        fs::path("gatherScatter"),
                        fs::path("3rd_party"),
                        fs::path("kernels")}) {
      
      fileBcast(srcPath / entry, NEKRS_HOME_NEW, comm.mpiComm, verbose);
    }
    setenv("NEKRS_HOME", std::string(NEKRS_HOME_NEW).c_str(), 1);
    setenv("NEKRS_KERNEL_DIR", std::string(NEKRS_HOME_NEW / "kernels").c_str(), 1);
    setenv("OGS_HOME", std::string(NEKRS_HOME_NEW / "gatherScatter").c_str(), 1);
  }

  {
    int rankLocal = rank;
    if(cacheLocal)
      MPI_Comm_rank(comm.mpiCommLocal, &rankLocal);
 
    if(rankLocal == 0) {
      std::string cache_dir;
      cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
      mkdir(cache_dir.c_str(), S_IRWXU);
    }
  }

  // Disables the automatic insertion of barriers
  // between separate OKL inner loop blocks.
  kernelInfo["okl/add_barriers"] = false;

  kernelInfo["defines/"
             "p_NVec"] = 3;
  kernelInfo["defines/"
             "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/"
             "dfloat"] = dfloatString;
  kernelInfo["defines/"
             "pfloat"] = pfloatString;
  kernelInfo["defines/"
             "dlong"] = dlongString;
  kernelInfo["defines/"
             "hlong"] = hlongString;

  if (device.mode() == "CUDA") {
  }

  if (device.mode() == "OpenCL") {
    kernelInfo["defines/"
               "hlong"] = "long";
  }

  if (device.mode() == "HIP") {
    warpSize = 64; // can be arch specific
  }

  serial = device.mode() == "Serial" || device.mode() == "OpenMP";

  const std::string extension = serial ? ".c" : ".okl";

  compileDummyKernel(*this);

  std::string kernelName, fileName;
  const auto oklpath = std::string(getenv("NEKRS_KERNEL_DIR"));
  kernelName = "copyDfloatToPfloat";
  fileName = oklpath + "/core/" + kernelName + extension;
  this->kernels.add(kernelName, fileName, this->kernelInfo);

  kernelName = "copyPfloatToDfloat";
  fileName = oklpath + "/core/" + kernelName + extension;
  this->kernels.add(kernelName, fileName, this->kernelInfo);

  occa::json properties;
  o_memPool = device.occaDevice().createMemoryPool(properties);
  o_memPool.setAlignment(ALIGN_SIZE);
}
