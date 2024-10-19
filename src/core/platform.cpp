#include <cstdlib>
#include <strings.h>

#include "platform.hpp"
#include "linAlg.hpp"
#include "flopCounter.hpp"
#include "fileUtils.hpp"

std::string createOptionsPrefix(std::string section)
{
  std::string prefix = section + std::string(" ");
  if (section.find("temperature") != std::string::npos) {
    prefix = std::string("scalar00 ");
  }
  std::transform(prefix.begin(), prefix.end(), prefix.begin(), [](unsigned char c) {
    return std::toupper(c);
  });
  return prefix;
}

namespace
{

static void compileDummyKernel(platform_t &plat)
{
  const bool buildNodeLocal = plat.cacheLocal;
  const std::string dummyKernelName = "myDummyKernelName";
  const std::string dummyKernelStr = std::string("@kernel void myDummyKernelName(int N) {"
                                                 "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {}"
                                                 "}");

  plat.device.occaDevice().buildKernelFromString(dummyKernelStr, dummyKernelName, plat.kernelInfo);
}

} // namespace

platform_t *platform_t::singleton = nullptr;

platform_t::platform_t(setupAide &_options, MPI_Comm _commg, MPI_Comm _comm)
    : options(_options), warpSize(32), comm(_commg, _comm), device(options, comm),
      timer(_comm, device.occaDevice(), 0, 0), kernelRequests(*this)
{
  int rank;
  MPI_Comm_rank(_comm, &rank);

  if (device.mode() == "Serial") {
    options.setArgs("ENABLE GS COMM OVERLAP", "FALSE");
  }

  if (comm.mpiCommSize == 1) {
    options.setArgs("ENABLE GS COMM OVERLAP", "FALSE");
  }

  if (comm.mpiRank == 0 && options.compareArgs("ENABLE GS COMM OVERLAP", "FALSE")) {
    std::cout << "ENABLE GS COMM OVERLAP disabled\n\n";
  }

  buildOnly = false;
  if (options.compareArgs("BUILD ONLY", "TRUE")) {
    buildOnly = true;
  }

  if (options.getArgs("CHECKPOINT ENGINE").empty()) {
    options.setArgs("CHECKPOINT ENGINE", "NEK");
  }

  exitValue = 0;

  // only relevant for SERIAL backend
  setenv("OCCA_MEM_BYTE_ALIGN", std::to_string(ALIGN_SIZE_BYTES).c_str(), 1);

  cacheLocal = 0;
  if (getenv("NEKRS_CACHE_LOCAL")) {
    cacheLocal = std::stoi(getenv("NEKRS_CACHE_LOCAL"));
  }

  cacheBcast = 0;
  if (getenv("NEKRS_CACHE_BCAST")) {
    cacheBcast = std::stoi(getenv("NEKRS_CACHE_BCAST"));
  }

  // build-only mode has to use cacheBcast as well otherwise include paths
  // change triggering a re-build
#if 0
  if (options.compareArgs("BUILD ONLY", "TRUE"))
    cacheBcast = 0;
#endif

  nekrsCheck(cacheLocal && cacheBcast,
             _comm,
             EXIT_FAILURE,
             "%s\n",
             "NEKRS_CACHE_LOCAL=1 and NEKRS_CACHE_BCAST=1 is incompatible!");

  oogs::gpu_mpi(std::stoi(getenv("NEKRS_GPU_MPI")));
  oogs::overlap(options.compareArgs("ENABLE GS COMM OVERLAP", "FALSE") ? 0 : 1);

  if (getenv("OOGS_SYNC_RECV")) {
    oogs::sync_recv(std::stoi(getenv("OOGS_SYNC_RECV")));
  }

  verbose = options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  timer.enableSync();
  if (options.compareArgs("ENABLE TIMER SYNC", "FALSE")) {
    timer.disableSync();
  }

  flopCounter = std::make_unique<flopCounter_t>();

  tmpDir = "";
  if (cacheBcast || cacheLocal) {
    if (getenv("NEKRS_LOCAL_TMP_DIR")) {
      tmpDir = getenv("NEKRS_LOCAL_TMP_DIR");
    } else {
      nekrsAbort(comm.mpiComm, EXIT_FAILURE, "%s\n", "NEKRS_LOCAL_TMP_DIR undefined!");
    }

    int rankLocal;
    MPI_Comm_rank(comm.mpiCommLocal, &rankLocal);

    if (rankLocal == 0) {
      nekrsCheck(!fs::exists(tmpDir),
                 MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "Cannot find NEKRS_LOCAL_TMP_DIR: %s\n",
                 tmpDir.c_str());
    }

    const auto multiSession = [&]() {
      int retVal;
      MPI_Comm_compare(comm.mpiComm, comm.mpiCommParent, &retVal);
      return (retVal == MPI_IDENT) ? false : true;
    }();

    if (multiSession) {
      int sessionID;
      options.getArgs("NEKNEK SESSION ID", sessionID);
      tmpDir = fs::path(tmpDir) / fs::path("sess" + std::to_string(sessionID));
      fs::create_directories(tmpDir);
    }
  }

  {
    int rankLocal = rank;
    if (cacheLocal) {
      MPI_Comm_rank(comm.mpiCommLocal, &rankLocal);
    }

    if (rankLocal == 0) {
      std::string cache_dir;
      cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
      mkdir(cache_dir.c_str(), S_IRWXU);
    }
  }

  kernelInfo["includes"].asArray();

  // Disables the automatic insertion of barriers
  // between separate OKL inner loop blocks.
  kernelInfo["okl/add_barriers"] = false;

  kernelInfo["defines/TOKEN_PASTE_(a,b)"] = std::string("a##b");
  kernelInfo["defines/TOKEN_PASTE(a,b)"] = std::string("TOKEN_PASTE_(a,b)");
  kernelInfo["defines/FUNC(a)"] = std::string("TOKEN_PASTE(a,SUFFIX)");

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
             "hlong"] = "long long"; // occa's parser doesn't know long long int

  kernelInfo["defines/"
             "p_PI"] = M_PI;

  if (device.mode() == "CUDA") {
    kernelInfo["defines/smXX"] = 1;
  }

  if (device.mode() == "OpenCL") {
    kernelInfo["defines/"
               "hlong"] = "long";
  }

  if (device.mode() == "HIP") {
    warpSize = 64; // can be arch specific
    kernelInfo["defines/gfxXX"] = 1;
  }

  if (device.mode() == "dpcpp") {
    kernelInfo["defines/XeHPC"] = 1;
  }

  serial = device.mode() == "Serial" || device.mode() == "OpenMP";

  if (serial) {
    kernelInfo["includes"] += "math.h";
  }

  const std::string floatingPointType = static_cast<std::string>(kernelInfo["defines/dfloat"]);
  if (floatingPointType.find("float") != std::string::npos) {
    kernelInfo["defines/FP32"] = 1;
  }

  const std::string extension = serial ? ".c" : ".okl";

  if (rank == 0) {
    compileDummyKernel(*this);
  }

  {
    occa::json properties;
    properties["resize_through_host"] = 1;
    deviceMemoryPool = device.occaDevice().createMemoryPool(properties);
    deviceMemoryPool.setAlignment(ALIGN_SIZE_BYTES);
  }

  {
    occa::json properties;
    properties["resize_through_host"] = 1;
    properties["host"] = true;
    memoryPool = device.occaDevice().createMemoryPool(properties);
    memoryPool.setAlignment(ALIGN_SIZE_BYTES);
  }
}

// input files required for JIT kernel compilation or load
void platform_t::bcastJITKernelSourceFiles()
{
  if (platform->verbose && comm.mpiRank == 0) {
    std::cout << "broadcast kernel sources to " << platform->tmpDir << std::endl;
  }

  const auto NEKRS_HOME_NEW = fs::path(tmpDir) / "nekrs";
  const auto srcPath = fs::path(getenv("NEKRS_HOME"));
  for (auto &entry : {fs::path("include"), fs::path("kernels")}) {

    fileBcast(srcPath / entry, NEKRS_HOME_NEW, comm.mpiComm, verbose);
  }

  setenv("NEKRS_HOME", std::string(NEKRS_HOME_NEW).c_str(), 1);
  setenv("NEKRS_KERNEL_DIR", std::string(NEKRS_HOME_NEW / "kernels").c_str(), 1);
}
