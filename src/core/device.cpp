#include <limits.h>
#include <unistd.h>
#include <regex>
#include "nekrsSys.hpp"
#include "device.hpp"
#include "platform.hpp"
#include "fileUtils.hpp"

namespace
{

void setOccaVars()
{
  std::string cache_dir;
  if (getenv("NEKRS_CACHE_DIR")) {
    cache_dir.assign(getenv("NEKRS_CACHE_DIR"));
  }

  if (!getenv("OCCA_CACHE_DIR")) {
    const std::string path = cache_dir + "/occa/";
    occa::env::OCCA_CACHE_DIR = path;
    setenv("OCCA_CACHE_DIR", path.c_str(), 1);
  }

  if (!getenv("OCCA_DIR")) {
    occa::env::OCCA_DIR = std::string(getenv("NEKRS_HOME")) + "/";
    setenv("OCCA_DIR", occa::env::OCCA_DIR.c_str(), 1);
  }

  occa::env::OCCA_INSTALL_DIR = occa::env::OCCA_DIR;
}

int compileDummyAtomicKernel(device_t &device)
{
  const std::string dummyKernelName = "simpleAtomicAdd";
  const std::string dummyKernelStr = std::string("@kernel void simpleAtomicAdd(int N, double * result) {"
                                                 "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {"
                                                 "    @atomic result[0] += 1;"
                                                 "  }"
                                                 "}");

  occa::properties noKernelInfo;

  auto simpleAtomicAddKernel =
      device.occaDevice().buildKernelFromString(dummyKernelStr, dummyKernelName, noKernelInfo);

  auto o_result = device.occaDevice().malloc(sizeof(double));
  double initialValue = 0.0;
  o_result.copyFrom(&initialValue, sizeof(double));

  constexpr int N = 1000;
  auto expectedValue = static_cast<double>(N);
  double actualValue = 0.0;

  auto eps = 10 * std::numeric_limits<double>::epsilon();

  simpleAtomicAddKernel(N, o_result);

  o_result.copyTo(&actualValue, sizeof(double));

  return std::abs(actualValue - expectedValue) < eps;
}

bool atomicsAvailable(device_t &device, MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  int atomicSupported = 1;

  if (rank == 0) {
    try {
      atomicSupported = compileDummyAtomicKernel(device);
    } catch (std::exception &e) {
      atomicSupported = 0;
    }
  }

  MPI_Bcast(&atomicSupported, 1, MPI_INT, 0, comm);

  return atomicSupported;
}

std::string extractKernelName(const std::string &fullPath)
{
  std::regex kernelNameRegex(R"((.+)\/(.+)\.)");
  std::smatch kernelNameMatch;
  const auto foundKernelName = std::regex_search(fullPath, kernelNameMatch, kernelNameRegex);

  // capture group
  // 0:   /path/to/install/nekrs/kernels/cds/advectMeshVelocityHex3D.okl
  // 1:   /path/to/install/nekrs/kernels/cds
  // 2:   advectMeshVelocityHex3D.okl

  return (foundKernelName && kernelNameMatch.size() == 3) ? kernelNameMatch[2].str() : "";
}

occa::properties adjustKernelProps(const std::string &fileName, const occa::properties &props_)
{
  occa::properties props = props_;
  if (fileName.find(".okl") != std::string::npos) {
    props["okl/enabled"] = true;
    props["defines/__okl__"] = 1;
  } else {
    props["okl/enabled"] = false;
  }
  props["defines/SUFFIX"] = ""; // not used anymore
  return props;
}

} // namespace

occa::kernel device_t::wrapperCompileKernel(const std::string &fileName,
                                            const occa::properties &props_,
                                            const std::string &suffix) const
{
  if (!_compilationEnabled) {
    nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s", "illegal call detected after 'finish' declaration\n");
  }

  if (fileName.empty()) {
    nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s", "Empty fileName\n");
  }

  auto props = props_;
  props["build/compile_only"] = true;
  props["build/load_only"] = false;

  return _device.buildKernel(fileName, "" /* dummy */, adjustKernelProps(fileName, props));
}

occa::kernel device_t::wrapperLoadKernel(const std::string &fileName,
                                         const std::string &kernelName,
                                         const occa::properties &props_,
                                         const std::string &suffix) const
{
  nekrsCheck(platform->options.compareArgs("BUILD ONLY", "TRUE"),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s",
             "illegal call during BUILD ONLY mode\n");

#if 0
  const std::string cacheDir0 = occa::env::OCCA_CACHE_DIR;

  // redirect
  if (platform->cacheBcast) {
    const std::string cacheDir = platform->tmpDir / fs::path("occa/");
    occa::env::OCCA_CACHE_DIR = cacheDir;
    setenv("OCCA_CACHE_DIR", cacheDir.c_str(), 1);
  }
#endif

  auto props = props_;
  props["build/load_only"] = true;
  props["build/compile_only"] = false;
  auto knl = _device.buildKernel(fileName, kernelName, adjustKernelProps(fileName, props));

  nekrsCheck(!knl.isInitialized(),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "Cannot load kernel <%s>\n",
             kernelName.c_str());

#if 0
  // restore
  if (platform->cacheBcast) {
    occa::env::OCCA_CACHE_DIR = cacheDir0;
    setenv("OCCA_CACHE_DIR", cacheDir0.c_str(), 1);
  }
#endif

  return knl;
}

occa::kernel device_t::compileKernel(const std::string &fileName,
                                     const occa::properties &props,
                                     const std::string &suffix,
                                     const MPI_Comm &commIn) const
{
  const auto collective = [&commIn]() {
    int tmp;
    MPI_Comm_compare(commIn, MPI_COMM_SELF, &tmp);
    return (tmp == MPI_UNEQUAL) ? true : false;
  }();

  MPI_Comm comm = commIn;
  if (collective) {
    if (platform->cacheLocal) {
      comm = _comm.mpiCommLocal;
    }
    if (platform->cacheBcast) {
      comm = _comm.mpiComm;
    }
  }

  const auto buildRank = [&comm]() {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return (rank == 0) ? true : false;
  }();

  occa::kernel knl;
  if (buildRank) {
    knl = this->wrapperCompileKernel(fileName, props, suffix);
  }
  MPI_Barrier(comm); // finish compilation

#if 0
  const std::string binaryFileName = ...;
  if (collective && platform->cacheBcast) {
    auto dstPath = fs::path(platform->tmpDir) / fs::path("occa/");
    fileBcast(fs::path(binaryFileName).parent_path(), dstPath, comm, platform->verbose);
  }
#endif

  return knl;
}

occa::kernel device_t::loadKernel(const std::string &fileName,
                                  const std::string &kernelName,
                                  const occa::properties &props,
                                  const std::string &suffix) const
{
  if (_compileWhenLoad) {
    this->compileKernel(fileName, props, suffix, platform->comm.mpiComm);
  }

  return this->wrapperLoadKernel(fileName, kernelName, props, suffix);
}

occa::kernel device_t::loadKernel(const std::string &fileName,
                                  const occa::properties &props,
                                  const std::string &suffix) const
{
  return this->loadKernel(fileName, extractKernelName(fileName), props, suffix);
}

occa::memory device_t::mallocHost(size_t Nbytes)
{
  occa::properties props;
  props["host"] = true;

  void *buffer = std::calloc(Nbytes, 1);
  occa::memory h_scratch = _device.malloc(Nbytes, props);
  h_scratch.copyFrom(buffer);
  std::free(buffer);

  return h_scratch;
}

occa::memory device_t::malloc(size_t Nbytes, const occa::properties &properties)
{
  void *buffer = std::calloc(Nbytes, 1);
  occa::memory o_returnValue = _device.malloc(Nbytes, properties);
  o_returnValue.copyFrom(buffer);
  std::free(buffer);
  return o_returnValue;
}

occa::memory device_t::malloc(size_t Nbytes, const void *src, const occa::properties &properties)
{
  auto props = properties;
  if (platform->serial) {
    props["use_host_pointer"] = true;
  }

  occa::memory o_returnValue = _device.malloc(Nbytes, src, props);
  return o_returnValue;
}

occa::memory device_t::malloc(size_t Nword, size_t wordSize, const occa::memory &src)
{
  occa::properties props;
  if (platform->serial) {
    props["use_host_pointer"] = true;
  }

  occa::memory o_returnValue = _device.malloc(Nword * wordSize, src, props);
  return o_returnValue;
}

occa::memory device_t::malloc(size_t Nword, size_t wordSize)
{
  void *buffer = std::calloc(Nword, wordSize);
  occa::memory o_returnValue = _device.malloc(Nword * wordSize);
  o_returnValue.copyFrom(buffer);
  std::free(buffer);
  return o_returnValue;
}

device_t::device_t(setupAide &options, comm_t &comm) : _comm(comm)
{
  // OCCA build stuff
  int deviceConfigSize = 4096;
  char deviceConfig[deviceConfigSize];
  int worldRank = _comm.mpiRank;

  int device_id = 0;

  if (options.compareArgs("DEVICE NUMBER", "LOCAL-RANK")) {
    device_id = _comm.localRank;
  } else {
    options.getArgs("DEVICE NUMBER", device_id);
  }

  occa::properties deviceProps;
  std::string requestedOccaMode;
  options.getArgs("THREAD MODEL", requestedOccaMode);

  if (strcasecmp(requestedOccaMode.c_str(), "CUDA") == 0) {
    if (!getenv("CUDA_CACHE_DISABLE")) {
      setenv("CUDA_CACHE_DISABLE", "1", 1);
    }
    snprintf(deviceConfig, deviceConfigSize, "{mode: 'CUDA', device_id: %d}", device_id);
  } else if (strcasecmp(requestedOccaMode.c_str(), "HIP") == 0) {
    snprintf(deviceConfig, deviceConfigSize, "{mode: 'HIP', device_id: %d}", device_id);
  } else if (strcasecmp(requestedOccaMode.c_str(), "DPCPP") == 0) {
    int plat = 0;
    options.getArgs("PLATFORM NUMBER", plat);
    snprintf(deviceConfig,
             deviceConfigSize,
             "{mode: 'dpcpp', device_id: %d, platform_id: %d}",
             device_id,
             plat);
  } else if (strcasecmp(requestedOccaMode.c_str(), "OPENCL") == 0) {
    int plat = 0;
    options.getArgs("PLATFORM NUMBER", plat);
    snprintf(deviceConfig,
             deviceConfigSize,
             "{mode: 'OpenCL', device_id: %d, platform_id: %d}",
             device_id,
             plat);
  } else if (strcasecmp(requestedOccaMode.c_str(), "OPENMP") == 0) {
    nekrsCheck(true, _comm.mpiComm, EXIT_FAILURE, "%s\n", "OpenMP backend currently not supported!");
    snprintf(deviceConfig, deviceConfigSize, "{mode: 'OpenMP'}");
  } else if (strcasecmp(requestedOccaMode.c_str(), "CPU") == 0 ||
             strcasecmp(requestedOccaMode.c_str(), "SERIAL") == 0) {
    snprintf(deviceConfig, deviceConfigSize, "{mode: 'Serial'}");
    options.setArgs("THREAD MODEL", "SERIAL");
    options.getArgs("THREAD MODEL", requestedOccaMode);
  } else {
    nekrsCheck(true, _comm.mpiComm, EXIT_FAILURE, "%s\n", "Invalid requested backend!");
  }

#if 1
  if (options.compareArgs("BUILD ONLY", "TRUE")) {
    if (!getenv("OCCA_VERBOSE")) {
      occa::settings()["device/verbose"] = true;
      occa::settings()["kernel/verbose"] = true;
      occa::settings()["memory/verbose"] = true;
    }
  }
#endif
  setOccaVars();

  if (worldRank == 0) {
    printf("Initializing device \n");
  }
  this->_device.setup(static_cast<std::string>(deviceConfig));

  if (worldRank == 0) {
    std::cout << "active occa mode: " << this->mode() << "\n\n";
  }

  nekrsCheck(strcasecmp(requestedOccaMode.c_str(), this->mode().c_str()) != 0,
             _comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "Active occa mode does not match selected backend!");

  // overwrite compiler settings to ensure
  // compatability of libocca and kernelLauchner
  if (this->mode() != "Serial") {
    std::string buf;
    buf.assign(getenv("NEKRS_MPI_UNDERLYING_COMPILER"));
    setenv("OCCA_CXX", buf.c_str(), 1);
    buf.assign(getenv("NEKRS_CXXFLAGS"));
    setenv("OCCA_CXXFLAGS", buf.c_str(), 1);
  }

  _compilationEnabled = true;
  _compileWhenLoad = false;
  _device_id = device_id;

  deviceAtomic = atomicsAvailable(*this, _comm.mpiComm);
}

size_t device_t::memoryUsage() const
{
  return platform->device.occaDevice().memoryAllocated();
}

void device_t::printMemoryUsage(MPI_Comm comm) const
{
  const auto maxMemSizes = [&]() {
    std::vector<uint64_t> work;
    work.push_back(platform->device.occaDevice().maxMemoryAllocated());
    work.push_back(platform->deviceMemoryPool.size());
    work.push_back(platform->memoryPool.size());

    MPI_Allreduce(MPI_IN_PLACE, work.data(), work.size(), MPI_UINT64_T, MPI_MAX, comm);
    return std::make_tuple(work[0], work[1], work[2]);
  }();

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank == 0) {
    int width = 12;
    std::cout << "occa max memory usage: " << std::setw(width) << std::right << std::get<0>(maxMemSizes)
              << " bytes" << std::endl;
    std::cout << "  deviceMemoryPool:    " << std::setw(width) << std::right << std::get<1>(maxMemSizes)
              << " bytes" << std::endl;
    std::cout << "  mempool:             " << std::setw(width) << std::right << std::get<2>(maxMemSizes)
              << " bytes" << std::endl
              << std::flush;
  }
}
