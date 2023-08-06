#include <unistd.h>
#include <regex>
#include "nrssys.hpp"
#include "device.hpp"
#include "platform.hpp"
#include "fileUtils.hpp"

namespace {

void setOccaVars()
{
  std::string cache_dir;
  if(getenv("NEKRS_CACHE_DIR"))
    cache_dir.assign(getenv("NEKRS_CACHE_DIR"));

  if (!getenv("OCCA_CACHE_DIR")) {
    const std::string path= cache_dir + "/occa/";
    occa::env::OCCA_CACHE_DIR = path;
    setenv("OCCA_CACHE_DIR", path.c_str(), 1);
  }

  if (!getenv("OCCA_DIR")){
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
    }
    catch (std::exception &e) {
      atomicSupported = 0;
    }
  }

  MPI_Bcast(&atomicSupported, 1, MPI_INT, 0, comm);

  return atomicSupported;
}
} // namespace

occa::kernel device_t::buildNativeKernel(const std::string &fileName,
                                         const std::string &kernelName,
                                         const occa::properties &props) const
{
  occa::properties nativeProperties = props;
  nativeProperties["okl/enabled"] = false;
  if (this->mode() == "OpenMP")
    nativeProperties["defines/__NEKRS__OMP__"] = 1;
  return _device.buildKernel(fileName, kernelName, nativeProperties);
}

occa::kernel device_t::buildKernel(const std::string &fullPath, const occa::properties &props) const
{
  const std::string noSuffix = std::string("");
  return this->buildKernel(fullPath, props, noSuffix);
}

occa::kernel device_t::buildKernel(const std::string &fullPath,
                                   const occa::properties &props,
                                   const std::string &suffix) const
{
  const std::string fileName = fullPath;
  std::string kernelName;
  std::regex kernelNameRegex(R"((.+)\/(.+)\.)");
  std::smatch kernelNameMatch;
  const bool foundKernelName = std::regex_search(fullPath, kernelNameMatch, kernelNameRegex);

  // e.g. /path/to/install/nekrs/kernels/cds/advectMeshVelocityHex3D.okl

  // Full string
  // 0:   /path/to/install/nekrs/kernels/cds/advectMeshVelocityHex3D.okl

  // First capture group
  // 1:   /path/to/install/nekrs/kernels/cds

  // Second capture group (kernel name)
  // 2:   advectMeshVelocityHex3D.okl
  if (foundKernelName) {
    if (kernelNameMatch.size() == 3) {
      kernelName = kernelNameMatch[2].str();
    }
  }

  return this->buildKernel(fileName, kernelName, props, suffix);
}

occa::kernel device_t::buildKernel(const std::string &fileName,
                                   const std::string &kernelName,
                                   const occa::properties &props,
                                   const std::string &suffix) const
{

  if (fileName.find(".okl") != std::string::npos) {
    occa::properties propsWithSuffix = props;
    propsWithSuffix["kernelNameSuffix"] = suffix;

    propsWithSuffix["defines/__okl__"] = 1;

    if (this->mode() == "CUDA")
      propsWithSuffix["defines/smXX"] = 1;
    if (this->mode() == "HIP")
      propsWithSuffix["defines/gfxXX"] = 1;

    const std::string floatingPointType = static_cast<std::string>(propsWithSuffix["defines/dfloat"]);

    if (floatingPointType.find("float") != std::string::npos) {
      propsWithSuffix["defines/FP32"] = 1;
    }

    std::string newKernelName = kernelName;
    if (props.has("defines/p_knl")) {
      const int kernelVariant = static_cast<int>(props["defines/p_knl"]);
      newKernelName += "_v" + std::to_string(kernelVariant);
    };

    return _device.buildKernel(fileName, newKernelName, propsWithSuffix);
  }
  else {
    std::string newKernelName = kernelName;
    if (props.has("defines/p_knl")) {
      const int kernelVariant = static_cast<int>(props["defines/p_knl"]);
      newKernelName += "_v" + std::to_string(kernelVariant);
    };

    occa::properties propsWithSuffix = props;
    propsWithSuffix["defines/SUFFIX"] = suffix;
    propsWithSuffix["defines/TOKEN_PASTE_(a,b)"] = std::string("a##b");
    propsWithSuffix["defines/TOKEN_PASTE(a,b)"] = std::string("TOKEN_PASTE_(a,b)");
    propsWithSuffix["defines/FUNC(a)"] = std::string("TOKEN_PASTE(a,SUFFIX)");
    newKernelName += suffix;

    return this->buildNativeKernel(fileName, newKernelName, propsWithSuffix);
  }
}

occa::kernel device_t::buildKernel(const std::string &fileName,
                                   const std::string &kernelName,
                                   const occa::properties &props) const
{
  const std::string suffix("");
  const int rank = platform->cacheLocal ? _comm.localRank : _comm.mpiRank;
  MPI_Comm localCommunicator = platform->cacheLocal ? _comm.mpiCommLocal : _comm.mpiComm;

  const auto OCCA_CACHE_DIR0 = std::string(occa::env::OCCA_CACHE_DIR);
  const auto OCCA_CACHE_DIR = platform->cacheBcast ? 
                              std::string(platform->tmpDir / fs::path("occa/")) :
                              OCCA_CACHE_DIR0;

  occa::kernel constructedKernel;

  // rank0 compiles + loads, then all other just load
  for (int pass = 0; pass < 2; ++pass) {
    occa::env::OCCA_CACHE_DIR = (pass == 0) ? OCCA_CACHE_DIR0 : OCCA_CACHE_DIR; 
    if ((pass == 0 && rank == 0) || (pass == 1 && rank != 0)) {
      constructedKernel = this->buildKernel(fileName, kernelName, props, suffix);
    }

    if(pass == 0) {
      if(platform->cacheBcast) {
        const auto srcPath = (fs::path(constructedKernel.binaryFilename()).parent_path());
        const auto dstPath = OCCA_CACHE_DIR / fs::path("cache/");
        fileBcast(srcPath, dstPath, _comm.mpiComm, platform->verbose);
      } else {
        MPI_Barrier(localCommunicator);
      }
    }
    occa::env::OCCA_CACHE_DIR = OCCA_CACHE_DIR0;
  }

  int isInitializedMin = constructedKernel.isInitialized();
  MPI_Allreduce(MPI_IN_PLACE, &isInitializedMin, 1, MPI_INT, MPI_MIN, _comm.mpiComm);
  int isInitializedMax = constructedKernel.isInitialized();
  MPI_Allreduce(MPI_IN_PLACE, &isInitializedMax, 1, MPI_INT, MPI_MAX, _comm.mpiComm);
  nrsCheck(isInitializedMin != isInitializedMax,  _comm.mpiComm, EXIT_FAILURE,
           "Kernel status of %s inconsistent across ranks\n", constructedKernel.name().c_str());

  return constructedKernel;
}

occa::kernel device_t::buildKernel(const std::string &fullPath,
                                   const occa::properties &props,
                                   const std::string &suffix,
                                   bool buildRank0) const
{
  occa::kernel constructedKernel;

  if (buildRank0) {
    const auto OCCA_CACHE_DIR0 = std::string(occa::env::OCCA_CACHE_DIR);
    const auto OCCA_CACHE_DIR = platform->cacheBcast ? 
                                std::string(platform->tmpDir / fs::path("occa/")) :
                                OCCA_CACHE_DIR0;

    const int rank = platform->cacheLocal ? _comm.localRank : _comm.mpiRank;
    MPI_Comm localCommunicator = platform->cacheLocal ? _comm.mpiCommLocal : _comm.mpiComm;

    // rank0 compiles + loads, then all other just load
    for (int pass = 0; pass < 2; ++pass) {
      occa::env::OCCA_CACHE_DIR = (pass == 0) ? OCCA_CACHE_DIR0 : OCCA_CACHE_DIR ;
      if ((pass == 0 && rank == 0) || (pass == 1 && rank != 0)) {
        constructedKernel = this->buildKernel(fullPath, props, suffix);
      }

      if(pass == 0) {
        if(platform->cacheBcast) {
          const auto srcPath = (fs::path(constructedKernel.binaryFilename()).parent_path());
          const auto dstPath = OCCA_CACHE_DIR / fs::path("cache/");
          fileBcast(srcPath, dstPath, _comm.mpiComm, platform->verbose);
        } else {
          MPI_Barrier(localCommunicator);
        }
      }
      occa::env::OCCA_CACHE_DIR = OCCA_CACHE_DIR0;
    }

    int isInitializedMin = constructedKernel.isInitialized();
    MPI_Allreduce(MPI_IN_PLACE, &isInitializedMin, 1, MPI_INT, MPI_MIN, _comm.mpiComm);
    int isInitializedMax = constructedKernel.isInitialized();
    MPI_Allreduce(MPI_IN_PLACE, &isInitializedMax, 1, MPI_INT, MPI_MAX, _comm.mpiComm);
    nrsCheck(isInitializedMin != isInitializedMax,  _comm.mpiComm, EXIT_FAILURE,
             "Kernel status of %s inconsistent across ranks\n", constructedKernel.name().c_str());

  } else {
    constructedKernel = this->buildKernel(fullPath, props, suffix);
  }

  return constructedKernel;
}

occa::kernel
device_t::buildKernel(const std::string &fullPath, const occa::properties &props, bool buildRank0) const
{
  std::string noSuffix = std::string("");
  return this->buildKernel(fullPath, props, noSuffix, buildRank0);
}

occa::memory device_t::mallocHost(size_t Nbytes)
{
  occa::properties props;
  props["host"] = true;

  void *buffer = std::calloc(Nbytes, 1);
  occa::memory h_scratch = _device.malloc(Nbytes, buffer, props);
  std::free(buffer);
  return h_scratch;
}

occa::memory device_t::malloc(size_t Nbytes, const occa::properties &properties)
{
  void *buffer = std::calloc(Nbytes, 1);
  occa::memory o_returnValue = _device.malloc(Nbytes, buffer, properties);
  std::free(buffer);
  return o_returnValue;
}

occa::memory device_t::malloc(size_t Nbytes, const void *src, const occa::properties &properties)
{
  void *buffer;
  buffer = std::calloc(Nbytes, 1);
  const void *init_ptr = (src) ? src : buffer;
  occa::memory o_returnValue = _device.malloc(Nbytes, init_ptr, properties);
  std::free(buffer);
  return o_returnValue;
}

occa::memory device_t::malloc(size_t Nword, size_t wordSize, occa::memory src)
{
  return _device.malloc(Nword * wordSize, src);
}

occa::memory device_t::malloc(size_t Nword, size_t wordSize)
{
  void *buffer = std::calloc(Nword, wordSize);
  occa::memory o_returnValue = _device.malloc(Nword * wordSize, buffer);
  std::free(buffer);
  return o_returnValue;
}

device_t::device_t(setupAide &options, comm_t &comm) : _comm(comm)
{
  // OCCA build stuff
  char deviceConfig[4096];
  int worldRank = _comm.mpiRank;

  int device_id = 0;

  if (options.compareArgs("DEVICE NUMBER", "LOCAL-RANK")) {
    device_id = _comm.localRank;
  }
  else {
    options.getArgs("DEVICE NUMBER", device_id);
  }

  occa::properties deviceProps;
  std::string requestedOccaMode;
  options.getArgs("THREAD MODEL", requestedOccaMode);

  if (strcasecmp(requestedOccaMode.c_str(), "CUDA") == 0) {
    if(!getenv("CUDA_CACHE_DISABLE"))
      setenv("CUDA_CACHE_DISABLE", "1", 1);
    sprintf(deviceConfig, "{mode: 'CUDA', device_id: %d}", device_id);
  }
  else if (strcasecmp(requestedOccaMode.c_str(), "HIP") == 0) {
    sprintf(deviceConfig, "{mode: 'HIP', device_id: %d}", device_id);
  }
  else if(strcasecmp(requestedOccaMode.c_str(), "DPCPP") == 0) {
    int plat = 0;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "{mode: 'dpcpp', device_id: %d, platform_id: %d}", device_id, plat);
  }
  else if (strcasecmp(requestedOccaMode.c_str(), "OPENCL") == 0) {
    int plat = 0;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "{mode: 'OpenCL', device_id: %d, platform_id: %d}", device_id, plat);
  }
  else if (strcasecmp(requestedOccaMode.c_str(), "OPENMP") == 0) {
    nrsCheck(true, _comm.mpiComm, EXIT_FAILURE,
             "%s\n", "OpenMP backend currently not supported!"); 
    sprintf(deviceConfig, "{mode: 'OpenMP'}");
  }
  else if (strcasecmp(requestedOccaMode.c_str(), "CPU") == 0 ||
           strcasecmp(requestedOccaMode.c_str(), "SERIAL") == 0) {
    sprintf(deviceConfig, "{mode: 'Serial'}");
    options.setArgs("THREAD MODEL", "SERIAL");
    options.getArgs("THREAD MODEL", requestedOccaMode);
  }
  else {
    nrsCheck(true, _comm.mpiComm, EXIT_FAILURE,
             "%s\n", "Invalid requested backend!"); 
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

  if (worldRank == 0)
    printf("Initializing device \n");
  this->_device.setup((std::string)deviceConfig);

  if (worldRank == 0)
    std::cout << "active occa mode: " << this->mode() << "\n\n";

  nrsCheck(strcasecmp(requestedOccaMode.c_str(), this->mode().c_str()) != 0, 
           _comm.mpiComm, EXIT_FAILURE,
           "%s\n", "Active occa mode does not match selected backend!"); 

  // overwrite compiler settings to ensure
  // compatability of libocca and kernelLauchner
  if (this->mode() != "Serial") {
    std::string buf;
    buf.assign(getenv("NEKRS_MPI_UNDERLYING_COMPILER"));
    setenv("OCCA_CXX", buf.c_str(), 1);
    buf.assign(getenv("NEKRS_CXXFLAGS"));
    setenv("OCCA_CXXFLAGS", buf.c_str(), 1);
  }

  _device_id = device_id;

  deviceAtomic = atomicsAvailable(*this, _comm.mpiComm);
}
