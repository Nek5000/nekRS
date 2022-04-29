#include "device.hpp" 
#include "platform.hpp"
#include <unistd.h>
#include <regex>

occa::kernel
device_t::buildNativeKernel(const std::string &fileName,
                         const std::string &kernelName,
                         const occa::properties &props) const
{
  occa::properties nativeProperties = props;
  nativeProperties["okl/enabled"] = false;
  if(_verbose)
    nativeProperties["verbose"] = true;
  if(this->mode() == "OpenMP")
    nativeProperties["defines/__NEKRS__OMP__"] = 1;
  return _device.buildKernel(fileName, kernelName, nativeProperties);
}

occa::kernel
device_t::buildKernel(const std::string &fullPath,
                      const occa::properties &props) const
{
  const std::string noSuffix = std::string("");
  return this->buildKernel(fullPath, props, noSuffix);
}

occa::kernel
device_t::buildKernel(const std::string &fullPath,
                      const occa::properties &props,
                      const std::string & suffix) const
{
  const std::string fileName = fullPath;
  std::string kernelName;
  std::regex kernelNameRegex(R"((.+)\/(.+)\.)");
  std::smatch kernelNameMatch;
  const bool foundKernelName = std::regex_search(fullPath, kernelNameMatch, kernelNameRegex);

  // e.g. /path/to/install/nekrs/okl/cds/advectMeshVelocityHex3D.okl

  // Full string
  // 0:   /path/to/install/nekrs/okl/cds/advectMeshVelocityHex3D.okl

  // First capture group
  // 1:   /path/to/install/nekrs/okl/cds

  // Second capture group (kernel name)
  // 2:   advectMeshVelocityHex3D.okl
  if(foundKernelName){
    if(kernelNameMatch.size() == 3){
      kernelName = kernelNameMatch[2].str();
    }
  }

  return this->buildKernel(fileName, kernelName, props, suffix);
}

occa::kernel
device_t::buildKernel(const std::string &fileName,
                             const std::string &kernelName,
                             const occa::properties &props,
                             const std::string& suffix) const
{

  if(fileName.find(".okl") != std::string::npos){
    occa::properties propsWithSuffix = props;
    propsWithSuffix["kernelNameSuffix"] = suffix;
    if(_verbose)
      propsWithSuffix["verbose"] = true;

    if (this->mode() == "CUDA")
      propsWithSuffix["defines/smXX"] = 1;
    if (this->mode() == "HIP")
      propsWithSuffix["defines/gfxXX"] = 1;

    const std::string floatingPointType = static_cast<std::string>(propsWithSuffix["defines/dfloat"]);

    if (floatingPointType.find("float") != std::string::npos) {
      propsWithSuffix["defines/FP32"] = 1;
    }

    // if p_knl is defined, add _v(p_knl) to the kernel name
    std::string newKernelName = kernelName;
    if (props.has("defines/p_knl")) {
      const int kernelVariant = static_cast<int>(props["defines/p_knl"]);
      newKernelName += "_v" + std::to_string(kernelVariant);
    };

    return _device.buildKernel(fileName, newKernelName, propsWithSuffix);
  }
  else{
    occa::properties propsWithSuffix = props;
    propsWithSuffix["defines/SUFFIX"] = suffix;
    propsWithSuffix["defines/TOKEN_PASTE_(a,b)"] = std::string("a##b");
    propsWithSuffix["defines/TOKEN_PASTE(a,b)"] = std::string("TOKEN_PASTE_(a,b)");
    propsWithSuffix["defines/FUNC(a)"] = std::string("TOKEN_PASTE(a,SUFFIX)");
    const std::string alteredName =  kernelName + suffix;
    return this->buildNativeKernel(fileName, alteredName, propsWithSuffix);
  }
}

occa::kernel
device_t::buildKernel(const std::string &fileName,
                             const std::string &kernelName,
                             const occa::properties &props) const
{

  const std::string suffix("");
  const bool buildNodeLocal = useNodeLocalCache();
  const int rank = buildNodeLocal ? _comm.localRank : _comm.mpiRank;
  MPI_Comm localCommunicator = buildNodeLocal ? _comm.mpiCommLocal : _comm.mpiComm;

  occa::kernel constructedKernel;
  for(int pass = 0; pass < 2; ++pass){
    if((pass == 0 && rank == 0) || (pass == 1 && rank != 0)){
      constructedKernel = this->buildKernel(fileName, kernelName, props, suffix);
    }
    MPI_Barrier(localCommunicator);
  }
  return constructedKernel;

}

occa::kernel
device_t::buildKernel(const std::string &fullPath,
                         const occa::properties &props,
                         const std::string & suffix,
                         bool buildRank0) const
{

  if(buildRank0){

    const bool buildNodeLocal = useNodeLocalCache();
    const int rank = buildNodeLocal ? _comm.localRank : _comm.mpiRank;
    MPI_Comm localCommunicator = buildNodeLocal ? _comm.mpiCommLocal : _comm.mpiComm;
    occa::kernel constructedKernel;
    for(int pass = 0; pass < 2; ++pass){
      if((pass == 0 && rank == 0) || (pass == 1 && rank != 0)){
        constructedKernel = this->buildKernel(fullPath, props, suffix);
      }
      MPI_Barrier(localCommunicator);
    }
    return constructedKernel;

  }

  return this->buildKernel(fullPath, props, suffix);

}

occa::kernel
device_t::buildKernel(const std::string &fullPath,
                         const occa::properties &props,
                         bool buildRank0) const
{
  std::string noSuffix = std::string("");
  return this->buildKernel(fullPath, props, noSuffix, buildRank0);
}

occa::memory
device_t::mallocHost(const size_t Nbytes)
{
  occa::properties props;
  props["host"] = true;
  
  void* buffer = std::calloc(Nbytes, 1);
  occa::memory h_scratch = _device.malloc(Nbytes, buffer, props);
  std::free(buffer);
  return h_scratch;
}

occa::memory
device_t::malloc(const size_t Nbytes, const occa::properties& properties)
{
  void* buffer = std::calloc(Nbytes, 1);
  occa::memory o_returnValue = _device.malloc(Nbytes, buffer, properties);
  std::free(buffer);
  return o_returnValue;
}

occa::memory
device_t::malloc(const size_t Nbytes, const void* src, const occa::properties& properties)
{
  void* buffer;
  buffer = std::calloc(Nbytes, 1);
  const void* init_ptr = (src) ? src : buffer;
  occa::memory o_returnValue = _device.malloc(Nbytes, init_ptr, properties);
  std::free(buffer);
  return o_returnValue;
}

occa::memory
device_t::malloc(const hlong Nword , const dlong wordSize, occa::memory src)
{
  return _device.malloc(Nword * wordSize, src);
}

occa::memory
device_t::malloc(const hlong Nword , const dlong wordSize)
{
  const size_t Nbytes = Nword * wordSize;
  void* buffer = std::calloc(Nword, wordSize);
  occa::memory o_returnValue = _device.malloc(Nword * wordSize, buffer);
  std::free(buffer);
  return o_returnValue;
}

device_t::device_t(setupAide& options, comm_t& comm)
:_comm(comm)
{
  _verbose = options.compareArgs("BUILD ONLY", "TRUE");

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int worldRank = _comm.mpiRank;

  int device_id = 0;

  if(options.compareArgs("DEVICE NUMBER", "LOCAL-RANK")) {
    device_id = _comm.localRank;
  } else {
    options.getArgs("DEVICE NUMBER",device_id);
  }

  occa::properties deviceProps;
  std::string requestedOccaMode; 
  options.getArgs("THREAD MODEL", requestedOccaMode);

  if(strcasecmp(requestedOccaMode.c_str(), "CUDA") == 0) {
    sprintf(deviceConfig, "{mode: 'CUDA', device_id: %d}", device_id);
  }else if(strcasecmp(requestedOccaMode.c_str(), "HIP") == 0) {
    sprintf(deviceConfig, "{mode: 'HIP', device_id: %d}",device_id);
  }else if(strcasecmp(requestedOccaMode.c_str(), "OPENCL") == 0) {
    int plat = 0;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "{mode: 'OpenCL', device_id: %d, platform_id: %d}", device_id, plat);
  }else if(strcasecmp(requestedOccaMode.c_str(), "OPENMP") == 0) {
    if(worldRank == 0) printf("OpenMP backend currently not supported!\n");
    ABORT(EXIT_FAILURE);
    sprintf(deviceConfig, "{mode: 'OpenMP'}");
  }else if(strcasecmp(requestedOccaMode.c_str(), "CPU") == 0 ||
           strcasecmp(requestedOccaMode.c_str(), "SERIAL") == 0) {
    sprintf(deviceConfig, "{mode: 'Serial'}");
    options.setArgs("THREAD MODEL", "SERIAL");
    options.getArgs("THREAD MODEL", requestedOccaMode);
  } else {
    if(worldRank == 0) printf("Invalid requested backend!\n");
    ABORT(EXIT_FAILURE);
  }

  if(worldRank == 0) printf("Initializing device \n");
  this->_device.setup((std::string)deviceConfig);
 
  if(worldRank == 0)
    std::cout << "active occa mode: " << this->mode() << "\n\n";

  if(strcasecmp(requestedOccaMode.c_str(), this->mode().c_str()) != 0) {
    if(worldRank == 0) printf("active occa mode does not match selected backend!\n");
    ABORT(EXIT_FAILURE);
  } 

  // overwrite compiler settings to ensure
  // compatability of libocca and kernelLauchner 
  if(this->mode() != "Serial") {
    std::string buf;
    buf.assign(getenv("NEKRS_MPI_UNDERLYING_COMPILER"));
    setenv("OCCA_CXX", buf.c_str(), 1);
    buf.assign(getenv("NEKRS_CXXFLAGS"));
    setenv("OCCA_CXXFLAGS", buf.c_str(), 1);
  }

  _device_id = device_id;

  deviceAtomic = this->mode() == "CUDA";
}