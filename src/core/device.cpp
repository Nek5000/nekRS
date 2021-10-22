#include "device.hpp" 
#include "platform.hpp"
#include <unistd.h>

occa::kernel
device_t::buildNativeKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props) const
{
  occa::properties nativeProperties = props;
  nativeProperties["okl/enabled"] = false;
  if(platform->options.compareArgs("BUILD ONLY", "TRUE"))
    nativeProperties["verbose"] = true;
  if(platform->device.mode() == "OpenMP")
    nativeProperties["defines/__NEKRS__OMP__"] = 1;
  return occa::device::buildKernel(filename, kernelName, nativeProperties);
}

occa::kernel
device_t::doBuildKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props,
                         std::string suffix) const
{
  if(filename.find(".okl") != std::string::npos){
    occa::properties propsWithSuffix = props;
    propsWithSuffix["kernelNameSuffix"] = suffix;
    if(platform->options.compareArgs("BUILD ONLY", "TRUE"))
      propsWithSuffix["verbose"] = true;
    return occa::device::buildKernel(filename, kernelName, propsWithSuffix);
  }
  else{
    occa::properties propsWithSuffix = props;
    propsWithSuffix["defines/SUFFIX"] = suffix;
    propsWithSuffix["defines/TOKEN_PASTE_(a,b)"] = std::string("a##b");
    propsWithSuffix["defines/TOKEN_PASTE(a,b)"] = std::string("TOKEN_PASTE_(a,b)");
    propsWithSuffix["defines/FUNC(a)"] = std::string("TOKEN_PASTE(a,SUFFIX)");
    const std::string alteredName =  kernelName + suffix;
    return this->buildNativeKernel(filename, alteredName, propsWithSuffix);
  }
}

occa::kernel
device_t::buildKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props,
                         std::string suffix,
                         bool buildRank0) const
{

  if(buildRank0){

    const bool buildNodeLocal = useNodeLocalCache();
    const int rank = buildNodeLocal ? platform->comm.localRank : platform->comm.mpiRank;
    MPI_Comm localCommunicator = buildNodeLocal ? platform->comm.mpiCommLocal : platform->comm.mpiComm;
    occa::kernel constructedKernel;
    for(int pass = 0; pass < 2; ++pass){
      if((pass == 0 && rank == 0) || (pass == 1 && rank != 0)){
        constructedKernel = doBuildKernel(filename, kernelName, props, suffix);
      }
      MPI_Barrier(localCommunicator);
    }
    return constructedKernel;

  }

  return doBuildKernel(filename, kernelName, props, suffix);

}

occa::memory
device_t::mallocHost(const size_t Nbytes)
{
  occa::properties props;
  props["host"] = true;
  
  void* buffer = std::calloc(Nbytes, 1);
  occa::memory h_scratch = occa::device::malloc(Nbytes, buffer, props);
  std::free(buffer);
  return h_scratch;
}

occa::memory
device_t::malloc(const size_t Nbytes, const occa::properties& properties)
{
  void* buffer = std::calloc(Nbytes, 1);
  occa::memory o_returnValue = occa::device::malloc(Nbytes, buffer, properties);
  std::free(buffer);
  return o_returnValue;
}

occa::memory
device_t::malloc(const size_t Nbytes, const void* src, const occa::properties& properties)
{
  void* buffer;
  buffer = std::calloc(Nbytes, 1);
  const void* init_ptr = (src) ? src : buffer;
  occa::memory o_returnValue = occa::device::malloc(Nbytes, init_ptr, properties);
  std::free(buffer);
  return o_returnValue;
}

occa::memory
device_t::malloc(const hlong Nword , const dlong wordSize, occa::memory src)
{
  return occa::device::malloc(Nword * wordSize, src);
}

occa::memory
device_t::malloc(const hlong Nword , const dlong wordSize)
{
  const size_t Nbytes = Nword * wordSize;
  void* buffer = std::calloc(Nword, wordSize);
  occa::memory o_returnValue = occa::device::malloc(Nword * wordSize, buffer);
  std::free(buffer);
  return o_returnValue;
}

device_t::device_t(setupAide& options, MPI_Comm commParent, MPI_Comm comm)
{
  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int device_id = 0;

  if(options.compareArgs("DEVICE NUMBER", "LOCAL-RANK")) {
    int worldRank;
    MPI_Comm_rank(commParent, &worldRank);

    MPI_Comm commLocal;
    MPI_Comm_split_type(commParent, MPI_COMM_TYPE_SHARED, worldRank, MPI_INFO_NULL, &commLocal);

    int localRank;
    MPI_Comm_rank(commLocal, &localRank);
    device_id = localRank;
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
    if(rank == 0) printf("OpenMP backend currently not supported!\n");
    ABORT(EXIT_FAILURE);
    sprintf(deviceConfig, "{mode: 'OpenMP'}");
  }else if(strcasecmp(requestedOccaMode.c_str(), "CPU") == 0 ||
           strcasecmp(requestedOccaMode.c_str(), "SERIAL") == 0) {
    sprintf(deviceConfig, "{mode: 'Serial'}");
    options.setArgs("THREAD MODEL", "SERIAL");
    options.getArgs("THREAD MODEL", requestedOccaMode);
  } else {
    if(rank == 0) printf("Invalid requested backend!\n");
    ABORT(EXIT_FAILURE);
  }

  if(rank == 0) printf("Initializing device \n");
  this->setup((std::string)deviceConfig);
  this->comm = comm;
 
  if(rank == 0)
    std::cout << "active occa mode: " << this->mode() << "\n\n";

  if(strcasecmp(requestedOccaMode.c_str(), this->mode().c_str()) != 0) {
    if(rank == 0) printf("active occa mode does not match selected backend!\n");
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