#include <cstdlib>
#include <strings.h>
#include "platform.hpp"
#include "nrs.hpp"
#include "linAlg.hpp"
#include "omp.h"

comm_t::comm_t(MPI_Comm _comm)
{
  mpiComm = _comm;
  MPI_Comm_rank(_comm, &mpiRank);
  MPI_Comm_size(_comm, &mpiCommSize);
}

deviceVector_t::deviceVector_t(const dlong _vectorSize, const dlong _nVectors, const dlong _wordSize, const std::string _vectorName)
: vectorSize(_vectorSize),
  nVectors(_nVectors),
  wordSize(_wordSize),
  vectorName(_vectorName)
{
  if(vectorSize <= 0 || nVectors <= 0 || wordSize <= 0) return; // bail
  o_vector = platform->device.malloc(vectorSize * nVectors, wordSize);
  // set slices
  for(int s = 0 ; s < nVectors; ++s){
    slices.push_back(o_vector + s * vectorSize * wordSize);
  }
}

occa::memory&
deviceVector_t::at(const int i)
{
  if(i >= nVectors){
    if(platform->comm.mpiRank == 0){
      printf("ERROR: deviceVector_t(%s) has %d size, but an attempt to access entry %i was made!\n",
        vectorName.c_str(),
        nVectors,
        i
      );
    }
    ABORT(EXIT_FAILURE);
    return o_vector;
  }
  occa::memory slice = o_vector + i * vectorSize * wordSize;
  return slices[i];
}




platform_t* platform_t::singleton = nullptr;
platform_t::platform_t(setupAide& _options, MPI_Comm _comm)
: options(_options),
  warpSize(32), // CUDA specific warp size
  device(options, _comm),
  timer(_comm, device, 0),
  comm(_comm)
{
  kernelInfo["defines/" "p_NVec"] = 3;
  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "dfloat"] = dfloatString;
  kernelInfo["defines/" "pfloat"] = pfloatString;
  kernelInfo["defines/" "dlong"] = dlongString;
  kernelInfo["defines/" "hlong"] = hlongString;

  if(device.mode() == "CUDA" && !getenv("OCCA_CUDA_COMPILER_FLAGS")) {
    kernelInfo["compiler_flags"] += "--ftz=true ";
    kernelInfo["compiler_flags"] += "--prec-div=false ";
    kernelInfo["compiler_flags"] += "--prec-sqrt=false ";
    kernelInfo["compiler_flags"] += "--use_fast_math ";
    kernelInfo["compiler_flags"] += "--fmad=true ";

    //kernelInfo["compiler_flags"] += "-Xptxas -dlcm=ca";
  }

  if(device.mode() == "OpenCL") {
    if(!getenv("OCCA_OPENCL_COMPILER_FLAGS")) {
      kernelInfo["compiler_flags"] += " -cl-std=CL2.0 ";
      kernelInfo["compiler_flags"] += " -cl-mad-enable ";
      kernelInfo["compiler_flags"] += " -cl-no-signed-zeros ";
      kernelInfo["compiler_flags"] += " -cl-unsafe-math-optimizations ";
      kernelInfo["compiler_flags"] += " -cl-fast-relaxed-math ";
    }
    kernelInfo["defines/" "hlong"]="long";
  }

  if(device.mode() == "HIP" && !getenv("OCCA_HIP_COMPILER_FLAGS")) {
    warpSize = 64;
    kernelInfo["compiler_flags"] += " -O3 ";
    kernelInfo["compiler_flags"] += " -ffp-contract=fast ";
    kernelInfo["compiler_flags"] += " -funsafe-math-optimizations ";
    kernelInfo["compiler_flags"] += " -ffast-math ";
  }
}
void memPool_t::allocate(const dlong offset, const dlong fields)
{
  ptr = (dfloat*) calloc(offset*fields, sizeof(dfloat));
  slice0 = ptr + 0 * offset;
  slice1 = ptr + 1 * offset;
  slice2 = ptr + 2 * offset;
  slice3 = ptr + 3 * offset;
  slice4 = ptr + 4 * offset;
  slice5 = ptr + 5 * offset;
  slice6 = ptr + 6 * offset;
  slice7 = ptr + 7 * offset;
  slice9 = ptr + 9 * offset;
  slice12 = ptr + 12 * offset;
  slice15 = ptr + 15 * offset;
  slice18 = ptr + 18 * offset;
  slice19 = ptr + 19 * offset;
}
void deviceMemPool_t::allocate(memPool_t& hostMemory, const dlong offset, const dlong fields)
{
  bytesAllocated = fields * offset * sizeof(dfloat);
  o_ptr = platform->device.malloc(offset*fields*sizeof(dfloat), hostMemory.slice0);
  slice0 = o_ptr.slice(0 * offset * sizeof(dfloat));
  slice1 = o_ptr.slice(1 * offset * sizeof(dfloat));
  slice2 = o_ptr.slice(2 * offset * sizeof(dfloat));
  slice3 = o_ptr.slice(3 * offset * sizeof(dfloat));
  slice4 = o_ptr.slice(4 * offset * sizeof(dfloat));
  slice5 = o_ptr.slice(5 * offset * sizeof(dfloat));
  slice6 = o_ptr.slice(6 * offset * sizeof(dfloat));
  slice7 = o_ptr.slice(7 * offset * sizeof(dfloat));
  slice9 = o_ptr.slice(9 * offset * sizeof(dfloat));
  slice12 = o_ptr.slice(12 * offset * sizeof(dfloat));
  slice15 = o_ptr.slice(15 * offset * sizeof(dfloat));
  slice18 = o_ptr.slice(18 * offset * sizeof(dfloat));
  slice19 = o_ptr.slice(19 * offset * sizeof(dfloat));
}
void
platform_t::create_mempool(const dlong offset, const dlong fields)
{
  mempool.allocate(offset, fields);
  o_mempool.allocate(mempool, offset, fields);
}

occa::kernel
device_t::buildNativeKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props) const
{
  occa::properties nativeProperties = props;
  nativeProperties["okl/enabled"] = false;
  return this->buildKernel(filename, kernelName, nativeProperties, comm);
}
occa::kernel
device_t::buildKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props) const
{
  return this->buildKernel(filename, kernelName, props, comm);
}
occa::kernel
device_t::buildKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props,
                         MPI_Comm comm) const
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  occa::kernel _kernel;
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
      _kernel = occa::device::buildKernel(filename, kernelName, props);
    }
    MPI_Barrier(comm);
  }
  return _kernel;
}
occa::memory
device_t::malloc(const dlong Nbytes, const occa::properties& properties)
{
  return occa::device::malloc(Nbytes, nullptr, properties);
}
occa::memory
device_t::malloc(const dlong Nbytes, const void* src, const occa::properties& properties)
{
  if(!src){
    if(Nbytes > bufferSize)
    {
      if(bufferSize > 0) std::free(_buffer);
      _buffer = std::calloc(Nbytes, 1);
      bufferSize = Nbytes;
    }
  }
  const void* init_ptr = (src) ? src : _buffer;
  return occa::device::malloc(Nbytes, init_ptr, properties);
}
occa::memory
device_t::malloc(const dlong Nword , const dlong wordSize, occa::memory src)
{
  return occa::device::malloc(Nword * wordSize, src);
}
occa::memory
device_t::malloc(const dlong Nword , const dlong wordSize)
{
  const dlong Nbytes = Nword * wordSize;
  if(Nbytes > bufferSize)
  {
    if(bufferSize > 0) std::free(_buffer);
    _buffer = std::calloc(Nword, wordSize);
    bufferSize = Nbytes;
  }
  return occa::device::malloc(Nword * wordSize, _buffer);
}
device_t::device_t(setupAide& options, MPI_Comm comm)
{
  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int device_id = 0;

  if(options.compareArgs("DEVICE NUMBER", "LOCAL-RANK")) {
    long int hostId = gethostid();

    long int* hostIds = (long int*) std::calloc(size,sizeof(long int));
    MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,comm);

    int totalDevices = 0;
    for (int r = 0; r < rank; r++)
      if (hostIds[r] == hostId) device_id++;
    for (int r = 0; r < size; r++)
      if (hostIds[r] == hostId) totalDevices++;
  } else {
    options.getArgs("DEVICE NUMBER",device_id);
  }

  occa::properties deviceProps;
  string requestedOccaMode; 
  options.getArgs("THREAD MODEL", requestedOccaMode);

  if(strcasecmp(requestedOccaMode.c_str(), "CUDA") == 0) {
    sprintf(deviceConfig, "{mode: 'CUDA', device_id: %d}", device_id);
  }else if(strcasecmp(requestedOccaMode.c_str(), "HIP") == 0) {
    sprintf(deviceConfig, "{mode: 'HIP', device_id: %d}",device_id);
  }else if(strcasecmp(requestedOccaMode.c_str(), "OPENCL") == 0) {
    int plat;
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
    buf.assign(getenv("NEKRS_CXX"));
    setenv("OCCA_CXX", buf.c_str(), 1);
    buf.assign(getenv("NEKRS_CXXFLAGS"));
    setenv("OCCA_CXXFLAGS", buf.c_str(), 1);
  }

  int Nthreads = 1;
  if(this->mode() != "OpenMP") omp_set_num_threads(Nthreads);

  bufferSize = 0;

  _device_id = device_id;
}
