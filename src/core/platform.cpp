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

platform_t* platform_t::singleton = nullptr;
platform_t::platform_t(setupAide& options, MPI_Comm _comm)
: device(options, _comm),
  timer(_comm, device, 0),
  comm(_comm)
{
  kernelInfo["defines/" "p_NVec"] = 3;
  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "dfloat4"] = dfloatString"4";
  kernelInfo["defines/" "dfloat8"] = dfloatString"8";
  kernelInfo["defines/" "dfloat"] = dfloatString;
  kernelInfo["defines/" "dlong"] = dlongString;
  kernelInfo["defines/" "hlong"] = hlongString;

  if(device.mode() == "CUDA") { // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += "--ftz=true ";
    kernelInfo["compiler_flags"] += "--prec-div=false ";
    kernelInfo["compiler_flags"] += "--prec-sqrt=false ";
    kernelInfo["compiler_flags"] += "--use_fast_math ";
    kernelInfo["compiler_flags"] += "--fmad=true "; // compiler option for cuda
    //kernelInfo["compiler_flags"] += "-Xptxas -dlcm=ca";
  }

  if(device.mode() == "OpenCL") { // add backend compiler optimization for OPENCL
    kernelInfo["compiler_flags"] += " -cl-std=CL2.0 ";
    kernelInfo["compiler_flags"] += " -cl-strict-aliasing ";
    kernelInfo["compiler_flags"] += " -cl-mad-enable ";
    kernelInfo["compiler_flags"] += " -cl-no-signed-zeros ";
    kernelInfo["compiler_flags"] += " -cl-unsafe-math-optimizations ";
    kernelInfo["compiler_flags"] += " -cl-fast-relaxed-math ";
  }

  if(device.mode() == "HIP") { // add backend compiler optimization for HIP
    kernelInfo["compiler_flags"] += " -O3 ";
    kernelInfo["compiler_flags"] += " -ffp-contract=fast ";
    // kernelInfo["compiler_flags"] += " -funsafe-math-optimizations ";
    // kernelInfo["compiler_flags"] += " -ffast-math ";
  }
}
void
platform_t::create_mempool(const dlong offset, const dlong fields)
{
  mempool = (dfloat*) calloc(offset*fields, sizeof(dfloat));
  o_mempool = device.malloc(offset*fields*sizeof(dfloat));
  o_slice0 = o_mempool.slice(0 * offset * sizeof(dfloat));
  o_slice1 = o_mempool.slice(1 * offset * sizeof(dfloat));
  o_slice2 = o_mempool.slice(2 * offset * sizeof(dfloat));
  o_slice3 = o_mempool.slice(3 * offset * sizeof(dfloat));
  o_slice4 = o_mempool.slice(4 * offset * sizeof(dfloat));
  o_slice5 = o_mempool.slice(5 * offset * sizeof(dfloat));
  o_slice6 = o_mempool.slice(6 * offset * sizeof(dfloat));
  o_slice7 = o_mempool.slice(7 * offset * sizeof(dfloat));
  o_slice9 = o_mempool.slice(9 * offset * sizeof(dfloat));
  o_slice12 = o_mempool.slice(12 * offset * sizeof(dfloat));
  o_slice15 = o_mempool.slice(15 * offset * sizeof(dfloat));
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

    long int* hostIds = (long int*) calloc(size,sizeof(long int));
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

  if(options.compareArgs("THREAD MODEL", "CUDA")) {
    sprintf(deviceConfig, "mode: 'CUDA', device_id: %d", device_id);
  }else if(options.compareArgs("THREAD MODEL", "HIP"))  {
    sprintf(deviceConfig, "mode: 'HIP', device_id: %d",device_id);
  }else if(options.compareArgs("THREAD MODEL", "OPENCL"))  {
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode: 'OpenCL', device_id: %d, platform_id: %d", device_id, plat);
  }else if(options.compareArgs("THREAD MODEL", "OPENMP"))  {
    sprintf(deviceConfig, "mode: 'OpenMP' ");
  }else  {
    sprintf(deviceConfig, "mode: 'Serial', memory: { use_host_pointer: true }");
    sprintf(deviceConfig, "mode: 'Serial'");
    options.setArgs("THREAD MODEL", "SERIAL");
  }

  if(rank == 0) printf("Initializing device\n");
  this->setup((std::string)deviceConfig);
  this->comm = comm;

  if (this->mode() == "Serial")
    options.setArgs("THREAD MODEL", "SERIAL");

  if(rank == 0)
    std::cout << "active occa mode: " << this->mode() << "\n\n";

#ifdef USE_OCCA_MEM_BYTE_ALIGN
  occa::env::OCCA_MEM_BYTE_ALIGN = USE_OCCA_MEM_BYTE_ALIGN;
#endif

  int Nthreads = 1;
  omp_set_num_threads(Nthreads);
}