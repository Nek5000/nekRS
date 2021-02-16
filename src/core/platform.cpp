#include "platform.hpp"
#include "nrs.hpp"
platform_t* platform_t::singleton = nullptr;
platform_t::platform_t(setupAide& options, MPI_Comm _comm)
: 
  device(occaDeviceConfig(options, _comm)),
  timer(_comm, device, 0),
  comm(_comm)
{
  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "p_BLOCKSIZE"] = BLOCKSIZE;
  if(sizeof(dfloat) == 4) {
    kernelInfo["defines/" "dfloat"] = "float";
    kernelInfo["defines/" "dfloat4"] = "float4";
    kernelInfo["defines/" "dfloat8"] = "float8";
  }
  if(sizeof(dfloat) == 8) {
    kernelInfo["defines/" "dfloat"] = "double";
    kernelInfo["defines/" "dfloat4"] = "double4";
    kernelInfo["defines/" "dfloat8"] = "double8";
  }

  if(sizeof(dlong) == 4)
    kernelInfo["defines/" "dlong"] = "int";
  if(sizeof(dlong) == 8)
    kernelInfo["defines/" "dlong"] = "long long int";
  if(sizeof(hlong) == 8)
    kernelInfo["defines/" "hlong"] = "long long int";

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

namespace occa{
occa::kernel
ParallelSafeDevice::buildKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props) const
{
  return this->buildKernel(filename, kernelName, props, comm);
}
occa::kernel
ParallelSafeDevice::buildKernel(const std::string &filename,
                         const std::string &kernelName,
                         const occa::properties &props,
                         MPI_Comm comm) const
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  occa::kernel _kernel;
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
      _kernel = device::buildKernel(filename, kernelName, props);
    }
    MPI_Barrier(comm);
  }
  return _kernel;
}
}