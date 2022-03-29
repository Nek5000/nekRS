#include <cstdlib>
#include <strings.h>
#include "platform.hpp"
#include "nrs.hpp"
#include "linAlg.hpp"
#include "omp.h"
#include <iostream>
#include "flopCounter.hpp"

namespace{

static void compileDummyKernel(const platform_t& plat)
{
  const bool buildNodeLocal = useNodeLocalCache();
  auto rank = buildNodeLocal ? plat.comm.localRank : plat.comm.mpiRank;
  const std::string dummyKernelName = "myDummyKernelName";
  const std::string dummyKernelStr = std::string(
      "@kernel void myDummyKernelName(int N) {"
      "  for (int i = 0; i < N; ++i; @tile(64, @outer, @inner)) {}"
      "}"
  );

  if(rank == 0){
    plat.device.occaDevice().buildKernelFromString(
      dummyKernelStr,
      dummyKernelName,
      plat.kernelInfo
    );
  }

}

}

comm_t::comm_t(MPI_Comm _commg, MPI_Comm _comm)
{

  mpiCommParent = _commg;
  mpiComm = _comm;
  MPI_Comm_rank(_comm, &mpiRank);
  MPI_Comm_size(_comm, &mpiCommSize);

  MPI_Comm_split_type(_comm, MPI_COMM_TYPE_SHARED, mpiRank, MPI_INFO_NULL, &mpiCommLocal);
  MPI_Comm_rank(mpiCommLocal, &localRank);
  MPI_Comm_size(mpiCommLocal, &mpiCommLocalSize);

}

deviceVector_t::deviceVector_t(const size_t _offset, const size_t _nVectors, const size_t _wordSize, const std::string _vectorName)
: 
  nVectors(_nVectors),
  wordSize(_wordSize),
  vectorName(_vectorName),
  offset(_offset)
{
  if(offset <= 0 || nVectors <= 0 || wordSize <= 0) {
    if(platform->comm.mpiRank == 0)
      printf("ERROR: deviceVector_t invalid input!\n");
    ABORT(EXIT_FAILURE);
  }

  o_vector = platform->device.malloc(nVectors * offset * wordSize);
  for(int s = 0; s < nVectors; ++s){
    slices.push_back(o_vector + s * offset * wordSize);
  }
}

occa::memory&
deviceVector_t::at(const int i)
{
  if(i >= nVectors){
    if(platform->comm.mpiRank == 0){
      printf("ERROR: deviceVector_t(%s) has %zu size, but an attempt to access entry %i was made!\n",
        vectorName.c_str(),
        nVectors,
        i
      );
    }
    ABORT(EXIT_FAILURE);
    return o_vector;
  }
  return slices[i];
}


platform_t* platform_t::singleton = nullptr;
platform_t::platform_t(setupAide& _options, MPI_Comm _commg, MPI_Comm _comm)
: options(_options),
  warpSize(32),
  comm(_commg, _comm),
  device(options, comm),
  timer(_comm, device.occaDevice(), 0),
  kernels(*this)
{

  flopCounter = std::make_unique<flopCounter_t>();

  kernelInfo["defines/" "p_NVec"] = 3;
  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "dfloat"] = dfloatString;
  kernelInfo["defines/" "pfloat"] = pfloatString;
  kernelInfo["defines/" "dlong"] = dlongString;
  kernelInfo["defines/" "hlong"] = hlongString;

  if(device.mode() == "CUDA" && !getenv("OCCA_CUDA_COMPILER_FLAGS")) {
    kernelInfo["compiler_flags"] += " -O3 ";
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
    warpSize = 64; // can be arch specific
    kernelInfo["compiler_flags"] += " -O3 ";
    kernelInfo["compiler_flags"] += " -ffp-contract=fast ";
    kernelInfo["compiler_flags"] += " -funsafe-math-optimizations ";
    kernelInfo["compiler_flags"] += " -ffast-math ";
  }

  serial = device.mode() == "Serial" ||
           device.mode() == "OpenMP";
  
  const std::string extension = serial ? ".c" : ".okl";
  
  compileDummyKernel(*this);

  std::string installDir, kernelName, fileName;
  installDir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = installDir + "/okl/";
  kernelName = "copyDfloatToPfloat";
  fileName = installDir + "/okl/core/" + kernelName + extension;
  this->kernels.add(kernelName, fileName, this->kernelInfo);

  kernelName = "copyPfloatToDfloat";
  fileName = installDir + "/okl/core/" + kernelName + extension;
  this->kernels.add(kernelName, fileName, this->kernelInfo);
}
void memPool_t::allocate(const dlong offset, const dlong fields)
{
  ptr = (dfloat*) calloc(offset*fields, sizeof(dfloat));
  if(fields > 0) slice0 = ptr + 0 * offset;
  if(fields > 1) slice1 = ptr + 1 * offset;
  if(fields > 2) slice2 = ptr + 2 * offset;
  if(fields > 3) slice3 = ptr + 3 * offset;
  if(fields > 4) slice4 = ptr + 4 * offset;
  if(fields > 5) slice5 = ptr + 5 * offset;
  if(fields > 6) slice6 = ptr + 6 * offset;
  if(fields > 7) slice7 = ptr + 7 * offset;
  if(fields > 9) slice9 = ptr + 9 * offset;
  if(fields > 12) slice12 = ptr + 12 * offset;
  if(fields > 15) slice15 = ptr + 15 * offset;
  if(fields > 18) slice18 = ptr + 18 * offset;
  if(fields > 19) slice19 = ptr + 19 * offset;
}
void deviceMemPool_t::allocate(memPool_t& hostMemory, const dlong offset, const dlong fields)
{
  bytesAllocated = fields * offset * sizeof(dfloat);
  o_ptr = platform->device.malloc(offset*fields*sizeof(dfloat), hostMemory.slice0);
  if(fields > 0) slice0 = o_ptr.slice(0 * offset * sizeof(dfloat));
  if(fields > 1) slice1 = o_ptr.slice(1 * offset * sizeof(dfloat));
  if(fields > 2) slice2 = o_ptr.slice(2 * offset * sizeof(dfloat));
  if(fields > 3) slice3 = o_ptr.slice(3 * offset * sizeof(dfloat));
  if(fields > 4) slice4 = o_ptr.slice(4 * offset * sizeof(dfloat));
  if(fields > 5) slice5 = o_ptr.slice(5 * offset * sizeof(dfloat));
  if(fields > 6) slice6 = o_ptr.slice(6 * offset * sizeof(dfloat));
  if(fields > 7) slice7 = o_ptr.slice(7 * offset * sizeof(dfloat));
  if(fields > 9) slice9 = o_ptr.slice(9 * offset * sizeof(dfloat));
  if(fields > 12) slice12 = o_ptr.slice(12 * offset * sizeof(dfloat));
  if(fields > 15) slice15 = o_ptr.slice(15 * offset * sizeof(dfloat));
  if(fields > 18) slice18 = o_ptr.slice(18 * offset * sizeof(dfloat));
  if(fields > 19) slice19 = o_ptr.slice(19 * offset * sizeof(dfloat));
}
void
platform_t::create_mempool(const dlong offset, const dlong fields)
{
  mempool.allocate(offset, fields);
  o_mempool.allocate(mempool, offset, fields);
}
