#include <cstdlib>
#include <strings.h>
#include "platform.hpp"
#include "nrs.hpp"
#include "linAlg.hpp"
#include "omp.h"
#include <iostream>

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

deviceVector_t::deviceVector_t(const dlong _vectorSize, const dlong _nVectors, const dlong _wordSize, const std::string _vectorName)
: vectorSize(_vectorSize),
  nVectors(_nVectors),
  wordSize(_wordSize),
  vectorName(_vectorName)
{
  if(vectorSize <= 0 || nVectors <= 0 || wordSize <= 0) {
    if(platform->comm.mpiRank == 0)
      printf("ERROR: deviceVector_t invalid input!\n");
    ABORT(EXIT_FAILURE);
  }

  const size_t offset = vectorSize * (size_t)wordSize;
  o_vector = platform->device.malloc(nVectors * offset);
  for(int s = 0; s < nVectors; ++s){
    slices.push_back(o_vector + s * offset);
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
  return slices[i];
}


platform_t* platform_t::singleton = nullptr;
platform_t::platform_t(setupAide& _options, MPI_Comm _commg, MPI_Comm _comm)
: options(_options),
  warpSize(32),
  device(options, _commg, _comm),
  timer(_comm, device, 0),
  comm(_commg, _comm),
  kernels(*this)
{
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

void
kernelRequestManager_t::add_kernel(const std::string& m_requestName,
                const std::string& m_fileName,
                const std::string& m_kernelName,
                const occa::properties& m_props,
                std::string m_suffix,
                bool checkUnique)
{
  this->add_kernel(kernelRequest_t{m_requestName, m_fileName, m_kernelName, m_props, m_suffix}, checkUnique);
}
void
kernelRequestManager_t::add_kernel(kernelRequest_t request, bool checkUnique)
{
  auto iterAndBoolPair = kernels.insert(request);
  if(checkUnique)
  {
    int unique = (iterAndBoolPair.second) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &unique, 1, MPI_INT, MPI_MIN, platformRef.comm.mpiComm);
    if(!unique){
      if(platformRef.comm.mpiRank == 0)
      {
        std::cout << "Error in kernelRequestManager_t::add_kernel\n";
        std::cout << "Request details:\n";
        std::cout << request.to_string();
      }
      ABORT(1);
    }
  }

  const std::string fileName = request.fileName;
  fileNameToRequestMap[fileName].insert(request);
}
occa::kernel
kernelRequestManager_t::getKernel(const std::string& request, bool checkValid) const
{
  if(checkValid){
    bool issueError = 0;
    issueError = !processed();
    issueError = (requestToKernelMap.count(request) == 0);

    int errorFlag = issueError ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &errorFlag, 1, MPI_INT, MPI_MAX, platformRef.comm.mpiComm);

    if(errorFlag){
      if(platformRef.comm.mpiRank == 0)
      {
        std::cout << "\n";
        std::cout << "Error in kernelRequestManager_t::getKernel():\n";
        std::cout << "Cannot find requested kernel " << request << "!\n";

        std::cout << "Available:\n";
        for(auto&& keyAndValue : requestToKernelMap)
        {
          std::cout << "\t" << keyAndValue.first << "\n";
        }
        std::cout << "===========================================================\n";
      }
      ABORT(1);
    }
  }

  return requestToKernelMap.at(request);
}



void
kernelRequestManager_t::compile()
{

  if(kernelsProcessed) return;
  kernelsProcessed = true;

  constexpr int maxCompilingRanks {100};

  int buildNodeLocal = 0;
  if(getenv("NEKRS_CACHE_LOCAL"))
    buildNodeLocal = std::stoi(getenv("NEKRS_CACHE_LOCAL"));

  const int rank = buildNodeLocal ? platformRef.comm.localRank : platformRef.comm.mpiRank;
  const int ranksCompiling =
    std::min(
      maxCompilingRanks,
      buildNodeLocal ?
        platformRef.comm.mpiCommLocalSize :
        platformRef.comm.mpiCommSize
    );

  std::vector<std::string> kernelFiles(fileNameToRequestMap.size());
  
  unsigned ctr = 0;
  for(auto&& fileNameAndRequests : fileNameToRequestMap)
  {
    kernelFiles[ctr] = fileNameAndRequests.first;
    ctr++;
  }

  const auto& device = platformRef.device;
  auto& requestToKernel = requestToKernelMap;
  auto& fileNameToRequest = fileNameToRequestMap;
  auto compileKernels = [&kernelFiles, &requestToKernel, &fileNameToRequest, &device, rank, ranksCompiling](){
    if(rank >= ranksCompiling) return;
    const unsigned nFiles = kernelFiles.size();
    for(unsigned fileId = 0; fileId < nFiles; ++fileId)
    {
      if(fileId % ranksCompiling == rank){
        const std::string fileName = kernelFiles[fileId];
        for(auto && kernelRequest : fileNameToRequest[fileName]){
          const std::string requestName = kernelRequest.requestName;
          const std::string fileName = kernelRequest.fileName;
          const std::string kernelName = kernelRequest.kernelName;
          const std::string suffix = kernelRequest.suffix;
          const occa::properties props = kernelRequest.props;
          auto kernel = device.buildKernel(fileName, kernelName, props, suffix);
          requestToKernel[requestName] = kernel;
        }
      }
    }
  };

  const auto& kernelRequests = this->kernels;
  auto loadKernels = [&requestToKernel, &kernelRequests,&device](){
    for(auto&& kernelRequest : kernelRequests)
    {
      const std::string requestName = kernelRequest.requestName;
      if(requestToKernel.count(requestName) == 0){
        const std::string fileName = kernelRequest.fileName;
        const std::string kernelName = kernelRequest.kernelName;
        const std::string suffix = kernelRequest.suffix;
        const occa::properties props = kernelRequest.props;
        auto kernel = device.buildKernel(fileName, kernelName, props, suffix);
        requestToKernel[requestName] = kernel;
      }
    }
  };

  MPI_Barrier(platform->comm.mpiComm);
  compileKernels();
  MPI_Barrier(platform->comm.mpiComm);
  loadKernels();
}
