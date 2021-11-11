#include <kernelRequestManager.hpp>
#include <platform.hpp>
kernelRequestManager_t::kernelRequestManager_t(const platform_t& m_platform)
: kernelsProcessed(false),
  platformRef(m_platform)
{}

void
kernelRequestManager_t::add(const std::string& m_requestName,
                const std::string& m_fileName,
                const occa::properties& m_props,
                std::string m_suffix,
                bool checkUnique)
{
  this->add(kernelRequest_t{m_requestName, m_fileName, m_props, m_suffix}, checkUnique);
}
void
kernelRequestManager_t::add(kernelRequest_t request, bool checkUnique)
{
  auto iterAndBoolPair = kernels.insert(request);
  if(checkUnique)
  {
    int unique = (iterAndBoolPair.second) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &unique, 1, MPI_INT, MPI_MIN, platformRef.comm.mpiComm);
    if(!unique){
      if(platformRef.comm.mpiRank == 0)
      {
        std::cout << "Error in kernelRequestManager_t::add\n";
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
kernelRequestManager_t::get(const std::string& request, bool checkValid) const
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

  const bool buildNodeLocal = useNodeLocalCache();

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
          const std::string suffix = kernelRequest.suffix;
          const occa::properties props = kernelRequest.props;

          // MPI staging already handled
          auto kernel = device.buildKernel(fileName, props, suffix, false);
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
        const std::string suffix = kernelRequest.suffix;
        const occa::properties props = kernelRequest.props;

        // MPI staging already handled
        auto kernel = device.buildKernel(fileName, props, suffix, false);
        requestToKernel[requestName] = kernel;
      }
    }
  };

  MPI_Barrier(platform->comm.mpiComm);
  compileKernels();
  MPI_Barrier(platform->comm.mpiComm);
  loadKernels();
}