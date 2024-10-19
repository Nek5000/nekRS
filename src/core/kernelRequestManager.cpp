#include "nekrsSys.hpp"
#include "kernelRequestManager.hpp"
#include "platform.hpp"
#include "fileUtils.hpp"
#include <unordered_set>
#include <regex>
#include "sha1.hpp"
#include "threadPool.hpp"

kernelRequestManager_t::kernelRequestManager_t(const platform_t &m_platform)
    : kernelsProcessed(false), platformRef(m_platform)
{
}

// add (autotuned) kernel for subsequent load
void kernelRequestManager_t::add(const std::string &requestName, occa::kernel kernel)
{
  if (!kernel.isInitialized()) {
    return;
  }

  kernelRequest_t req(requestName, kernel.sourceFilename(), kernel.properties(), "");
  req.kernel = kernel;
  this->add(req, false);
}

void kernelRequestManager_t::add(const std::string &m_requestName,
                                 const std::string &m_fileName,
                                 const occa::properties &m_props,
                                 std::string m_suffix,
                                 bool checkUnique)
{
  this->add(kernelRequest_t{m_requestName, m_fileName, m_props, m_suffix}, checkUnique);
}

void kernelRequestManager_t::add(kernelRequest_t request, bool checkUnique)
{
  auto [iter, inserted] = requests.insert(request);

  // checkUnique flag is typically set to false because we may add the same request
  // (source file + properties) multiple times.
  if (checkUnique) {
    int unique = (inserted) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &unique, 1, MPI_INT, MPI_MIN, platformRef.comm.mpiComm);

    nekrsCheck(!unique,
               platformRef.comm.mpiComm,
               EXIT_FAILURE,
               "request already exists:\n%s",
               request.to_string().c_str());
  }

  // if the request already exists, it's important to verify that it is indeed the same,
  // as inadvertently overwriting the existing entry could occur otherwise.
  if (!inserted) {
    auto exisitingProps = (requestMap.find(request.requestName)->second).props;
    nekrsCheck(request.props.hash() != exisitingProps.hash(),
               platformRef.comm.mpiComm,
               EXIT_FAILURE,
               "detected different kernel hash for same request\n%s",
               request.to_string().c_str());

    auto exisitingFileName = (requestMap.find(request.requestName)->second).fileName;
    nekrsCheck(request.fileName != exisitingFileName,
               platformRef.comm.mpiComm,
               EXIT_FAILURE,
               "detected different kernel hash for same request\n%s",
               request.to_string().c_str());

    return;
  }

  requestMap.insert({request.requestName, request});
}

occa::kernel kernelRequestManager_t::load(const std::string &requestName, const std::string &_kernelName)
{
  auto errTxt = [&]() {
    const auto valid = processed() && (requestMap.find(requestName) != requestMap.end());

    if (valid) {
      return std::string();
    }

    std::stringstream txt;
    txt << "\n";
    txt << "Cannot find request " << "<" << requestName << ">" << "\n";
    txt << "Available:\n";
    for (auto &keyAndValue : requestMap) {
      txt << "\t" << "<" << keyAndValue.second.requestName << ">" << "\n";
    }

    txt << "===========================================================\n";
    auto retVal = txt.str();

    return retVal;
  }();

  nekrsCheck(errTxt.size(), platformRef.comm.mpiComm, EXIT_FAILURE, "%s\n", errTxt.c_str());

  auto kernel = [&]() {
    const auto &req = requestMap.find(requestName)->second;

    auto reqKnl = req.kernel;
    if (reqKnl.isInitialized()) {
      return reqKnl; // request is mapped to a already loaded kernel
    }

    const auto kernelName = [&]() {
      if (_kernelName.empty()) {
        auto fullPath = req.fileName;
        std::regex kernelNameRegex(R"((.+)\/(.+)\.)");
        std::smatch kernelNameMatch;
        const auto foundKernelName = std::regex_search(fullPath, kernelNameMatch, kernelNameRegex);

        // capture group
        // 0:   /path/to/install/nekrs/kernels/cds/advectMeshVelocityHex3D.okl
        // 1:   /path/to/install/nekrs/kernels/cds
        // 2:   advectMeshVelocityHex3D.okl

        return (foundKernelName && kernelNameMatch.size() == 3) ? kernelNameMatch[2].str() : "";
      } else {
        return _kernelName;
      }
    }();

    if (kernelMap.find({req, kernelName}) != kernelMap.end()) {
      return kernelMap[{req, kernelName}];
    } else {
      return kernelMap[{req, kernelName}] =
                 platformRef.device.loadKernel(req.fileName, kernelName, req.props, req.suffix);
    }
  }();

  nekrsCheck(!kernel.isInitialized(),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "kernel <%s> for request <%s> could not be initialized!\n",
             _kernelName.c_str(),
             requestName.c_str());

  return kernel;
}

void kernelRequestManager_t::compile()
{
  if (kernelsProcessed) {
    return;
  } else {
    kernelsProcessed = true;
  }

  const auto &device = platformRef.device;

  constexpr int maxCompilingRanks{
      32}; // large enough to speed things up, small enough to control pressure on filesystem
  const int rank = platform->cacheLocal ? platformRef.comm.localRank : platformRef.comm.mpiRank;
  const int ranksCompiling =
      std::min(maxCompilingRanks,
               platform->cacheLocal ? platformRef.comm.mpiCommLocalSize : platformRef.comm.mpiCommSize);

  auto Nthreads = 1;
  if (getenv("NEKRS_JITC_NTHREADS")) {
    Nthreads = std::stoi(getenv("NEKRS_JITC_NTHREADS"));
  }

  if (platformRef.comm.mpiRank == 0 && (platform->verbose || platform->buildOnly)) {
    std::cout << "requests.size(): " << requests.size() << std::endl;
    std::cout << "Nthreads: " << Nthreads << std::endl;
  }

  {
    std::map<std::string, kernelRequest_t> map;
    for (auto &&req : requests) {
      const auto fileName = (requestMap.find(req.requestName)->second).fileName;
      const auto props = (requestMap.find(req.requestName)->second).props;
      const auto hash = SHA1::from_string(fileName + props.hash().getFullString());
      auto [iter, inserted] = map.insert({hash, req});
      const std::string txt =
          "request collision between <" + req.requestName + "> and <" + (iter->second).requestName + ">!";
      nekrsCheck(!inserted, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", txt.c_str());
    }
  }

  // compile requests (assumed to have a unique occa hash) on build ranks
  constexpr int hashLength = 16 + 1; // null-terminated
  auto hashes = (char *)std::calloc(requests.size() * hashLength, sizeof(char));

  ThreadPool pool(Nthreads);

  if (rank < ranksCompiling) {
    for (auto &&req : requests) {
      auto retVal = pool.enqueue([&]() {
        const auto reqId = std::distance(requests.begin(), requests.find(req));
        if (reqId % ranksCompiling != rank) {
          return;
        }
        try {
          if (platform->verbose || platform->buildOnly) {
            std::cout << "Compiling request <" << req.requestName << ">";
          }

          auto knl = device.compileKernel(req.fileName, req.props, req.suffix, MPI_COMM_SELF);
          const auto hash = knl.hash().getString();
          nekrsCheck(hash.size() != hashLength - 1, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "Invalid hash!");

          std::strncpy(hashes + reqId * hashLength, hash.c_str(), hashLength);

          if (platform->verbose || platform->buildOnly) {
            std::cout << " (" << hash << ") on rank " << rank << std::flush << std::endl;
          }

        } catch (const std::exception &e) {
          std::cerr << "Caught exception: " << e.what() << std::endl;
        }
      });
    }
  }

  // finish compilation
  pool.finish();
  MPI_Barrier(platform->comm.mpiComm);

  // a-posteriori check for duplicated hash causing a potential race condition
  // no parallel version available yet
  if (platform->comm.mpiCommSize == 1) {
    const auto duplicateHashFound = [&]() {
      if (platform->comm.mpiRank == 0) {
        std::unordered_set<std::string> encounteredHashes;
        for (const auto &req : requests) {
          const auto reqId = distance(requests.begin(), requests.find(req));
          char hash[hashLength];
          std::strncpy(hash, hashes + reqId * hashLength, hashLength);
          if (!encounteredHashes.insert(hash).second) {
            std::cerr << "duplicate hash <" << hash << "> found for request: " << req.requestName
                      << std::endl;
            return true;
          }
        }
        return false;
      }
      return false;
    }();
    nekrsCheck(duplicateHashFound,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "More than one compile request is using the same hash!");
  }

  free(hashes);

  // after this point it is illegal to compile kernels
  platform->device.compilationFinished();

  if (platform->cacheBcast && !platform->buildOnly) {
    const auto srcPath = fs::path(getenv("OCCA_CACHE_DIR"));
    const std::string cacheDir = platform->tmpDir / fs::path("occa/");
    fileBcast(srcPath, fs::path(cacheDir) / "..", platform->comm.mpiComm, platform->verbose);

    // redirect
    occa::env::OCCA_CACHE_DIR = cacheDir;
    setenv("OCCA_CACHE_DIR", cacheDir.c_str(), 1);
  }
}
