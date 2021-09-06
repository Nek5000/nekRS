#ifndef platform_hpp_
#define platform_hpp_
#include <occa.hpp>
#include <vector>
#include <mpi.h>
#include "nrssys.hpp"
#include "timer.hpp"
#include "inipp.hpp"
#include <set>
#include <map>
#include <vector>
class setupAide;
class linAlg_t;
class kernelRequestManager_t;

class deviceVector_t{
public:
// allow implicit conversion between this and the underlying occa::memory object
  operator occa::memory&(){ return o_vector; }
// allow implicit conversion between this and kernelArg (for passing to kernels)
  operator occa::kernelArg(){ return o_vector; }
  deviceVector_t(const dlong _vectorSize, const dlong _nVectors, const dlong _wordSize, std::string _vectorName = "");
  occa::memory& at(const int);
private:
  occa::memory o_vector;
  std::vector<occa::memory> slices;
  const dlong vectorSize;
  const dlong nVectors;
  const dlong wordSize;
  const std::string vectorName;
};

struct memPool_t{
  void allocate(const dlong offset, const dlong fields);
  dfloat* slice0, *slice1, *slice2, *slice3, *slice4, *slice5, *slice6, *slice7;
  dfloat* slice9, *slice12, *slice15, *slice18, *slice19;
  dfloat* ptr;
};
struct deviceMemPool_t{

  void allocate(memPool_t& hostMemory, const dlong offset, const dlong fields);
  occa::memory slice0, slice1, slice2, slice3, slice4, slice5, slice6, slice7;
  occa::memory slice9, slice12, slice15, slice18, slice19;
  occa::memory o_ptr;
  long long bytesAllocated;
};

class kernelRequestManager_t
{

  struct kernelRequest_t
  {
    inline bool operator==(const kernelRequest_t& other) const
    {
      return requestName == other.requestName;
    }
    inline bool operator<(const kernelRequest_t& other) const
    {
      return requestName < other.requestName;
    }
    inline bool operator> (const kernelRequest_t& other) const { return *this < other; }
    inline bool operator<=(const kernelRequest_t& other) const { return !(*this > other); }
    inline bool operator>=(const kernelRequest_t& other) const { return !(*this < other); }
    inline bool operator!=(const kernelRequest_t& other) const { return !(*this == other); }

    kernelRequest_t(const std::string& m_requestName,
                    const std::string& m_fileName,
                    const std::string& m_kernelName,
                    const occa::properties& m_props,
                    std::string m_suffix = std::string())
    :
    requestName(m_requestName),
    fileName(m_fileName),
    kernelName(m_kernelName),
    suffix(m_suffix),
    props(m_props)
    {}
    std::string requestName;
    std::string fileName;
    std::string kernelName;
    std::string suffix;
    occa::properties props;

    std::string to_string() const {
      std::ostringstream ss;
      ss << "requestName : " << requestName << "\n";
      ss << "fileName : " << fileName << "\n";
      ss << "kernelName : " << kernelName << "\n";
      ss << "suffix : " << suffix << "\n";
      ss << "props : " << props << "\n";;
      return ss.str();
    }
  };
public:
  kernelRequestManager_t(const platform_t& m_platform)
  : kernelsProcessed(false),
    platformRef(m_platform)
  {}
  void add_kernel(const std::string& m_requestName,
                  const std::string& m_fileName,
                  const std::string& m_kernelName,
                  const occa::properties& m_props,
                  std::string m_suffix = std::string(),
                  bool assertUnique = false);
  
  void compile();

  occa::kernel
  get(const std::string& request, bool checkValid = true) const;

  bool
  processed() const { return kernelsProcessed; }

private:
  const platform_t& platformRef;
  bool kernelsProcessed;
  std::set<kernelRequest_t> kernels;
  std::map<std::string, occa::kernel> requestToKernelMap;

  void add_kernel(kernelRequest_t request, bool assertUnique = true);

};

class device_t : public occa::device{
  public:
    device_t(setupAide& options, MPI_Comm comm);
    MPI_Comm comm;
    occa::memory malloc(const dlong Nbytes, const void* src = nullptr, const occa::properties& properties = occa::properties());
    occa::memory malloc(const dlong Nbytes, const occa::properties& properties);
    occa::memory malloc(const dlong Nwords, const dlong wordSize, occa::memory src);
    occa::memory malloc(const dlong Nwords, const dlong wordSize);

    occa::memory mallocHost(const dlong Nbytes);

    int id() const { return _device_id; }
    occa::kernel buildNativeKernel(const std::string &filename,
                             const std::string &kernelName,
                             const occa::properties &props) const;
    occa::kernel buildKernel(const std::string &filename,
                             const std::string &kernelName,
                             const occa::properties &props,
                             std::string suffix = std::string()) const;
  private:
    dlong bufferSize;
    int _device_id;
    void* _buffer;
};
struct comm_t{
  comm_t(MPI_Comm);
  MPI_Comm mpiComm;
  int mpiRank;
  int mpiCommSize;

  MPI_Comm localComm;
  int localCommSize;
  int localRank;
};
struct platform_t{
  setupAide& options;
  int warpSize;
  device_t device;
  occa::properties kernelInfo;
  timer::timer_t timer;
  comm_t comm;
  linAlg_t* linAlg;
  memPool_t mempool;
  deviceMemPool_t o_mempool;
  kernelRequestManager_t kernels;
  void create_mempool(const dlong offset, const dlong fields);
  platform_t(setupAide& _options, MPI_Comm _comm);
  inipp::Ini *par;

  static platform_t* getInstance(setupAide& _options, MPI_Comm _comm){
    if(!singleton)
      singleton = new platform_t(_options, _comm);
    return singleton;
  }
  static platform_t* getInstance(){
    return singleton;
  }
  private:
  static platform_t * singleton;
};
#endif
