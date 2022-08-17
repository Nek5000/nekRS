#ifndef platform_hpp_
#define platform_hpp_
#include <occa.hpp>
#include <mpi.h>
#include "flopCounter.hpp"
#include "nrssys.hpp"
#include "timer.hpp"
#include "comm.hpp"
#include "inipp.hpp"
#include "device.hpp"
#include "kernelRequestManager.hpp"
#include <set>
#include <map>
#include <vector>
#include <memory>
class setupAide;
class linAlg_t;
class flopCounter_t;

class deviceVector_t{
public:
// allow implicit conversion between this and the underlying occa::memory object
  operator occa::memory&(){ return o_vector; }
// allow implicit conversion between this and kernelArg (for passing to kernels)
  operator occa::kernelArg(){ return o_vector; }
  deviceVector_t(const size_t _offset, const size_t _nVectors, const size_t _wordSize, std::string _vectorName = "");
  occa::memory& at(const int);
  const size_t offset;
private:
  occa::memory o_vector;
  std::vector<occa::memory> slices;
  const size_t nVectors;
  const size_t wordSize;
  const std::string vectorName;
};

struct memPool_t{
  void allocate(const dlong offset, const dlong fields);
  dfloat *slice0 = nullptr;
  dfloat *slice1 = nullptr; 
  dfloat *slice2 = nullptr; 
  dfloat *slice3 = nullptr; 
  dfloat *slice4 = nullptr; 
  dfloat *slice5 = nullptr; 
  dfloat *slice6 = nullptr; 
  dfloat *slice7 = nullptr;
  dfloat *slice9 = nullptr;
  dfloat *slice12 = nullptr; 
  dfloat *slice15 = nullptr; 
  dfloat *slice18 = nullptr; 
  dfloat *slice19 = nullptr;
  dfloat *ptr = nullptr;
};
struct deviceMemPool_t{
  void allocate(memPool_t& hostMemory, const dlong offset, const dlong fields);
  occa::memory slice0;
  occa::memory slice1; 
  occa::memory slice2; 
  occa::memory slice3; 
  occa::memory slice4; 
  occa::memory slice5; 
  occa::memory slice6; 
  occa::memory slice7;
  occa::memory slice9;
  occa::memory slice12; 
  occa::memory slice15; 
  occa::memory slice18; 
  occa::memory slice19;
  occa::memory o_ptr;
  size_t bytesAllocated;
};

struct platform_t{
public:
  void create_mempool(const dlong offset, const dlong fields);
  platform_t(setupAide& _options, MPI_Comm _commg, MPI_Comm _comm);

  static platform_t* getInstance(setupAide& _options, MPI_Comm _commg, MPI_Comm _comm){
    if(!singleton)
      singleton = new platform_t(_options, _commg, _comm);
    return singleton;
  }
  static platform_t* getInstance(){
    return singleton;
  }
private:
  static platform_t * singleton;

public:
  setupAide& options;
  int warpSize;
  comm_t comm;
  device_t device;
  occa::properties kernelInfo;
  timer::timer_t timer;
  deviceMemPool_t o_mempool;
  kernelRequestManager_t kernels;
  inipp::Ini *par;
  bool serial;
  linAlg_t* linAlg;
  std::unique_ptr<flopCounter_t> flopCounter;
  memPool_t mempool;

  occa::kernel copyDfloatToPfloatKernel;
  occa::kernel copyPfloatToDfloatKernel;
};
#endif
