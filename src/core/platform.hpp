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
#include "device.tpp"
#include "kernelRequestManager.hpp"
#include <set>
#include <map>
#include <vector>
#include <memory>
class setupAide;
class linAlg_t;
class flopCounter_t;

struct platform_t{
public:
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
  occa::experimental::memoryPool o_memPool;
  kernelRequestManager_t kernels;
  inipp::Ini *par;
  bool serial;
  linAlg_t* linAlg;
  std::unique_ptr<flopCounter_t> flopCounter;
  int exitValue;
  std::string tmpDir;
  int verbose;
  bool cacheLocal;
  bool cacheBcast; 

  occa::kernel copyDfloatToPfloatKernel;
  occa::kernel copyPfloatToDfloatKernel;
};
#endif
