#ifndef platform_hpp_
#define platform_hpp_
#include <set>
#include "nekrsSys.hpp"
#include "QQt.hpp"
#include "flopCounter.hpp"
#include "timer.hpp"
#include "comm.hpp"
#include "par.hpp"
#include "device.hpp"
#include "device.tpp"
#include "solver.hpp"
#include "kernelRequestManager.hpp"

class setupAide;
class linAlg_t;
class flopCounter_t;

struct platform_t {
public:
  platform_t(setupAide &_options, MPI_Comm _commg, MPI_Comm _comm);
  void bcastJITKernelSourceFiles();

  static platform_t *getInstance(setupAide &_options, MPI_Comm _commg, MPI_Comm _comm)
  {
    if (!singleton) {
      singleton = new platform_t(_options, _commg, _comm);
    }
    return singleton;
  }

  static platform_t *getInstance()
  {
    return singleton;
  }

private:
  static platform_t *singleton;

public:
  setupAide &options;
  int warpSize;
  comm_t comm;
  device_t device;
  occa::properties kernelInfo;
  timer::timer_t timer;
  occa::memoryPool deviceMemoryPool;
  occa::memoryPool memoryPool;
  kernelRequestManager_t kernelRequests;
  Par *par;
  solver_t *solver;
  bool serial;
  linAlg_t *linAlg;
  std::unique_ptr<flopCounter_t> flopCounter;
  int exitValue;
  std::string tmpDir;
  int verbose;
  bool cacheLocal;
  bool cacheBcast;
  bool buildOnly;

  occa::kernel copyDfloatToPfloatKernel;
  occa::kernel copyPfloatToDfloatKernel;
  occa::kernel copyDfloatToDoubleKernel;
  occa::kernel copyDfloatToFloatKernel;
  occa::kernel copyDoubleToDfloatKernel;
  occa::kernel copyFloatToDfloatKernel;
};
#endif

#include "occaWrapper.hpp"
