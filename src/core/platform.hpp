#ifndef platform_hpp_
#define platform_hpp_
#include <occa.hpp>
#include <mpi.h>
#include "nrssys.hpp"
#include "timer.hpp"
class setupAide;
class linAlg_t;
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
};
class device_t : public occa::device{
  public:
    device_t(setupAide& options, MPI_Comm comm);
    MPI_Comm comm;
    occa::kernel buildNativeKernel(const std::string &filename,
                             const std::string &kernelName,
                             const occa::properties &props) const;
    occa::kernel buildKernel(const std::string &filename,
                             const std::string &kernelName,
                             const occa::properties &props) const;
    occa::kernel buildKernel(const std::string &filename,
                             const std::string &kernelName,
                             const occa::properties &props,
                             MPI_Comm comm) const;
    occa::memory malloc(const dlong Nbytes, const void* src = nullptr, const occa::properties& properties = occa::properties());
    occa::memory malloc(const dlong Nbytes, const occa::properties& properties);
    occa::memory malloc(const dlong Nwords, const dlong wordSize, occa::memory src);
    occa::memory malloc(const dlong Nwords, const dlong wordSize);
  private:
    dlong bufferSize;
    void* _buffer;
};
struct comm_t{
  comm_t(MPI_Comm);
  MPI_Comm mpiComm;
  int mpiRank;
  int mpiCommSize;
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
  void create_mempool(const dlong offset, const dlong fields);
  platform_t(setupAide& _options, MPI_Comm _comm);

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