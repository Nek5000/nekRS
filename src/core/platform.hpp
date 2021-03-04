#ifndef platform_hpp_
#define platform_hpp_
#include <occa.hpp>
#include <mpi.h>
#include "nrssys.hpp"
#include "timer.hpp"
class setupAide;
class linAlg_t;
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
};
struct comm_t{
  comm_t(MPI_Comm);
  MPI_Comm mpiComm;
  int mpiRank;
  int mpiCommSize;
};
struct platform_t{
  int warpSize;
  device_t device;
  occa::properties kernelInfo;
  timer::timer_t timer;
  comm_t comm;
  linAlg_t* linAlg;
  dfloat* mempool;
  occa::memory o_mempool;
  occa::memory o_slice0, o_slice1, o_slice2, o_slice3, o_slice4, o_slice5, o_slice6, o_slice7;
  occa::memory o_slice9, o_slice12, o_slice15, o_slice18;
  void create_mempool(const dlong offset, const dlong fields);
  platform_t(setupAide& options, MPI_Comm _comm);

  static platform_t* getInstance(setupAide& options, MPI_Comm _comm){
    if(!singleton)
      singleton = new platform_t(options, _comm);
    return singleton;
  }
  static platform_t* getInstance(){
    return singleton;
  }
  private:
  static platform_t * singleton;
};
#endif