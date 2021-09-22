#ifndef device_hpp_
#define device_hpp_
#include <string>
#include <occa.hpp>
#include <mpi.h>
#include "nrssys.hpp"

class setupAide;

class device_t : public occa::device{
  public:
    device_t(setupAide& options, MPI_Comm comm);
    MPI_Comm comm;
    occa::memory malloc(const hlong Nbytes, const void* src = nullptr, const occa::properties& properties = occa::properties());
    occa::memory malloc(const hlong Nbytes, const occa::properties& properties);
    occa::memory malloc(const hlong Nwords, const dlong wordSize, occa::memory src);
    occa::memory malloc(const hlong Nwords, const dlong wordSize);

    occa::memory mallocHost(const hlong Nbytes);

    int id() const { return _device_id; }
    occa::kernel buildNativeKernel(const std::string &filename,
                             const std::string &kernelName,
                             const occa::properties &props) const;
    occa::kernel buildKernel(const std::string &filename,
                             const std::string &kernelName,
                             const occa::properties &props,
                             std::string suffix = std::string()) const;
  private:
    int _device_id;
};
#endif
