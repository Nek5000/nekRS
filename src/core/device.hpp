#ifndef device_hpp_
#define device_hpp_
#include <string>
#include <occa.hpp>
#include <mpi.h>
#include "nrssys.hpp"

class setupAide;
class comm_t;

class device_t {
  public:
    device_t(setupAide& options, comm_t& comm);
    occa::memory
    malloc(size_t Nbytes, const void *src = nullptr, const occa::properties &properties = occa::properties());
    occa::memory malloc(size_t Nbytes, const occa::properties &properties);
    occa::memory malloc(size_t Nwords, size_t wordSize, occa::memory src);
    occa::memory malloc(size_t Nwords, size_t wordSize);

    occa::memory mallocHost(size_t Nbytes);

    int id() const { return _device_id; }
    const occa::device& occaDevice() const { return _device; }
    std::string mode() const { return _device.mode(); }
    occa::device& occaDevice() { return _device; }
    void finish() { _device.finish(); }

    occa::kernel buildKernel(const std::string &fullPath,
                             const occa::properties &props,
                             const std::string& suffix,
                             bool buildRank0) const;
    occa::kernel buildKernel(const std::string &fullPath,
                             const occa::properties &props,
                             bool buildRank0) const;

    // collective
    occa::kernel buildKernel(const std::string &fileName,
                             const std::string &kernelName,
                             const occa::properties &props) const;

    bool deviceAtomic;

  private:

    // non-collective
    occa::kernel buildKernel(const std::string &fullPath,
                             const occa::properties &props) const;
    occa::kernel buildKernel(const std::string &fullPath,
                             const occa::properties &props,
                             const std::string& suffix) const;
    occa::kernel buildKernel(const std::string &fileName,
                             const std::string &kernelName,
                             const occa::properties &props,
                             const std::string& suffix) const;

    occa::kernel buildNativeKernel(const std::string &fileName,
                             const std::string &kernelName,
                             const occa::properties &props) const;
    comm_t& _comm;
    occa::device _device;
    int _device_id;
};
#endif
