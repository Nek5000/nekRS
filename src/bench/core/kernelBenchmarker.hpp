#include "nrs.hpp"
#include <utility>
#include <functional>

std::pair<occa::kernel, double>
benchmarkKernel(std::function<occa::kernel(int kernelVariant)> kernelBuilder,
                std::function<void(occa::kernel &)> kernelRunner,
                std::function<void(int kernelVariant, double tKernel, int Ntests)> printCallback,
                const std::vector<int> &kernelVariants,
                int Ntests);

std::pair<occa::kernel, double>
benchmarkKernel(std::function<occa::kernel(int kernelVariant)> kernelBuilder,
                std::function<void(occa::kernel &)> kernelRunner,
                std::function<void(int kernelVariant, double tKernel, int Ntests)> printCallback,
                const std::vector<int> &kernelVariants,
                double targetTime);

template <typename T>
T maxRelErr(const std::vector<T>& uRef, const std::vector<T>& u, MPI_Comm comm, T absTol = 0)
{
  double err = 0;

  if (absTol > 0) {
    for (int i = 0; i < uRef.size(); ++i) {
      if (std::abs(uRef[i]) > absTol) {
        err = std::max(err, (double) std::abs((uRef[i] - u[i])/uRef[i]));
      }
    }
  } else {
    for (int i = 0; i < uRef.size(); ++i) {
      if(std::abs(uRef[i] - u[i]) > 500*std::numeric_limits<T>::epsilon()) {
        err = std::max(err, (double) std::abs((uRef[i] - u[i])/uRef[i]));
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, comm);
  return static_cast<T>(err);
}
