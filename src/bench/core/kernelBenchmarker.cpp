#include "kernelBenchmarker.hpp"
#include <limits>
#include "nrs.hpp"

namespace {
double run(int Nsamples, std::function<void(occa::kernel &)> kernelRunner, occa::kernel &kernel)
{
  platform->device.finish();
  MPI_Barrier(platform->comm.mpiComm);
  const double start = MPI_Wtime();

  for (int test = 0; test < Nsamples; ++test) {
    kernelRunner(kernel);
  }

  platform->device.finish();
  return (MPI_Wtime() - start) / Nsamples;
}
} // namespace
std::pair<occa::kernel, double>
benchmarkKernel(std::function<occa::kernel(int kernelVariant)> kernelBuilder,
                std::function<void(occa::kernel &)> kernelRunner,
                std::function<void(int kernelVariant, double tKernel, int Ntests)> printCallback,
                const std::vector<int> &kernelVariants,
                int Ntests)
{
  occa::kernel fastestKernel;
  double fastestTime = std::numeric_limits<double>::max();
  for (auto &&kernelVariant : kernelVariants) {

    auto candidateKernel = kernelBuilder(kernelVariant);

    if(platform->options.compareArgs("BUILD ONLY", "FALSE")){
      // warmup
      double elapsed = run(10, kernelRunner, candidateKernel);

      double candidateKernelTiming = run(Ntests, kernelRunner, candidateKernel);
      double tMax;
      double tMin;
      MPI_Allreduce(&candidateKernelTiming, &tMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
      MPI_Allreduce(&candidateKernelTiming, &tMin, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);

      const double tRatio = tMax/tMin;
      if (platform->comm.mpiRank == 0 && tRatio > 1.1)
        printf("WARNING: kernel timings differ by up to %.2 across ranks!\n", tRatio);

      candidateKernelTiming = tMax;

      if (candidateKernelTiming < fastestTime) {
        fastestTime = candidateKernelTiming;
        fastestKernel = candidateKernel;
      }

      printCallback(kernelVariant, candidateKernelTiming, Ntests);
    } else {
      fastestKernel = candidateKernel;
    }
  }

  return std::make_pair(fastestKernel, fastestTime);
}

std::pair<occa::kernel, double>
benchmarkKernel(std::function<occa::kernel(int kernelVariant)> kernelBuilder,
                std::function<void(occa::kernel &)> kernelRunner,
                std::function<void(int kernelVariant, double tKernel, int Ntests)> printCallback,
                const std::vector<int> &kernelVariants,
                double targetTime)
{
  occa::kernel fastestKernel;
  double fastestTime = std::numeric_limits<double>::max();
  for (auto &&kernelVariant : kernelVariants) {

    auto candidateKernel = kernelBuilder(kernelVariant);
    if(platform->options.compareArgs("BUILD ONLY", "FALSE")){

      // warmup
      double elapsed = run(10, kernelRunner, candidateKernel);

      // evaluation
      int Ntests = static_cast<int>(targetTime / elapsed);
      MPI_Allreduce(MPI_IN_PLACE, &Ntests, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);

      double candidateKernelTiming = run(Ntests, kernelRunner, candidateKernel);
      double tMax;
      double tMin;
      MPI_Allreduce(&candidateKernelTiming, &tMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
      MPI_Allreduce(&candidateKernelTiming, &tMin, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);

      const double tRatio = tMax/tMin;
      if (platform->comm.mpiRank == 0 && tRatio > 1.1)
        printf("WARNING: kernel timings differ by up to %.2 across ranks!\n", tRatio);

      candidateKernelTiming = tMax;

      if (candidateKernelTiming < fastestTime) {
        fastestTime = candidateKernelTiming;
        fastestKernel = candidateKernel;
      }

      printCallback(kernelVariant, candidateKernelTiming, Ntests);
    } else {
      fastestKernel = candidateKernel;
    }
  }

  return std::make_pair(fastestKernel, fastestTime);
}
