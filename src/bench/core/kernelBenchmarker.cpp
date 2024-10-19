#include "kernelBenchmarker.hpp"
#include <limits>

namespace
{
constexpr int Nbaseline{100};
constexpr int Nwarmup{500};

double run(int Nsamples, std::function<void(occa::kernel &)> kernelRunner, occa::kernel &kernel)
{
  platform->device.finish();
  MPI_Barrier(platform->comm.mpiComm);
  const double start = MPI_Wtime();

  for (int test = 0; test < Nsamples; ++test) {
    kernelRunner(kernel);
  }

  // reset in case something weired happened because of an invalid kernel
  std::feclearexcept(FE_ALL_EXCEPT);

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

    MPI_Barrier(platform->comm.mpiComm);
    auto candidateKernel = kernelBuilder(kernelVariant);

    if (!candidateKernel.isInitialized()) {
      continue;
    }

    if (platform->options.compareArgs("BUILD ONLY", "FALSE")) {
      // warmup
      run(Nwarmup, kernelRunner, candidateKernel);
      double candidateKernelTiming = run(Ntests, kernelRunner, candidateKernel);
      
      double tMax;
      double tMin;
      MPI_Allreduce(&candidateKernelTiming, &tMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
      MPI_Allreduce(&candidateKernelTiming, &tMin, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);

      const double tRatio = tMax / tMin;
      if (platform->comm.mpiRank == 0 && tRatio > 1.1) {
        printf("WARNING: kernel[%d] timings differ by up to %.2f across ranks!\n", kernelVariant, tRatio);
      }

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

  const auto check = !platform->options.compareArgs("REGISTER ONLY", "TRUE") && !platform->options.compareArgs("BUILD ONLY", "TRUE");

  for (auto &&kernelVariant : kernelVariants) {

    MPI_Barrier(platform->comm.mpiComm);
    auto candidateKernel = kernelBuilder(kernelVariant);

    if (!candidateKernel.isInitialized()) {
      continue; // no need to proceed
    }

    if (check) {
      double elapsed = run(Nbaseline, kernelRunner, candidateKernel);
      int Ntests = std::max(1, static_cast<int>(targetTime / elapsed));
      MPI_Allreduce(MPI_IN_PLACE, &Ntests, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);

      double candidateKernelTiming = run(Ntests, kernelRunner, candidateKernel);

      double tMax;
      double tMin;
      MPI_Allreduce(&candidateKernelTiming, &tMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
      MPI_Allreduce(&candidateKernelTiming, &tMin, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);

      const double tRatio = tMax / tMin;
      if (platform->comm.mpiRank == 0 && tRatio > 1.1) {
        printf("WARNING: kernel[%d] timings differ by up to %.2f across ranks!\n", kernelVariant, tRatio);
      }

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

  nekrsCheck(!fastestKernel.isInitialized() && check,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "Cannot find valid kernel variant!");

  return std::make_pair(fastestKernel, fastestTime);
}
