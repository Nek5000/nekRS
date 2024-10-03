#include "benchmarkFDM.hpp"
#include <vector>
#include <numeric>
#include <iostream>

#include "randomVector.hpp"
#include "kernelBenchmarker.hpp"
#include <tuple>
#include <map>

namespace {
struct CallParameters {
  int Nelements;
  int Nq_e;
  size_t wordSize;
  bool useRAS;
  std::string suffix;
};
} // namespace

namespace std {
template <> struct less<CallParameters> {
  bool operator()(const CallParameters &lhs, const CallParameters &rhs) const
  {
    auto tier = [](const CallParameters &v) {
      return std::tie(v.Nelements, v.Nq_e, v.wordSize, v.useRAS, v.suffix);
    };
    return tier(lhs) < tier(rhs);
  }
};
} // namespace std

namespace {
std::map<CallParameters, occa::kernel> cachedResults;
}

template <typename T>
occa::kernel benchmarkFDM(int Nelements,
                          int Nq_e,
                          size_t wordSize,
                          bool useRAS,
                          int verbosity,
                          T NtestsOrTargetTime,
                          bool runAutotuner,
                          std::string suffix)
{
  if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
    Nelements = 1;
  }

  CallParameters params{Nelements, Nq_e, wordSize, useRAS, suffix};

  if (cachedResults.count(params) > 0) {
    return cachedResults.at(params);
  }

  const auto Nq = Nq_e - 2;
  const auto Np = std::pow(Nq, 3);
  const auto N_e = Nq_e - 1;
  const auto N = Nq - 1;
  const auto Np_e = Nq_e * Nq_e * Nq_e;

  occa::properties props = platform->kernelInfo + meshKernelProperties(N); // regular, non-extended mesh
  if (wordSize == 4)
    props["defines/pfloat"] = "float";
  else
    props["defines/pfloat"] = "dfloat";

  props["defines/p_Nq_e"] = Nq_e;
  props["defines/p_Np_e"] = Np_e;

  if (useRAS) {
    props["defines/p_restrict"] = 1;
  }
  else {
    props["defines/p_restrict"] = 0;
  }

  const std::string oklpath(getenv("NEKRS_KERNEL_DIR"));
  const std::string kernelName = "fusedFDM";
  const std::string fileName = oklpath + "/elliptic/" + kernelName;
  const std::string ext = platform->serial ? ".c" : ".okl";

  auto benchmarkFDMWithPrecision = [&](auto sampleWord) {
    using FPType = decltype(sampleWord);
    const auto wordSize = sizeof(FPType);

    std::vector<int> kernelVariants;
    if (platform->serial) {
     kernelVariants = {0};
    } else {
      kernelVariants = {0,1,2,3,4,11};
    }

    auto buildKernel = [&props, &fileName, &ext, &kernelName, &suffix](int ver)
    {
      auto newProps = props;
      newProps["defines/p_knl"] = ver;
      const auto verSuffix = "_v" + std::to_string(ver);

      if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
        const auto reqName = std::string(fileName) + "::" + std::string(newProps.hash().getString());
        platform->kernelRequests.add(reqName, fileName + ext, newProps, suffix);
        return occa::kernel();
      } else {
        return platform->device.loadKernel(fileName + ext, kernelName + verSuffix, newProps, suffix);
      }
    };

    auto referenceKernel = buildKernel(kernelVariants.front());

    if (!runAutotuner) {
      return std::make_pair(referenceKernel, -1.0);
    }

    auto Sx = randomVector<FPType>(Nelements * Nq_e * Nq_e, 0.1, 0.2, true);
    auto Sy = randomVector<FPType>(Nelements * Nq_e * Nq_e, 0.2, 0.3, true);
    auto Sz = randomVector<FPType>(Nelements * Nq_e * Nq_e, 0.3, 0.4, true);
    auto invL = randomVector<FPType>(Nelements * Np_e, 0, 1, true);
    auto Su = std::vector<FPType>(Nelements * ((useRAS) ? Np : Np_e), 0);
    auto u = randomVector<FPType>(Nelements * Np_e, 0, 1, true);
    auto invDegree = randomVector<FPType>(Nelements * Np_e, 0, 1, true);

    // elementList[e] = e
    std::vector<int> elementList(Nelements);
    std::iota(elementList.begin(), elementList.end(), 0);
    auto o_elementList = platform->device.malloc(Nelements * sizeof(int), elementList.data());

    auto o_Sx = platform->device.malloc(Sx.size() * wordSize, Sx.data());
    auto o_Sy = platform->device.malloc(Sy.size() * wordSize, Sy.data());
    auto o_Sz = platform->device.malloc(Sz.size() * wordSize, Sz.data());
    auto o_invL = platform->device.malloc(invL.size() * wordSize, invL.data());
    auto o_Su = platform->device.malloc(Su.size() * wordSize);
    o_Su.copyFrom(Su.data());
    auto o_u = platform->device.malloc(u.size() * wordSize);
    o_u.copyFrom(u.data());
    auto o_invDegree = platform->device.malloc(invDegree.size() * wordSize, invDegree.data());

    auto resetFields = [&]() {
      o_u.copyFrom(u.data()); // kernel reads and writes o_u 
    };

    auto kernelRunner = [&](occa::kernel &kernel) {
      // resetFields(); // disabling otherwise it would be included in the timer

      if (useRAS)
        kernel(Nelements, o_elementList, o_Su, o_Sx, o_Sy, o_Sz, o_invL, o_invDegree, o_u);
      else
        kernel(Nelements, o_elementList, o_Su, o_Sx, o_Sy, o_Sz, o_invL, o_u);
    };

    auto fdmKernelBuilder = [&](int kernelVariant) {
      auto kernel = buildKernel(kernelVariant);
      if (!kernel.isInitialized()) return occa::kernel();

      auto dumpResult = [&]() {
        std::vector<FPType> res;
        if (useRAS) {
          res.resize(Nelements * Np);
        }
        else {
          res.resize(Nelements * Np_e);
        }
        o_Su.copyTo(res.data());

        return res;
      };

      resetFields();
      kernelRunner(referenceKernel);
      auto referenceResult = dumpResult();

      resetFields();
      kernelRunner(kernel);
      auto result = dumpResult();

      const auto absTol = 1e-2;
      const auto err = maxRelErr<FPType>(referenceResult, result, platform->comm.mpiComm, absTol);
      const auto scale = 10 * range<FPType>(referenceResult, absTol);

      if (err > scale * std::numeric_limits<FPType>::epsilon() || std::isnan(err)) {
        if (platform->comm.mpiRank == 0 && verbosity > 1) {
          std::cout << "fdm: Ignore version " << kernelVariant
                    << " as correctness check failed with " << err << std::endl;
        }

        // pass un-initialized kernel to skip this kernel variant
        kernel = occa::kernel();
      }

      return kernel;
    };

    auto printPerformanceInfo = [&](int kernelVariant, double elapsed, int Ntests, bool skipPrint) {
      // print statistics
      const double GDOFPerSecond = (Nelements * (N_e * N_e * N_e) / elapsed) / 1.e9;

      size_t bytesPerElem = (3 * Np_e + 3 * Nq_e * Nq_e) * wordSize;
      const double bw = (Nelements * bytesPerElem / elapsed) / 1.e9;

      double flopsPerElem = 12 * Nq_e * Np_e + Np_e;
      const double gflops = (Nelements * flopsPerElem / elapsed) / 1.e9;

#ifdef _OPENMP
      const int Nthreads = omp_get_max_threads();
#else
      const int Nthreads = 1;
#endif

      if (platform->comm.mpiRank == 0 && !skipPrint) {
        if (verbosity > 0) {
          std::cout << "fdm:";
        }
        if (verbosity > 1) {
          std::cout << "MPItasks=" << platform->comm.mpiCommSize << " OMPthreads=" << Nthreads
                    << " NRepetitions=" << Ntests;
        }
        if (verbosity > 0) {
          std::cout << " N=" << N_e;

          if (verbosity > 1)
            std::cout << " Nelements=" << Nelements;

          if (verbosity > 1)
            std::cout << " elapsed time=" << elapsed;

          std::cout << " wordSize=" << 8 * wordSize << " GDOF/s=" << GDOFPerSecond << " GB/s=" << bw
                    << " GFLOPS/s=" << gflops << " kernelVer=" << kernelVariant << "\n";
        }
      }
    };

    auto printCallBack = [&](int kernelVariant, double elapsed, int Ntests) {
      printPerformanceInfo(kernelVariant, elapsed, Ntests, verbosity < 2);
    };

    auto kernelAndTime =
        benchmarkKernel(fdmKernelBuilder, kernelRunner, printCallBack, kernelVariants, NtestsOrTargetTime);

    if (kernelAndTime.first.properties().has("defines/p_knl") &&
        !platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      int bestKernelVariant = static_cast<int>(kernelAndTime.first.properties()["defines/p_knl"]);

      // print only the fastest kernel
      if (verbosity == 1) {
        printPerformanceInfo(bestKernelVariant, kernelAndTime.second, 0, false);
      }
    }
 
    return kernelAndTime;
  };

  occa::kernel kernel;

  if (wordSize == sizeof(float)) {
    float p = 0.0;
    auto kernelAndTime = benchmarkFDMWithPrecision(p);
    kernel = kernelAndTime.first;
  }
  else {
    double p = 0.0;
    auto kernelAndTime = benchmarkFDMWithPrecision(p);
    kernel = kernelAndTime.first;
  }

  cachedResults[params] = kernel;

  return kernel;
}

template occa::kernel benchmarkFDM<int>(int Nelements,
                                        int Nq_e,
                                        size_t wordSize,
                                        bool useRAS,
                                        int verbosity,
                                        int Ntests,
                                        bool runAutotuner,
                                        std::string suffix);

template occa::kernel benchmarkFDM<double>(int Nelements,
                                           int Nq_e,
                                           size_t wordSize,
                                           bool useRAS,
                                           int verbosity,
                                           double targetTime,
                                           bool runAutotuner,
                                           std::string suffix);
