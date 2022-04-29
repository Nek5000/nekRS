#include "benchmarkAx.hpp"
#include <vector>
#include <iostream>
#include <numeric>
#include "nrs.hpp"

#include "kernelBenchmarker.hpp"
#include "randomVector.hpp"
#include "omp.h"
#include <tuple>
#include <map>

namespace{
struct CallParameters{
  int Nelements;
  int Nq;
  int Ng;
  bool constCoeff;
  bool poisson;
  bool computeGeom;
  int wordSize;
  int Ndim;
  bool stressForm;
  std::string suffix;
};
}

namespace std
{
  template<> struct less<CallParameters>
  {
    bool operator() (const CallParameters& lhs, const CallParameters& rhs) const
    {
      auto tier = [](const CallParameters &v) {
        return std::tie(v.Nelements,
                        v.Nq,
                        v.Ng,
                        v.constCoeff,
                        v.poisson,
                        v.computeGeom,
                        v.wordSize,
                        v.Ndim,
                        v.stressForm,
                        v.suffix);
      };
      return tier(lhs) < tier(rhs);
    }
  };
}

namespace{
std::map<CallParameters, occa::kernel> cachedResults;
}

template <typename T>
occa::kernel benchmarkAx(int Nelements,
                         int Nq,
                         int Ng,
                         bool constCoeff,
                         bool poisson,
                         bool computeGeom,
                         int wordSize,
                         int Ndim,
                         bool stressForm,
                         int verbosity,
                         T NtestsOrTargetTime,
                         bool requiresBenchmark,
                         std::string suffix)
{
  CallParameters params{Nelements, Nq, Ng, constCoeff, poisson, computeGeom, wordSize, Ndim, stressForm, suffix};

  if(cachedResults.count(params) > 0){
    return cachedResults.at(params);
  }

  const auto N = Nq-1;

  const auto Np = Nq * Nq * Nq;
  const auto Nq_g = Ng + 1;
  const int Np_g = Nq_g * Nq_g * Nq_g;

  occa::properties props = platform->kernelInfo + meshKernelProperties(N);
  if (wordSize == 4)
    props["defines/dfloat"] = "float";
  if (Ng != N) {
    props["defines/p_Nq_g"] = Nq_g;
    props["defines/p_Np_g"] = Np_g;
  }
  if (poisson)
    props["defines/p_poisson"] = 1;

  std::string kernelName = "elliptic";
  if (Ndim > 1){
    kernelName += stressForm ? "Stress" : "Block";
  }
  kernelName += "PartialAx";
  if (!constCoeff)
    kernelName += "Coeff";
  if (Ng != N) {
    if (computeGeom) {
      if (Ng == 1) {
        kernelName += "Trilinear";
      }
      else {
        printf("Unsupported g-order=%d\n", Ng);
        exit(1);
      }
    }
    else {
      printf("for now g-order != p-order requires --computeGeom!\n");
      exit(1);
      kernelName += "Ngeom";
    }
  }
  kernelName += "Hex3D";
  if (Ndim > 1 && !stressForm)
    kernelName += "_N" + std::to_string(Ndim);

  auto benchmarkAxWithPrecision = [&](auto sampleWord) {
    using FPType = decltype(sampleWord);
    const auto wordSize = sizeof(FPType);
    constexpr int p_Nggeo{7};

    int Nkernels = 1;
    if (kernelName == "ellipticPartialAxHex3D")
      Nkernels = 7;
    std::vector<int> kernelVariants;
    if (platform->serial) {
      kernelVariants.push_back(0);
    }
    else {
      for (int knl = 0; knl < Nkernels; ++knl) {

#if 0
        // v3 requires Nq^3 < 1024 (max threads/thread block on CUDA/HIP)
        if (knl == 3 && Np > 1024)
          continue;
#else
        // disable v3 for now, since correctness check is off
        if (knl == 3)
          continue;
#endif
        kernelVariants.push_back(knl);
      }
    }

    const std::string installDir(getenv("NEKRS_HOME"));

    // only a single choice, no need to run benchmark
    if (kernelVariants.size() == 1 && !requiresBenchmark) {

      auto newProps = props;
      if (kernelName == "ellipticPartialAxHex3D" && !platform->serial) {
        newProps["defines/p_knl"] = kernelVariants.back();
      }

      const std::string ext = platform->serial ? ".c" : ".okl";
      const std::string fileName = installDir + "/okl/elliptic/" + kernelName + ext;

      return std::make_pair(platform->device.buildKernel(fileName, newProps, suffix, true), -1.0);
    }

    auto DrV = randomVector<FPType>(Nq * Nq);
    auto ggeo = randomVector<FPType>(Np_g * Nelements * p_Nggeo);
    auto q = randomVector<FPType>((Ndim * Np) * Nelements);
    auto Aq = randomVector<FPType>((Ndim * Np) * Nelements);
    auto exyz = randomVector<FPType>((3 * Np_g) * Nelements);
    auto gllwz = randomVector<FPType>(2 * Nq_g);
    auto lambda = randomVector<FPType>(2 * Np * Nelements);

    // elementList[e] = e
    std::vector<dlong> elementList(Nelements);
    std::iota(elementList.begin(), elementList.end(), 0);
    auto o_elementList = platform->device.malloc(Nelements * sizeof(dlong), elementList.data());

    auto o_D = platform->device.malloc(Nq * Nq * wordSize, DrV.data());
    auto o_S = o_D;
    auto o_ggeo = platform->device.malloc(Np_g * Nelements * p_Nggeo * wordSize, ggeo.data());
    auto o_q = platform->device.malloc((Ndim * Np) * Nelements * wordSize, q.data());
    auto o_Aq = platform->device.malloc((Ndim * Np) * Nelements * wordSize, Aq.data());
    auto o_exyz = platform->device.malloc((3 * Np_g) * Nelements * wordSize, exyz.data());
    auto o_gllwz = platform->device.malloc(2 * Nq_g * wordSize, gllwz.data());

    auto o_lambda = platform->device.malloc(2 * Np * Nelements * wordSize, lambda.data());

    occa::kernel referenceKernel;
    {
      auto newProps = props;
      if (!platform->serial)
        newProps["defines/p_knl"] = kernelVariants.front();

      const std::string ext = platform->serial ? ".c" : ".okl";
      const std::string fileName = installDir + "/okl/elliptic/" + kernelName + ext;

      referenceKernel = platform->device.buildKernel(fileName, newProps, suffix, true);
    }

    auto kernelRunner = [&](occa::kernel &kernel) {
      const int loffset = 0;
      const int offset = Nelements * Np;
      if (computeGeom) {
        kernel(Nelements, offset, loffset, o_elementList, o_exyz, o_gllwz, o_D, o_S, o_lambda, o_q, o_Aq);
      }
      else {
        kernel(Nelements, offset, loffset, o_elementList, o_ggeo, o_D, o_S, o_lambda, o_q, o_Aq);
      }
    };

    auto axKernelBuilder = [&](int kernelVariant) {
      auto newProps = props;
      if (!platform->serial)
        newProps["defines/p_knl"] = kernelVariant;

      const std::string ext = platform->serial ? ".c" : ".okl";
      const std::string fileName = installDir + "/okl/elliptic/" + kernelName + ext;

      auto kernel = platform->device.buildKernel(fileName, newProps, suffix, true);

      if(platform->options.compareArgs("BUILD ONLY", "TRUE")) return kernel;

      std::vector<FPType> refResults((Ndim * Np) * Nelements);
      std::vector<FPType> results((Ndim * Np) * Nelements);

      kernelRunner(referenceKernel);
      o_Aq.copyTo(refResults.data(), refResults.size() * sizeof(FPType));

      kernelRunner(kernel);
      o_Aq.copyTo(results.data(), results.size() * sizeof(FPType));

      FPType err = 0.0;
      for (int i = 0; i < refResults.size(); ++i) {
        err = std::max(err, std::abs(refResults[i] - results[i]));
      }

      if (platform->comm.mpiRank == 0 && verbosity > 1) {
        std::cout << "Error in kernel compared to reference implementation " << kernelVariant << ": " << err
                  << std::endl;
      }

      return kernel;
    };

    auto printPerformanceInfo = [&](int kernelVariant, double elapsed, int Ntests, bool skipPrint) {
      const bool BKmode = constCoeff && poisson;

      // print statistics
      const dfloat GDOFPerSecond = (Nelements * Ndim * (N * N * N) / elapsed) / 1.e9;

      size_t bytesMoved = Ndim * 2 * Np * wordSize; // x, Ax
      bytesMoved += 6 * Np_g * wordSize;            // geo
      if (!constCoeff)
        bytesMoved += 3 * Np * wordSize; // lambda1, lambda2, Jw
      const double bw = (Nelements * bytesMoved / elapsed) / 1.e9;

      double flopCount = Np * 12 * Nq + 15 * Np;
      if (!constCoeff)
        flopCount += 5 * Np;
      const double gflops = Ndim * (flopCount * Nelements / elapsed) / 1.e9;
      const int Nthreads = omp_get_max_threads();

      if (platform->comm.mpiRank == 0 && !skipPrint) {
        if (verbosity > 0) {
          std::cout << "Ax:";
        }
        if (verbosity > 1) {
          std::cout << " MPItasks=" << platform->comm.mpiCommSize << " OMPthreads=" << Nthreads
                    << " NRepetitions=" << Ntests;
        }
        if (verbosity > 0) {
          if (Ndim > 1)
            std::cout << " Ndim=" << Ndim;

          std::cout << " N=" << N;
          if (Ng != N)
            std::cout << " Ng=" << Ng;

          if (verbosity > 1)
            std::cout << " Nelements=" << Nelements;

          if (verbosity > 1)
            std::cout << " elapsed time=" << elapsed;

          std::cout << " wordSize=" << 8 * wordSize << " GDOF/s=" << GDOFPerSecond << " GB/s=" << bw
                    << " GFLOPS/s=" << gflops << " bkMode=" << BKmode << " kernelVer=" << kernelVariant
                    << "\n";
        }
      }
    };

    auto printCallBack = [&](int kernelVariant, double elapsed, int Ntests) {
      printPerformanceInfo(kernelVariant, elapsed, Ntests, verbosity < 2);
    };

    auto kernelAndTime =
        benchmarkKernel(axKernelBuilder, kernelRunner, printCallBack, kernelVariants, NtestsOrTargetTime);

    if (kernelAndTime.first.properties().has("defines/p_knl") && platform->options.compareArgs("BUILD ONLY","FALSE")) {
      int bestKernelVariant = static_cast<int>(kernelAndTime.first.properties()["defines/p_knl"]);

      // print only the fastest kernel
      if (verbosity == 1) {
        printPerformanceInfo(bestKernelVariant, kernelAndTime.second, 0, false);
      }
    }

    free(o_D);
    free(o_S);
    free(o_ggeo);
    free(o_q);
    free(o_Aq);
    free(o_exyz);
    free(o_gllwz);
    free(o_lambda);
    free(o_elementList);

    return kernelAndTime;
  };

  occa::kernel kernel;

  if (wordSize == sizeof(float)) {
    float p = 0.0;
    auto kernelAndTime = benchmarkAxWithPrecision(p);
    kernel = kernelAndTime.first;
  }
  else {
    double p = 0.0;
    auto kernelAndTime = benchmarkAxWithPrecision(p);
    kernel = kernelAndTime.first;
  }
  
  cachedResults[params] = kernel;

  return kernel;
}

template occa::kernel benchmarkAx<int>(int Nelements,
                                       int Nq,
                                       int Ng,
                                       bool constCoeff,
                                       bool poisson,
                                       bool computeGeom,
                                       int wordSize,
                                       int Ndim,
                                       bool stressForm,
                                       int verbosity,
                                       int Ntests,
                                       bool requiresBenchmark,
                                       std::string suffix);

template occa::kernel benchmarkAx<double>(int Nelements,
                                          int Nq,
                                          int Ng,
                                          bool constCoeff,
                                          bool poisson,
                                          bool computeGeom,
                                          int wordSize,
                                          int Ndim,
                                          bool stressForm,
                                          int verbosity,
                                          double targetTime,
                                          bool requiresBenchmark,
                                          std::string suffix);
