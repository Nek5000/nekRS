#include "benchmarkAx.hpp"
#include <vector>
#include <iostream>
#include <numeric>
#include "nrs.hpp"

#include "kernelBenchmarker.hpp"
#include "randomVector.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif
#include <tuple>
#include <map>

namespace {
struct CallParameters {
  int Nelements;
  int Nq;
  int Ng;
  bool constCoeff;
  bool poisson;
  bool computeGeom;
  size_t wordSize;
  int Ndim;
  bool stressForm;
  std::string suffix;
};
} // namespace

namespace std {
template <> struct less<CallParameters> {
  bool operator()(const CallParameters &lhs, const CallParameters &rhs) const
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
} // namespace std

namespace {
std::map<CallParameters, occa::kernel> cachedResults;
}

template <typename T>
occa::kernel benchmarkAx(int Nelements,
                         int Nq,
                         int Ng,
                         bool constCoeff,
                         bool poisson,
                         bool computeGeom,
                         size_t wordSize,
                         int Ndim,
                         bool stressForm,
                         int verbosity,
                         T NtestsOrTargetTime,
                         bool requiresBenchmark,
                         std::string suffix)
{
  if (platform->options.compareArgs("BUILD ONLY", "TRUE")) {
    Nelements = 1;
  }

  CallParameters
      params{Nelements, Nq, Ng, constCoeff, poisson, computeGeom, wordSize, Ndim, stressForm, suffix};

  if (cachedResults.count(params) > 0) {
    return cachedResults.at(params);
  }

  const auto N = Nq - 1;

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

  if (constCoeff)
    props["defines/p_lambda"] = 0;
  else
    props["defines/p_lambda"] = 1;

  std::string kernelName = "elliptic";
  if (Ndim > 1) {
    kernelName += stressForm ? "Stress" : "Block";
  }
  kernelName += "PartialAx";
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

  auto benchmarkAxWithPrecision = [&](auto sampleWord) {
    using FPType = decltype(sampleWord);
    const auto wordSize = sizeof(FPType);
    constexpr int p_Nggeo{7}, p_Nvgeo{12};

    std::vector<int> kernelVariants;

    if (platform->serial) {
      const int Nkernels = 1;
      for (int knl = 0; knl < Nkernels; ++knl)
        kernelVariants.push_back(knl);
    }
    else {
      if (kernelName == "ellipticPartialAxCoeffHex3D") {
        const int Nkernels = 8;
        for (int knl = 0; knl < Nkernels; ++knl)
          kernelVariants.push_back(knl);

        kernelVariants.erase(kernelVariants.begin() + 3); // correctness check is off
      }
      if (kernelName == "ellipticStressPartialAxCoeffHex3D") {
        const int Nkernels = 2;
        for (int knl = 0; knl < Nkernels; ++knl)
          kernelVariants.push_back(knl);
        int n_plane = 1;
        switch (Nq) {
        case 4: n_plane = 2; break;
        case 5: n_plane = 1; break;
        case 6: n_plane = 2; break;
        case 7: n_plane = 1; break;
        case 8: n_plane = 1; break;
        case 9: n_plane = 3; break;
        case 10: n_plane = 1; break;
        case 11: n_plane = 1; break;
        }
        props["defines/n_plane"] = n_plane;
        props["defines/pts_per_thread"] = Nq/n_plane;              
      }
      if (kernelName == "ellipticBlockPartialAxCoeffHex3D") {
        const int Nkernels = 2;
        for (int knl = 0; knl < Nkernels; ++knl)
          kernelVariants.push_back(knl);

        int n_plane = 1;
        switch (Nq) {
        case 4:
          n_plane = 2;
          break;
        case 5:
          n_plane = 1;
          break;
        case 6:
          n_plane = 3;
          break;
        case 7:
          n_plane = 1;
          break;
        case 8:
          n_plane = 2;
          break;
        case 9:
          n_plane = 3;
          break;
        case 10:
          n_plane = 2;
          break;
        case 11:
          n_plane = 1;
          break;
        }
        props["defines/n_plane"] = n_plane;
        props["defines/pts_per_thread"] = Nq / n_plane;
      }
    }
    const std::string oklpath(getenv("NEKRS_KERNEL_DIR"));

    // only a single choice, no need to run benchmark
    if (kernelVariants.size() == 1 || !requiresBenchmark) {

      auto newProps = props;
      newProps["defines/p_knl"] = kernelVariants.front();

      const std::string ext = platform->serial ? ".c" : ".okl";
      const std::string fileName = oklpath + "/elliptic/" + kernelName + ext;

      return std::make_pair(platform->device.buildKernel(fileName, newProps, suffix, true), -1.0);
    }

    auto DrV = randomVector<FPType>(Nq * Nq);
    auto ggeo = randomVector<FPType>(Np_g * Nelements * p_Nggeo);
    auto vgeo = randomVector<FPType>(Np * Nelements * p_Nvgeo);
    auto q = randomVector<FPType>((Ndim * Np) * Nelements);
    auto Aq = randomVector<FPType>((Ndim * Np) * Nelements);
    auto exyz = randomVector<FPType>((3 * Np_g) * Nelements);
    auto gllwz = randomVector<FPType>(2 * Nq_g);
    auto lambda0 = randomVector<FPType>(Np * Nelements);
    auto lambda1 = randomVector<FPType>(Np * Nelements);

    // elementList[e] = e
    std::vector<dlong> elementList(Nelements);
    std::iota(elementList.begin(), elementList.end(), 0);
    auto o_elementList = platform->device.malloc(Nelements * sizeof(dlong), elementList.data());

    auto o_D = platform->device.malloc(Nq * Nq * wordSize, DrV.data());
    auto o_S = o_D;
    auto o_ggeo = platform->device.malloc(Np_g * Nelements * p_Nggeo * wordSize, ggeo.data());
    auto o_vgeo = platform->device.malloc(Np * Nelements * p_Nvgeo * wordSize, vgeo.data());    
    auto o_q = platform->device.malloc((Ndim * Np) * Nelements * wordSize, q.data());
    auto o_Aq = platform->device.malloc((Ndim * Np) * Nelements * wordSize, Aq.data());
    auto o_exyz = platform->device.malloc((3 * Np_g) * Nelements * wordSize, exyz.data());
    auto o_gllwz = platform->device.malloc(2 * Nq_g * wordSize, gllwz.data());

    auto o_lambda0 = platform->device.malloc(Np * Nelements * wordSize, lambda0.data());
    auto o_lambda1 = platform->device.malloc(Np * Nelements * wordSize, lambda1.data());

    occa::kernel referenceKernel;
    {
      auto newProps = props;
      newProps["defines/p_knl"] = kernelVariants.front();

      const std::string ext = platform->serial ? ".c" : ".okl";
      const std::string fileName = oklpath + "/elliptic/" + kernelName + ext;

      referenceKernel = platform->device.buildKernel(fileName, newProps, suffix, true);
    }

    auto kernelRunner = [&](occa::kernel &kernel) {
      const int loffset = 0;
      const int offset = Nelements * Np;
      if (computeGeom) {
        kernel(Nelements, offset, loffset, o_elementList, o_exyz, o_gllwz, o_D, o_S, o_lambda0, o_lambda1, o_q, o_Aq);
      }
      else {
        if (!stressForm) {
          kernel(Nelements, offset, loffset, o_elementList, o_ggeo, o_D, o_S, o_lambda0, o_lambda1, o_q, o_Aq);
        } else {
          kernel(Nelements, offset, loffset, o_elementList, o_vgeo, o_D, o_S, o_lambda0, o_lambda1, o_q, o_Aq);
        }
      }
    };

    auto axKernelBuilder = [&](int kernelVariant) {
      auto newProps = props;
      newProps["defines/p_knl"] = kernelVariant;

      const std::string ext = platform->serial ? ".c" : ".okl";
      const std::string fileName = oklpath + "/elliptic/" + kernelName + ext;

      auto kernel = platform->device.buildKernel(fileName, newProps, suffix, true);
      
      if (platform->options.compareArgs("BUILD ONLY", "TRUE"))
        return kernel;

      std::vector<FPType> refResults((Ndim * Np) * Nelements);
      std::vector<FPType> results((Ndim * Np) * Nelements);

      kernelRunner(referenceKernel);
      o_Aq.copyTo(refResults.data(), refResults.size() * sizeof(FPType));

      kernelRunner(kernel);
      o_Aq.copyTo(results.data(), results.size() * sizeof(FPType));

      double err = 0.0;
      for (int i = 0; i < refResults.size(); ++i) {
        err = std::max(err, (double) std::abs((refResults[i] - results[i])/refResults[i]));
      }
      MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
  
      const auto tol = 400. * std::numeric_limits<FPType>::epsilon();
      if (err > tol) {
        if (platform->comm.mpiRank == 0 && verbosity > 1) {
          std::cout << "Ax: Ignore kernel " << kernelVariant
                    << " because error of " << err
                    << " is too large compared to reference\n";
        }
  
        // pass un-initialized kernel to skip this kernel variant
        kernel = occa::kernel();
      }

      return kernel;
    };

    auto printPerformanceInfo = [&](int kernelVariant, double elapsed, int Ntests, bool skipPrint) {
      // print statistics
      const dfloat GDOFPerSecond = (Nelements * Ndim * (N * N * N) / elapsed) / 1.e9;

      size_t bytesMoved = Ndim * 2 * Np * wordSize; // x, Ax
      bytesMoved += 6 * Np_g * wordSize;            // geo

      if(!poisson || stressForm)
        bytesMoved += 1 * Np * wordSize; // Jw

      if (!constCoeff) {
        bytesMoved += 1 * Np * wordSize; // lambda1
        if(!poisson) bytesMoved += 1 * Np * wordSize; // lambda2
      }

      if (stressForm)
        bytesMoved += 3 * Np_g * wordSize; 

      const double bw = (Nelements * bytesMoved / elapsed) / 1.e9;

      double flopCount = Np * 12 * Nq + 15 * Np;
      if(constCoeff)
        flopCount += 1 * Np;
       else 
        flopCount += 3 * Np;

      if(!poisson)
        flopCount += 3 * Np;

      if (stressForm)
        flopCount += 21 * Np;

      const double gflops = Ndim * (flopCount * Nelements / elapsed) / 1.e9;
#ifdef _OPENMP
      const int Nthreads = omp_get_max_threads();
#else
      const int Nthreads = 1;
#endif

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
                    << " GFLOPS/s=" << gflops << " constCoeff=" << constCoeff << " poisson=" << poisson
                    << " kernelVer=" << kernelVariant << "\n";
        }
      }
    };

    auto printCallBack = [&](int kernelVariant, double elapsed, int Ntests) {
      printPerformanceInfo(kernelVariant, elapsed, Ntests, verbosity < 2);
    };

    auto kernelAndTime =
        benchmarkKernel(axKernelBuilder, kernelRunner, printCallBack, kernelVariants, NtestsOrTargetTime);

    if (kernelAndTime.first.properties().has("defines/p_knl") &&
        platform->options.compareArgs("BUILD ONLY", "FALSE")) {
      int bestKernelVariant = static_cast<int>(kernelAndTime.first.properties()["defines/p_knl"]);

      // print only the fastest kernel
      if (verbosity == 1) {
        printPerformanceInfo(bestKernelVariant, kernelAndTime.second, 0, false);
      }
    }

    free(o_D);
    free(o_S);
    free(o_ggeo);
    free(o_vgeo);
    free(o_q);
    free(o_Aq);
    free(o_exyz);
    free(o_gllwz);
    free(o_lambda0);
    free(o_lambda1);
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
                                       size_t wordSize,
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
                                          size_t wordSize,
                                          int Ndim,
                                          bool stressForm,
                                          int verbosity,
                                          double targetTime,
                                          bool requiresBenchmark,
                                          std::string suffix);
