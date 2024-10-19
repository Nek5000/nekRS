#include "benchmarkAx.hpp"
#include <vector>
#include <iostream>
#include <numeric>

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
                         bool runAutotuner,
                         std::string suffix)
{
  const std::string oklpath(getenv("NEKRS_KERNEL_DIR"));
  const std::string ext = platform->serial ? ".c" : ".okl";

  if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
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
  
  props["defines/p_cubNq"] = Nq; // Needed for const differentiation matrices
  std::string diffDataFile = oklpath + "/mesh/constantGLLDifferentiationMatrices.h";
  props["includes"].asArray();
  props["includes"] += diffDataFile.c_str();

  if (wordSize == 4) {
    props["defines/dfloat"] = "float";
    props["defines/FP32"] = 1;
  }
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
  const std::string fileName = oklpath + "/elliptic/" + kernelName;

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
        const int Nkernels = 10;
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
        const int Nkernels = 5;
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

    auto buildConstMatrixReader = [&props, &oklpath](const std::string& _fileName) {
      const auto filePath = oklpath + _fileName;

      if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
        const auto reqName = std::string(fs::path(filePath).filename()) + "::" +  std::string(props.hash().getString());
        platform->kernelRequests.add(reqName, filePath, props);
        return occa::kernel();
      } else {
        return platform->device.loadKernel(filePath, props);
      }
    };
    auto readConstDMatrixKernel = buildConstMatrixReader("/nrs/readCubDMatrix.okl");

    auto ggeo = randomVector<FPType>(Np_g * Nelements * p_Nggeo, 0, 1, true);
    auto vgeo = randomVector<FPType>(Np * Nelements * p_Nvgeo, 0, 1, true);
    auto q = randomVector<FPType>((Ndim * Np) * Nelements, 0, 1, true);
    auto Aq = randomVector<FPType>((Ndim * Np) * Nelements, 0, 1, true);
    auto exyz = randomVector<FPType>((3 * Np_g) * Nelements, 0, 1, true);
    auto gllwz = randomVector<FPType>(2 * Nq_g, 0, 1, true);
    auto lambda0 = randomVector<FPType>(Np * Nelements, 0.01, 0.02, true);
    auto lambda1 = randomVector<FPType>(Np * Nelements, 0.2, 0.3, true);

    // elementList[e] = e
    std::vector<dlong> elementList(Nelements);
    std::iota(elementList.begin(), elementList.end(), 0);
    auto o_elementList = platform->device.malloc(Nelements * sizeof(dlong), elementList.data());

    auto o_D = platform->device.malloc(Nq * Nq * wordSize);
    if(readConstDMatrixKernel.isInitialized()) {readConstDMatrixKernel(o_D);}

    auto o_S = o_D;
    auto o_ggeo = platform->device.malloc(Np_g * Nelements * p_Nggeo * wordSize, ggeo.data());
    auto o_vgeo = platform->device.malloc(Np * Nelements * p_Nvgeo * wordSize, vgeo.data());    
    auto o_q = platform->device.malloc((Ndim * Np) * Nelements * wordSize, q.data());
    auto o_exyz = platform->device.malloc((3 * Np_g) * Nelements * wordSize, exyz.data());
    auto o_gllwz = platform->device.malloc(2 * Nq_g * wordSize, gllwz.data());

    auto o_Aq = platform->device.malloc((Ndim * Np) * Nelements * wordSize);
    o_Aq.copyFrom(Aq.data());

    auto o_lambda0 = platform->device.malloc(Np * Nelements * wordSize, lambda0.data());
    auto o_lambda1 = platform->device.malloc(Np * Nelements * wordSize, lambda1.data());

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
      auto kernel = buildKernel(kernelVariant);
      if (!kernel.isInitialized()) return occa::kernel();

      std::vector<FPType> refResults((Ndim * Np) * Nelements);
      std::vector<FPType> results((Ndim * Np) * Nelements);
      
      // Reset o_Aq for each kernel variant to avoid buggy/no-op kernels
      // from falsely passing verification.
      o_Aq.copyFrom(refResults.data());
      kernelRunner(referenceKernel);
      o_Aq.copyTo(refResults.data());
      
      o_Aq.copyFrom(results.data());
      kernelRunner(kernel);
      o_Aq.copyTo(results.data());

      const auto absTol = 1e-2;
      const auto err = maxRelErr<FPType>(refResults, results, platform->comm.mpiComm, absTol);
      const auto scale = 10 * range<FPType>(refResults, absTol);

      if (err > scale * std::numeric_limits<FPType>::epsilon() || std::isnan(err)) {
        if (platform->comm.mpiRank == 0 && verbosity > 1) {
          std::cout << "Ax: Ignore version " << kernelVariant
                    << " as correctness check failed with " << err << std::endl;
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
        !platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
 
      // print only the fastest kernel
      if (verbosity == 1) {
        int bestKernelVariant = static_cast<int>(kernelAndTime.first.properties()["defines/p_knl"]);
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
                                       bool runAutotuner,
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
                                          bool runAutotuner,
                                          std::string suffix);
