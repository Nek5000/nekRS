#include "benchmarkAdvsub.hpp"
#include <vector>
#include <numeric>
#include <iostream>
#include "nrs.hpp"

#include "randomVector.hpp"
#include "kernelBenchmarker.hpp"
#ifdef _OPENMP
#include "omp.h"
#endif

namespace {

// for a given Nq, return the largest cubNq
const std::map<int, int> maximumCubaturePoints = {
    {2, 3},
    {3, 5},
    {4, 6},
    {5, 8},
    {6, 9},
    {7, 11},
    {8, 12},
    {9, 14},
    {10, 15},
    {11, 17},
    {12, 18},
    {13, 20},
    {14, 21},
};
struct CallParameters {
  int Nfields;
  int Nelements;
  int Nq;
  int cubNq;
  int nEXT;
  bool dealias;
  bool isScalar;
};
} // namespace

namespace std {
template <> struct less<CallParameters> {
  bool operator()(const CallParameters &lhs, const CallParameters &rhs) const
  {
    auto tier = [](const CallParameters &v) {
      return std::tie(v.Nfields, v.Nelements, v.Nq, v.cubNq, v.nEXT, v.dealias, v.isScalar);
    };
    return tier(lhs) < tier(rhs);
  }
};
} // namespace std

namespace {
std::map<CallParameters, occa::kernel> cachedResults;
}

template <typename T>
occa::kernel benchmarkAdvsub(int Nfields,
                             int Nelements,
                             int Nq,
                             int cubNq,
                             int nEXT,
                             bool dealias,
                             bool isScalar,
                             int verbosity,
                             T NtestsOrTargetTime,
                             bool requiresBenchmark)
{
  if (platform->options.compareArgs("BUILD ONLY", "TRUE")) {
    Nelements = 1;
  }

  CallParameters params{
      Nfields,
      Nelements,
      Nq,
      cubNq,
      nEXT,
      dealias,
      isScalar,
  };

  if (cachedResults.count(params) > 0) {
    return cachedResults.at(params);
  }

  nrsCheck(Nq > 14, platform->comm.mpiComm, EXIT_FAILURE, 
           "%s\n", "Nq > 14 is unsupported");

  const auto largestCubNq = maximumCubaturePoints.at(Nq);

  nrsCheck(cubNq > largestCubNq, platform->comm.mpiComm, EXIT_FAILURE, 
           "cubNq > %d is unsupported\n", largestCubNq); 

  if (!dealias || cubNq < Nq) {
    cubNq = Nq;
  }

  const int N = Nq - 1;
  const int cubN = cubNq - 1;
  const int Np = Nq * Nq * Nq;
  const int cubNp = cubNq * cubNq * cubNq;
  int fieldOffset = Np * Nelements;
  const int pageW = ALIGN_SIZE / sizeof(dfloat);
  if (fieldOffset % pageW)
    fieldOffset = (fieldOffset / pageW + 1) * pageW;
  int cubatureOffset = std::max(fieldOffset, Nelements * cubNp);
  if (cubatureOffset % pageW)
    cubatureOffset = (cubatureOffset / pageW + 1) * pageW;

  occa::properties props = platform->kernelInfo + meshKernelProperties(N);
  props["defines"].asObject();
  props["includes"].asArray();
  props["header"].asArray();
  props["flags"].asObject();
  props["include_paths"].asArray();

  constexpr int NVfields{3};

  props["defines/p_cubNq"] = cubNq;
  props["defines/p_cubNp"] = cubNp;
  props["defines/p_nEXT"] = nEXT;
  props["defines/p_NVfields"] = NVfields;
  props["defines/p_MovingMesh"] = platform->options.compareArgs("MOVING MESH", "TRUE");

  const std::string oklpath(getenv("NEKRS_KERNEL_DIR"));

  std::string diffDataFile = oklpath + "/mesh/constantDifferentiationMatrices.h";
  std::string interpDataFile = oklpath + "/mesh/constantInterpolationMatrices.h";
  std::string diffInterpDataFile =
      oklpath + "/mesh/constantDifferentiationInterpolationMatrices.h";

  props["includes"] += diffDataFile.c_str();
  props["includes"] += interpDataFile.c_str();
  props["includes"] += diffInterpDataFile.c_str();

  std::string fileName = oklpath + "/bench/advsub/readCubDMatrix.okl";
  auto readCubDMatrixKernel = platform->device.buildKernel(fileName, props, true);

  fileName = oklpath + "/bench/advsub/readIMatrix.okl";
  auto readIMatrixKernel = platform->device.buildKernel(fileName, props, true);

  std::string kernelName;
  if (dealias) {
    kernelName = "subCycleStrongCubatureVolumeHex3D";
  }
  else {
    kernelName = "subCycleStrongVolumeHex3D";
  }

  const std::string ext = (platform->device.mode() == "Serial") ? ".c" : ".okl";
  fileName = oklpath + "/nrs/" + kernelName + ext;

  if (isScalar) {
    fileName = oklpath + "/cds/" + kernelName + ext;
  }

  // currently lacking a native implementation of the non-dealiased kernel
  if (!dealias) {
    fileName = oklpath + "/nrs/" + kernelName + ".okl";
    if (isScalar) {
      fileName = oklpath + "/cds/" + kernelName + ".okl";
    }
  }

  std::vector<int> kernelVariants = {0};
  if (!platform->serial && dealias && !isScalar) {

   std::vector<int> kernelSearchSpace = { 6, 7, 8, 9, 16 };
    for (auto i : kernelSearchSpace) {
      // v12 requires cubNq <=13
      if (i == 11 && cubNq > 13)
        continue;

      // v14 requires cubNq <=12
      if (i == 14 && cubNq > 12)
        continue;

      // v14 requires cubNq <=12
      if (i == 16 && cubNq > 14)
        continue;

      kernelVariants.push_back(i);
    }
  }
  else if (!platform->serial && dealias && isScalar) {
    kernelVariants.push_back(8);
  }

  if (kernelVariants.size() == 1 || !requiresBenchmark) {
    auto newProps = props;
    if (!platform->serial && dealias)
      newProps["defines/p_knl"] = kernelVariants.front();
    return platform->device.buildKernel(fileName, newProps, true);
  }

  occa::kernel referenceKernel;
  {
    auto newProps = props;
    if (!platform->serial && dealias)
      newProps["defines/p_knl"] = kernelVariants.front();
    referenceKernel = platform->device.buildKernel(fileName, newProps, true);
  }

  const auto wordSize = sizeof(dfloat);

  auto invLMM = randomVector<dfloat>(fieldOffset * nEXT);
  auto cubD = randomVector<dfloat>(cubNq * cubNq);
  auto NU = randomVector<dfloat>(Nfields * fieldOffset);
  auto conv = randomVector<dfloat>(NVfields * cubatureOffset * nEXT);
  auto cubInterpT = randomVector<dfloat>(Nq * cubNq);
  auto Ud = randomVector<dfloat>(Nfields * fieldOffset);
  auto BdivW = randomVector<dfloat>(fieldOffset * nEXT);

  // elementList[e] = e
  std::vector<dlong> elementList(Nelements);
  std::iota(elementList.begin(), elementList.end(), 0);
  auto o_elementList = platform->device.malloc(Nelements * sizeof(dlong), elementList.data());

  auto o_invLMM = platform->device.malloc((nEXT * wordSize) * fieldOffset, invLMM.data());
  auto o_cubD = platform->device.malloc(cubNq * cubNq * wordSize, cubD.data());
  auto o_NU = platform->device.malloc((Nfields * wordSize) * fieldOffset, NU.data());
  auto o_conv = platform->device.malloc((NVfields * nEXT * wordSize) * cubatureOffset, conv.data());
  auto o_cubInterpT = platform->device.malloc(Nq * cubNq * wordSize, cubInterpT.data());
  auto o_Ud = platform->device.malloc((Nfields * wordSize) * fieldOffset, Ud.data());
  auto o_BdivW = platform->device.malloc((nEXT * wordSize) * fieldOffset, BdivW.data());

  // popular cubD, cubInterpT with correct data
  readCubDMatrixKernel(o_cubD);
  readIMatrixKernel(o_cubInterpT);

  auto kernelRunner = [&](occa::kernel &subcyclingKernel) {
    const auto c0 = 0.1;
    const auto c1 = 0.2;
    const auto c2 = 0.3;
    if (!dealias) {
      subcyclingKernel(Nelements,
                       o_elementList,
                       o_cubD,
                       fieldOffset,
                       0,
                       o_invLMM,
                       o_BdivW,
                       c0,
                       c1,
                       c2,
                       o_conv,
                       o_Ud,
                       o_NU);
    }
    else {
      subcyclingKernel(Nelements,
                       o_elementList,
                       o_cubD,
                       o_cubInterpT,
                       fieldOffset,
                       cubatureOffset,
                       0,
                       o_invLMM,
                       o_BdivW,
                       c0,
                       c1,
                       c2,
                       o_conv,
                       o_Ud,
                       o_NU);
    }
  };

  auto advSubKernelBuilder = [&](int kernelVariant) {
    auto newProps = props;
    if (!platform->serial && dealias)
      newProps["defines/p_knl"] = kernelVariant;

    auto kernel = platform->device.buildKernel(fileName, newProps, true);
    if (platform->options.compareArgs("BUILD ONLY", "TRUE"))
      return kernel;

    // perform correctness check
    std::vector<dfloat> referenceResults(Nfields * fieldOffset);
    std::vector<dfloat> results(Nfields * fieldOffset);

    kernelRunner(referenceKernel);
    o_NU.copyTo(referenceResults.data(), referenceResults.size() * sizeof(dfloat));

    kernelRunner(kernel);
    o_NU.copyTo(results.data(), results.size() * sizeof(dfloat));

    const auto tol = 1000. * std::numeric_limits<dfloat>::epsilon();
    double err = 0;
    for (auto i = 0; i < results.size(); ++i) {
      const auto refValue = std::abs(referenceResults[i]);
      const auto denom = refValue > tol ? refValue : 1.0;
      const auto absDiff = std::abs(results[i] - referenceResults[i]);

      // ignore values that are already near machine epsilon
      if(absDiff > tol){
        const auto relDiff = absDiff / denom;
        err = std::max(err, (double) relDiff);
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);

    if (err > tol) {
      if (platform->comm.mpiRank == 0 && verbosity > 1) {
        std::cout << "advSub: Ignore kernel " << kernelVariant 
                  << " because error of " << err
                  << " is too large compared to reference\n";
      }

      // pass un-initialized kernel to skip this kernel variant
      kernel = occa::kernel();
    }

    return kernel;
  };

  auto printPerformanceInfo = [&](int kernelVariant, double elapsed, int Ntests, bool skipPrint) {
    const dfloat GDOFPerSecond = Nfields * (Nelements * (N * N * N) / elapsed) / 1.e9;

    size_t bytesPerElem = 2 * Nfields * Np; // Ud, NU
    bytesPerElem += Np;                     // inv mass matrix
    bytesPerElem += Nfields * cubNp * nEXT; // U(r,s,t)

    size_t otherBytes = cubNq * cubNq; // D
    if (cubNq > Nq) {
      otherBytes += Nq * cubNq; // interpolator
    }
    otherBytes *= wordSize;
    bytesPerElem *= wordSize;
    const double bw = ((Nelements * bytesPerElem + otherBytes) / elapsed) / 1.e9;

    double flopCount = 0.0; // per elem basis
    if (dealias) {
      flopCount += 6. * cubNp * nEXT;            // extrapolate U(r,s,t) to current time
      flopCount += 6. * cubNp * cubNq * Nfields; // apply Dcub
      flopCount += 3. * Np * Nfields;            // compute NU
      flopCount += 4. * Nq * (cubNp + cubNq * cubNq * Nq + cubNq * Nq * Nq) * Nfields; // interpolation
    }
    else {
      flopCount = Nq * Nq * Nq * (6. * Nq + 6. * nEXT + 8.) * Nfields;
    }
    const double gflops = (flopCount * Nelements / elapsed) / 1.e9;
#ifdef _OPENMP
    const int Nthreads = omp_get_max_threads();
#else
    const int Nthreads = 1;
#endif

    if (platform->comm.mpiRank == 0 && !skipPrint) {

      if (verbosity > 0) {
        std::cout << "advSub:";
      }

      if (verbosity > 1) {
        std::cout << " MPItasks=" << platform->comm.mpiCommSize << " OMPthreads=" << Nthreads
                  << " NRepetitions=" << Ntests;
      }
      if (verbosity > 0) {
        std::cout << " N=" << N;
        if (dealias) {
          std::cout << " cubN=" << cubN;
        }

        if (verbosity > 1) {
          std::cout << " nEXT=" << nEXT;
        }

        std::cout << " Nfields=" << Nfields;
        if (verbosity > 1) {
          std::cout << " Nelements=" << Nelements;
          std::cout << " elapsed time=" << elapsed;
        }
        std::cout << " wordSize=" << 8 * wordSize << " GDOF/s=" << GDOFPerSecond << " GB/s=" << bw
                  << " GFLOPS/s=" << gflops << " kernelVer=" << kernelVariant << "\n";
      }
    }
  };

  auto printCallBack = [&](int kernelVariant, double elapsed, int Ntests) {
    printPerformanceInfo(kernelVariant, elapsed, Ntests, verbosity < 2);
  };

  auto kernelAndTime =
      benchmarkKernel(advSubKernelBuilder, kernelRunner, printCallBack, kernelVariants, NtestsOrTargetTime);

  if (kernelAndTime.first.properties().has("defines/p_knl") &&
      platform->options.compareArgs("BUILD ONLY", "FALSE")) {
    int bestKernelVariant = static_cast<int>(kernelAndTime.first.properties()["defines/p_knl"]);

    // print only the fastest kernel
    if (verbosity == 1) {
      printPerformanceInfo(bestKernelVariant, kernelAndTime.second, 0, false);
    }
  }

  free(o_elementList);
  free(o_invLMM);
  free(o_cubD);
  free(o_NU);
  free(o_conv);
  free(o_cubInterpT);
  free(o_Ud);

  cachedResults[params] = kernelAndTime.first;

  return kernelAndTime.first;
}

template occa::kernel benchmarkAdvsub<int>(int Nfields,
                                           int Nelements,
                                           int Nq,
                                           int cubNq,
                                           int nEXT,
                                           bool dealias,
                                           bool isScalar,
                                           int verbosity,
                                           int Ntests,
                                           bool requiresBenchmark);

template occa::kernel benchmarkAdvsub<double>(int Nfields,
                                              int Nelements,
                                              int Nq,
                                              int cubNq,
                                              int nEXT,
                                              bool dealias,
                                              bool isScalar,
                                              int verbosity,
                                              double targetTime,
                                              bool requiresBenchmark);
