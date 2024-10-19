#include "benchmarkAdvsub.hpp"
#include <vector>
#include <numeric>
#include <iostream>

#include "randomVector.hpp"
#include "kernelBenchmarker.hpp"
#ifdef _OPENMP
#include "omp.h"
#endif

namespace
{

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

namespace std
{
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

namespace
{
std::map<CallParameters, occa::kernel> cachedResults;
}

template <typename T>
occa::kernel benchmarkAdvsub(int Nfields,
                             dlong Nelements,
                             int Nq,
                             int cubNq,
                             int nEXT,
                             bool dealias,
                             bool isScalar,
                             int verbosity,
                             T NtestsOrTargetTime,
                             bool runAutotuner)
{
  if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
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

  nekrsCheck(Nq > 14, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "Nq > 14 is unsupported");

  const auto largestCubNq = maximumCubaturePoints.at(Nq);

  nekrsCheck(cubNq > largestCubNq,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "cubNq > %d is unsupported\n",
             largestCubNq);

  if (!dealias || cubNq < Nq) {
    cubNq = Nq;
  }

  const int N = Nq - 1;
  const int Np = Nq * Nq * Nq;
  const int Ntotal = Np * Nelements;
  const dlong fieldOffset = alignStride<dfloat>(Ntotal);

  const int cubN = cubNq - 1;
  const int cubNp = cubNq * cubNq * cubNq;
  const dlong cubatureOffset = alignStride<dfloat>(cubNp * Nelements); 

  const int NVfields = 3;

  occa::properties props = platform->kernelInfo + meshKernelProperties(N);
  props["defines"].asObject();
  props["includes"].asArray();
  props["header"].asArray();
  props["flags"].asObject();
  props["include_paths"].asArray();

  props["defines/p_cubNq"] = cubNq;
  props["defines/p_cubNp"] = cubNp;
  props["defines/p_nEXT"] = nEXT;

  props["defines/p_NVfields"] = NVfields;
  props["defines/p_MovingMesh"] = platform->options.compareArgs("MOVING MESH", "TRUE");

  const std::string oklpath(getenv("NEKRS_KERNEL_DIR"));

  std::string diffDataFile = oklpath + "/mesh/constantDifferentiationMatrices.h";
  std::string interpDataFile = oklpath + "/mesh/constantInterpolationMatrices.h";
  std::string diffInterpDataFile = oklpath + "/mesh/constantDifferentiationInterpolationMatrices.h";

  props["includes"] += diffDataFile.c_str();
  props["includes"] += interpDataFile.c_str();
  props["includes"] += diffInterpDataFile.c_str();

  if (platform->device.mode() == "dpcpp") {
    props["simd_length"] = 16;
  }

  const std::string ext = (platform->device.mode() == "Serial" && dealias) ? ".c" : ".okl";
  const std::string kernelName = std::string("subCycleStrong") + (dealias ? "Cubature" : "") + std::string("VolumeHex3D");
  const std::string fileName = oklpath + (isScalar ? "/nrs/cds/" : "/nrs/") +  kernelName + ext;

  std::vector<int> kernelVariants = {0};
  if (!platform->serial && dealias) {

    if (!isScalar) {

      std::vector<int> kernelSearchSpace = {6, 7, 8, 9, 16};
      for (auto i : kernelSearchSpace) {
        // v12 requires cubNq <=13
        if (i == 11 && cubNq > 13) {
          continue;
        }

        // v14 requires cubNq <=12
        if (i == 14 && cubNq > 12) {
          continue;
        }

        // v14 requires cubNq <=12
        if (i == 16 && cubNq > 14) {
          continue;
        }

        kernelVariants.push_back(i);
      }
    } else {
      kernelVariants.push_back(8);
    }
  }

  auto buildKernel = [&props, &fileName, &kernelName, &isScalar](int ver)
  {
    auto newProps = props;
    newProps["defines/p_knl"] = ver;
    const auto verSuffix = "_v" + std::to_string(ver);

    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      const auto reqName = std::string(fileName) + "::" + std::string(newProps.hash().getString());
      platform->kernelRequests.add(reqName, fileName, newProps);
      return occa::kernel();
    } else {
      auto knl = platform->device.loadKernel(fileName, kernelName + verSuffix, newProps);
      return knl;
    }
  };

  auto referenceKernel = buildKernel(kernelVariants.front());

  if (!runAutotuner) {
    return referenceKernel;
  }

  const auto wordSize = sizeof(dfloat);

#if 1
  auto lambda = randomVector<dfloat>(4, 1, 2, true);
#else
  std::vector<dfloat> lambda(4, 1.0);
#endif

  auto generateField = [&](const int Nfields, const int Nq, const dlong offset, const dfloat lambda)
  {
    const dlong Np = Nq * Nq * Nq;
    const dlong Nlocal = Nelements * Np; 
    std::vector<dfloat> out((Nfields > 1) ? Nfields * offset : Nlocal, 0.0);

   // convert to [-1, 1] 
    auto convertToRange =  [&](int n) {
      return 2.0 * (static_cast<double>(n) / Nq) - 1;
    };

    for (int f = 0; f < Nfields; f++) {
      for (int e = 0; e < Nelements; e++) {
        for (int i = 0; i < Nq; i++) {
          for (int j = 0; j < Nq; j++) {
            for (int k = 0; k < Nq; k++) {
              const auto x = convertToRange(i);
              const auto y = convertToRange(j);
              const auto z = convertToRange(k);

              const auto id = i * j * k + e * Np + f * offset;
              out[id] = (f+1) * lambda * sin(M_PI * x + e) * sin(M_PI * y) * sin(M_PI * z);
            }   
          }
        }
      }
    }
    return out;
  };

  std::vector<dlong> elementList(Nelements);
  std::iota(elementList.begin(), elementList.end(), 0);
  auto o_elementList = platform->device.malloc(elementList.size() * sizeof(dlong), elementList.data());

#if 1
  auto cubD = randomVector<dfloat>(cubNq * cubNq, 1, 2, true);
#else
  std::vector<dfloat> cubD(cubNq * cubNq, 1.0);
#endif
  auto o_cubD = platform->device.malloc(cubD.size() * wordSize, cubD.data());

#if 1
  auto cubInterpT = randomVector<dfloat>(Nq * cubNq, 1, 2, true);
#else
  std::vector<dfloat> cubInterpT(Nq * cubNq, 1.0);
#endif

  auto o_cubInterpT = platform->device.malloc(cubInterpT.size() * wordSize, cubInterpT.data());

  auto conv = generateField(NVfields * nEXT, cubNq, cubatureOffset, lambda[0]);
  auto o_conv = platform->device.malloc(conv.size() * wordSize, conv.data());

  auto invLMM = generateField(nEXT, Nq, fieldOffset, lambda[1]);
  auto o_invLMM = platform->device.malloc(invLMM.size() * wordSize, invLMM.data());

  auto BdivW = generateField(nEXT, Nq, fieldOffset, lambda[2]);
  auto o_BdivW = platform->device.malloc(BdivW.size() * wordSize, BdivW.data());

  auto Ud = generateField(Nfields, Nq, fieldOffset, lambda[3]);
  auto o_Ud = platform->device.malloc(Ud.size() * wordSize, Ud.data());

  std::vector<dfloat> NU(Nfields * fieldOffset, 0);
  auto o_NU = platform->device.malloc(NU.size() * wordSize, NU.data());

  auto buildKernel2 = [&props, &oklpath](const std::string& _fileName)
  {
    const auto fileName = oklpath + _fileName;

    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      const auto reqName = std::string(fs::path(fileName).filename()) + "::" + std::string(props.hash().getString());
      platform->kernelRequests.add(reqName, fileName, props);
      return occa::kernel();
    } else {
      return platform->device.loadKernel(fileName, props);
    }
  };

  auto readCubDMatrixKernel = buildKernel2("/nrs/readCubDMatrix.okl");
  auto readIMatrixKernel = buildKernel2("/nrs/readIMatrix.okl");
  if(readCubDMatrixKernel.isInitialized() && readIMatrixKernel.isInitialized()) {
    readCubDMatrixKernel(o_cubD);
    readIMatrixKernel(o_cubInterpT);
  }

  auto kernelRunner = [&](occa::kernel &subcyclingKernel) {
    const dfloat c0 = 0.1;
    const dfloat c1 = 0.2;
    const dfloat c2 = 0.3;

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
    } else {
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
    auto kernel = buildKernel(kernelVariant);
    if (!kernel.isInitialized()) return occa::kernel();

    kernelRunner(referenceKernel);
    auto o_NUref = platform->device.malloc(o_NU.size());
    o_NU.copyTo(o_NUref);

    kernelRunner(kernel);

    for(int i = 0; i < Nfields; i++) {
      std::vector<dfloat> referenceResults(Ntotal);
      std::vector<dfloat> results(referenceResults.size());

      o_NUref.copyTo(referenceResults.data(), referenceResults.size() * wordSize, i*fieldOffset * wordSize);
      o_NU.copyTo(results.data(), results.size() * wordSize, i*fieldOffset * wordSize);

      const auto absTol = 1.0;
      const auto err = maxRelErr<dfloat>(referenceResults, results, platform->comm.mpiComm, absTol);
      const auto scale = 100 * range<dfloat>(referenceResults, absTol);
      const auto eps = scale * std::numeric_limits<dfloat>::epsilon();
 
      if (err > eps || std::isnan(err)) {
        if (platform->comm.mpiRank == 0 && verbosity > 1) {
          std::cout << "advSub: Ignore version " << kernelVariant << " as correctness check failed with err=" << err
                    << std::endl;
        }
 
        // pass un-initialized kernel to skip this kernel variant
        return occa::kernel();
      } else {
        if (platform->comm.mpiRank == 0 && verbosity > 1) {
          std::cout << "advSub: kernel version " << kernelVariant << " passed correctness check with err=" << err
                    << std::endl;
        }
      }
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
    } else {
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
      !platform->options.compareArgs("REGISTER ONLY", "TRUE")) {

    // print only the fastest kernel
    if (verbosity == 1) {
      int bestKernelVariant = static_cast<int>(kernelAndTime.first.properties()["defines/p_knl"]);
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
                                           bool runAutotuner);

template occa::kernel benchmarkAdvsub<double>(int Nfields,
                                              int Nelements,
                                              int Nq,
                                              int cubNq,
                                              int nEXT,
                                              bool dealias,
                                              bool isScalar,
                                              int verbosity,
                                              double targetTime,
                                              bool runAutotuner);
