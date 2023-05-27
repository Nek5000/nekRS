#if !defined(nekrs_udf_hpp_)
#define nekrs_udf_hpp_

#define ins_t nrs_t

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "parReader.hpp"
#include "bcMap.hpp"
#include "constantFlowRate.hpp"
#include "postProcessing.hpp"
#include <functional>

#define CIPASS                                                                                               \
{                                                                                                            \
if (platform->comm.mpiRank == 0)                                                                             \
printf("TESTS passed \n");                                                                                   \
platform->exitValue = 0;                                                                                     \
}
#define CIFAIL                                                                                               \
{                                                                                                            \
if (platform->comm.mpiRank == 0)                                                                             \
printf("TESTS failed!\n");                                                                                   \
platform->exitValue += 1;                                                                                    \
}

extern "C" {
void UDF_Setup0(MPI_Comm comm, setupAide &options);
void UDF_Setup(nrs_t *nrs);
void UDF_LoadKernels(occa::properties &kernelInfo);
void UDF_AutoLoadKernels(occa::properties &kernelInfo);
void UDF_AutoLoadPlugins(occa::properties &kernelInfo);
void UDF_ExecuteStep(nrs_t *nrs, dfloat time, int tstep);
}

using udfsetup0 = void (*)(MPI_Comm, setupAide &);
using udfsetup = void (*)(nrs_t *);
using udfloadKernels = void (*)(occa::properties &);
using udfautoloadKernels = void (*)(occa::properties &);
using udfautoloadPlugins = void (*)(occa::properties &);
using udfexecuteStep = void (*)(nrs_t *, dfloat, int);

using udfuEqnSource = std::function<void(nrs_t *, dfloat, occa::memory, occa::memory)>;
using udfsEqnSource = std::function<void(nrs_t *, dfloat, occa::memory, occa::memory)>;
using udfproperties =
    std::function<void(nrs_t *, dfloat, occa::memory, occa::memory, occa::memory, occa::memory)>;
using udfdiv = std::function<void(nrs_t *, dfloat, occa::memory)>;
using udfconv = std::function<int(nrs_t *, int)>;
using udfPreFluid = std::function<void(nrs_t*, dfloat, int)>;
using udfPostScalar = std::function<void(nrs_t*, dfloat, int)>;

struct UDF {
  udfsetup0 setup0;
  udfsetup setup;
  udfloadKernels loadKernels;
  udfautoloadKernels autoloadKernels;
  udfautoloadPlugins autoloadPlugins;
  udfexecuteStep executeStep;
  udfuEqnSource uEqnSource;
  udfsEqnSource sEqnSource;
  udfproperties properties;
  udfdiv div;
  udfconv timeStepConverged;
  udfPreFluid preFluid;
  udfPostScalar postScalar;
};

extern UDF udf;

void oudfFindDirichlet(std::string &field);
void oudfFindNeumann(std::string &field);
void oudfInit(setupAide &options);
void udfBuild(const std::string& udfFile, setupAide &options);
void udfLoad();
void udfEcho();
void *udfLoadFunction(const char *fname, int errchk);
occa::kernel oudfBuildKernel(occa::properties kernelInfo, const char *function);

#endif
