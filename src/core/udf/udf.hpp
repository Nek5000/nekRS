#if !defined(nekrs_udf_hpp_)
#define nekrs_udf_hpp_

#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "iofldFactory.hpp"
#include "opSEM.hpp" 
#include "tinyexpr.h"

extern "C" {
void UDF_Setup0(MPI_Comm comm, setupAide &options);
void UDF_Setup();
void UDF_LoadKernels(deviceKernelProperties& kernelInfo);
void UDF_AutoLoadKernels(occa::properties &kernelInfo);
void UDF_AutoLoadPlugins(occa::properties &kernelInfo);
void UDF_ExecuteStep(double time, int tstep);
void UDF_Finalize();
}

using udfsetup0 = void (*)(MPI_Comm, setupAide &);
using udfsetup = void (*)();
using udfloadKernels = void (*)(deviceKernelProperties &);
using udfautoloadKernels = void (*)(occa::properties &);
using udfautoloadPlugins = void (*)(occa::properties &);
using udfexecuteStep = void (*)(double, int);
using udffinalize = void (*)();

struct UDF {
  udfsetup0 setup0;
  udfsetup setup;
  udfloadKernels loadKernels;
  udfautoloadKernels autoloadKernels;
  udfautoloadPlugins autoloadPlugins;
  udfexecuteStep executeStep;
  udffinalize finalize;
};

extern UDF udf;

void oudfFindDirichlet(const std::string &field);
void oudfFindNeumann(const std::string &field);
void oudfFindRobin(const std::string &field);
void oudfInit(setupAide &options);
void udfBuild(setupAide &options);
void udfLoad();
int udfMake(setupAide &options, const std::string &solverName, int rank);
void udfEcho();
void *udfLoadFunction(const char *fname, int errchk);
void udfUnload();
occa::kernel oudfBuildKernel(occa::properties kernelInfo, const std::string& kernelName);

#endif
