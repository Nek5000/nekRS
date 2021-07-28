#if !defined(nekrs_udf_hpp_)
#define nekrs_udf_hpp_

#define ins_t nrs_t

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "parReader.hpp"

extern "C" {
void UDF_Setup0(MPI_Comm comm, setupAide &options);
void UDF_Setup(nrs_t* nrs);
void UDF_LoadKernels(nrs_t* nrs);
void UDF_ExecuteStep(nrs_t* nrs, dfloat time, int tstep);
}

typedef void (* udfsetup0)(MPI_Comm comm, setupAide &options);
typedef void (* udfsetup)(nrs_t* nrs);
typedef void (* udfloadKernels)(nrs_t* nrs);
typedef void (* udfexecuteStep)(nrs_t* nrs, dfloat time, int tstep);

typedef void (* udfuEqnSource)(nrs_t* nrs, dfloat time, occa::memory o_U, occa::memory o_FU);
typedef void (* udfsEqnSource)(nrs_t* nrs, dfloat time, occa::memory o_S, occa::memory o_SU);
typedef void (* udfproperties)(nrs_t* nrs, dfloat time, occa::memory o_U,
                               occa::memory o_S, occa::memory o_UProp,
                               occa::memory o_SProp);
typedef void (* udfdiv)(nrs_t* nrs, dfloat time, occa::memory o_div);
typedef int (* udfconv)(nrs_t* nrs, int stage);

struct UDF
{
  udfsetup0 setup0;
  udfsetup setup;
  udfloadKernels loadKernels;
  udfexecuteStep executeStep;
  udfuEqnSource uEqnSource;
  udfsEqnSource sEqnSource;
  udfproperties properties;
  udfdiv div;
  udfconv converged;
};

extern UDF udf;

int udfBuild(const char* casename, const char* udfFile, int buildOnly);
void udfLoad(const char* casename);
void* udfLoadFunction(const char* casename, const char* fname, int errchk);
occa::kernel udfBuildKernel(nrs_t* nrs, const char* function);

#endif
