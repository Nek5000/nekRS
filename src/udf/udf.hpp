#if !defined(nekrs_udf_hpp_)
#define nekrs_udf_hpp_

#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"

extern int ciMode;

extern "C" {
  void UDF_Setup0(MPI_Comm comm, setupAide &options);
  void UDF_Setup(ins_t *ins);
  void UDF_LoadKernels(ins_t *ins);
  void UDF_ExecuteStep(ins_t *ins, dfloat time, int tstep);
};


typedef void (*udfsetup0)(MPI_Comm comm, setupAide &options);
typedef void (*udfsetup)(ins_t *ins);
typedef void (*udfloadKernels)(ins_t *ins);
typedef void (*udfexecuteStep)(ins_t *ins, dfloat time, int tstep);

typedef void (*udfuEqnSource)(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_FU);

typedef void (*udfsEqnSource)(ins_t *ins, dfloat time, occa::memory o_S, occa::memory o_SU);

typedef struct
{
  udfsetup0 setup0;
  udfsetup setup;
  udfloadKernels loadKernels;
  udfexecuteStep executeStep;
  udfuEqnSource uEqnSource;
  udfsEqnSource sEqnSource;
} UDF;

extern UDF udf;

void udfBuild(const char *udfFile);
void udfLoad(void);
void *udfLoadFunction(const char *fname, int errchk);
occa::kernel udfBuildKernel(ins_t *ins, const char *function);

#endif
