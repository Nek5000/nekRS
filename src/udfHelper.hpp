#if !defined(nekrs_udfhelper_hpp_)
#define nekrs_udfhelper_hpp_

#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"

typedef void (*udfsetup0)(MPI_Comm comm, setupAide &options);
typedef void (*udfsetup)(ins_t *ins);
typedef void (*udfloadKernels)(ins_t *ins);
typedef void (*udfexecuteStep)(ins_t *ins, dfloat time, int tstep);
typedef void (*udfinit)(ins_t *ins);

typedef void (*udfuEqnSource)(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_FU);
typedef void (*udfsEqnSource)(ins_t *ins, dfloat time, occa::memory o_S, occa::memory o_SU);

typedef struct
{
  udfinit init;
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
void udfInit(ins_t *ins);
void *udfLoadFunction(const char *fname, int errchk);
occa::kernel udfBuildKernel(ins_t *ins, const char *function);

#endif
