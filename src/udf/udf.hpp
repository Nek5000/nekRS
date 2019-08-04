#if !defined(nekrs_udf_hpp_)
#define nekrs_udf_hpp_

#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"

extern int ciMode;

extern "C" {
  void UDF_Setup0(MPI_Comm comm, setupAide &options);
  void UDF_Setup(ins_t *ins);
  void UDF_LoadKernels(ins_t *ins, const occa::properties kernelInfo);
  void UDF_ExecuteStep(ins_t *ins, dfloat time, int tstep);
};


typedef void (*udfsetup0)(MPI_Comm comm, setupAide &options);
typedef void (*udfsetup)(ins_t *ins);
typedef void (*udfloadKernels)(ins_t *ins, const occa::properties kernelInfo);
typedef void (*udfexecuteStep)(ins_t *ins, dfloat time, int tstep);

typedef void (*udfvelocityForce)(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_FU);

typedef struct
{
  udfsetup0 setup0;
  udfsetup setup;
  udfloadKernels loadKernels;
  udfexecuteStep executeStep;
  udfvelocityForce velocityForce;
} UDF;

extern UDF udf;

void udfBuild(const char *udfFile);
void udfLoad(void);
void *udfLoadFunction(const char *fname, int errchk);
occa::kernel udfBuildKernel(mesh_t *mesh, const char *function, const occa::properties kernelInfo);

#endif
