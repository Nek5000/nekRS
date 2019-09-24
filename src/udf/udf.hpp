#if !defined(nekrs_udf_hpp_)
#define nekrs_udf_hpp_

#include "udfHelper.hpp"
#include "libParanumal.hpp"
#include "nekInterfaceAdapter.hpp"

extern "C" {
  int  UDF_Init(ins_t *ins);
  void UDF_Setup0(MPI_Comm comm, setupAide &options);
  void UDF_Setup(ins_t *ins);
  void UDF_LoadKernels(ins_t *ins);
  void UDF_ExecuteStep(ins_t *ins, dfloat time, int tstep);
};

namespace udf_{
  int init(ins_t *ins);

  dfloat dt();
  int    dt(double &dt);

  dfloat startTime();

  int isOutputStep();
  int readRestartFile();

  dlong nElements();
  int Np();
}

#endif
