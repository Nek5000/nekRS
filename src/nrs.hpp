#if !defined(nekrs_nrs_hpp_)
#define nekrs_nrs_hpp_

#include <string>
#include <mpi.h>
#include "setupAide.hpp"
#include "environment.hpp"

namespace nrs {
  setupAide setup(MPI_Comm comm, int buildOnly, int sizeTarget,
                  int ciMode, string setupFile);
  void runStep(double time, int tstep);
  void copyToNek(double time, int tstep);
  void udfExecuteStep(double time, int tstep, int isOutputStep); 
  void nekOutfld(void);
}

#endif
