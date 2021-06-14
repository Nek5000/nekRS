
#include "ogstypes.h"
#include "ogs_FINDPTS.hpp"
#include "ogsKernels_FINDPTS.hpp"

occa::kernel ogs::initFindptsKernel(MPI_Comm comm, occa::device device,
                                    dlong D, const dlong *n) {
  std::string oklDir;
  oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
  oklDir += "/okl/ogsFindpts/";

  occa::properties kernelInfo;

  kernelInfo["defines/ " "p_NR"] = n[0];
  kernelInfo["defines/ " "p_NS"] = n[1];
  if (D == 3){
    kernelInfo["defines/ " "p_NT"] = n[2];
  } else {
    kernelInfo["defines/ " "p_NT"] = 1;
  }
  kernelInfo["defines/ " "dlong"] = dlongString;
  kernelInfo["defines/ " "hlong"] = hlongString;
  kernelInfo["defines/ " "dfloat"] = dfloatString;

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  occa::kernel findpts_local_eval;

  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
      if (D == 2) {
        findpts_local_eval = device.buildKernel(oklDir+"findpts_local_eval.okl", "findpts_local_eval_2", kernelInfo);
      } else {
        findpts_local_eval = device.buildKernel(oklDir+"findpts_local_eval.okl", "findpts_local_eval_3", kernelInfo);
      }
      findpts_local_eval.dontUseRefs();
    }
    MPI_Barrier(comm);
  }

  return findpts_local_eval;
}
