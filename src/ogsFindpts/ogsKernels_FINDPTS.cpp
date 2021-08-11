
#include "ogstypes.h"
#include "ogs_FINDPTS.hpp"
#include "ogsKernels_FINDPTS.hpp"
#include <cfloat>

std::pair<occa::kernel, occa::kernel> ogs::initFindptsKernel(MPI_Comm comm, occa::device device,
                                                             dlong D, const dlong *n) {

  static occa::kernel findpts_local;
  static occa::kernel findpts_local_eval;
  static occa::device device_cached;
  static dlong D_cached = -1;
  static dlong n_cached[3];

  if (D_cached != D || device_cached != device
      || n_cached[0] != n[0] || n_cached[1] != n[1] || (D == 3 && n_cached[2] != n[2])) {

    D_cached = D;
    device_cached = device;
    n_cached[0] = n[0];
    n_cached[1] = n[1];
    if(D == 3) n_cached[2] = n[2];

    std::string oklDir;
    oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
    oklDir += "/okl/ogsFindpts/";

    occa::properties kernelInfo;
    kernelInfo["defines"].asObject();
    kernelInfo["includes"].asArray();
    kernelInfo["header"].asArray();
    kernelInfo["flags"].asObject();

    kernelInfo["defines/ " "p_D"]  = D;
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
    kernelInfo["defines/ " "DBL_MAX"] = DBL_MAX;

    kernelInfo["includes"] += oklDir+"findpts.okl.hpp";
    kernelInfo["includes"] += oklDir+"poly.okl.hpp";

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    for (int r=0;r<2;r++){
      if ((r==0 && rank==0) || (r==1 && rank>0)) {
        findpts_local = device.buildKernel(oklDir+"findpts_local.okl", "ogs_findpts_local", kernelInfo);
        if (D == 2) {
          findpts_local_eval = device.buildKernel(oklDir+"findpts_local_eval.okl", "findpts_local_eval_2", kernelInfo);
        } else {
          findpts_local_eval = device.buildKernel(oklDir+"findpts_local_eval.okl", "findpts_local_eval_3", kernelInfo);
        }
      }
      MPI_Barrier(comm);
    }
  }

  return std::pair<occa::kernel, occa::kernel>(findpts_local_eval, findpts_local);
}
