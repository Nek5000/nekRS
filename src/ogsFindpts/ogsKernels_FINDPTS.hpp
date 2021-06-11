
#ifndef OGS_KERNELS_FINDPTS_HPP
#define OGS_KERNELS_FINDPTS_HPP

#include "ogsKernels.hpp"

namespace ogs {
      occa::kernel initFindptsKernel(MPI_Comm comm, occa::device device, dlong D, const dlong *n);
}

#endif
