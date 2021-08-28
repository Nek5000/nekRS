/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */

#include "elliptic.h"
#include <string>
#include "platform.hpp"
#include "linAlg.hpp"

void ellipticBuildPreconditionerKernels(elliptic_t* elliptic, occa::properties kernelInfo)
{
  
  mesh_t* mesh      = elliptic->mesh;

  std::string prefix = "Hex3D";
  std::string filename, kernelName;

  kernelInfo["defines/pfloat"] = pfloatString;
  kernelInfo["defines/" "p_eNfields"] = elliptic->Nfields;

  occa::properties pfloatKernelInfo = kernelInfo;
  pfloatKernelInfo["defines/dfloat"] = pfloatString;
  pfloatKernelInfo["defines/pfloat"] = pfloatString;

  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0) printf("loading elliptic preconditioner kernels ... ");
  fflush(stdout);

  const std::string orderSuffix = std::string("_") + std::to_string(mesh->N);

  {
      const std::string oklpath = install_dir + "/okl/core/";
      std::string filename;

      filename = oklpath + "mask.okl";
      mesh->maskKernel =
        platform->device.buildKernel(filename,
                                 "mask",
                                 kernelInfo,
                                 orderSuffix);

      filename = oklpath + "mask.okl";
      mesh->maskPfloatKernel =
        platform->device.buildKernel(filename,
                                 "mask",
                                 pfloatKernelInfo,
                                 orderSuffix);
        filename = install_dir + "/okl/elliptic/ellipticLinAlg.okl";
        elliptic->fusedCopyDfloatToPfloatKernel =
          platform->device.buildKernel(filename,
                                   "fusedCopyDfloatToPfloat",
                                   kernelInfo,
                                 orderSuffix);
        elliptic->copyDfloatToPfloatKernel =
          platform->device.buildKernel(filename,
                                   "copyDfloatToPfloat",
                                   kernelInfo,
                                 orderSuffix);

        elliptic->copyPfloatToDPfloatKernel =
          platform->device.buildKernel(filename,
                                   "copyPfloatToDfloat",
                                   kernelInfo,
                                 orderSuffix);

        elliptic->scaledAddPfloatKernel =
          platform->device.buildKernel(filename,
                                   "scaledAdd",
                                   kernelInfo,
                                 orderSuffix);
        elliptic->dotMultiplyPfloatKernel =
          platform->device.buildKernel(filename,
                                   "dotMultiply",
                                   kernelInfo,
                                 orderSuffix);
        filename = install_dir + "/okl/elliptic/chebyshev.okl";
        elliptic->updateSmoothedSolutionVecKernel =
          platform->device.buildKernel(filename,
                                   "updateSmoothedSolutionVec",
                                   kernelInfo,
                                 orderSuffix);
        elliptic->updateChebyshevSolutionVecKernel =
          platform->device.buildKernel(filename,
                                   "updateChebyshevSolutionVec",
                                   kernelInfo,
                                 orderSuffix);

        elliptic->updateIntermediateSolutionVecKernel =
          platform->device.buildKernel(filename,
                                   "updateIntermediateSolutionVec",
                                   kernelInfo,
                                 orderSuffix);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel);
  fflush(stdout);

}
