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
#include "platform.hpp"

namespace{

std::string gen_suffix(const elliptic_t * elliptic, const char * floatString)
{
  const std::string precision = std::string(floatString);
  if(precision.find(pfloatString) != std::string::npos){
    return std::string("_") + std::to_string(elliptic->mesh->N) + std::string("pfloat");
  }
  else{
    return std::string("_") + std::to_string(elliptic->mesh->N);
  }
  
}

}

elliptic_t* ellipticBuildMultigridLevel(elliptic_t* baseElliptic, int Nc, int Nf)
{
  
  elliptic_t* elliptic = new elliptic_t();
  memcpy(elliptic,baseElliptic,sizeof(elliptic_t));

  mesh_t* mesh = createMeshMG(baseElliptic->mesh, Nc);
  elliptic->mesh = mesh;

  { // setup an unmasked gs handle
    ogs_t *ogs = NULL;
    ellipticOgs(mesh,
                mesh->Nlocal,
                /* nFields */ 1,
                /* offset */ 0,
                elliptic->EToB,
                elliptic->Nmasked,
                elliptic->o_maskIds,
                elliptic->NmaskedLocal,
                elliptic->o_maskIdsLocal,
                elliptic->NmaskedGlobal,
                elliptic->o_maskIdsGlobal,
                &ogs);
    elliptic->ogs = ogs;
    elliptic->o_invDegree = elliptic->ogs->o_invDegree;
  }

  const std::string suffix = "Hex3D";

  std::string kernelName;

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();

  ellipticBuildPreconditionerKernels(elliptic);

  const std::string poissonPrefix = elliptic->poisson ? "poisson-" : "";

  if(Nc > 1 || elliptic->options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE"))
  {
      const std::string AxSuffix = elliptic->coeffFieldPreco ? "CoeffHex3D" : "Hex3D";
      // check for trilinear
      if(elliptic->elementType != HEXAHEDRA) {
        kernelName = "ellipticPartialAx" + AxSuffix;
      }else {
        if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
          kernelName = "ellipticPartialAxTrilinear" + AxSuffix;
        else
          kernelName = "ellipticPartialAx" + AxSuffix;
      }

      {
        const std::string kernelSuffix = gen_suffix(elliptic, dfloatString);
        elliptic->AxKernel = platform->kernels.get(poissonPrefix + kernelName + kernelSuffix);
      }
      {
        const std::string kernelSuffix = gen_suffix(elliptic, pfloatString);
        elliptic->AxPfloatKernel =
          platform->kernels.get(poissonPrefix + kernelName + kernelSuffix);
      }
  }

  elliptic->precon = new precon_t();

  {
    const std::string kernelSuffix =
        std::string("_Nf_") + std::to_string(Nf) + std::string("_Nc_") + std::to_string(Nc);

    kernelName = "ellipticPreconCoarsen" + suffix;
    elliptic->precon->coarsenKernel = platform->kernels.get(kernelName + kernelSuffix);
    kernelName = "ellipticPreconProlongate" + suffix;
    elliptic->precon->prolongateKernel = platform->kernels.get(kernelName + kernelSuffix);

  }

  elliptic->o_lambdaPfloat = platform->device.malloc(2 * mesh->Nelements * mesh->Np, sizeof(pfloat));
  elliptic->o_lambda = platform->device.malloc(2 * mesh->Nelements * mesh->Np, sizeof(dfloat));

  const int Nfq = Nf+1;
  const int Ncq = Nc+1;
  dfloat* fToCInterp = (dfloat*) calloc(Nfq * Ncq, sizeof(dfloat));
  InterpolationMatrix1D(Nf, Nfq, baseElliptic->mesh->r, Ncq, mesh->r, fToCInterp);
  elliptic->o_interp = platform->device.malloc(Nfq * Ncq * sizeof(dfloat), fToCInterp);

  elliptic->precon->coarsenKernel(2 * mesh->Nelements, elliptic->o_interp, baseElliptic->o_lambda, elliptic->o_lambda);

  elliptic->copyDfloatToPfloatKernel(2 * mesh->Nelements * mesh->Np,
    elliptic->o_lambda,
    elliptic->o_lambdaPfloat);
  
  free(fToCInterp);


  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel);
  fflush(stdout);

  return elliptic;
}
