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

#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"
#include "platform.hpp"
#include "linAlg.hpp"

void meshGeometricFactorsHex3D(mesh3D* mesh, int ifcub)
{
  double tStart = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("computing geometric factors ... "); fflush(stdout);

  mesh->vgeo = (dfloat*) calloc(mesh->Nelements * mesh->Nvgeo * mesh->Np, sizeof(dfloat));

  mesh->cubvgeo = (dfloat*) calloc(mesh->Nelements * mesh->Nvgeo * mesh->cubNp, sizeof(dfloat));

  mesh->ggeo    = (dfloat*) calloc(mesh->Nelements * mesh->Nggeo * mesh->Np,    sizeof(dfloat));

  mesh->o_vgeo = platform->device.malloc(mesh->Nelements * mesh->Nvgeo * mesh->Np * sizeof(dfloat), mesh->vgeo);

  mesh->o_cubvgeo =  platform->device.malloc(mesh->Nelements * mesh->Nvgeo * mesh->cubNp * sizeof(dfloat), mesh->cubvgeo);

  mesh->o_ggeo    =  platform->device.malloc(mesh->Nelements * mesh->Nggeo * mesh->Np * sizeof(dfloat), mesh->ggeo);

  mesh->LMM = (dfloat*) calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));
  mesh->invLMM = (dfloat*) calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));
  mesh->o_LMM =
    platform->device.calloc(mesh->Nelements * mesh->Np ,  sizeof(dfloat));
  mesh->o_invLMM =
    platform->device.calloc(mesh->Nelements * mesh->Np ,  sizeof(dfloat));
  mesh->o_D = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);
  mesh->o_x =
    platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->x);
  mesh->o_y =
    platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->y);
  mesh->o_z =
    platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->z);
  mesh->o_gllw =
    platform->device.malloc( mesh->Nq * sizeof(dfloat), mesh->gllw);
  mesh->o_cubw =
    platform->device.malloc( mesh->cubNq * sizeof(dfloat), mesh->cubw);
  dfloat* cubInterpT = (dfloat*) calloc(mesh->cubNq * mesh->Nq, sizeof(dfloat));
  for(int n = 0; n < mesh->Nq; ++n)
    for(int m = 0; m < mesh->cubNq; ++m) {
      cubInterpT[m + n * mesh->cubNq] = mesh->cubInterp[m * mesh->Nq + n];
    }
  mesh->o_cubInterpT =
    platform->device.malloc(mesh->Nq * mesh->cubNq * sizeof(dfloat),
                        cubInterpT);

  occa::memory o_scratch = platform->device.calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));

  mesh->geometricFactorsKernel(
        mesh->Nelements,
        ifcub,
        mesh->o_D,
        mesh->o_gllw,
        mesh->o_x,
        mesh->o_y,
        mesh->o_z,
        mesh->o_cubInterpT,
        mesh->o_cubw,
        mesh->o_LMM,
        mesh->o_vgeo,
        mesh->o_ggeo,
        mesh->o_cubvgeo,
        o_scratch
    );


  // copy into host buffers
  mesh->o_LMM.copyTo(mesh->LMM, mesh->Nelements * mesh->Np * sizeof(dfloat));
  mesh->o_vgeo.copyTo(mesh->vgeo, mesh->Nelements * mesh->Nvgeo * mesh->Np * sizeof(dfloat));
  mesh->o_cubvgeo.copyTo(mesh->cubvgeo, mesh->Nelements * mesh->Nvgeo * mesh->cubNp * sizeof(dfloat));
  mesh->o_ggeo.copyTo(mesh->ggeo, mesh->Nelements * mesh->Nggeo * mesh->Np * sizeof(dfloat));

  // compute mesh quality metrics
  const dfloat minJ = platform->linAlg->min(mesh->Nelements * mesh->Np, o_scratch, platform->comm.mpiComm);
  const dfloat maxJ = platform->linAlg->max(mesh->Nelements * mesh->Np, o_scratch, platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0){
    printf("J [%g,%g]\n", minJ, maxJ);
  }
  mesh->volume = platform->linAlg->sum(mesh->Nelements * mesh->Np, mesh->o_LMM, platform->comm.mpiComm);


  o_scratch.free();
  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
}
