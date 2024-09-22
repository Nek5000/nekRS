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
#include <string.h>
#include  "mpi.h"

#include "mesh3D.h"
#include "platform.hpp"
#include "bcMap.hpp"

void meshOccaPopulateDeviceHex3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo)
{
  mesh->o_elementInfo = platform->device.malloc<dlong>(mesh->Nelements, 
		                            mesh->elementInfo);

  if(mesh->cubNq > 1) { 
    dfloat* cubDWT = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
    dfloat* cubProjectT = (dfloat*) calloc(mesh->cubNq * mesh->Nq, sizeof(dfloat));
    dfloat* cubInterpT = (dfloat*) calloc(mesh->cubNq * mesh->Nq, sizeof(dfloat));
    for(int n = 0; n < mesh->Nq; ++n)
      for(int m = 0; m < mesh->cubNq; ++m) {
        cubProjectT[n + m * mesh->Nq] = mesh->cubProject[n * mesh->cubNq + m];
        cubInterpT[m + n * mesh->cubNq] = mesh->cubInterp[m * mesh->Nq + n];
      }
    for(int n = 0; n < mesh->cubNq; ++n)
      for(int m = 0; m < mesh->cubNq; ++m)
        cubDWT[n + m * mesh->cubNq] = mesh->cubDW[n * mesh->cubNq + m];
 
    mesh->o_cubInterpT =
      platform->device.malloc<dfloat>(mesh->Nq * mesh->cubNq,
                          cubInterpT);
    mesh->o_cubProjectT =
      platform->device.malloc<dfloat>(mesh->Nq * mesh->cubNq,
                          cubProjectT);
    mesh->o_cubDWT =
      platform->device.malloc<dfloat>(mesh->cubNq * mesh->cubNq,
                          cubDWT);
    mesh->o_cubD =
      platform->device.malloc<dfloat>(mesh->cubNq * mesh->cubNq,
                          mesh->cubD);
    mesh->o_cubDWmatrices = 
      platform->device.malloc<dfloat>(mesh->cubNq * mesh->cubNq, cubDWT);

    mesh->o_cubw =
      platform->device.malloc<dfloat>( mesh->cubNq, mesh->cubw);
 
    mesh->o_cubDiffInterpT = mesh->o_cubDWmatrices;

    mesh->o_cubvgeo =
        platform->device.malloc<dfloat>(mesh->Nelements * mesh->Nvgeo * mesh->cubNp);
 
    free(cubProjectT);
    free(cubInterpT);
    free(cubDWT);
  }

  mesh->o_D = platform->device.malloc<dfloat>(mesh->Nq * mesh->Nq, mesh->D);
  mesh->o_DW = platform->device.malloc<dfloat>(mesh->Nq * mesh->Nq, mesh->DW);
  dfloat* DT = (dfloat*) calloc(mesh->Nq * mesh->Nq, sizeof(dfloat));
  for(int j = 0; j < mesh->Nq; ++j)
    for(int i = 0; i < mesh->Nq; ++i)
      DT[i * mesh->Nq + j] = mesh->D[j * mesh->Nq + i];
  mesh->o_DT = platform->device.malloc<dfloat>(mesh->Nq * mesh->Nq);
  mesh->o_DT.copyFrom(DT);
  free(DT);

  mesh->o_gllw =
    platform->device.malloc<dfloat>(mesh->Nq, mesh->gllw);

  mesh->o_gllz = platform->device.malloc<dfloat>(mesh->Nq, mesh->gllz);

  mesh->o_LMM =
    platform->device.malloc<dfloat>(mesh->Nlocal);

  mesh->o_invLMM =
    platform->device.malloc<dfloat>(mesh->Nlocal);

  mesh->o_vgeo =
    platform->device.malloc<dfloat>(mesh->Nlocal * mesh->Nvgeo);

  mesh->o_sgeo =
      platform->device.malloc<dfloat>(mesh->Nelements * mesh->Nfaces * mesh->Nfp * mesh->Nsgeo);

  mesh->o_ggeo =
    platform->device.malloc<dfloat>(mesh->Nlocal * mesh->Nggeo);

  mesh->o_vmapM =
      platform->device.malloc<dlong>(mesh->Nelements * mesh->Nfp * mesh->Nfaces, mesh->vmapM);

  mesh->o_EToB =
    platform->device.malloc<int>(mesh->Nelements * mesh->Nfaces, mesh->EToB);

  mesh->o_faceNodes =
    platform->device.malloc<int>(mesh->Nfaces * mesh->Nfp, mesh->faceNodes);

  kernelInfo += meshKernelProperties(mesh->N);
}

void meshOccaSetup3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo)
{
  meshOccaPopulateDeviceHex3D(mesh, newOptions, kernelInfo);
}
