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

void reportMemoryUsage(occa::device &device, const char* mess)
{
  size_t bytes = device.memoryAllocated();

  printf("%s: bytes allocated = %lu\n", mess, bytes);
}

void meshOccaPopulateDeviceHex3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo)
{
  mesh->o_elementInfo = platform->device.malloc(mesh->Nelements * sizeof(dlong), 
		                            mesh->elementInfo);
 
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
    platform->device.malloc(mesh->Nq * mesh->cubNq * sizeof(dfloat),
                        cubInterpT);
  mesh->o_cubProjectT =
    platform->device.malloc(mesh->Nq * mesh->cubNq * sizeof(dfloat),
                        cubProjectT);
  mesh->o_cubDWT =
    platform->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat),
                        cubDWT);
  mesh->o_cubD =
    platform->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat),
                        mesh->cubD);
  mesh->o_cubDWmatrices = 
    platform->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat), cubDWT);

  mesh->o_cubDiffInterpT = mesh->o_cubDWmatrices;

  free(cubProjectT);
  free(cubInterpT);
  free(cubDWT);

  mesh->LMM = (dfloat*) calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));
  mesh->invLMM = (dfloat*) calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));
  mesh->o_LMM =
    platform->device.malloc(mesh->Nelements * mesh->Np ,  sizeof(dfloat));
  mesh->o_invLMM =
    platform->device.malloc(mesh->Nelements * mesh->Np ,  sizeof(dfloat));

  mesh->o_D = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);
  mesh->o_DW = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->DW);
  dfloat* DT = (dfloat*) calloc(mesh->Nq * mesh->Nq,sizeof(dfloat));
  for(int j = 0; j < mesh->Nq; ++j)
    for(int i = 0; i < mesh->Nq; ++i)
      DT[i * mesh->Nq + j] = mesh->D[j * mesh->Nq + i];
  mesh->o_DT = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), DT); //dummy
  free(DT);

  mesh->o_vgeo =
    platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nvgeo * sizeof(dfloat),
                        mesh->vgeo);
  mesh->o_sgeo =
      platform->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->Nfp * mesh->Nsgeo * sizeof(dfloat),
                              mesh->sgeo);

  mesh->o_ggeo =
    platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(dfloat),
                        mesh->ggeo);
  if (mesh->cubNq - 1) {
    mesh->o_cubvgeo =
        platform->device.malloc(mesh->Nelements * mesh->Nvgeo * mesh->cubNp * sizeof(dfloat), mesh->cubvgeo);
  }

  mesh->o_vmapM =
    platform->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                        mesh->vmapM);
  mesh->o_vmapP =
    platform->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                        mesh->vmapP);
  mesh->o_EToB =
    platform->device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int),
                        mesh->EToB);

  mesh->o_x =
    platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->x);
  mesh->o_y =
    platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->y);
  mesh->o_z =
    platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->z);
  mesh->o_gllw =
    platform->device.malloc( mesh->Nq * sizeof(dfloat), mesh->gllw);
  mesh->o_gllz = platform->device.malloc(mesh->Nq * sizeof(dfloat), mesh->gllz);
  mesh->o_cubw =
    platform->device.malloc( mesh->cubNq * sizeof(dfloat), mesh->cubw);
  mesh->o_faceNodes =
    platform->device.malloc(mesh->Nfaces * mesh->Nfp * sizeof(int), mesh->faceNodes);

  if(mesh->totalHaloPairs > 0) {
    mesh->o_haloElementList =
      platform->device.malloc(mesh->totalHaloPairs * sizeof(dlong), mesh->haloElementList);

    mesh->o_haloBuffer =
      platform->device.malloc(mesh->totalHaloPairs * mesh->Np * mesh->Nfields ,  sizeof(dfloat));

    mesh->o_haloGetNodeIds =
      platform->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloGetNodeIds);

    mesh->o_haloPutNodeIds =
      platform->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloPutNodeIds);
  }

  kernelInfo += meshKernelProperties(mesh->N);
}

void meshOccaSetup3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo)
{
  meshOccaPopulateDeviceHex3D(mesh, newOptions, kernelInfo);
}
