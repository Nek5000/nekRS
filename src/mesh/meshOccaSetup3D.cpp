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

void reportMemoryUsage(occa::device &device, const char* mess)
{
  size_t bytes = device.memoryAllocated();

  printf("%s: bytes allocated = %lu\n", mess, bytes);
}

void meshOccaPopulateDevice3D(mesh3D* mesh, setupAide &newOptions, occa::properties &kernelInfo)
{
  
  // find elements that have all neighbors on this process
  dlong* internalElementIds = (dlong*) calloc(mesh->Nelements, sizeof(dlong));
  dlong* notInternalElementIds = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  dlong Ninterior = 0, NnotInterior = 0;
  for(dlong e = 0; e < mesh->Nelements; ++e) {
    int flag = 0;
    for(int f = 0; f < mesh->Nfaces; ++f)
      if(mesh->EToP[e * mesh->Nfaces + f] != -1)
        flag = 1;
    if(!flag)
      internalElementIds[Ninterior++] = e;
    else
      notInternalElementIds[NnotInterior++] = e;
  }

  mesh->o_elementInfo = platform->device.malloc(mesh->Nelements * sizeof(dlong), 
		                            mesh->elementInfo);
 
  mesh->NinternalElements = Ninterior;
  mesh->NnotInternalElements = NnotInterior;
  if(Ninterior)
    mesh->o_internalElementIds    = platform->device.malloc(Ninterior * sizeof(dlong),
                                                        internalElementIds);

  if(NnotInterior > 0)
    mesh->o_notInternalElementIds = platform->device.malloc(NnotInterior * sizeof(dlong),
                                                        notInternalElementIds);

  if (mesh->Nverts == 8) {    // hardcoded for hexes
    //lumped mass matrix
    mesh->MM = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    for (int k = 0; k < mesh->Nq; k++)
      for (int j = 0; j < mesh->Nq; j++)
        for (int i = 0; i < mesh->Nq; i++) {
          int n = i + j * mesh->Nq + k * mesh->Nq * mesh->Nq;
          mesh->MM[n + n * mesh->Np] = mesh->gllw[i] * mesh->gllw[j] * mesh->gllw[k];
        }

    mesh->LIFT = (dfloat*) calloc(mesh->Np * mesh->Nfaces * mesh->Nfp, sizeof(dfloat));

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

    dfloat* LIFTT = (dfloat*) calloc(mesh->Np * mesh->Nfaces * mesh->Nfp, sizeof(dfloat));

    mesh->o_LIFTT =
      platform->device.malloc(1 * sizeof(dfloat)); // dummy

    mesh->intx = (dfloat*) calloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp, sizeof(dfloat));
    mesh->inty = (dfloat*) calloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp, sizeof(dfloat));
    mesh->intz = (dfloat*) calloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp, sizeof(dfloat));

    dfloat* ix = (dfloat*) calloc(mesh->cubNq * mesh->Nq,sizeof(dfloat));
    dfloat* iy = (dfloat*) calloc(mesh->cubNq * mesh->Nq,sizeof(dfloat));
    dfloat* iz = (dfloat*) calloc(mesh->cubNq * mesh->Nq,sizeof(dfloat));
    for(dlong e = 0; e < mesh->Nelements; ++e)
      for(int f = 0; f < mesh->Nfaces; ++f) {
        //interpolate in i
        for(int ny = 0; ny < mesh->Nq; ++ny)
          for(int nx = 0; nx < mesh->cubNq; ++nx) {
            ix[nx + mesh->cubNq * ny] = 0;
            iy[nx + mesh->cubNq * ny] = 0;
            iz[nx + mesh->cubNq * ny] = 0;

            for(int m = 0; m < mesh->Nq; ++m) {
              dlong vid = m + ny * mesh->Nq + f * mesh->Nfp + e * mesh->Nfp * mesh->Nfaces;
              dlong idM = mesh->vmapM[vid];

              dfloat xm = mesh->x[idM];
              dfloat ym = mesh->y[idM];
              dfloat zm = mesh->z[idM];

              dfloat Inm = mesh->cubInterp[m + nx * mesh->Nq];
              ix[nx + mesh->cubNq * ny] += Inm * xm;
              iy[nx + mesh->cubNq * ny] += Inm * ym;
              iz[nx + mesh->cubNq * ny] += Inm * zm;
            }
          }

        //interpolate in j and store
        for(int ny = 0; ny < mesh->cubNq; ++ny)
          for(int nx = 0; nx < mesh->cubNq; ++nx) {
            dfloat x = 0.0, y = 0.0, z = 0.0;

            for(int m = 0; m < mesh->Nq; ++m) {
              dfloat xm = ix[nx + m * mesh->cubNq];
              dfloat ym = iy[nx + m * mesh->cubNq];
              dfloat zm = iz[nx + m * mesh->cubNq];

              dfloat Inm = mesh->cubInterp[m + ny * mesh->Nq];
              x += Inm * xm;
              y += Inm * ym;
              z += Inm * zm;
            }

            dlong id = nx + ny * mesh->cubNq + f * mesh->cubNfp + e * mesh->Nfaces * mesh->cubNfp;
            mesh->intx[id] = x;
            mesh->inty[id] = y;
            mesh->intz[id] = z;
          }
      }
    free(ix);
    free(iy);
    free(iz);

    mesh->LMM = (dfloat*) calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));
    mesh->o_LMM =
      platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat));
    mesh->invLMM = (dfloat*) calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));
    mesh->o_invLMM =
      platform->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat));

    mesh->o_MM =
      platform->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                          mesh->MM); //dummy

    mesh->o_D = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);

    mesh->o_DW = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->DW);

    mesh->o_Dmatrices = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);

    dfloat* DT = (dfloat*) calloc(mesh->Nq * mesh->Nq,sizeof(dfloat));
    for(int j = 0; j < mesh->Nq; ++j)
      for(int i = 0; i < mesh->Nq; ++i)
        DT[i * mesh->Nq + j] = mesh->D[j * mesh->Nq + i];

    mesh->o_Smatrices = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), DT); //dummy

    mesh->o_vgeo =
      platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nvgeo * sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      platform->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->Nfp * mesh->Nsgeo * sizeof(dfloat),
                          mesh->sgeo);

    mesh->o_ggeo =
      platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(dfloat),
                          mesh->ggeo);

    mesh->o_cubvgeo =
      platform->device.malloc(mesh->Nelements * mesh->Nvgeo * mesh->cubNp * sizeof(dfloat),
                          mesh->cubvgeo);

    /*
    mesh->o_cubsgeo =
      platform->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * mesh->Nsgeo *
                          sizeof(dfloat),
                          mesh->cubsgeo);
    mesh->o_cubggeo =
      platform->device.malloc(mesh->Nelements * mesh->Nggeo * mesh->cubNp * sizeof(dfloat),
                          mesh->cubggeo);
    */		  

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

    mesh->o_cubDWmatrices = platform->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat), cubDWT);

    // just neeeded to combine quad and hex cub kernels
    mesh->o_cubDiffInterpT = mesh->o_cubDWmatrices;

    mesh->o_intx =
      platform->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * sizeof(dfloat),
                          mesh->intx);

    mesh->o_inty =
      platform->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * sizeof(dfloat),
                          mesh->inty);

    mesh->o_intz =
      platform->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * sizeof(dfloat),
                          mesh->intz);

    mesh->o_intInterpT = platform->device.malloc(mesh->cubNq * mesh->Nq * sizeof(dfloat));
    mesh->o_intInterpT.copyFrom(mesh->o_cubInterpT);

    mesh->o_intLIFTT = platform->device.malloc(mesh->cubNq * mesh->Nq * sizeof(dfloat));
    mesh->o_intLIFTT.copyFrom(mesh->o_cubProjectT);

  } else {
    printf("Nverts = %d: unknown element type!\n",mesh->Nverts);
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

  if(mesh->totalHaloPairs > 0) {
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      platform->device.malloc(mesh->totalHaloPairs * sizeof(dlong), mesh->haloElementList);

    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    //printf("mesh->Nfields = %d\n", mesh->Nfields);
    mesh->o_haloBuffer =
      platform->device.malloc(mesh->totalHaloPairs * mesh->Np * mesh->Nfields * sizeof(dfloat));

    mesh->o_haloGetNodeIds =
      platform->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloGetNodeIds);

    mesh->o_haloPutNodeIds =
      platform->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloPutNodeIds);
  }

  kernelInfo["defines/" "p_dim"] = 3;
  kernelInfo["defines/" "p_Nfields"] = mesh->Nfields;
  kernelInfo["defines/" "p_N"] = mesh->N;
  kernelInfo["defines/" "p_Nq"] = mesh->N + 1;
  kernelInfo["defines/" "p_Np"] = mesh->Np;
  kernelInfo["defines/" "p_Nfp"] = mesh->Nfp;
  kernelInfo["defines/" "p_Nfaces"] = mesh->Nfaces;
  kernelInfo["defines/" "p_NfacesNfp"] = mesh->Nfp * mesh->Nfaces;
  kernelInfo["defines/" "p_Nvgeo"] = mesh->Nvgeo;
  kernelInfo["defines/" "p_Nsgeo"] = mesh->Nsgeo;
  kernelInfo["defines/" "p_Nggeo"] = mesh->Nggeo;

  kernelInfo["defines/" "p_NXID"] = NXID;
  kernelInfo["defines/" "p_NYID"] = NYID;
  kernelInfo["defines/" "p_NZID"] = NZID;
  kernelInfo["defines/" "p_SJID"] = SJID;
  kernelInfo["defines/" "p_IJID"] = IJID;
  kernelInfo["defines/" "p_IHID"] = IHID;
  kernelInfo["defines/" "p_WSJID"] = WSJID;
  kernelInfo["defines/" "p_WIJID"] = WIJID;
  kernelInfo["defines/" "p_STXID"] = STXID;
  kernelInfo["defines/" "p_STYID"] = STYID;
  kernelInfo["defines/" "p_STZID"] = STZID;
  kernelInfo["defines/" "p_SBXID"] = SBXID;
  kernelInfo["defines/" "p_SBYID"] = SBYID;
  kernelInfo["defines/" "p_SBZID"] = SBZID;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp * mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"] = maxNodes;

  kernelInfo["defines/" "p_Lambda2"] = 0.5f;

  kernelInfo["defines/" "p_cubNq"] = mesh->cubNq;
  kernelInfo["defines/" "p_cubNfp"] = mesh->cubNfp;
  kernelInfo["defines/" "p_cubNp"] = mesh->cubNp;
  kernelInfo["defines/" "p_intNfp"] = mesh->intNfp;
  kernelInfo["defines/" "p_intNfpNfaces"] = mesh->intNfp * mesh->Nfaces;

  kernelInfo["defines/" "p_G00ID"] = G00ID;
  kernelInfo["defines/" "p_G01ID"] = G01ID;
  kernelInfo["defines/" "p_G02ID"] = G02ID;
  kernelInfo["defines/" "p_G11ID"] = G11ID;
  kernelInfo["defines/" "p_G12ID"] = G12ID;
  kernelInfo["defines/" "p_G22ID"] = G22ID;
  kernelInfo["defines/" "p_GWJID"] = GWJID;

  kernelInfo["defines/" "p_RXID"] = RXID;
  kernelInfo["defines/" "p_SXID"] = SXID;
  kernelInfo["defines/" "p_TXID"] = TXID;

  kernelInfo["defines/" "p_RYID"] = RYID;
  kernelInfo["defines/" "p_SYID"] = SYID;
  kernelInfo["defines/" "p_TYID"] = TYID;

  kernelInfo["defines/" "p_RZID"] = RZID;
  kernelInfo["defines/" "p_SZID"] = SZID;
  kernelInfo["defines/" "p_TZID"] = TZID;

  kernelInfo["defines/" "p_JID"] = JID;
  kernelInfo["defines/" "p_JWID"] = JWID;
  kernelInfo["defines/" "p_IJWID"] = IJWID;
}

void meshOccaSetup3D(mesh3D* mesh, setupAide &newOptions, occa::properties &kernelInfo)
{
  meshOccaPopulateDevice3D(mesh, newOptions, kernelInfo);
}