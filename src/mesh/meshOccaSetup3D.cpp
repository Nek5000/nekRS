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

  //  printf("NinteriorElements = %d, NnotInternalElements = %d\n", Ninterior, NnotInterior);

  mesh->NinternalElements = Ninterior;
  mesh->NnotInternalElements = NnotInterior;
  if(Ninterior)
    mesh->o_internalElementIds    = mesh->device.malloc(Ninterior * sizeof(dlong),
                                                        internalElementIds);

  if(NnotInterior > 0)
    mesh->o_notInternalElementIds = mesh->device.malloc(NnotInterior * sizeof(dlong),
                                                        notInternalElementIds);

  // // OCCA allocate device memory (remember to go back for halo)
  // mesh->o_q =
  //   mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  // mesh->o_rhsq =
  //   mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  // mesh->o_resq =
  //   mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);

  //  reportMemoryUsage(mesh->device, "meshOccaSetup3D: before operators ");

  if(mesh->Nfaces == 4) {
    // build Dr, Ds, LIFT transposes
    dfloat* DrT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    dfloat* DsT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    dfloat* DtT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        DrT[n + m * mesh->Np] = mesh->Dr[n * mesh->Np + m];
        DsT[n + m * mesh->Np] = mesh->Ds[n * mesh->Np + m];
        DtT[n + m * mesh->Np] = mesh->Dt[n * mesh->Np + m];
      }

    // build Dr, Ds transposes
    dfloat* DrstT = (dfloat*) calloc(3 * mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        DrstT[n + m * mesh->Np] = mesh->Dr[n * mesh->Np + m];
        DrstT[n + m * mesh->Np + mesh->Np * mesh->Np] = mesh->Ds[n * mesh->Np + m];
        DrstT[n + m * mesh->Np + 2 * mesh->Np * mesh->Np] = mesh->Dt[n * mesh->Np + m];
      }

    dfloat* LIFTT = (dfloat*) calloc(mesh->Np * mesh->Nfaces * mesh->Nfp, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Nfaces * mesh->Nfp; ++m)
        LIFTT[n + m * mesh->Np] = mesh->LIFT[n * mesh->Nfp * mesh->Nfaces + m];

    // =============== BB operators [added by NC] ===============

    // deriv operators: transpose from row major to column major
    int* D0ids = (int*) calloc(mesh->Np * 4,sizeof(int));
    int* D1ids = (int*) calloc(mesh->Np * 4,sizeof(int));
    int* D2ids = (int*) calloc(mesh->Np * 4,sizeof(int));
    int* D3ids = (int*) calloc(mesh->Np * 4,sizeof(int));
    dfloat* Dvals = (dfloat*) calloc(mesh->Np * 4,sizeof(dfloat));

    int* L0ids = (int*) calloc(mesh->Nfp * 7,sizeof(int));
    dfloat* L0vals = (dfloat*) calloc(mesh->Nfp * 7,sizeof(dfloat)); // tridiag
    int* ELids = (int*) calloc(mesh->Np * mesh->max_EL_nnz,sizeof(int));
    dfloat* ELvals = (dfloat*) calloc(mesh->Np * mesh->max_EL_nnz,sizeof(dfloat));

    for (int i = 0; i < mesh->Np; ++i)
      for (int j = 0; j < 4; ++j) {
        D0ids[i + j * mesh->Np] = mesh->D0ids[j + i * 4];
        D1ids[i + j * mesh->Np] = mesh->D1ids[j + i * 4];
        D2ids[i + j * mesh->Np] = mesh->D2ids[j + i * 4];
        D3ids[i + j * mesh->Np] = mesh->D3ids[j + i * 4];
        Dvals[i + j * mesh->Np] = mesh->Dvals[j + i * 4];
      }

    for (int i = 0; i < mesh->Nfp; ++i)
      for (int j = 0; j < 7; ++j) {
        L0ids [i + j * mesh->Nfp] = mesh->L0ids [j + i * 7];
        L0vals[i + j * mesh->Nfp] = mesh->L0vals[j + i * 7];
      }

    for (int i = 0; i < mesh->Np; ++i)
      for (int j = 0; j < mesh->max_EL_nnz; ++j) {
        ELids [i + j * mesh->Np] = mesh->ELids [j + i * mesh->max_EL_nnz];
        ELvals[i + j * mesh->Np] = mesh->ELvals[j + i * mesh->max_EL_nnz];
      }
    // =============== end BB stuff =============================

    if(mesh->cubNp) {
      dfloat* cubDrWT = (dfloat*) calloc(mesh->cubNp * mesh->Np, sizeof(dfloat));
      dfloat* cubDsWT = (dfloat*) calloc(mesh->cubNp * mesh->Np, sizeof(dfloat));
      dfloat* cubDtWT = (dfloat*) calloc(mesh->cubNp * mesh->Np, sizeof(dfloat));
      dfloat* cubDrstWT = (dfloat*) calloc(3 * mesh->cubNp * mesh->Np, sizeof(dfloat));
      dfloat* cubProjectT = (dfloat*) calloc(mesh->cubNp * mesh->Np, sizeof(dfloat));
      dfloat* cubInterpT = (dfloat*) calloc(mesh->cubNp * mesh->Np, sizeof(dfloat));
      for(int n = 0; n < mesh->Np; ++n)
        for(int m = 0; m < mesh->cubNp; ++m) {
          cubDrWT[n + m * mesh->Np] = mesh->cubDrW[n * mesh->cubNp + m];
          cubDsWT[n + m * mesh->Np] = mesh->cubDsW[n * mesh->cubNp + m];
          cubDtWT[n + m * mesh->Np] = mesh->cubDtW[n * mesh->cubNp + m];

          cubDrstWT[n + m * mesh->Np] = mesh->cubDrW[n * mesh->cubNp + m];
          cubDrstWT[n + m * mesh->Np + mesh->cubNp * mesh->Np] = mesh->cubDsW[n * mesh->cubNp + m];
          cubDrstWT[n + m * mesh->Np + 2 * mesh->cubNp *
                    mesh->Np] = mesh->cubDtW[n * mesh->cubNp + m];

          cubProjectT[n + m * mesh->Np] = mesh->cubProject[n * mesh->cubNp + m];
          cubInterpT[m + n * mesh->cubNp] = mesh->cubInterp[m * mesh->Np + n];
        }

      mesh->o_cubInterpT =
        mesh->device.malloc(mesh->Np * mesh->cubNp * sizeof(dfloat),
                            cubInterpT);

      mesh->o_cubProjectT =
        mesh->device.malloc(mesh->Np * mesh->cubNp * sizeof(dfloat),
                            cubProjectT);

      mesh->o_cubDrWT =
        mesh->device.malloc(mesh->Np * mesh->cubNp * sizeof(dfloat),
                            cubDrWT);

      mesh->o_cubDsWT =
        mesh->device.malloc(mesh->Np * mesh->cubNp * sizeof(dfloat),
                            cubDsWT);

      mesh->o_cubDtWT =
        mesh->device.malloc(mesh->Np * mesh->cubNp * sizeof(dfloat),
                            cubDtWT);

      mesh->o_cubDWmatrices = mesh->device.malloc(3 * mesh->cubNp * mesh->Np * sizeof(dfloat),
                                                  cubDrstWT);
    }

    if(mesh->intNfp) {
      // build surface integration matrix transposes
      dfloat* intLIFTT = (dfloat*) calloc(mesh->Np * mesh->Nfaces * mesh->intNfp, sizeof(dfloat));
      dfloat* intInterpT =
        (dfloat*) calloc(mesh->Nfp * mesh->Nfaces * mesh->intNfp, sizeof(dfloat));
      for(int n = 0; n < mesh->Np; ++n)
        for(int m = 0; m < mesh->Nfaces * mesh->intNfp; ++m)
          intLIFTT[n + m * mesh->Np] = mesh->intLIFT[n * mesh->intNfp * mesh->Nfaces + m];
      for(int n = 0; n < mesh->intNfp * mesh->Nfaces; ++n)
        for(int m = 0; m < mesh->Nfp; ++m)
          intInterpT[n + m * mesh->Nfaces * mesh->intNfp] = mesh->intInterp[n * mesh->Nfp + m];

      mesh->o_intInterpT =
        mesh->device.malloc(mesh->Nfp * mesh->Nfaces * mesh->intNfp * sizeof(dfloat),
                            intInterpT);

      mesh->o_intLIFTT =
        mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->intNfp * sizeof(dfloat),
                            intLIFTT);

      // printf("Integration number of points: %d \n",mesh->intNfp);
      mesh->intx = (dfloat*) calloc(mesh->Nelements * mesh->Nfaces * mesh->intNfp, sizeof(dfloat));
      mesh->inty = (dfloat*) calloc(mesh->Nelements * mesh->Nfaces * mesh->intNfp, sizeof(dfloat));
      mesh->intz = (dfloat*) calloc(mesh->Nelements * mesh->Nfaces * mesh->intNfp, sizeof(dfloat));

      for(dlong e = 0; e < mesh->Nelements; ++e)
        for(int f = 0; f < mesh->Nfaces; ++f)
          for(int n = 0; n < mesh->intNfp; ++n) {
            dfloat ix = 0, iy = 0, iz = 0;
            for(int m = 0; m < mesh->Nfp; ++m) {
              dlong vid = mesh->vmapM[m + f * mesh->Nfp + e * mesh->Nfp * mesh->Nfaces];
              dfloat xm = mesh->x[vid];
              dfloat ym = mesh->y[vid];
              dfloat zm = mesh->z[vid];
              dfloat Inm = mesh->intInterp[m + n * mesh->Nfp + f * mesh->intNfp * mesh->Nfp]; // Fixed
              ix += Inm * xm;
              iy += Inm * ym;
              iz += Inm * zm;
            }
            dlong id = n + f * mesh->intNfp + e * mesh->Nfaces * mesh->intNfp;
            mesh->intx[id] = ix;
            mesh->inty[id] = iy;
            mesh->intz[id] = iz;
          }

      mesh->o_intx =
        mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->intNfp * sizeof(dfloat),
                            mesh->intx);

      mesh->o_inty =
        mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->intNfp * sizeof(dfloat),
                            mesh->inty);

      mesh->o_intz =
        mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->intNfp * sizeof(dfloat),
                            mesh->intz);
    }

    //    reportMemoryUsage(mesh->device, "meshOccaSetup3D: after operators and integration grids ");

    // =============== Bernstein-Bezier allocations [added by NC] ============

    // create packed BB indexes
    mesh->o_D0ids = mesh->device.malloc(mesh->Np * 4 * sizeof(int),D0ids);
    mesh->o_D1ids = mesh->device.malloc(mesh->Np * 4 * sizeof(int),D1ids);
    mesh->o_D2ids = mesh->device.malloc(mesh->Np * 4 * sizeof(int),D2ids);
    mesh->o_D3ids = mesh->device.malloc(mesh->Np * 4 * sizeof(int),D3ids);
    mesh->o_Dvals = mesh->device.malloc(mesh->Np * 4 * sizeof(dfloat),Dvals);

    unsigned char* packedDids = (unsigned char*) malloc(mesh->Np * 3 * 4 * sizeof(unsigned char));

    for(int n = 0; n < 4 * mesh->Np; ++n) {
      if(D1ids[n] < D0ids[n]) printf("bugger: D0id > D1id\n");
      if(D2ids[n] < D0ids[n]) printf("bugger: D0id > D2id\n");
      if(D3ids[n] < D0ids[n]) printf("bugger: D0id > D3id\n");
    }

    for(int n = 0; n < mesh->Np; ++n) {
      packedDids[n * 4 + 0] = D1ids[n + 0 * mesh->Np] - D0ids[n + 0 * mesh->Np];
      packedDids[n * 4 + 1] = D1ids[n + 1 * mesh->Np] - D0ids[n + 1 * mesh->Np];
      packedDids[n * 4 + 2] = D1ids[n + 2 * mesh->Np] - D0ids[n + 2 * mesh->Np];
      packedDids[n * 4 + 3] = D1ids[n + 3 * mesh->Np] - D0ids[n + 3 * mesh->Np];

      packedDids[4 * mesh->Np + n * 4 + 0] = D2ids[n + 0 * mesh->Np] - D0ids[n + 0 * mesh->Np];
      packedDids[4 * mesh->Np + n * 4 + 1] = D2ids[n + 1 * mesh->Np] - D0ids[n + 1 * mesh->Np];
      packedDids[4 * mesh->Np + n * 4 + 2] = D2ids[n + 2 * mesh->Np] - D0ids[n + 2 * mesh->Np];
      packedDids[4 * mesh->Np + n * 4 + 3] = D2ids[n + 3 * mesh->Np] - D0ids[n + 3 * mesh->Np];

      packedDids[8 * mesh->Np + n * 4 + 0] = D3ids[n + 0 * mesh->Np] - D0ids[n + 0 * mesh->Np];
      packedDids[8 * mesh->Np + n * 4 + 1] = D3ids[n + 1 * mesh->Np] - D0ids[n + 1 * mesh->Np];
      packedDids[8 * mesh->Np + n * 4 + 2] = D3ids[n + 2 * mesh->Np] - D0ids[n + 2 * mesh->Np];
      packedDids[8 * mesh->Np + n * 4 + 3] = D3ids[n + 3 * mesh->Np] - D0ids[n + 3 * mesh->Np];
    }

    mesh->o_packedDids = mesh->device.malloc(mesh->Np * 3 * 4 * sizeof(unsigned char),packedDids);

    mesh->o_L0ids  = mesh->device.malloc(mesh->Nfp * 7 * sizeof(int),L0ids);
    mesh->o_L0vals = mesh->device.malloc(mesh->Nfp * 7 * sizeof(dfloat),L0vals);
    mesh->o_ELids  = mesh->device.malloc(mesh->Np * mesh->max_EL_nnz * sizeof(int),ELids);
    mesh->o_ELvals = mesh->device.malloc(mesh->Np * mesh->max_EL_nnz * sizeof(dfloat),ELvals);
    // =============== end Bernstein-Bezier section [added by NC] ============

    //build element stiffness matrices
    dfloat* SrrT, * SrsT, * SrtT;
    dfloat* SsrT, * SssT, * SstT;
    dfloat* StrT, * StsT, * SttT;

    mesh->Srr = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Srs = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Srt = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Ssr = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Sss = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Sst = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Str = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Sts = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Stt = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    for (int n = 0; n < mesh->Np; n++)
      for (int m = 0; m < mesh->Np; m++)
        for (int k = 0; k < mesh->Np; k++)
          for (int l = 0; l < mesh->Np; l++) {
            mesh->Srr[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Srs[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
            mesh->Srt[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dt[m + k * mesh->Np];
            mesh->Ssr[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Sss[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
            mesh->Sst[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dt[m + k * mesh->Np];
            mesh->Str[m + n * mesh->Np] += mesh->Dt[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Sts[m + n * mesh->Np] += mesh->Dt[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
            mesh->Stt[m + n * mesh->Np] += mesh->Dt[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dt[m + k * mesh->Np];
          }
    SrrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SrsT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SrtT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SsrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SssT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SstT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    StrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    StsT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SttT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    for (int n = 0; n < mesh->Np; n++) {
      for (int m = 0; m < mesh->Np; m++) {
#if 0
        SrrT[m + n * mesh->Np] = mesh->Srr[n + m * mesh->Np];
        SrsT[m + n * mesh->Np] = mesh->Srs[n + m * mesh->Np];
        SrtT[m + n * mesh->Np] = mesh->Srt[n + m * mesh->Np];
        SsrT[m + n * mesh->Np] = mesh->Ssr[n + m * mesh->Np];
        SssT[m + n * mesh->Np] = mesh->Sss[n + m * mesh->Np];
        SstT[m + n * mesh->Np] = mesh->Sst[n + m * mesh->Np];
        StrT[m + n * mesh->Np] = mesh->Str[n + m * mesh->Np];
        StsT[m + n * mesh->Np] = mesh->Sts[n + m * mesh->Np];
        SttT[m + n * mesh->Np] = mesh->Stt[n + m * mesh->Np];
#else
        SrrT[m + n * mesh->Np] = mesh->Srr[n + m * mesh->Np];
        SrsT[m + n * mesh->Np] = mesh->Srs[n + m * mesh->Np] + mesh->Ssr[n + m * mesh->Np];
        SrtT[m + n * mesh->Np] = mesh->Srt[n + m * mesh->Np] + mesh->Str[n + m * mesh->Np];
        SssT[m + n * mesh->Np] = mesh->Sss[n + m * mesh->Np];
        SstT[m + n * mesh->Np] = mesh->Sst[n + m * mesh->Np] + mesh->Sts[n + m * mesh->Np];
        SttT[m + n * mesh->Np] = mesh->Stt[n + m * mesh->Np];
#endif
      }
    }

    dfloat* ST = (dfloat*) calloc(6 * mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        ST[n + m * mesh->Np + 0 * mesh->Np * mesh->Np] = mesh->Srr[n * mesh->Np + m];
        ST[n + m * mesh->Np + 1 * mesh->Np * mesh->Np] = mesh->Srs[n * mesh->Np + m] +
                                                         mesh->Ssr[n * mesh->Np + m];
        ST[n + m * mesh->Np + 2 * mesh->Np * mesh->Np] = mesh->Srt[n * mesh->Np + m] +
                                                         mesh->Str[n * mesh->Np + m];
        ST[n + m * mesh->Np + 3 * mesh->Np * mesh->Np] = mesh->Sss[n * mesh->Np + m];
        ST[n + m * mesh->Np + 4 * mesh->Np * mesh->Np] = mesh->Sst[n * mesh->Np + m] +
                                                         mesh->Sts[n * mesh->Np + m];
        ST[n + m * mesh->Np + 5 * mesh->Np * mesh->Np] = mesh->Stt[n * mesh->Np + m];
      }


    mesh->o_Dr = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), mesh->Dr);
    mesh->o_Ds = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), mesh->Ds);
    mesh->o_Dt = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), mesh->Dt);

    mesh->o_DrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), DrT);
    mesh->o_DsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), DsT);
    mesh->o_DtT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), DtT);

    mesh->o_Dmatrices = mesh->device.malloc(3 * mesh->Np * mesh->Np * sizeof(dfloat), DrstT);

    mesh->o_LIFT =
      mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->Nfp * sizeof(dfloat),
                          mesh->LIFT);

    mesh->o_LIFTT =
      mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->Nfp * sizeof(dfloat),
                          LIFTT);

    mesh->o_MM =
      mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                          mesh->MM);

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nvgeo * sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->Nsgeo * sizeof(dfloat),
                          mesh->sgeo);

    mesh->o_cubsgeo = mesh->o_sgeo; //dummy cubature geo factors

    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nggeo * sizeof(dfloat),
                          mesh->ggeo);

    mesh->o_cubvgeo =   mesh->device.malloc(sizeof(dfloat));// dummy

    mesh->o_SrrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrrT);
    mesh->o_SrsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrsT);
    mesh->o_SrtT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrtT);
    mesh->o_SsrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SsrT);
    mesh->o_SssT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SssT);
    mesh->o_SstT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SstT);
    mesh->o_StrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), StrT);
    mesh->o_StsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), StsT);
    mesh->o_SttT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SttT);

    mesh->o_Smatrices = mesh->device.malloc(6 * mesh->Np * mesh->Np * sizeof(dfloat), ST);

    free(DrstT);
    free(ST);
  } else if (mesh->Nverts == 8) {    // hardcoded for hexes
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
      mesh->device.malloc(1 * sizeof(dfloat)); // dummy

    //    reportMemoryUsage(mesh->device, "meshOccaSetup3D: before intX ");

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
      mesh->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat));
    mesh->invLMM = (dfloat*) calloc(mesh->Nelements * mesh->Np, sizeof(dfloat));
    mesh->o_invLMM =
      mesh->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat));

    mesh->o_MM =
      mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                          mesh->MM); //dummy

    mesh->o_D = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);

    mesh->o_DW = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->DW);

    mesh->o_Dmatrices = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);

    dfloat* DT = (dfloat*) calloc(mesh->Nq * mesh->Nq,sizeof(dfloat));
    for(int j = 0; j < mesh->Nq; ++j)
      for(int i = 0; i < mesh->Nq; ++i)
        DT[i * mesh->Nq + j] = mesh->D[j * mesh->Nq + i];

    mesh->o_Smatrices = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), DT); //dummy

    //    reportMemoryUsage(mesh->device, "meshOccaSetup3D: before geofactors ");

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Np * mesh->Nvgeo * sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->Nfp * mesh->Nsgeo * sizeof(dfloat),
                          mesh->sgeo);

    //    reportMemoryUsage(mesh->device, "meshOccaSetup3D: before vgeo,sgeo ");

    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(dfloat),
                          mesh->ggeo);

    mesh->o_cubvgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nvgeo * mesh->cubNp * sizeof(dfloat),
                          mesh->cubvgeo);

    mesh->o_cubsgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * mesh->Nsgeo *
                          sizeof(dfloat),
                          mesh->cubsgeo);

    mesh->o_cubggeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nggeo * mesh->cubNp * sizeof(dfloat),
                          mesh->cubggeo);

    mesh->o_cubInterpT =
      mesh->device.malloc(mesh->Nq * mesh->cubNq * sizeof(dfloat),
                          cubInterpT);

    mesh->o_cubProjectT =
      mesh->device.malloc(mesh->Nq * mesh->cubNq * sizeof(dfloat),
                          cubProjectT);

    mesh->o_cubDWT =
      mesh->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat),
                          cubDWT);

    mesh->o_cubD =
      mesh->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat),
                          mesh->cubD);

    mesh->o_cubDWmatrices = mesh->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat), cubDWT);

    // just neeeded to combine quad and hex cub kernels
    mesh->o_cubDiffInterpT = mesh->o_cubDWmatrices;

    //    reportMemoryUsage(mesh->device, "meshOccaSetup3D: after geofactors ");

    mesh->o_intx =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * sizeof(dfloat),
                          mesh->intx);

    mesh->o_inty =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * sizeof(dfloat),
                          mesh->inty);

    mesh->o_intz =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->cubNfp * sizeof(dfloat),
                          mesh->intz);

    mesh->o_intInterpT = mesh->device.malloc(mesh->cubNq * mesh->Nq * sizeof(dfloat));
    mesh->o_intInterpT.copyFrom(mesh->o_cubInterpT);

    mesh->o_intLIFTT = mesh->device.malloc(mesh->cubNq * mesh->Nq * sizeof(dfloat));
    mesh->o_intLIFTT.copyFrom(mesh->o_cubProjectT);

    //    reportMemoryUsage(mesh->device, "meshOccaSetup3D: after intX ");
  } else {
    printf("Nverts = %d: unknown element type!\n",mesh->Nverts);
  }

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                        mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                        mesh->vmapP);

  mesh->o_EToB =
    mesh->device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int),
                        mesh->EToB);

  mesh->o_x =
    mesh->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->x);

  mesh->o_y =
    mesh->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->y);

  mesh->o_z =
    mesh->device.malloc(mesh->Nelements * mesh->Np * sizeof(dfloat), mesh->z);

  if(mesh->totalHaloPairs > 0) {
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs * sizeof(dlong), mesh->haloElementList);

    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    //printf("mesh->Nfields = %d\n", mesh->Nfields);
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs * mesh->Np * mesh->Nfields * sizeof(dfloat));

    mesh->o_haloGetNodeIds =
      mesh->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloGetNodeIds);

    mesh->o_haloPutNodeIds =
      mesh->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloPutNodeIds);
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

  kernelInfo["defines/" "p_max_EL_nnz"] = mesh->max_EL_nnz; // for Bernstein Bezier lift

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

#if 0
  // TW: these should be defined at the solver setup
  int NblockV = 256 / mesh->Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"] = NblockV;

  int NblockS = 256 / maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"] = NblockS;
#endif

  kernelInfo["defines/" "p_Lambda2"] = 0.5f;

  kernelInfo["defines/" "p_cubNq"] = mesh->cubNq;
  kernelInfo["defines/" "p_cubNfp"] = mesh->cubNfp;
  kernelInfo["defines/" "p_cubNp"] = mesh->cubNp;
  kernelInfo["defines/" "p_intNfp"] = mesh->intNfp;
  kernelInfo["defines/" "p_intNfpNfaces"] = mesh->intNfp * mesh->Nfaces;

  if(sizeof(dfloat) == 4) {
    kernelInfo["defines/" "dfloat"] = "float";
    kernelInfo["defines/" "dfloat4"] = "float4";
    kernelInfo["defines/" "dfloat8"] = "float8";
  }
  if(sizeof(dfloat) == 8) {
    kernelInfo["defines/" "dfloat"] = "double";
    kernelInfo["defines/" "dfloat4"] = "double4";
    kernelInfo["defines/" "dfloat8"] = "double8";
  }

  if(sizeof(dlong) == 4)
    kernelInfo["defines/" "dlong"] = "int";
  if(sizeof(dlong) == 8)
    kernelInfo["defines/" "dlong"] = "long long int";

  if(mesh->device.mode() == "CUDA") { // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += "--ftz=true ";
    kernelInfo["compiler_flags"] += "--prec-div=false ";
    kernelInfo["compiler_flags"] += "--prec-sqrt=false ";
    kernelInfo["compiler_flags"] += "--use_fast_math ";
    kernelInfo["compiler_flags"] += "--fmad=true "; // compiler option for cuda
    //kernelInfo["compiler_flags"] += "-Xptxas -dlcm=ca";
  }

  if(mesh->device.mode() == "OpenCL") { // add backend compiler optimization for OPENCL
    kernelInfo["compiler_flags"] += " -cl-std=CL2.0 ";
    kernelInfo["compiler_flags"] += " -cl-strict-aliasing ";
    kernelInfo["compiler_flags"] += " -cl-mad-enable ";
    kernelInfo["compiler_flags"] += " -cl-no-signed-zeros ";
    kernelInfo["compiler_flags"] += " -cl-unsafe-math-optimizations ";
    kernelInfo["compiler_flags"] += " -cl-fast-relaxed-math ";
  }

  if(mesh->device.mode() == "HIP") { // add backend compiler optimization for HIP
    kernelInfo["compiler_flags"] += " -O3 ";
    kernelInfo["compiler_flags"] += " -ffp-contract=fast ";
    // kernelInfo["compiler_flags"] += " -funsafe-math-optimizations ";
    // kernelInfo["compiler_flags"] += " -ffast-math ";
  }

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
  //make seperate stream for halo exchange
  mesh->defaultStream = mesh->device.getStream();
  mesh->dataStream = mesh->device.createStream();
  mesh->computeStream = mesh->device.createStream();
  mesh->device.setStream(mesh->defaultStream);

  meshOccaPopulateDevice3D(mesh, newOptions, kernelInfo);
}

void meshOccaCloneDevice(mesh_t* donorMesh, mesh_t* mesh)
{
  mesh->device = donorMesh->device;

  mesh->defaultStream = donorMesh->defaultStream;
  mesh->dataStream = donorMesh->dataStream;
  mesh->computeStream = donorMesh->computeStream;
}
