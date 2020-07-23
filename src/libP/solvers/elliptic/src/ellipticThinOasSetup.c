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

typedef struct
{
  dlong element, elementN;
  int face, faceN, rankN;
}oasFacePair_t;

/* comparison function that orders halo element/face
   based on their indexes */
int compareOasHaloFaces(const void* a,
                        const void* b)
{
  oasFacePair_t* fa = (oasFacePair_t*) a;
  oasFacePair_t* fb = (oasFacePair_t*) b;

  if(fa->rankN < fb->rankN) return -1;
  if(fa->rankN > fb->rankN) return +1;

  if(fa->elementN < fb->elementN) return -1;
  if(fa->elementN > fb->elementN) return +1;

  if(fa->faceN < fb->faceN) return -1;
  if(fa->faceN > fb->faceN) return +1;

  return 0;
}

void ellipticThinOasSetup(elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;
  const dfloat lambda = elliptic->lambda[0];

  int size = mesh->size;
  int rank = mesh->rank;

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;

  // --------------------------------------------------------------------------------------------------------------
  // 1. construct array of unique discontinuous node numbers including halo elements
  // --------------------------------------------------------------------------------------------------------------

  // number of degrees of freedom on this rank
  hlong Nnum = Np * Nelements;

  // create a global numbering system
  hlong* globalStarts = (hlong*) calloc(mesh->size + 1, sizeof(hlong));
  hlong* globalIds = (hlong*) calloc((Nelements + mesh->totalHaloPairs) * Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts + 1, 1, MPI_HLONG, mesh->comm);
  for(int r = 0; r < mesh->size; ++r)
    globalStarts[r + 1] = globalStarts[r] + globalStarts[r + 1];

  /* so find number of elements on each rank */
  dlong* rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong* rankStarts = (hlong*) calloc(mesh->size + 1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG,
                rankNelements, 1, MPI_DLONG, mesh->comm);

  //find offsets
  for(int r = 0; r < mesh->size; ++r)
    rankStarts[r + 1] = rankStarts[r] + rankNelements[r];
  //use the offsets to set a global id
  for (dlong e = 0; e < Nelements; e++)
    for (int n = 0; n < Np; n++)
      globalIds[e * Np + n] = n + (e + rankStarts[rank]) * Np;

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong* idSendBuffer = (hlong*) calloc(Np * mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, Np * sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements * Np);
    free(idSendBuffer);
  }

  // --------------------------------------------------------------------------------------------------------------
  // 2. construct Nelements*(N+3)^3 array
  // --------------------------------------------------------------------------------------------------------------

  int oasNq = mesh->Nq + 2;
  int oasNp = oasNq * oasNq * oasNq;

  hlong NtotalOgs = oasNp * mesh->Nelements;
  hlong* nodeIds = (hlong*) calloc(NtotalOgs, sizeof(hlong));

  int* offsets = (int*) calloc(mesh->Nfaces, sizeof(int));
  offsets[0] = +mesh->Nq * mesh->Nq;
  offsets[1] = +mesh->Nq;
  offsets[2] = -1;
  offsets[3] = -mesh->Nq;
  offsets[4] = +1;
  offsets[5] = -mesh->Nq * mesh->Nq;

  hlong cnt = 0;

  elliptic->oasMapP = (dlong*) calloc(mesh->Nelements * mesh->Nfaces * mesh->Nfp, sizeof(dfloat));

  for(hlong e = 0; e < mesh->Nelements; ++e) {
    for(int n = 0; n < oasNp; ++n) {
      nodeIds[e * oasNp + n] = -1; // gs will ignore this ?
    }
    for(int k = 0; k < mesh->Nq; ++k)
      for(int j = 0; j < mesh->Nq; ++j)
        for(int i = 0; i < mesh->Nq; ++i) {
          hlong id    = i + j * mesh->Nq + k * mesh->Nq * mesh->Nq;
          hlong idOas = i + 1 + (j + 1) * oasNq + (k + 1) * oasNq * oasNq;

          nodeIds[e * oasNp + idOas] = globalIds[e * mesh->Np + id];
        }

    for(int fM = 0; fM < mesh->Nfaces; ++fM) {
      int fP = mesh->EToF[e * mesh->Nfaces + fM];
      if(fP < 0) fP = fM;
      for(int n = 0; n < mesh->Nfp; ++n) {
        // usual trace node locations
        hlong idP = mesh->vmapP[e * mesh->Nfp * mesh->Nfaces + n];

        elliptic->oasMapP[cnt++] = idP + offsets[fP]; // double check for bcs
      }
    }
  }

  // --------------------------------------------------------------------------------------------------------------
  // 3. create thin halo exchange [ to use regular halo get and put ]
  // --------------------------------------------------------------------------------------------------------------

  // create a list of element/faces with halo neighbor
  oasFacePair_t* haloElements =
    (oasFacePair_t*) calloc(mesh->totalHaloPairs, sizeof(oasFacePair_t));

  cnt = 0;
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int f = 0; f < mesh->Nfaces; ++f) {
      dlong ef = e * mesh->Nfaces + f;
      if(mesh->EToP[ef] != -1) {
        haloElements[cnt].element  = e;
        haloElements[cnt].face     = f;
        haloElements[cnt].elementN = mesh->EToE[ef];
        haloElements[cnt].faceN    = mesh->EToF[ef];
        haloElements[cnt].rankN    = mesh->EToP[ef];
        ++cnt;
      }
    }

  // sort the face pairs in order the destination requires
  qsort(haloElements, mesh->totalHaloPairs, sizeof(oasFacePair_t), compareOasHaloFaces);

  // record the outgoing order for elements
  elliptic->oasHaloElementList = (dlong*) calloc(mesh->totalHaloPairs, sizeof(dlong));
  for(dlong i = 0; i < mesh->totalHaloPairs; ++i) {
    dlong e = haloElements[i].element;
    elliptic->oasHaloElementList[i] = e;
  }

  // record the outgoing node ids for trace nodes
  elliptic->oasHaloGetNodeIds = (dlong*) calloc(mesh->totalHaloPairs * mesh->Nfp, sizeof(dlong));
  elliptic->oasHaloPutNodeIds = (dlong*) calloc(mesh->totalHaloPairs * mesh->Nfp, sizeof(dlong));

  cnt = 0;
  for(dlong i = 0; i < mesh->totalHaloPairs; ++i) {
    dlong eM = haloElements[i].element;
    int fM = haloElements[i].face;
    int fP = haloElements[i].faceN;
    for(int n = 0; n < mesh->Nfp; ++n) {
      elliptic->oasHaloGetNodeIds[cnt] = eM * mesh->Np + mesh->faceNodes[fM * mesh->Nfp + n] +
                                         offsets[fM];
      ++cnt;
    }
  }

  // now arrange for incoming nodes
  cnt = mesh->Nelements;
  dlong ncnt = 0;
  for(int r = 0; r < size; ++r)
    for(dlong e = 0; e < mesh->Nelements; ++e)
      for(int f = 0; f < mesh->Nfaces; ++f) {
        dlong ef = e * mesh->Nfaces + f;
        if(mesh->EToP[ef] == r) {
          mesh->EToE[ef] = cnt;
          int fP = mesh->EToF[ef];
          for(int n = 0; n < mesh->Nfp; ++n) {
            elliptic->oasHaloPutNodeIds[ncnt] = cnt * mesh->Np +
                                                mesh->faceNodes[fP * mesh->Nfp + n] + offsets[fP];
            ++ncnt;
          }
          ++cnt; // next halo element
        }
      }


  // --------------------------------------------------------------------------------------------------------------
  // 4. construct diagonal scaling for fast diagonal inverse
  // --------------------------------------------------------------------------------------------------------------

  // hack estimate for Jacobian scaling
  dfloat* diagInvOp = (dfloat*) calloc(oasNp * mesh->Nelements, sizeof(dfloat));

  for(dlong e = 0; e < mesh->Nelements; ++e) {
    // S = Jabc*(wa*wb*wc*lambda + wb*wc*Da'*wa*Da + wa*wc*Db'*wb*Db + wa*wb*Dc'*wc*Dc)
    // S = Jabc*wa*wb*wc*(lambda*I+1/wa*Da'*wa*Da + 1/wb*Db'*wb*Db + 1/wc*Dc'*wc*Dc)

    dfloat Jhrinv2 = 0, Jhsinv2 = 0, Jhtinv2 = 0, J = 0;

    for(int n = 0; n < mesh->Np; ++n) {
      dlong base = mesh->Nvgeo * mesh->Np * e + n;

      dfloat Jn  = mesh->vgeo[base + JID * mesh->Np];

      dfloat rx = mesh->vgeo[base + RXID * mesh->Np];
      dfloat ry = mesh->vgeo[base + RYID * mesh->Np];
      dfloat rz = mesh->vgeo[base + RZID * mesh->Np];

      dfloat sx = mesh->vgeo[base + SXID * mesh->Np];
      dfloat sy = mesh->vgeo[base + SYID * mesh->Np];
      dfloat sz = mesh->vgeo[base + SZID * mesh->Np];

      dfloat tx = mesh->vgeo[base + TXID * mesh->Np];
      dfloat ty = mesh->vgeo[base + TYID * mesh->Np];
      dfloat tz = mesh->vgeo[base + TZID * mesh->Np];

      dfloat gradr2 = rx * rx + ry * ry + rz * rz;
      dfloat grads2 = sx * sx + sy * sy + sz * sz;
      dfloat gradt2 = tx * tx + ty * ty + tz * tz;

      J = mymax(J, Jn);

      Jhrinv2 = mymax(Jhrinv2, gradr2);
      Jhsinv2 = mymax(Jhsinv2, grads2);
      Jhtinv2 = mymax(Jhtinv2, gradt2);
    }

    for(int k = 0; k < oasNq; ++k)
      for(int j = 0; j < oasNq; ++j)
        for(int i = 0; i < oasNq; ++i) {
          dlong pid = i + j * oasNq + k * oasNq * oasNq + e * oasNp;

          diagInvOp[pid] =
            1. / (J * lambda +
                  Jhrinv2 * mesh->oasDiagOp[i] +
                  Jhsinv2 * mesh->oasDiagOp[j] +
                  Jhtinv2 * mesh->oasDiagOp[k]);

          //	  printf("diagInvOp[%d] = %lf\n", pid, diagInvOp[pid]);
        }
  }

  // --------------------------------------------------------------------------------------------------------------
  // 5. construct ogs for the Oas patches end game gather-scatter
  // --------------------------------------------------------------------------------------------------------------
  elliptic->oasOgs = ogsSetup(NtotalOgs, nodeIds, mesh->comm, verbose, mesh->device);

  // --------------------------------------------------------------------------------------------------------------
  // 6. populate device arrays
  // --------------------------------------------------------------------------------------------------------------

  elliptic->o_oasMapP =
    mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->Nfp * sizeof(dfloat),
                        elliptic->oasMapP);
  elliptic->o_oasHaloElementList = mesh->device.malloc(mesh->totalHaloPairs * sizeof(dlong),
                                                       elliptic->o_oasHaloElementList);
  elliptic->o_oasHaloGetNodeIds =
    mesh->device.malloc(mesh->totalHaloPairs * mesh->Nfp * sizeof(dlong),
                        elliptic->oasHaloGetNodeIds);
  elliptic->o_oasHaloPutNodeIds =
    mesh->device.malloc(mesh->totalHaloPairs * mesh->Nfp * sizeof(dlong),
                        elliptic->oasHaloPutNodeIds);

  // operators
  elliptic->o_oasForward  = mesh->device.malloc(oasNq * oasNq * sizeof(dfloat), mesh->oasForward);
  elliptic->o_oasBack     = mesh->device.malloc(oasNq * oasNq * sizeof(dfloat), mesh->oasBack);

  elliptic->o_oasDiagInvOp =
    mesh->device.malloc(oasNp * mesh->Nelements * sizeof(dfloat), diagInvOp);

  elliptic->oasNq = oasNq;
  elliptic->oasNp = oasNp;

  // --------------------------------------------------------------------------------------------------------------
  // 7. free local stuff
  // --------------------------------------------------------------------------------------------------------------
  free(globalStarts);
  free(globalIds);
  free(rankNelements);
  free(rankStarts);
  free(nodeIds);
  free(offsets);
  free(haloElements);
  free(diagInvOp);
}
