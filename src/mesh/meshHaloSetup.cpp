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
#include "platform.hpp"
#include "mesh.h"

struct facePair_t
{
  dlong element, elementN;
  int face, faceN, rankN;
};

/* comparison function that orders halo element/face
   based on their indexes */
int compareHaloFaces(const void* a,
                     const void* b)
{
  facePair_t* fa = (facePair_t*) a;
  facePair_t* fb = (facePair_t*) b;

  if(fa->rankN < fb->rankN) return -1;
  if(fa->rankN > fb->rankN) return +1;

  if(fa->elementN < fb->elementN) return -1;
  if(fa->elementN > fb->elementN) return +1;

  if(fa->faceN < fb->faceN) return -1;
  if(fa->faceN > fb->faceN) return +1;

  return 0;
}

// set up halo infomation for inter-processor MPI
// exchange of trace nodes
void meshHaloSetup(mesh_t* mesh)
{
  // MPI info
  int rank, size;
  rank = platform->comm.mpiRank;
  size = platform->comm.mpiCommSize;

  mesh->totalHaloPairs = 0;
  return;

  // non-blocking MPI isend/irecv requests (used in meshHaloExchange)
  mesh->haloSendRequests = calloc(size, sizeof(MPI_Request));
  mesh->haloRecvRequests = calloc(size, sizeof(MPI_Request));

  // count number of halo element nodes to swap
  mesh->NhaloPairs = (int*) calloc(size, sizeof(int));
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int f = 0; f < mesh->Nfaces; ++f) {
      int r = mesh->EToP[e * mesh->Nfaces + f]; // rank of neighbor
      if(r != -1) {
        mesh->totalHaloPairs += 1;
        mesh->NhaloPairs[r] += 1;
      }
    }

  // count number of MPI messages in halo exchange
  mesh->NhaloMessages = 0;
  for(int r = 0; r < size; ++r)
    if(mesh->NhaloPairs[r])
      ++mesh->NhaloMessages;

  // create a list of element/faces with halo neighbor
  facePair_t* haloElements =
    (facePair_t*) calloc(mesh->totalHaloPairs, sizeof(facePair_t));

  dlong cnt = 0;
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
  qsort(haloElements, mesh->totalHaloPairs, sizeof(facePair_t), compareHaloFaces);

  // record the outgoing order for elements
  mesh->haloElementList = (dlong*) calloc(mesh->totalHaloPairs, sizeof(dlong));
  for(dlong i = 0; i < mesh->totalHaloPairs; ++i) {
    dlong e = haloElements[i].element;
    mesh->haloElementList[i] = e;
  }

  // record the outgoing node ids for trace nodes
  mesh->haloGetNodeIds = (dlong*) calloc(mesh->totalHaloPairs * mesh->Nfp, sizeof(dlong));
  mesh->haloPutNodeIds = (dlong*) calloc(mesh->totalHaloPairs * mesh->Nfp, sizeof(dlong));

  cnt = 0;
  for(dlong i = 0; i < mesh->totalHaloPairs; ++i) {
    dlong eM = haloElements[i].element;
    int fM = haloElements[i].face;
    int fP = haloElements[i].faceN;
    for(int n = 0; n < mesh->Nfp; ++n) {
      mesh->haloGetNodeIds[cnt] = eM * mesh->Np + mesh->faceNodes[fM * mesh->Nfp + n];
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
            mesh->haloPutNodeIds[ncnt] = cnt * mesh->Np + mesh->faceNodes[fP * mesh->Nfp + n];
            ++ncnt;
          }
          ++cnt; // next halo element
        }
      }

  free(haloElements);
}

void meshHaloPhysicalNodes(mesh_t* mesh)
{
  // create halo extension for x,y arrays
  dlong totalHaloNodes = mesh->totalHaloPairs * mesh->Np;
  dlong localNodes     = mesh->Nelements * mesh->Np;

  if(totalHaloNodes) {
    // temporary send buffer
    dfloat* sendBuffer = (dfloat*) calloc(totalHaloNodes, sizeof(dfloat));
 
    // send halo data and recv into extended part of arrays
    meshHaloExchange(mesh, mesh->Np * sizeof(dfloat), mesh->x, sendBuffer, mesh->x + localNodes);
    meshHaloExchange(mesh, mesh->Np * sizeof(dfloat), mesh->y, sendBuffer, mesh->y + localNodes);
    if(mesh->dim == 3)
      meshHaloExchange(mesh, mesh->Np * sizeof(dfloat), mesh->z, sendBuffer, mesh->z + localNodes);
 
    free(sendBuffer);
  }
}
