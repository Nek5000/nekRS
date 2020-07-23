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
#include <stddef.h>

#include "mesh.h"

void meshParallelGatherScatterSetup(mesh_t* mesh,
                                    dlong N,
                                    hlong* globalIds,
                                    MPI_Comm &comm,
                                    int verbose)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  mesh->ogs = ogsSetup(N, globalIds, comm, verbose, mesh->device);

  //use the gs to find what nodes are local to this rank
  int* minRank = (int*) calloc(N,sizeof(int));
  int* maxRank = (int*) calloc(N,sizeof(int));
  for (dlong i = 0; i < N; i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
  }

  ogsGatherScatter(minRank, ogsInt, ogsMin, mesh->ogs); //minRank[n] contains the smallest rank taking part in the gather of node n
  ogsGatherScatter(maxRank, ogsInt, ogsMax, mesh->ogs); //maxRank[n] contains the largest rank taking part in the gather of node n

  // count elements that contribute to global C0 gather-scatter
  dlong globalCount = 0;
  dlong localCount = 0;
  for(dlong e = 0; e < mesh->Nelements; ++e) {
    int isHalo = 0;
    for(int n = 0; n < mesh->Np; ++n) {
      dlong id = e * mesh->Np + n;
      if ((minRank[id] != rank) || (maxRank[id] != rank)) {
        isHalo = 1;
        break;
      }
    }
    globalCount += isHalo;
    localCount += 1 - isHalo;
  }

  mesh->globalGatherElementList = (dlong*) calloc(globalCount, sizeof(dlong));
  mesh->localGatherElementList  = (dlong*) calloc(localCount, sizeof(dlong));

  globalCount = 0;
  localCount = 0;

  for(dlong e = 0; e < mesh->Nelements; ++e) {
    int isHalo = 0;
    for(int n = 0; n < mesh->Np; ++n) {
      dlong id = e * mesh->Np + n;
      if ((minRank[id] != rank) || (maxRank[id] != rank)) {
        isHalo = 1;
        break;
      }
    }
    if(isHalo)
      mesh->globalGatherElementList[globalCount++] = e;
    else
      mesh->localGatherElementList[localCount++] = e;
  }
  //printf("local = %d, global = %d\n", localCount, globalCount);

  mesh->NglobalGatherElements = globalCount;
  mesh->NlocalGatherElements = localCount;

  if(globalCount)
    mesh->o_globalGatherElementList =
      mesh->device.malloc(globalCount * sizeof(dlong), mesh->globalGatherElementList);

  if(localCount)
    mesh->o_localGatherElementList =
      mesh->device.malloc(localCount * sizeof(dlong), mesh->localGatherElementList);
}
