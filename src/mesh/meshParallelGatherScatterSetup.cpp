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
#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"

void meshParallelGatherScatterSetup(mesh_t *mesh,
                                    dlong N,
                                    hlong *globalIds,
                                    MPI_Comm &comm,
                                    oogs_mode gsMode,
                                    int verbose)
{

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (platform->comm.mpiRank == 0)
    std::cout << "meshParallelGatherScatterSetup N=" << mesh->N << "\n";
  mesh->ogs = ogsSetup(N, globalIds, comm, verbose, platform->device.occaDevice());

  // use the gs to find what nodes are local to this rank
  int *minRank = (int *)calloc(N, sizeof(int));
  int *maxRank = (int *)calloc(N, sizeof(int));
  for (dlong i = 0; i < N; i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
  }

  ogsGatherScatter(minRank,
                   ogsInt,
                   ogsMin,
                   mesh->ogs); // minRank[n] contains the smallest rank taking part in the gather of node n
  ogsGatherScatter(maxRank,
                   ogsInt,
                   ogsMax,
                   mesh->ogs); // maxRank[n] contains the largest rank taking part in the gather of node n

  int overlap = 1;
  if (platform->options.compareArgs("ENABLE GS COMM OVERLAP", "FALSE"))
    overlap = 0;

  // count elements that contribute to global C0 gather-scatter
  dlong globalCount = 0;
  dlong localCount = 0;
  for (dlong e = 0; e < mesh->Nelements; ++e) {
    int isHalo = 1;
    if (overlap) {
      isHalo = 0;
      for (int n = 0; n < mesh->Np; ++n) {
        dlong id = e * mesh->Np + n;
        if ((minRank[id] != rank) || (maxRank[id] != rank)) {
          isHalo = 1;
          break;
        }
      }
    }
    globalCount += isHalo;
    localCount += 1 - isHalo;
  }

  mesh->elementList = (dlong *)calloc(mesh->Nelements, sizeof(dlong));
  mesh->globalGatherElementList = (dlong *)calloc(globalCount, sizeof(dlong));
  mesh->localGatherElementList = (dlong *)calloc(localCount, sizeof(dlong));

  for (dlong e = 0; e < mesh->Nelements; ++e) {
    mesh->elementList[e] = e;
  }

  globalCount = 0;
  localCount = 0;

  for (dlong e = 0; e < mesh->Nelements; ++e) {
    int isHalo = 1;
    if (overlap) {
      isHalo = 0;
      for (int n = 0; n < mesh->Np; ++n) {
        dlong id = e * mesh->Np + n;
        if ((minRank[id] != rank) || (maxRank[id] != rank)) {
          isHalo = 1;
          break;
        }
      }
    }
    if (isHalo)
      mesh->globalGatherElementList[globalCount++] = e;
    else
      mesh->localGatherElementList[localCount++] = e;
  }

  free(minRank);
  free(maxRank);

  mesh->NglobalGatherElements = globalCount;
  mesh->NlocalGatherElements = localCount;

  mesh->o_elementList = platform->device.malloc<dlong>(mesh->Nelements, mesh->elementList);

  if (globalCount)
    mesh->o_globalGatherElementList =
        platform->device.malloc<dlong>(globalCount, mesh->globalGatherElementList);

  if (localCount)
    mesh->o_localGatherElementList =
        platform->device.malloc<dlong>(localCount, mesh->localGatherElementList);

  { // sanity check
    int err = 0;
    hlong NpGlobal = static_cast<hlong>(mesh->Nelements) * mesh->Np;
    MPI_Allreduce(MPI_IN_PLACE, &NpGlobal, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);

    occa::memory o_tmp = platform->device.malloc<dfloat>(mesh->Nlocal);
    platform->linAlg->fill(mesh->Nlocal, 1.0, o_tmp);
    ogsGatherScatter(o_tmp, ogsDfloat, ogsAdd, mesh->ogs);

    std::vector<dfloat> tmp(o_tmp.length());
    o_tmp.copyTo(tmp.data());

    hlong sum = 0;
    const dfloat eps = 0.1;
    for (int i = 0; i < mesh->Nlocal; i++) {
      dfloat val = tmp[i]*mesh->ogs->invDegree[i];
      hlong valInt = val/(abs(val)) * (abs(val)+0.5); // rounded to nearest integer
      sum += valInt;
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);

    if (sum - NpGlobal != 0) {
      err++;
    }
    nrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "invDegree sanity check failed");
  }
}
