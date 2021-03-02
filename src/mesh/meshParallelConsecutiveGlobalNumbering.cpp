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

struct parallelNode2_t
{
  dlong localId;
  hlong globalId;
  dlong recvId;
  hlong newGlobalId;
  int originalRank;
  int ownerRank;
};

// compare on global indices
int parallelCompareGlobalIndices(const void* a, const void* b)
{
  parallelNode2_t* fa = (parallelNode2_t*) a;
  parallelNode2_t* fb = (parallelNode2_t*) b;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  return 0;
}

// compare on global indices
int parallelCompareSourceIndices(const void* a, const void* b)
{
  parallelNode2_t* fa = (parallelNode2_t*) a;
  parallelNode2_t* fb = (parallelNode2_t*) b;

  if(fa->originalRank < fb->originalRank) return -1;
  if(fa->originalRank > fb->originalRank) return +1;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

// compare on global indices
int parallelCompareOwners2(const void* a, const void* b)
{
  parallelNode2_t* fa = (parallelNode2_t*) a;
  parallelNode2_t* fb = (parallelNode2_t*) b;

  if(fa->ownerRank < fb->ownerRank) return -1;
  if(fa->ownerRank > fb->ownerRank) return +1;

  return 0;
}

// squeeze gaps out of a globalNumbering of local nodes (arranged in NpNum blocks
void meshParallelConsecutiveGlobalNumbering(mesh_t* mesh,
                                            dlong Nnum,
                                            hlong* globalNumbering,
                                            int* globalOwners,
                                            hlong* globalStarts)
{
  
  int rank = platform->comm.mpiRank;
  int size = platform->comm.mpiCommSize;

  // count how many nodes to send to each process
  dlong* allCounts   = (dlong*) calloc(size, sizeof(dlong));

  int* sendCounts = (int*) calloc(size,sizeof(int));
  int* recvCounts = (int*) calloc(size,sizeof(int));
  int* sendOffsets = (int*) calloc(size + 1,sizeof(int));
  int* recvOffsets = (int*) calloc(size + 1,sizeof(int));

  dlong cnt = 0;
  for(dlong n = 0; n < Nnum; ++n) {
    if (globalNumbering[n] < 0) continue; //skip negative ids
    sendCounts[globalOwners[n]]++;
    cnt++;
  }

  dlong Nlocal = cnt; //number of unmasked nodes

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, platform->comm.mpiComm);

  // find send and recv offsets for gather
  dlong recvNtotal = 0;
  for(int r = 0; r < size; ++r) {
    sendOffsets[r + 1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r + 1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }

  printf("Nlocal = %d\n", Nlocal);

  // populate parallel nodes to send
  parallelNode2_t* sendNodes = NULL;

  sendNodes = (parallelNode2_t*) calloc(Nlocal + 1, sizeof(parallelNode2_t));

  // Make the MPI_PARALLELFACE_T data type
  MPI_Datatype MPI_PARALLELFACE_T;
  MPI_Datatype dtype[6] = {MPI_DLONG, MPI_HLONG, MPI_DLONG, MPI_HLONG, MPI_INT, MPI_INT};
  int blength[6] = {1, 1, 1, 1, 1, 1};
  MPI_Aint addr[6], displ[6];
  MPI_Get_address ( &(sendNodes[0]             ), addr + 0);
  MPI_Get_address ( &(sendNodes[0].globalId    ), addr + 1);
  MPI_Get_address ( &(sendNodes[0].recvId      ), addr + 2);
  MPI_Get_address ( &(sendNodes[0].newGlobalId ), addr + 3);
  MPI_Get_address ( &(sendNodes[0].originalRank), addr + 4);
  MPI_Get_address ( &(sendNodes[0].ownerRank   ), addr + 5);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  displ[5] = addr[5] - addr[0];
  MPI_Type_create_struct (6, blength, displ, dtype, &MPI_PARALLELFACE_T);
  MPI_Type_commit (&MPI_PARALLELFACE_T);

  cnt = 0;
  for(dlong n = 0; n < Nnum; ++n) {
    if (globalNumbering[n] < 0) continue; //skip negative ids
    sendNodes[cnt].localId = n;
    sendNodes[cnt].globalId = globalNumbering[n];
    sendNodes[cnt].newGlobalId = -1;
    sendNodes[cnt].originalRank = rank;
    sendNodes[cnt].ownerRank = globalOwners[n];
    cnt++;
  }

  // sort by global index
  qsort(sendNodes, Nlocal, sizeof(parallelNode2_t), parallelCompareOwners2);

  parallelNode2_t* recvNodes = NULL;
  recvNodes = (parallelNode2_t*) calloc(recvNtotal + 1, sizeof(parallelNode2_t));

  // load up node data to send (NEED TO SCALE sendCounts, sendOffsets etc by sizeof(parallelNode2_t)
  MPI_Alltoallv(sendNodes, sendCounts, sendOffsets, MPI_PARALLELFACE_T,
                recvNodes, recvCounts, recvOffsets, MPI_PARALLELFACE_T,
                platform->comm.mpiComm);

  for (dlong n = 0; n < recvNtotal; n++) recvNodes[n].recvId = n;

  // sort by global index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode2_t), parallelCompareGlobalIndices);

  // renumber unique nodes starting from 0 (need to be careful about zeros)
  cnt = 0;
  if (recvNtotal) recvNodes[0].newGlobalId = cnt;
  for(dlong n = 1; n < recvNtotal; ++n) {
    if(recvNodes[n].globalId != recvNodes[n - 1].globalId) // new node
      ++cnt;
    recvNodes[n].newGlobalId = cnt;
  }
  if (recvNtotal) ++cnt; // increment to actual number of unique nodes on this rank

  // collect unique node counts from all processes
  MPI_Allgather(&cnt, 1, MPI_DLONG, allCounts, 1, MPI_DLONG, platform->comm.mpiComm);

  // cumulative sum of unique node counts => starting node index for each process
  for(int r = 0; r < size; ++r)
    globalStarts[r + 1] = globalStarts[r] + allCounts[r];

  // shift numbering
  for(dlong n = 0; n < recvNtotal; ++n)
    recvNodes[n].newGlobalId += globalStarts[rank];

  // sort by rank, local index
  qsort(recvNodes, recvNtotal, sizeof(parallelNode2_t), parallelCompareSourceIndices);

  // reverse all to all to reclaim nodes
  MPI_Alltoallv(recvNodes, recvCounts, recvOffsets, MPI_PARALLELFACE_T,
                sendNodes, sendCounts, sendOffsets, MPI_PARALLELFACE_T,
                platform->comm.mpiComm);

  // extract new global indices and push back to original numbering array
  for(dlong n = 0; n < Nlocal; ++n) {
    // shuffle incoming nodes based on local id
    dlong id = sendNodes[n].localId;
    globalNumbering[id] = sendNodes[n].newGlobalId;
  }

  MPI_Barrier(platform->comm.mpiComm);
  MPI_Type_free(&MPI_PARALLELFACE_T);

  free(sendNodes);
  free(recvNodes);

  free(sendCounts);
  free(recvCounts);
  free(sendOffsets);
  free(recvOffsets);
  free(allCounts);
}
