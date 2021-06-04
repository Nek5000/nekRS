/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"

namespace parAlmond {

void adjustPartition(agmgLevel *level, hlong* FineToCoarse, setupAide options) {
  // MPI info
  int rank, size;
  rank = agmg::rank;
  size = agmg::size;

  dlong N = level->A->Nrows;

  //Need to establish 'ownership' of aggregates

  //Keep the current partitioning for STRONGNODES.
  // The rank that had the strong node for each aggregate owns the aggregate
  if (options.compareArgs("PARALMOND PARTITION", "STRONGNODES")) return;

  //populate aggregate array
  hlong gNumAggs = level->globalAggStarts[size]; //total number of aggregates

  parallelAggregate_t *sendAggs;
  if (N)
    sendAggs = (parallelAggregate_t *) calloc(N,sizeof(parallelAggregate_t));
  else
    sendAggs = (parallelAggregate_t *) calloc(1,sizeof(parallelAggregate_t));

  for (dlong i=0;i<N;i++) {
    sendAggs[i].fineId = i;
    sendAggs[i].originRank = rank;

    sendAggs[i].coarseId = FineToCoarse[i];

    //set a temporary owner. Evenly distibute aggregates amoungst ranks
    sendAggs[i].ownerRank = (int) (FineToCoarse[i]*size)/gNumAggs;
  }

  // Make the MPI_PARALLEL_AGGREGATE data type
  MPI_Datatype MPI_PARALLEL_AGGREGATE;
  MPI_Datatype dtype[5] = {MPI_DLONG, MPI_HLONG, MPI_HLONG, MPI_INT, MPI_INT};
  int blength[5] = {1, 1, 1, 1, 1};
  MPI_Aint addr[5], displ[5];
  MPI_Get_address ( &(sendAggs[0]            ), addr+0);
  MPI_Get_address ( &(sendAggs[0].coarseId   ), addr+1);
  MPI_Get_address ( &(sendAggs[0].newCoarseId), addr+2);
  MPI_Get_address ( &(sendAggs[0].originRank ), addr+3);
  MPI_Get_address ( &(sendAggs[0].ownerRank  ), addr+4);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  MPI_Type_create_struct (5, blength, displ, dtype, &MPI_PARALLEL_AGGREGATE);
  MPI_Type_commit (&MPI_PARALLEL_AGGREGATE);

  //sort by owning rank for all_reduce
  qsort(sendAggs, N, sizeof(parallelAggregate_t), compareOwner);

  int *sendCounts = (int *) calloc(size,sizeof(int));
  int *recvCounts = (int *) calloc(size,sizeof(int));
  int *sendOffsets = (int *) calloc(size+1,sizeof(int));
  int *recvOffsets = (int *) calloc(size+1,sizeof(int));

  for(dlong i=0;i<N;++i)
    sendCounts[sendAggs[i].ownerRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, agmg::comm);

  // find send and recv offsets for gather
  dlong recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }
  parallelAggregate_t *recvAggs = (parallelAggregate_t *) calloc(recvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(sendAggs, sendCounts, sendOffsets, MPI_PARALLEL_AGGREGATE,
                recvAggs, recvCounts, recvOffsets, MPI_PARALLEL_AGGREGATE,
                agmg::comm);

  //sort by coarse aggregate number, and then by original rank
  qsort(recvAggs, recvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates here
  dlong NumUniqueAggs =0;
  if (recvNtotal) NumUniqueAggs++;
  for (dlong i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) NumUniqueAggs++;

  //get their locations in the array
  dlong *aggStarts;
  if (NumUniqueAggs)
    aggStarts = (dlong *) calloc(NumUniqueAggs+1,sizeof(dlong));
  dlong cnt = 1;
  for (dlong i=1;i<recvNtotal;i++)
    if(recvAggs[i].coarseId!=recvAggs[i-1].coarseId) aggStarts[cnt++] = i;
  aggStarts[NumUniqueAggs] = recvNtotal;


  if (options.compareArgs("PARALMOND PARTITION", "DISTRIBUTED")) { //rank that contributes most to the aggregate ownes it
    //use a random dfloat for each rank to break ties.
    dfloat rand = (dfloat) drand48();
    dfloat *gRands = (dfloat *) calloc(size,sizeof(dfloat));
    MPI_Allgather(&rand, 1, MPI_DFLOAT, gRands, 1, MPI_DFLOAT, agmg::comm);

    //determine the aggregates majority owner
    int *rankCounts = (int *) calloc(size,sizeof(int));
    for (dlong n=0;n<NumUniqueAggs;n++) {
      //populate randomizer
      for (int r=0;r<size;r++)
        rankCounts[r] = gRands[r];

      //count the number of contributions to the aggregate from the separate ranks
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++)
        rankCounts[recvAggs[i].originRank]++;

      //find which rank is contributing the most to this aggregate
      int ownerRank = 0;
      dfloat maxEntries = rankCounts[0];
      for (int r=1;r<size;r++) {
        if (rankCounts[r]>maxEntries) {
          ownerRank = r;
          maxEntries = rankCounts[r];
        }
      }

      //set this aggregate's owner
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++)
        recvAggs[i].ownerRank = ownerRank;
    }
    free(gRands); free(rankCounts);
  } else { //default SATURATE: always choose the lowest rank to own the aggregate
    for (dlong n=0;n<NumUniqueAggs;n++) {

      int minrank = size;

      //count the number of contributions to the aggregate from the separate ranks
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++){

        minrank = (recvAggs[i].originRank<minrank) ? recvAggs[i].originRank : minrank;
      }

      //set this aggregate's owner
      for (dlong i=aggStarts[n];i<aggStarts[n+1];i++)
        recvAggs[i].ownerRank = minrank;
    }
  }
  free(aggStarts);

  //sort by owning rank
  qsort(recvAggs, recvNtotal, sizeof(parallelAggregate_t), compareOwner);

  int *newSendCounts = (int *) calloc(size,sizeof(int));
  int *newRecvCounts = (int *) calloc(size,sizeof(int));
  int *newSendOffsets = (int *) calloc(size+1,sizeof(int));
  int *newRecvOffsets = (int *) calloc(size+1,sizeof(int));

  for(dlong i=0;i<recvNtotal;++i)
    newSendCounts[recvAggs[i].ownerRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(newSendCounts, 1, MPI_INT, newRecvCounts, 1, MPI_INT, agmg::comm);

  // find send and recv offsets for gather
  dlong newRecvNtotal = 0;
  for(int r=0;r<size;++r){
    newSendOffsets[r+1] = newSendOffsets[r] + newSendCounts[r];
    newRecvOffsets[r+1] = newRecvOffsets[r] + newRecvCounts[r];
    newRecvNtotal += newRecvCounts[r];
  }
  parallelAggregate_t *newRecvAggs = (parallelAggregate_t *) calloc(newRecvNtotal,sizeof(parallelAggregate_t));

  MPI_Alltoallv(   recvAggs, newSendCounts, newSendOffsets, MPI_PARALLEL_AGGREGATE,
                newRecvAggs, newRecvCounts, newRecvOffsets, MPI_PARALLEL_AGGREGATE,
                agmg::comm);

  //sort by coarse aggregate number, and then by original rank
  qsort(newRecvAggs, newRecvNtotal, sizeof(parallelAggregate_t), compareAgg);

  //count the number of unique aggregates this rank owns
  dlong numAggs = 0;
  if (newRecvNtotal) numAggs++;
  for (dlong i=1;i<newRecvNtotal;i++)
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) numAggs++;

  //determine a global numbering of the aggregates
  dlong *lNumAggs = (dlong*) calloc(size,sizeof(dlong));
  MPI_Allgather(&numAggs, 1, MPI_DLONG, lNumAggs, 1, MPI_INT, agmg::comm);

  level->globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    level->globalAggStarts[r+1] = level->globalAggStarts[r] + lNumAggs[r];

  //set the new global coarse index
  cnt = level->globalAggStarts[rank];
  if (newRecvNtotal) newRecvAggs[0].newCoarseId = cnt;
  for (dlong i=1;i<newRecvNtotal;i++) {
    if(newRecvAggs[i].coarseId!=newRecvAggs[i-1].coarseId) cnt++;

    newRecvAggs[i].newCoarseId = cnt;
  }

  //sort by owning rank
  qsort(newRecvAggs, newRecvNtotal, sizeof(parallelAggregate_t), compareOrigin);

  for(int r=0;r<size;r++) sendCounts[r] = 0;
  for(int r=0;r<=size;r++) {
    sendOffsets[r] = 0;
    recvOffsets[r] = 0;
  }

  for(dlong i=0;i<newRecvNtotal;++i)
    sendCounts[newRecvAggs[i].originRank]++;

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, agmg::comm);

  // find send and recv offsets for gather
  recvNtotal = 0;
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r];
  }

  //send the aggregate data back
  MPI_Alltoallv(newRecvAggs, sendCounts, sendOffsets, MPI_PARALLEL_AGGREGATE,
                   sendAggs, recvCounts, recvOffsets, MPI_PARALLEL_AGGREGATE,
                agmg::comm);

  //clean up
  MPI_Barrier(agmg::comm);
  MPI_Type_free(&MPI_PARALLEL_AGGREGATE);

  free(recvAggs);
  free(sendCounts);  free(recvCounts);
  free(sendOffsets); free(recvOffsets);
  free(newRecvAggs);
  free(newSendCounts);  free(newRecvCounts);
  free(newSendOffsets); free(newRecvOffsets);

  //record the new FineToCoarse map
  for (dlong i=0;i<N;i++)
    FineToCoarse[sendAggs[i].fineId] = sendAggs[i].newCoarseId;

  free(sendAggs);
}

} //namespace parAlmond