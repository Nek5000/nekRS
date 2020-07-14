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
#include "mpi.h"
#include "mesh2D.h"

#define bitRange 10

typedef struct {
  int id;
  int level;
  dfloat weight;

  // 4 for maximum number of vertices per element in 2D
  int v[4];
  dfloat EX[4], EY[4];

  int cRank;
  int cId;
  int type;
} cElement_t;

typedef struct {
  int Nelements;
  int offSet;
} cluster_t;

typedef struct {
  int Nelements;
  int offSet;
  int rank;

  int destId;
  int destOffset;
  int destRank;

  dfloat weight;
  unsigned int index; //hilbert index
} parallelCluster_t;

//This is linked form meshGeometricPartition2D.c
unsigned int hilbert2D(unsigned int n, unsigned int index1, unsigned int index2);
void bogusMatch(void *a, void *b);

dfloat improveClusteredPartition2D(int rank, int size, MPI_Comm comm,
				   int *Nclusters, parallelCluster_t **parallelClusters);

// compare the Morton indices for two clusters
int compareIndex2D(const void *a, const void *b){

  parallelCluster_t *ca = (parallelCluster_t*) a;
  parallelCluster_t *cb = (parallelCluster_t*) b;
  
  if(ca->index < cb->index) return -1;
  if(ca->index > cb->index) return  1;
  
  return 0;
}

// compare the Morton indices for two clusters
int compareRank2D(const void *a, const void *b){

  parallelCluster_t *ca = (parallelCluster_t*) a;
  parallelCluster_t *cb = (parallelCluster_t*) b;
  
  if(ca->rank < cb->rank) return -1;
  if(ca->rank > cb->rank) return  1;

  if(ca->offSet < cb->offSet) return -1;
  if(ca->offSet > cb->offSet) return  1;
  
  return 0;
}

// geometric partition of clusters of elements in 2D mesh using Morton ordering + parallelSort
dfloat meshClusteredGeometricPartition2D(mesh2D *mesh, int Nclusters, cluster_t *clusters, 
					 int *Nelements, cElement_t **elements){

  int rank, size;
  rank = mesh->rank;
  size = mesh->size;

  int maxNclusters;
  MPI_Allreduce(&Nclusters, &maxNclusters, 1, MPI_INT, MPI_MAX,mesh->comm);
  maxNclusters = 2*((maxNclusters+1)/2);
  
  // fix maxNclusters
  parallelCluster_t *parallelClusters 
    = (parallelCluster_t*) calloc(maxNclusters, sizeof(parallelCluster_t));

  // local bounding box of element centers
  dfloat mincx = 1e9, maxcx = -1e9;
  dfloat mincy = 1e9, maxcy = -1e9;

  // compute cluster centers on this process
  for(int cnt=0;cnt<Nclusters;++cnt){
    int id = clusters[cnt].offSet;
    dfloat cx = 0, cy = 0;
    for (int e=0;e<clusters[cnt].Nelements;e++) {
      for(int n=0;n<mesh->Nverts;++n){
        cx += (*elements)[id+e].EX[n];
        cy += (*elements)[id+e].EY[n];
      }
    }
    cx /= (mesh->Nverts*clusters[cnt].Nelements);
    cy /= (mesh->Nverts*clusters[cnt].Nelements);
    
    mincx = mymin(mincx, cx);
    maxcx = mymax(maxcx, cx);
    mincy = mymin(mincy, cy);
    maxcy = mymax(maxcy, cy);
  }

  dfloat delta = 1e-4;
  mincx -= delta;
  mincy -= delta;
  maxcx += delta;
  maxcy += delta;

  // find global bounding box of cluster centers
  dfloat gmincx, gmincy, gmaxcx, gmaxcy;
  MPI_Allreduce(&mincx, &gmincx, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
  MPI_Allreduce(&mincy, &gmincy, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
  MPI_Allreduce(&maxcx, &gmaxcx, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
  MPI_Allreduce(&maxcy, &gmaxcy, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

  // choose sub-range of Morton lattice coordinates to embed cluster centers in
  unsigned int Nboxes = (((unsigned int)1)<<(bitRange-1));
  
  // compute Morton index for each cluster
  for(int cnt=0;cnt<Nclusters;++cnt){
    // cluster center coordinates
    dfloat cx = 0, cy = 0;
    parallelClusters[cnt].weight = 0.;
    int id = clusters[cnt].offSet;
    for (int e=0;e<clusters[cnt].Nelements;e++) {
      for(int n=0;n<mesh->Nverts;++n){
        cx += (*elements)[id+e].EX[n];
        cy += (*elements)[id+e].EY[n];
      }
      parallelClusters[cnt].weight += (*elements)[id+e].weight;
    }
    cx /= (mesh->Nverts*clusters[cnt].Nelements);
    cy /= (mesh->Nverts*clusters[cnt].Nelements);

    unsigned int ix = (cx-gmincx)*Nboxes/(gmaxcx-gmincx);
    unsigned int iy = (cy-gmincy)*Nboxes/(gmaxcy-gmincy);

    //fill the parallel cluster struct
    parallelClusters[cnt].index =  hilbert2D(Nboxes, ix, iy);
    parallelClusters[cnt].Nelements = clusters[cnt].Nelements;
    parallelClusters[cnt].offSet = clusters[cnt].offSet;
    parallelClusters[cnt].rank = rank;
  }

  // pad cluster array with dummy clusters
  for(int n=Nclusters;n<maxNclusters;++n){
    parallelClusters[n].Nelements = -1;
    parallelClusters[n].index = hilbert2D(Nboxes, Nboxes-1, Nboxes-1);
  }

  // odd-even parallel sort of cluster capsules based on their Morton index
  parallelSort(mesh->size, mesh->rank, mesh->comm,
	       maxNclusters, parallelClusters, sizeof(parallelCluster_t),
	       compareIndex2D, bogusMatch);

  int newNclusters =0;
  for (int n=0;n<maxNclusters;n++)
    newNclusters += (parallelClusters[n].Nelements != -1);

  //Do an initial partitioning
  dfloat localTotalWeight = 0.;
  for (int n=0; n<newNclusters; n++) 
    localTotalWeight += parallelClusters[n].weight;

  dfloat *totalWeights = (dfloat *) calloc(size,sizeof(dfloat));
  dfloat *weightOffsets = (dfloat *) calloc(size+1,sizeof(dfloat));
  
  MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, mesh->comm);

  for (int r=0; r<size; r++)
    weightOffsets[r+1] = weightOffsets[r] + totalWeights[r];

  dfloat globalTotalWeight = weightOffsets[size];
  dfloat chunkSize = globalTotalWeight/((dfloat)size);

  int *Nsend = (int *) calloc(size, sizeof(int));
  int *Nrecv = (int *) calloc(size, sizeof(int));
  int *Ncount = (int *) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size, sizeof(int));
  int *recvOffsets = (int*) calloc(size, sizeof(int));

  //determine the destination rank based on which chunk the cluster is in
  localTotalWeight = weightOffsets[rank];
  for (int n=0; n<newNclusters; n++) {
    int destRank = (int) (localTotalWeight/chunkSize);
    Nsend[destRank]++; 
    localTotalWeight += parallelClusters[n].weight;
  }

  // find send offsets
  for(int r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];
  
  // exchange byte counts 
  MPI_Alltoall(Nsend, 1, MPI_INT, Nrecv, 1, MPI_INT, mesh->comm);
  
  // count incoming clusters
  newNclusters = 0;
  for(int r=0;r<size;++r){
    newNclusters += Nrecv[r];
    Nrecv[r] *= sizeof(parallelCluster_t);
    Nsend[r] *= sizeof(parallelCluster_t);
    sendOffsets[r] *= sizeof(parallelCluster_t);
  }
  for(int r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  parallelCluster_t *tmpParallelClusters = (parallelCluster_t *) calloc(newNclusters, sizeof(parallelCluster_t));
  
  // exchange parallel clusters
  MPI_Alltoallv(parallelClusters, Nsend, sendOffsets, MPI_CHAR,
                tmpParallelClusters, Nrecv, recvOffsets, MPI_CHAR, mesh->comm);

  if (parallelClusters) free(parallelClusters);
  parallelClusters = tmpParallelClusters;

  //improve the partitioning by exchanging elements between neighboring prcesses
  dfloat partQuality = improveClusteredPartition2D(mesh->rank, mesh->size, mesh->comm,
						   &newNclusters, &parallelClusters);

  //now that we're partitioned and (hopefully) balance2Dd, send the elements

  // count number of elements that should end up on this process
  int newNelements = 0;
  for(int n=0;n<newNclusters;n++)
    newNelements += parallelClusters[n].Nelements;

  //record the destination info
  if (newNclusters) {
    parallelClusters[0].destId = 0;
    parallelClusters[0].destOffset = 0;
    parallelClusters[0].destRank = rank;
  }
  for (int n=1; n<newNclusters; n++) {
    parallelClusters[n].destId = n;
    parallelClusters[n].destOffset = parallelClusters[n-1].destOffset + parallelClusters[n-1].Nelements;
    parallelClusters[n].destRank = rank;
  }

  //sort by original rank and offset
  qsort(parallelClusters, newNclusters, sizeof(parallelCluster_t), compareRank2D);

  //reset counters
  for(int r=0;r<size;++r)
    Nsend[r] =0;

  for (int n=0;n<newNclusters;n++) 
    Nsend[parallelClusters[n].rank]++;

  for(int r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];

  // exchange byte counts 
  MPI_Alltoall(Nsend, 1, MPI_INT, Nrecv, 1, MPI_INT, mesh->comm);

  for(int r=0;r<size;++r){
    Nrecv[r] *= sizeof(parallelCluster_t);
    Nsend[r] *= sizeof(parallelCluster_t);
    sendOffsets[r] *= sizeof(parallelCluster_t);
  }
  for(int r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  parallelCluster_t *recvParallelClusters;
  if (Nclusters) 
    recvParallelClusters = (parallelCluster_t*) calloc(Nclusters, sizeof(parallelCluster_t));


  MPI_Alltoallv(parallelClusters, Nsend, sendOffsets, MPI_CHAR,
              recvParallelClusters, Nrecv, recvOffsets, MPI_CHAR, mesh->comm);

  //build the array of elements to send
  cElement_t *sendElements = (cElement_t *) calloc(1,sizeof(cElement_t));
  cElement_t *recvElements = (cElement_t *) calloc(1,sizeof(cElement_t));

  if (*Nelements) sendElements = (cElement_t *) calloc(*Nelements,sizeof(cElement_t));
  if (newNelements) recvElements = (cElement_t *) calloc(newNelements,sizeof(cElement_t));

  //reset send counts
  for (int r=0; r<size; r++)
    Nsend[r] = 0;

  for (int n=0;n<Nclusters;n++) {
    Nsend[recvParallelClusters[n].destRank] += recvParallelClusters[n].Nelements;
  }

  // find send offsets
  for(int r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];

  //build the array of elements to send
  for (int n=0;n<Nclusters;n++) {
    int destRank = recvParallelClusters[n].destRank;
    int cnt = recvParallelClusters[n].Nelements;

    int sendId = sendOffsets[destRank] + Ncount[destRank];
    int id = recvParallelClusters[n].offSet;
    memcpy(sendElements+sendId, *elements+id, cnt*sizeof(cElement_t)); 
    Ncount[destRank] += cnt;
  }
  free(recvParallelClusters);

  // exchange element counts 
  MPI_Alltoall(Nsend, 1, MPI_INT, Nrecv, 1, MPI_INT, mesh->comm);
  
  for(int r=0;r<size;++r){
    Nrecv[r] *= sizeof(cElement_t);
    Nsend[r] *= sizeof(cElement_t);
    sendOffsets[r] *= sizeof(cElement_t);
  }
  for(int r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  MPI_Alltoallv(sendElements, Nsend, sendOffsets, MPI_CHAR,
                recvElements, Nrecv, recvOffsets, MPI_CHAR, mesh->comm);
  free(sendElements);

  //write the clusters in the proper order
  cluster_t *newClusters = (cluster_t *) calloc(newNclusters,sizeof(cluster_t));
  cElement_t *newElements = (cElement_t *) calloc(newNelements,sizeof(cElement_t));
  int cnt =0;
  for (int n=0;n<newNclusters;n++) {
    int id = parallelClusters[n].destId;
    newClusters[id].Nelements = parallelClusters[n].Nelements;
    newClusters[id].offSet = parallelClusters[n].destOffset;
    for (int e=0;e<parallelClusters[n].Nelements;e++) {
      memcpy(newElements + newClusters[id].offSet+e, recvElements+cnt++, sizeof(cElement_t));
    }
  }
  free(recvElements);
  free(parallelClusters);

  *Nelements = newNelements;

  if (*elements) free(*elements);
  *elements = newElements;

  return partQuality;
};



//swap clusters between neighboring processes to try and improve the partitioning
void balance2D(int rank, int size, MPI_Comm comm,
	       int rankL, int rankR, dfloat *weightL, dfloat *weightR, 
              int *Nclusters, parallelCluster_t **parallelClusters) {
  
  
  int tag = 999;
  MPI_Request recv, send;
  MPI_Status status;

  if (rank==rankL) {
    if ( *weightL > *weightR) {
      //count number of clusters to send to proc
      int Nsend = 0;
      for (int cnt=*Nclusters-1;cnt>-1;cnt--) {
        dfloat w = (*parallelClusters)[cnt].weight;
        if ((*weightL-w)>=(*weightR+w)) {
          //sending this cluster improves the balance2D
          *weightL -= w;
          *weightR += w;
          Nsend++; 
        } else if((*weightL-w) > *weightR) { 
          //sending makes the neighbor have a higher weight, but it improves the balance2D
          *weightL -= w;
          *weightR += w;
          Nsend++; 
          break;
        } else {
          break;
        }
      }

      MPI_Isend(&Nsend, 1, MPI_INT,  rankR, tag, comm, &send);
      MPI_Wait(&send, &status);

      if (Nsend) {
        *Nclusters -= Nsend;

        MPI_Isend((*parallelClusters) + *Nclusters, Nsend*sizeof(parallelCluster_t), MPI_CHAR,  rankR, tag, comm, &send);
        MPI_Wait(&send, &status);
      }
    } else if ( *weightL < *weightR) {
      int Nrecv;
      MPI_Irecv(&Nrecv, 1, MPI_INT,  rankR, tag, comm, &recv);
      MPI_Wait(&recv, &status);

      if (Nrecv) {
        parallelCluster_t *newParallelClusters = (parallelCluster_t *) calloc(*Nclusters+Nrecv,sizeof(parallelCluster_t));
        memcpy(newParallelClusters,*parallelClusters,*Nclusters*sizeof(parallelCluster_t));

        MPI_Irecv(newParallelClusters+*Nclusters, Nrecv*sizeof(parallelCluster_t), MPI_CHAR,  rankR, tag, comm, &recv);
        MPI_Wait(&recv, &status);
        
        for (int n=*Nclusters;n<*Nclusters+Nrecv;n++) {
          dfloat w = newParallelClusters[n].weight;
          *weightL += w;
          *weightR -= w;
        }

        *Nclusters += Nrecv;
        free(*parallelClusters);
        *parallelClusters = newParallelClusters;
      }
    }
  } else if (rank==rankR) {
    if (*weightL < *weightR) {
      //count number of clusters to send to proc
      int Nsend = 0;
      for (int cnt=0;cnt<*Nclusters;cnt++) {
        dfloat w = (*parallelClusters)[cnt].weight;
        if ((*weightR-w)>=(*weightL+w)) {
          //sending this cluster improves the balance2D
          *weightR -= w;
          *weightL += w;
          Nsend++; 
        } else if((*weightR-w) > *weightL) { 
          //sending makes the neighbor have a higher weight, but it improves the balance2D
          *weightR -= w;
          *weightL += w;
          Nsend++; 
          break;
        } else {
          break;
        }
      }

      MPI_Isend(&Nsend, 1, MPI_INT,  rankL, tag, comm, &send);
      MPI_Wait(&send, &status);

      if (Nsend) {
        *Nclusters -= Nsend;
        parallelCluster_t *newParallelClusters = (parallelCluster_t *) calloc(*Nclusters,sizeof(parallelCluster_t));
        memcpy(newParallelClusters,(*parallelClusters) + Nsend,*Nclusters*sizeof(parallelCluster_t));

        MPI_Isend(*parallelClusters, Nsend*sizeof(parallelCluster_t), MPI_CHAR,  rankL, tag, comm, &send);
        MPI_Wait(&send, &status);

        free(*parallelClusters);
        *parallelClusters = newParallelClusters;
      }
    } else if (*weightL > *weightR) {
      int Nrecv;
      MPI_Irecv(&Nrecv, 1, MPI_INT,  rankL, tag, comm, &recv);
      MPI_Wait(&recv, &status);

      if (Nrecv) {
        parallelCluster_t *tmpParallelClusters = (parallelCluster_t *) calloc(Nrecv,sizeof(parallelCluster_t));

        MPI_Irecv(tmpParallelClusters, Nrecv*sizeof(parallelCluster_t), MPI_CHAR, rankL, tag, comm, &recv);
        MPI_Wait(&recv, &status);

        for (int n=0;n<Nrecv;n++) {
          dfloat w = tmpParallelClusters[n].weight;
          *weightR += w;
          *weightL -= w;
        }

        *Nclusters += Nrecv;
        parallelCluster_t *newParallelClusters = (parallelCluster_t *) calloc(*Nclusters,sizeof(parallelCluster_t));
        memcpy(newParallelClusters,tmpParallelClusters,Nrecv*sizeof(parallelCluster_t));
        memcpy(newParallelClusters+Nrecv,*parallelClusters,(*Nclusters-Nrecv)*sizeof(parallelCluster_t));

        free(tmpParallelClusters);
        free(*parallelClusters);
        *parallelClusters = newParallelClusters;
      }
    }
  }
}


dfloat improveClusteredPartition2D(int rank, int size, MPI_Comm comm,
				   int *Nclusters, parallelCluster_t **parallelClusters){

  int tag = 999;

  MPI_Request recv, send;
  MPI_Status status;

  dfloat *totalWeights = (dfloat *) calloc(size,sizeof(dfloat));
  dfloat quality;

  while (true) {

    dfloat localTotalWeight = 0.;
    for (int n=0; n<*Nclusters; n++) 
      localTotalWeight += (*parallelClusters)[n].weight;
    
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, comm);

    dfloat maxTotalWeight, minTotalWeight;
    MPI_Allreduce(&localTotalWeight, &minTotalWeight, 1, MPI_DFLOAT, MPI_MIN, comm);
    MPI_Allreduce(&localTotalWeight, &maxTotalWeight, 1, MPI_DFLOAT, MPI_MAX, comm);

    quality = minTotalWeight/maxTotalWeight;

    //ends
    if ((rank==0)||(rank==size-1)) 
      balance2D(rank, size, comm,size-1,0,totalWeights+size-1, totalWeights+0, Nclusters,parallelClusters);

    //resync
    localTotalWeight = totalWeights[rank];
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, comm);

    //evens
    if (( (rank%2) == 0)&&(rank+1 < size))
      balance2D(rank, size, comm,rank,rank+1,totalWeights+rank, totalWeights+rank+1, Nclusters,parallelClusters);
    if (( (rank%2) == 1)&&(rank-1 > -1))
      balance2D(rank, size, comm,rank-1,rank,totalWeights+rank-1, totalWeights+rank, Nclusters,parallelClusters);

    //resync
    localTotalWeight = totalWeights[rank];
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, comm);

    //odds
    if (((rank%2) == 0)&&(rank-1 > -1))
      balance2D(rank, size, comm,rank-1,rank,totalWeights+rank-1, totalWeights+rank, Nclusters,parallelClusters);
    if (((rank%2) == 1)&&(rank+1 < size))
      balance2D(rank, size, comm,rank,rank+1,totalWeights+rank, totalWeights+rank+1, Nclusters,parallelClusters);

    //resync
    localTotalWeight = totalWeights[rank];
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, comm);
    MPI_Allreduce(&localTotalWeight, &minTotalWeight, 1, MPI_DFLOAT, MPI_MIN, comm);
    MPI_Allreduce(&localTotalWeight, &maxTotalWeight, 1, MPI_DFLOAT, MPI_MAX, comm);

    dfloat newQuality = minTotalWeight/maxTotalWeight;

    if (newQuality == quality) break; //no change
  }

  return quality;
} 
