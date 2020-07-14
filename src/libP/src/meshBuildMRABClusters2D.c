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

int compareCluster2D(const void *a, const void *b) {
  cElement_t *na = (cElement_t *) a;
  cElement_t *nb = (cElement_t *) b;

  if (na->cRank < nb->cRank) return -1;
  if (nb->cRank < na->cRank) return +1;

  if (na->cId < nb->cId) return -1;
  if (nb->cId < na->cId) return +1;  

  return 0;
}


void meshBuildMRABClusters2D(mesh2D *mesh, int lev, dfloat *weights, int *levels,
            int *Nclusters, cluster_t **clusters, int *Nelements, cElement_t **elements) {

  int rank, size;

  rank = mesh->rank;
  size = mesh->size;

  // minimum {vertex id % size}
  int *Nsend = (int*) calloc(size, sizeof(int));
  int *Nrecv = (int*) calloc(size, sizeof(int));
  int *Ncount = (int*) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size, sizeof(int));
  int *recvOffsets = (int*) calloc(size, sizeof(int));
  int *sendCounts = (int*) calloc(size, sizeof(int));


  //build element struct
  *elements = (cElement_t *) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(cElement_t));
  for (int e=0;e<mesh->Nelements;e++) {
    (*elements)[e].id = e;
    (*elements)[e].level = 0.;
    if (levels) (*elements)[e].level = levels[e];
    
    (*elements)[e].weight = 1.;
    if (weights) (*elements)[e].weight = weights[e];

    for(int n=0;n<mesh->Nverts;++n){
      (*elements)[e].v[n] = mesh->EToV[e*mesh->Nverts+n];
      (*elements)[e].EX[n] = mesh->EX[e*mesh->Nverts+n];
      (*elements)[e].EY[n] = mesh->EY[e*mesh->Nverts+n];
    }
    (*elements)[e].type = mesh->elementInfo[e];

    //initialize the clustering numbering
    (*elements)[e].cId = e; 
    (*elements)[e].cRank = rank;
  }

  cElement_t *sendBuffer = (cElement_t *) calloc(mesh->totalHaloPairs,sizeof(cElement_t));

  //propagate clusters 
  int allDone = 0;
  int rankDone, done;
  while(!allDone) {
    meshHaloExchange(mesh, sizeof(cElement_t), *elements, sendBuffer, *elements + mesh->Nelements);

    rankDone = 1;
    //local clustering
    done = 0;
    while(!done) {
      done = 1;
      for (int e=0;e<mesh->Nelements;e++) {
        for (int f=0;f<mesh->Nfaces;f++) {
          int eP = mesh->EToE[e*mesh->Nfaces +f];
          if (eP>-1) {
            if (((*elements)[eP].level<lev+1)||((*elements)[e].level<lev+1)){
              if (compareCluster2D(*elements+eP,*elements+e)<0) {
                (*elements)[e].cRank = (*elements)[eP].cRank;
                (*elements)[e].cId   = (*elements)[eP].cId;
                done = 0;
                rankDone = 0;
              }
            }
          }
        }
      }    
    }

    MPI_Allreduce(&rankDone, &allDone, 1, MPI_INT, MPI_SUM, mesh->comm);
    allDone /= size;
  }

  //clusters have been built
  //transfer them to their owning rank

  qsort((*elements), mesh->Nelements, sizeof(cElement_t), compareCluster2D);

  //set up exchange along MPI interfaces
  for (int r=0;r<size;r++)
    Nsend[r] = 0;

  for(int e=0;e<mesh->Nelements;++e)
    ++Nsend[(*elements)[e].cRank];

  // find send offsets
  sendOffsets[0] = 0;
  for(int r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];
  
  // exchange byte counts 
  MPI_Alltoall(Nsend, 1, MPI_INT,
         Nrecv, 1, MPI_INT,
         mesh->comm);
  
  // count incoming faces
  int allNrecv = 0;
  for(int r=0;r<size;++r){
    allNrecv += Nrecv[r];
    Nrecv[r] *= sizeof(cElement_t);
    Nsend[r] *= sizeof(cElement_t);
    sendOffsets[r] *= sizeof(cElement_t);
  }
  for(int r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  // buffer for recvied elements
  cElement_t *recvElements = (cElement_t*) calloc(allNrecv, sizeof(cElement_t));  

  // exchange parallel faces
  MPI_Alltoallv(*elements, Nsend, sendOffsets, MPI_CHAR,
      recvElements, Nrecv, recvOffsets, MPI_CHAR,
      mesh->comm);

  free(*elements);
  *elements = recvElements;
  *Nelements = allNrecv;

  qsort((*elements), *Nelements, sizeof(cElement_t), compareCluster2D);

  //build cluster lists
  // the lists are already sorted by cluster, so we just scan for different indices
  *Nclusters = 0;
  if (*Nelements) {
    (*Nclusters)++;
    for (int e=1;e<*Nelements;e++) {
      if ((*elements)[e].cId != (*elements)[e-1].cId) (*Nclusters)++;
    }

    *clusters = (cluster_t *) calloc(*Nclusters,sizeof(cluster_t));

    int cnt  = 0;
    int ecnt = 1;
    (*clusters)[0].Nelements = 1;
    (*clusters)[0].offSet = 0;
    for (int e=1;e<*Nelements;e++) {
      if ((*elements)[e].cId != (*elements)[e-1].cId) {
        cnt++;
        (*clusters)[cnt].offSet = e;
        (*clusters)[cnt].Nelements = 1;
      } else {
        (*clusters)[cnt].Nelements++;
      }
    }
  }
}
