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


typedef struct{

  int baseRank;
  hlong baseId;

}parallelNode_t;


// uniquely label each node with a global index, used for gatherScatter
void meshParallelConnectNodes(mesh_t *mesh){

  int rank, size;
  rank = mesh->rank; 
  size = mesh->size; 

  dlong localNodeCount = mesh->Np*mesh->Nelements;
  dlong *allLocalNodeCounts = (dlong*) calloc(size, sizeof(dlong));

  MPI_Allgather(&localNodeCount,    1, MPI_DLONG,
                allLocalNodeCounts, 1, MPI_DLONG,
                mesh->comm);
  
  hlong gatherNodeStart = 0;
  for(int r=0;r<rank;++r)
    gatherNodeStart += allLocalNodeCounts[r];
  
  free(allLocalNodeCounts);

  // form continuous node numbering (local=>virtual gather)
  parallelNode_t *localNodes =
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,
                             sizeof(parallelNode_t));

  // use local numbering
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dlong id = e*mesh->Np+n;

      localNodes[id].baseRank = rank;
      localNodes[id].baseId = 1 + id + mesh->Nnodes + gatherNodeStart;

    }

    // use vertex ids for vertex nodes to reduce iterations
    for(int v=0;v<mesh->Nverts;++v){
      dlong id = e*mesh->Np + mesh->vertexNodes[v];
      hlong gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      localNodes[id].baseId = gid; 
    }
  }

  dlong localChange = 0, gatherChange = 1;

  parallelNode_t *sendBuffer =
    (parallelNode_t*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(parallelNode_t));

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    meshHaloExchange(mesh, mesh->Np*sizeof(parallelNode_t),
                     localNodes, sendBuffer, localNodes+localNodeCount);

    // compare trace nodes
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
        dlong id  = e*mesh->Nfp*mesh->Nfaces + n;
        dlong idM = mesh->vmapM[id];
        dlong idP = mesh->vmapP[id];
        hlong gidM = localNodes[idM].baseId;
        hlong gidP = localNodes[idP].baseId;

        int baseRankM = localNodes[idM].baseRank;
        int baseRankP = localNodes[idP].baseRank;
        
        if(gidM<gidP || (gidP==gidM && baseRankM<baseRankP)){
          ++localChange;
          localNodes[idP].baseRank    = localNodes[idM].baseRank;
          localNodes[idP].baseId      = localNodes[idM].baseId;
        }
        
        if(gidP<gidM || (gidP==gidM && baseRankP<baseRankM)){
          ++localChange;
          localNodes[idM].baseRank    = localNodes[idP].baseRank;
          localNodes[idM].baseId      = localNodes[idP].baseId;
        }
      }
    }

    // sum up changes
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_DLONG, MPI_SUM, mesh->comm);
  }

  //make a locally-ordered version
  mesh->globalIds = (hlong*) calloc(localNodeCount, sizeof(hlong));
  for(dlong id=0;id<localNodeCount;++id){
    mesh->globalIds[id] = localNodes[id].baseId;    
  }
  
  free(localNodes);
  free(sendBuffer);
}
