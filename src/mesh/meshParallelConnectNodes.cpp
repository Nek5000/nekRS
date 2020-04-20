#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"
#include "nekInterfaceAdapter.hpp"

extern int nrsBuildOnly;
extern int isTmesh;

typedef struct{

  int baseRank;
  hlong baseId;

}parallelNode_t;


// uniquely label each node with a global index, used for gatherScatter
void meshNekParallelConnectNodes(mesh_t *mesh)
{
  int rank, size;
  rank = mesh->rank; 
  size = mesh->size; 

  dlong localNodeCount = mesh->Np*mesh->Nelements;

  mesh->globalIds = (hlong*) calloc(localNodeCount, sizeof(hlong));
  hlong ngv = nek_set_glo_num(mesh->N+1, isTmesh);
  for(dlong id=0;id<localNodeCount;++id){
    mesh->globalIds[id] = nekData.glo_num[id];    
  }
}

void meshParallelConnectNodes(mesh_t *mesh){

  if(!nrsBuildOnly) {
    // hotfix as libP version seems to be broken
    meshNekParallelConnectNodes(mesh);
    return;
  }

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
