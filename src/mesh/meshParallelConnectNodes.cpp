#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"
#include "nekInterfaceAdapter.hpp"

// uniquely label each node with a global index, used for gatherScatter
void meshParallelConnectNodes(mesh_t *mesh)
{
  int rank, size;
  rank = mesh->rank; 
  size = mesh->size; 

  dlong localNodeCount = mesh->Np*mesh->Nelements;

  mesh->globalIds = (hlong*) calloc(localNodeCount, sizeof(hlong));
  hlong ngv = nek_set_glo_num(mesh->N+1, 0);
  for(dlong id=0;id<localNodeCount;++id){
    mesh->globalIds[id] = nekData.glo_num[id];    
  }
}
