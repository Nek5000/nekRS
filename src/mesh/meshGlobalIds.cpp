#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"
#include "nekInterfaceAdapter.hpp"

struct parallelNode_t
{
  int baseRank;
  hlong baseId;
};

// uniquely label each node with a global index, used for gatherScatter
void meshNekParallelConnectNodes(mesh_t* mesh)
{
  int rank, size;
  rank = platform->comm.mpiRank;
  size = platform->comm.mpiCommSize;

  dlong localNodeCount = mesh->Np * mesh->Nelements;

  mesh->globalIds = (hlong*) calloc(localNodeCount, sizeof(hlong));
  hlong ngv = nek::set_glo_num(mesh->N + 1, mesh->cht);
  for(dlong id = 0; id < localNodeCount; ++id)
    mesh->globalIds[id] = nekData.glo_num[id];
}

void meshGlobalIds(mesh_t* mesh)
{
  meshNekParallelConnectNodes(mesh);
  return;
}
