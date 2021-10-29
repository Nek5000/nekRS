#include "elliptic.h"
#include "platform.hpp"

void ellipticOgs(mesh_t *mesh, dlong _Nlocal, int nFields, dlong offset, int *BCType, int BCTypeOffset,
                 dlong& Nmasked, occa::memory& o_mapB, occa::memory& o_maskIds, ogs_t **ogs)
{
  const int Nlocal = (nFields == 1) ? _Nlocal : nFields*offset;
  const int largeNumber = 1 << 20;

  int *mapB = (int*) calloc(Nlocal, sizeof(int));
  for(int fld = 0; fld < nFields; fld++) {
    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int n = 0; n < mesh->Np; n++)
        mapB[n + e * mesh->Np + fld * offset] = largeNumber;
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bc = mesh->EToB[f + e * mesh->Nfaces];
        if (bc > 0) {
          int BCFlag = BCType[bc + BCTypeOffset * fld];
          for (int n = 0; n < mesh->Nfp; n++) {
            int fid = mesh->faceNodes[n + f * mesh->Nfp];
            mapB[fid + e * mesh->Np + fld * offset] =
              mymin(BCFlag, mapB[fid + e * mesh->Np + fld * offset]);
          }
        }
      }
    }
  }
  ogsGatherScatterMany(mapB,
                       nFields,
                       offset,
                       ogsInt,
                       ogsMin,
                       mesh->ogs);

  Nmasked = 0;
  for(int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == largeNumber) {
        mapB[n + fld * offset] = 0;
      } else if (mapB[n + fld * offset] == 1) { //Dirichlet boundary
        Nmasked++;
      }
    }
  }
  o_mapB = platform->device.malloc(Nlocal * sizeof(int), mapB);
  dlong *maskIds = (dlong*) calloc(Nmasked, sizeof(dlong));

  Nmasked = 0;
  for(int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == 1) maskIds[Nmasked++] = n + fld * offset;
    }
  }
  if(Nmasked) o_maskIds = platform->device.malloc(Nmasked * sizeof(dlong), maskIds);
  free(mapB);

  if(! *ogs) {
    if(nFields > 1) {
      if(platform->comm.mpiRank == 0)
        printf("Creating a masked gs handle for nFields > 1 is currently not supported!\n");
      ABORT(EXIT_FAILURE);
    }

    hlong* maskedGlobalIds = (hlong*) calloc(mesh->Nlocal,sizeof(hlong));
    memcpy(maskedGlobalIds, mesh->globalIds, mesh->Nlocal * sizeof(hlong));
    for (dlong n = 0; n < Nmasked; n++) maskedGlobalIds[maskIds[n]] = 0;
    *ogs = ogsSetup(mesh->Nlocal, maskedGlobalIds, platform->comm.mpiComm, 1, platform->device.occaDevice());
    free(maskedGlobalIds);
  }
  free(maskIds);
}
