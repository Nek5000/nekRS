#include "elliptic.h"
#include "platform.hpp"

void ellipticOgs(mesh_t *mesh,
                 dlong mNlocal,
                 int nFields,
                 dlong offset,
                 int *EToB,
                 dlong &Nmasked,
                 occa::memory &o_maskIds,
                 dlong &NmaskedLocal,
                 occa::memory &o_maskIdsLocal,
                 dlong &NmaskedGlobal,
                 occa::memory &o_maskIdsGlobal,
                 ogs_t **ogs)
{
  const int Nlocal = (nFields == 1) ? mNlocal : nFields * offset;
  const int largeNumber = 1 << 20;

  int *mapB = (int*) calloc(Nlocal, sizeof(int));
  for(int fld = 0; fld < nFields; fld++) {
    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int n = 0; n < mesh->Np; n++)
        mapB[n + e * mesh->Np + fld * offset] = largeNumber;
      for (int f = 0; f < mesh->Nfaces; f++) {
        const int fOffset = fld * mesh->Nelements * mesh->Nfaces;
        int bc = EToB[f + e * mesh->Nfaces + fOffset];
        if (bc > 0) {
          for (int n = 0; n < mesh->Nfp; n++) {
            int fid = mesh->faceNodes[n + f * mesh->Nfp];
            mapB[fid + e * mesh->Np + fld * offset] =
                std::min(bc, mapB[fid + e * mesh->Np + fld * offset]); // DIRICHLET wins
          }
        }
      }
    }
  }
  ogsGatherScatterMany(mapB,
                       nFields,
                       offset,
                       ogsInt,
                       ogsMin, // DIRICHLET wins
                       mesh->ogs);

  Nmasked = 0;
  for(int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == largeNumber) {
        mapB[n + fld * offset] = 0;
      }
      else if (mapB[n + fld * offset] == DIRICHLET) {
        Nmasked++;
      }
    }
  }
  dlong *maskIds = (dlong*) calloc(Nmasked, sizeof(dlong));

  Nmasked = 0;
  for(int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == DIRICHLET)
        maskIds[Nmasked++] = n + fld * offset;
    }
  }
  if(Nmasked) o_maskIds = platform->device.malloc(Nmasked * sizeof(dlong), maskIds);

  NmaskedLocal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong el = 0; el < mesh->NlocalGatherElements; ++el) {
      const dlong elemOffset = mesh->localGatherElementList[el] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == DIRICHLET)
          NmaskedLocal++;
      }
    }
  }
  dlong *localMaskIds = (dlong *)calloc(NmaskedLocal, sizeof(dlong));
  NmaskedLocal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong el = 0; el < mesh->NlocalGatherElements; ++el) {
      const dlong elemOffset = mesh->localGatherElementList[el] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == DIRICHLET)
          localMaskIds[NmaskedLocal++] = n + fld * offset;
      }
    }
  }
  if (NmaskedLocal)
    o_maskIdsLocal = platform->device.malloc(NmaskedLocal * sizeof(dlong), localMaskIds);
  free(localMaskIds);

  NmaskedGlobal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong eg = 0; eg < mesh->NglobalGatherElements; ++eg) {
      const dlong elemOffset = mesh->globalGatherElementList[eg] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == DIRICHLET)
          NmaskedGlobal++;
      }
    }
  }
  dlong *globalMaskIds = (dlong *)calloc(NmaskedGlobal, sizeof(dlong));
  NmaskedGlobal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong eg = 0; eg < mesh->NglobalGatherElements; ++eg) {
      const dlong elemOffset = mesh->globalGatherElementList[eg] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == DIRICHLET)
          globalMaskIds[NmaskedGlobal++] = n + fld * offset;
      }
    }
  }
  if (NmaskedGlobal)
    o_maskIdsGlobal = platform->device.malloc(NmaskedGlobal * sizeof(dlong), globalMaskIds);
  free(globalMaskIds);

  free(mapB);

  if(! *ogs) {
    nrsCheck(nFields > 1, platform->comm.mpiComm, EXIT_FAILURE,
             "%s\n", "Creating a masked gs handle for nFields > 1 is currently not supported!");

    hlong* maskedGlobalIds = (hlong*) calloc(mesh->Nlocal,sizeof(hlong));
    memcpy(maskedGlobalIds, mesh->globalIds, mesh->Nlocal * sizeof(hlong));
    for (dlong n = 0; n < Nmasked; n++) maskedGlobalIds[maskIds[n]] = 0;
    *ogs = ogsSetup(mesh->Nlocal, maskedGlobalIds, platform->comm.mpiComm, 1, platform->device.occaDevice());
    free(maskedGlobalIds);
  }
  free(maskIds);
}
