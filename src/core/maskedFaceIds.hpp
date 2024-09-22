#ifndef maskedFaceIds_H
#define maskedFaceIds_H 1

#include "platform.hpp"
#include "mesh.h"

static std::tuple<dlong, occa::memory, dlong, occa::memory, dlong, occa::memory>
maskedFaceIds(mesh_t *mesh,
              dlong mNlocal,
              int nFields,
              dlong offset,
              int *EToB,
              int dirichletBcTypeId)
{
  dlong Nmasked;
  occa::memory o_maskIds;
  dlong NmaskedLocal;
  occa::memory o_maskIdsLocal;
  dlong NmaskedGlobal;
  occa::memory o_maskIdsGlobal;

  const int Nlocal = (nFields == 1) ? mNlocal : nFields * offset;
  const int largeNumber = 1 << 20;

  std::vector<int> mapB(Nlocal);
 
  // assumes mapB is identical for all nFields
  for (int fld = 0; fld < nFields; fld++) {

    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int n = 0; n < mesh->Np; n++) {
        mapB[n + e * mesh->Np + fld * offset] = largeNumber;
      }
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
  ogsGatherScatterMany(mapB.data(),
                       nFields,
                       offset,
                       ogsInt,
                       ogsMin, // DIRICHLET wins
                       mesh->ogs);

  Nmasked = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == largeNumber) {
        mapB[n + fld * offset] = 0;
      } else if (mapB[n + fld * offset] == dirichletBcTypeId) {
        Nmasked++;
      }
    }
  }
  std::vector<dlong> maskIds(Nmasked);

  Nmasked = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == dirichletBcTypeId) {
        maskIds[Nmasked++] = n + fld * offset;
      }
    }
  }
  if (Nmasked) {
    o_maskIds = platform->device.malloc<dlong>(Nmasked);
    o_maskIds.copyFrom(maskIds.data());
  }

  NmaskedLocal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong el = 0; el < mesh->NlocalGatherElements; ++el) {
      const dlong elemOffset = mesh->localGatherElementList[el] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          NmaskedLocal++;
        }
      }
    }
  }

  std::vector<dlong> localMaskIds(NmaskedLocal);
  NmaskedLocal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong el = 0; el < mesh->NlocalGatherElements; ++el) {
      const dlong elemOffset = mesh->localGatherElementList[el] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          localMaskIds[NmaskedLocal++] = n + fld * offset;
        }
      }
    }
  }
  if (NmaskedLocal) {
    o_maskIdsLocal = platform->device.malloc<dlong>(NmaskedLocal);
    o_maskIdsLocal.copyFrom(localMaskIds.data());
  }

  NmaskedGlobal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong eg = 0; eg < mesh->NglobalGatherElements; ++eg) {
      const dlong elemOffset = mesh->globalGatherElementList[eg] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          NmaskedGlobal++;
        }
      }
    }
  }
  std::vector<dlong> globalMaskIds(NmaskedGlobal);
  NmaskedGlobal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong eg = 0; eg < mesh->NglobalGatherElements; ++eg) {
      const dlong elemOffset = mesh->globalGatherElementList[eg] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          globalMaskIds[NmaskedGlobal++] = n + fld * offset;
        }
      }
    }
  }
  if (NmaskedGlobal) {
    o_maskIdsGlobal = platform->device.malloc<dlong>(NmaskedGlobal);
    o_maskIdsGlobal.copyFrom(globalMaskIds.data());
  }


  return {Nmasked, o_maskIds, NmaskedLocal, o_maskIdsLocal, NmaskedGlobal, o_maskIdsGlobal};  
}

#endif
