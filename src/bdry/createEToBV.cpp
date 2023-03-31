#include <vector>
#include "elliptic.h"

void createEToBV(const mesh_t* mesh, const int* EToB, occa::memory& o_EToBV)
{
  const int largeNumber = 1 << 20;

  std::vector<int> EToBV(mesh->Nlocal, largeNumber);
  for (dlong e = 0; e < mesh->Nelements; e++) {
    for (int f = 0; f < mesh->Nfaces; f++) {
      int bc = EToB[f + e * mesh->Nfaces];
      if (bc > 0) {
        for (int n = 0; n < mesh->Nfp; n++) {
          int fid = mesh->faceNodes[n + f * mesh->Nfp];
          EToBV[fid + e * mesh->Np] = std::min(bc, EToBV[fid + e * mesh->Np]);
        }
      }
    }
  }

  ogsGatherScatter(EToBV.data(), ogsInt, ogsMin, mesh->ogs);

  for (dlong n = 0; n < mesh->Nlocal; n++) {
    if (EToBV[n] == largeNumber) {
      EToBV[n] = 0;
    }
  }

  o_EToBV.copyFrom(EToBV.data(), EToBV.size() * sizeof(int));
}
