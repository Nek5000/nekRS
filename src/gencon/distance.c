#include <gencon-impl.h>

int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

int findMinNeighborDistance(Mesh mesh) {
  Point p = mesh->elements.ptr;
  Point e = p + mesh->nVertex * mesh->nelt;

  int nDim = mesh->nDim;
  int nVertex = mesh->nVertex;

  uint i, j, k, neighbor;
  GenmapScalar d;
  for (i = 0; i < mesh->elements.n; i += nVertex) {
    for (j = 0; j < nVertex; j++) {
      p[i + j].dx = GENMAP_SCALAR_MAX;
      for (k = 0; k < mesh->nNeighbors; k++) {
        neighbor = NEIGHBOR_MAP[j][k];
        if (nDim == 3)
          d = distance3D(p[i + j], p[i + neighbor]);
        else
          d = distance2D(p[i + j], p[i + neighbor]);
        p[i + j].dx = min(p[i + j].dx, d);
      }
    }
  }

  return 0;
}
