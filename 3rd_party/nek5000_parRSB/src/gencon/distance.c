#include <gencon-impl.h>

int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

double diff_sqr(double x, double y) { return (x - y) * (x - y); }

double distance_2d(struct Point_private *a, struct Point_private *b) {
  return diff_sqr(a->x[0], b->x[0]) + diff_sqr(a->x[1], b->x[1]);
}

double distance_3d(struct Point_private *a, struct Point_private *b) {
  return distance_2d(a, b) + diff_sqr(a->x[2], b->x[2]);
}

int findMinNeighborDistance(Mesh mesh) {
  Point p = mesh->elements.ptr;
  Point e = p + mesh->nVertex * mesh->nelt;

  int nDim = mesh->nDim;
  int nVertex = mesh->nVertex;

  uint i, j, k;
  int neighbor;
  GenmapScalar d;

  if (nDim == 3) {
    for (i = 0; i < mesh->elements.n; i += nVertex) {
      for (j = 0; j < nVertex; j++) {
        p[i + j].dx = GENMAP_SCALAR_MAX;
        for (k = 0; k < mesh->nNeighbors; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          d = distance_3d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = min(p[i + j].dx, d);
        }
      }
    }
  } else if (nDim == 2) {
    for (i = 0; i < mesh->elements.n; i += nVertex) {
      for (j = 0; j < nVertex; j++) {
        p[i + j].dx = GENMAP_SCALAR_MAX;
        for (k = 0; k < mesh->nNeighbors; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          d = distance_2d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = min(p[i + j].dx, d);
        }
      }
    }
  } else {
    return 1;
  }

  return 0;
}
