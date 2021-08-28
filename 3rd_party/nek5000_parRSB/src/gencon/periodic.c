#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gencon-impl.h>
#include <genmap-impl.h>
#include <sort.h>

int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{1, 5, 7, 3}, {2, 4, 8, 6},
                                                   {1, 2, 6, 5}, {3, 7, 8, 4},
                                                   {1, 3, 4, 2}, {5, 6, 8, 7}};

int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{3, 1, 0, 0}, {2, 4, 0, 0},
                                                   {1, 2, 0, 0}, {4, 3, 0, 0},
                                                   {0, 0, 0, 0}, {0, 0, 0, 0}};

struct minPair_private {
  uint proc;
  ulong orig;
  ulong min;
};
typedef struct minPair_private *minPair;

int compressPeriodicVertices(Mesh mesh, struct comm *c, buffer *bfr) {
  parallel_sort(struct Point_private, &mesh->elements, globalId, gs_long, 0, 0,
                c, bfr);

  Point points = mesh->elements.ptr;
  uint npoints = mesh->elements.n;

  sint i, nunique = 0;
  if (npoints > 0) {
    slong current = points[0].globalId;
    points[0].globalId = nunique;
    for (i = 1; i < npoints; i++)
      if (points[i].globalId == current)
        points[i].globalId = nunique;
      else {
        current = points[i].globalId, ++nunique;
        points[i].globalId = nunique;
      }
  }

  slong out[2][1], buf[2][1], in[1];
  if (npoints > 0)
    in[0] = nunique + 1;
  else
    in[0] = 0;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];

  for (i = 0; i < npoints; i++)
    points[i].globalId += start;

  return 0;
}

ulong findMinBelowI(ulong min, uint I, struct array *arr) {
  minPair ptr = arr->ptr;

  uint i;
  for (i = 0; i < I; i++)
    if (ptr[i].orig == min)
      return ptr[i].min;
  return min;
}

int renumberPeriodicVertices(Mesh mesh, struct comm *c, struct array *matched,
                             buffer *buf) {
  minPair ptr = matched->ptr;
  uint size = matched->n;

  slong *ids;
  GenmapMalloc(size, &ids);
  sint i;
  for (i = 0; i < size; i++)
    ids[i] = ptr[i].orig;

  struct gs_data *t = gs_setup(ids, size, c, 0, gs_pairwise, 0);

  for (i = 0; i < size; i++)
    ids[i] = ptr[i].min;

  gs(ids, gs_long, gs_min, 0, t, buf);

  for (i = 0; i < size; i++)
    ptr[i].min = ids[i];

  GenmapFree(ids);
  gs_free(t);

  sarray_sort_2(struct minPair_private, ptr, size, orig, 1, min, 1, buf);

  struct array compressed;
  array_init(struct minPair_private, &compressed, 10);
  compressed.n = 0;
  if (size > 0)
    array_cat(struct minPair_private, &compressed, ptr, 1);

  for (i = 1; i < size; i++)
    if (ptr[i].orig != ptr[i - 1].orig)
      array_cat(struct minPair_private, &compressed, &ptr[i], 1);

  ptr = compressed.ptr;
  size = compressed.n;
  for (i = 0; i < size; i++)
    ptr[i].proc = 0;

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct minPair_private, &compressed, proc, 1, &cr);
  crystal_free(&cr);

  sint rank = c->id;
  ptr = compressed.ptr;
  size = compressed.n;
  if (rank == 0) {
    sarray_sort_2(struct minPair_private, ptr, size, orig, 1, min, 1, buf);
    for (i = 0; i < size; i++)
      ptr[i].min = findMinBelowI(ptr[i].min, i, &compressed);
  }

  uint sizec = 0;
  if (rank == 0)
    sizec = size;
  size = mesh->elements.n;

  GenmapCalloc(size + sizec, &ids);

  Point pnt = mesh->elements.ptr;
  for (i = 0; i < size; i++)
    ids[i] = pnt[i].globalId;
  for (i = 0; i < sizec; i++)
    ids[size + i] = ptr[i].orig;

  t = gs_setup(ids, size + sizec, c, 0, gs_pairwise, 0);

  for (i = 0; i < size; i++)
    ids[i] = pnt[i].globalId;
  for (i = 0; i < sizec; i++)
    ids[size + i] = ptr[i].min;

  gs(ids, gs_long, gs_min, 0, t, buf);

  for (i = 0; i < size; i++)
    pnt[i].globalId = ids[i];

  gs_free(t);
  GenmapFree(ids);
  array_free(&compressed);
}

int findConnectedPeriodicPairs(Mesh mesh, BoundaryFace f_, BoundaryFace g_,
                               struct array *matched) {
  struct Boundary_private f = *f_, g = *g_;

  int nvf = mesh->nVertex / 2;
  int nDim = mesh->nDim;

  int i, j;
  GenmapScalar fMax = 0.0, gMax = 0.0;

  for (i = 0; i < nDim; i++) {
    GenmapScalar meanF = 0.0, meanG = 0.0;

    for (j = 0; j < nvf; j++) {
      fMax = max(fMax, fabs(f.face.vertex[j].x[i]));
      gMax = max(gMax, fabs(g.face.vertex[j].x[i]));
      meanF += f.face.vertex[j].x[i];
      meanG += g.face.vertex[j].x[i];
    }

    for (j = 0; j < nvf; j++) {
      f.face.vertex[j].x[i] -= (meanF / nvf);
      g.face.vertex[j].x[i] -= (meanG / nvf);
    }
  }

  int shift = 0, k;
  GenmapScalar d2Min = 1.e20, d2;
  for (i = 0; i < nvf; i++) {
    d2 = 0.0;
    for (j = 0; j < nvf; j++) {
      k = (j + i) % nvf;
      k = nvf - 1 - k;
      if (nDim == 3)
        d2 += distance3D(f.face.vertex[j], g.face.vertex[k]);
      else if (nDim == 2)
        d2 += distance2D(f.face.vertex[j], g.face.vertex[k]);
    }
    if (d2 < d2Min) {
      d2Min = d2;
      shift = i;
    }
  }
  d2Min = sqrt(d2Min);

  GenmapScalar fgMax = max(fMax, gMax);
  GenmapScalar tol = (1e-3) * fgMax;
  if (d2Min > tol) {
    fprintf(stderr,
            "Faces did not match: (d2Min,tol,face1,face2): "
            "%lf %lf %lld %lld\n",
            d2Min, tol, f.faceId, g.faceId);
    exit(1);
  }
#if defined(GENMAP_DEBUG)
  printf("Periodic face match (elementId,faceId): (%lld %lld) "
         "and (%lld %lld)\n",
         f.elementId, f.faceId, g.elementId, g.faceId);
#endif

  struct minPair_private m;
  for (i = 0; i < nvf; i++) {
    k = (i + shift) % nvf;
    k = nvf - 1 - k;
    m.min = min(f.face.vertex[i].globalId, g.face.vertex[k].globalId);
    m.orig = max(f.face.vertex[i].globalId, g.face.vertex[k].globalId);
    array_cat(struct minPair_private, matched, &m, 1);
  }
}

int findConnectedPeriodicFaces(Mesh mesh, struct array *matched) {
  sint bSize = mesh->boundary.n;
  BoundaryFace ptr = mesh->boundary.ptr;
  sint i, j;

  for (i = 0; i < bSize - 1; i++)
    for (j = i + 1; j < bSize; j++)
      if (ptr[j].bc[0] == ptr[i].elementId && ptr[j].bc[1] == ptr[i].faceId) {
        findConnectedPeriodicPairs(mesh, &ptr[i], &ptr[j], matched);
      }
}

int gatherMatchingPeriodicFaces(Mesh mesh, struct comm *c) {
  int size = c->np, rank = c->id;

  BoundaryFace bPtr = mesh->boundary.ptr;
  int nFaces = mesh->boundary.n;

  slong nelgt = mesh->nelgt;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  sint i;
  slong eid;
  for (i = 0; i < nFaces; i++) {
    eid = max(bPtr[i].bc[0], bPtr[i].elementId);
#if defined(GENMAP_DEBUG)
    printf("Send matching (%lld,%lld) to (%ld,%ld).\n", bPtr[i].elementId,
           bPtr[i].faceId, bPtr[i].bc[0], bPtr[i].bc[1]);
#endif
    if (eid < N)
      bPtr[i].proc = eid / nelt;
    else
      bPtr[i].proc = (eid - N) / (nelt + 1) + size - nrem;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Boundary_private, &mesh->boundary, proc, 1, &cr);
  crystal_free(&cr);
}

int setPeriodicFaceCoordinates(Mesh mesh, struct comm *c, buffer *buf) {
  BoundaryFace bPtr = mesh->boundary.ptr;
  sint bSize = mesh->boundary.n;
  if (bSize == 0)
    return 0;

  Point ePtr = mesh->elements.ptr;
  sint eSize = mesh->elements.n;
  if (eSize == 0)
    return 0;

  /* Need boundary array to be sorted by elementId */
  sarray_sort(struct Boundary_private, bPtr, bSize, elementId, 1, buf);

  /* Need element array to be sorted by sequenceId */
  sarray_sort(struct Point_private, ePtr, eSize, sequenceId, 1, buf);

  int faces[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
  if (mesh->nDim == 3)
    memcpy(faces, faces3D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));
  else
    memcpy(faces, faces2D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));

  sint i = 0, k = 0;
  int nv = mesh->nVertex, nvf = mesh->nVertex / 2, j;
  while (i < bSize) {
    while (k < eSize && ePtr[k].elementId < bPtr[i].elementId)
      k += nv;
    // copy vertices to boundary face
    if (k < eSize && ePtr[k].elementId == bPtr[i].elementId) {
      int faceId = bPtr[i].faceId;
      for (j = 0; j < nvf; j++)
        bPtr[i].face.vertex[j] = ePtr[k + faces[faceId][j] - 1];

#if defined(GENMAP_DEBUG)
      printf("Periodic BC (element,face):(%ld,%ld)\n", bPtr[i].bc[0],
             bPtr[i].bc[1]);
#endif
    }
    i++;
  }
}

int matchPeriodicFaces(Mesh mesh, struct comm *c, buffer *bfr) {
  setPeriodicFaceCoordinates(mesh, c, bfr);
  gatherMatchingPeriodicFaces(mesh, c);

  struct array matched;
  array_init(struct minPair_private, &matched, 10);
  matched.n = 0;

  findConnectedPeriodicFaces(mesh, &matched);
  renumberPeriodicVertices(mesh, c, &matched, bfr);
  array_free(&matched);

  compressPeriodicVertices(mesh, c, bfr);
  sendBack(mesh, c, bfr);

  return 0;
}
