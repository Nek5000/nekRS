#include "con-impl.h"

typedef struct {
  ulong sequenceId;
  int nnbrs;
  ulong elementId;
  ulong vertexId;
  uint workProc;
} vertex;

typedef struct {
  uint *offsets;
  ulong *elements;
  ulong *globalIds;
  uint size;
} VToEMap;

typedef struct {
  ulong id;
} LongID;

typedef struct {
  ulong procId;
} ProcID;

static VToEMap *getVToEMap(Mesh m, struct comm *c, buffer *bfr) {
  sint nelt = m->nelt;
  sint nv = m->nv;

  slong out[2][1], buf[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  ulong elemId = out[0][0];
  ulong sequenceId = elemId * nv;

  size_t size = nelt * (size_t)nv;
  struct array vertices;
  array_init(vertex, &vertices, size);

  // Create (globalId, elementId) pairs and send them to globalId % np
  Point ptr = m->elements.ptr;
  sint i, j;
  for (i = 0; i < nelt; i++) {
    for (j = 0; j < nv; j++) {
      ulong globalId = ptr[i * nv + j].globalId + 1;
      vertex t = {.elementId = elemId,
                  .sequenceId = sequenceId,
                  .vertexId = globalId,
                  .workProc = globalId % c->np};
      array_cat(vertex, &vertices, &t, 1);
      sequenceId++;
    }

    elemId++;
  }

  sarray_sort_2(vertex, vertices.ptr, vertices.n, vertexId, 1, elementId, 1,
                bfr);

  struct array vtcsCmpct;
  array_init(vertex, &vtcsCmpct, 10);
  vertex *vPtr = (vertex *)vertices.ptr;

  if (vertices.n > 0) {
    vertex prev = vPtr[0];
    array_cat(vertex, &vtcsCmpct, &prev, 1);

    for (i = 1; i < vertices.n; i++) {
      if ((vPtr[i].elementId != prev.elementId) ||
          (vPtr[i].vertexId != prev.vertexId)) {
        prev = vPtr[i];
        array_cat(vertex, &vtcsCmpct, &prev, 1);
      }
    }
  }
  array_free(&vertices);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(vertex, &vtcsCmpct, workProc, 1, &cr);

  // Find all the elements which share globalId and send the union
  // back to all the processors which has globalId
  // FIXME: Assumes quads or hexes
  vPtr = vtcsCmpct.ptr;
  sarray_sort_2(vertex, vPtr, vtcsCmpct.n, vertexId, 1, workProc, 0, bfr);

  struct array a;
  array_init(vertex, &a, 10);
  struct array procs;
  array_init(ProcID, &procs, 10);

  vPtr = vtcsCmpct.ptr;
  sint s = 0, e;
  vertex t;
  ProcID p;
  while (s < vtcsCmpct.n) {
    procs.n = 0;

    p.procId = vPtr[s].workProc;
    array_cat(ProcID, &procs, &p, 1);
    for (e = s + 1; e < vtcsCmpct.n && vPtr[s].vertexId == vPtr[e].vertexId;
         e++) {
      if (vPtr[e].workProc != p.procId) {
        p.procId = vPtr[e].workProc;
        array_cat(ProcID, &procs, &p, 1);
      }
    }

    ProcID *pPtr = procs.ptr;
    e = MIN(e, vtcsCmpct.n);
    for (i = 0; i < procs.n; i++) {
      t.workProc = pPtr[i].procId;
      for (j = s; j < e; j++) {
        t.vertexId = vPtr[j].vertexId;
        t.sequenceId = vPtr[j].sequenceId;
        t.elementId = vPtr[j].elementId;
        array_cat(vertex, &a, &t, 1);
      }
    }
    s = e;
  }
  array_free(&vtcsCmpct);
  array_free(&procs);

  sarray_transfer(vertex, &a, workProc, 1, &cr);
  sarray_sort_2(vertex, a.ptr, a.n, vertexId, 1, elementId, 1, bfr);
  crystal_free(&cr);

  // create the map
  if (a.n == 0)
    return NULL;

  VToEMap *map = calloc(1, sizeof(VToEMap));
  map->elements = calloc(a.n, sizeof(ulong));

  uint nGIds = 1, prev = 0;
  vertex *aPtr = (vertex *)a.ptr;
  for (i = 1; i < a.n; i++) {
    if (aPtr[i].vertexId != aPtr[prev].vertexId)
      nGIds++;
    prev = i;
  }

  map->size = nGIds;
  map->globalIds = calloc(nGIds, sizeof(ulong));
  map->offsets = calloc(nGIds + 1, sizeof(ulong));

  map->elements[0] = aPtr[0].elementId;
  map->globalIds[0] = aPtr[0].vertexId;
  map->offsets[0] = 0;

  prev = 0;
  uint nOffsets = 0;
  for (i = 1; i < a.n; i++) {
    if (aPtr[i].vertexId != aPtr[prev].vertexId) {
      nOffsets++;
      map->globalIds[nOffsets] = aPtr[i].vertexId;
      map->offsets[nOffsets] = prev = i;
    }
    map->elements[i] = aPtr[i].elementId;
  }
  map->offsets[++nOffsets] = a.n;
  assert(nOffsets == nGIds);

  array_free(&a);

  return map;
}

// key must be present in globalIds
static int getPosition(VToEMap *map, ulong key) {
  ulong *globalIds = map->globalIds;

  int begin = 0;
  int end = map->size;
  int mid = 0;
  while (begin < end) {
    mid = (begin + end) / 2;

    if (key == globalIds[mid])
      return mid;
    else if (key < globalIds[mid])
      end = mid;
    else
      begin = mid;
  };

  if (globalIds[mid] != key)
    return -1;
  return mid;
}

static void freeVToEMap(VToEMap *map) {
  free(map->globalIds);
  free(map->offsets);
  free(map->elements);
  free(map);
}

int faceCheck(Mesh mesh, struct comm *c, buffer *bfr) {
  VToEMap *map = getVToEMap(mesh, c, bfr);

  sint nelt = mesh->nelt;
  sint ndim = mesh->ndim;

  int faces[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
  if (ndim == 3)
    memcpy(faces, faces3D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));
  else
    memcpy(faces, faces2D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));

  Point ptr = mesh->elements.ptr;
  int nf = (ndim == 3) ? 6 : 4;
  int nfv = (ndim == 3) ? 4 : 2;
  int nv = (ndim == 3) ? 8 : 4;

  struct array shared;
  array_init(LongID, &shared, 200);

  int err = 0;

  int i, j, k, l;
  for (i = 0; i < nelt && err == 0; i++) {
    for (j = 0; j < nf && err == 0; j++) {
      shared.n = 0;

      for (k = 0; k < nfv; k++) {
        ulong globalId = ptr[i * nv + faces[j][k] - 1].globalId + 1;
        int indx = getPosition(map, globalId);
        assert(indx >= 0);
        LongID elemId;
        for (l = map->offsets[indx]; l < map->offsets[indx + 1]; l++) {
          elemId.id = map->elements[l];
          array_cat(LongID, &shared, &elemId, 1);
        }
      }

      sarray_sort(LongID, shared.ptr, shared.n, id, 1, bfr);

      ulong prev = 0;
      int ncount = 1;
      LongID *sptr = shared.ptr;
      for (l = 1; l < shared.n; l++) {
        if (sptr[l].id != sptr[prev].id) {
          if (ncount == 3) {
            err = 1;
            break;
          }
          prev = l;
          ncount = 1;
        } else
          ncount++;
      }

      if (ncount == 3) {
        err = 1;
        break;
      }
    }
  }

  array_free(&shared);
  freeVToEMap(map);

  return err;
}

int elementCheck(Mesh mesh, struct comm *c, buffer *bfr) {
  uint nelt = mesh->nelt;
  uint ndim = mesh->ndim;
  int nv = (ndim == 3) ? 8 : 4;

  LongID globalIds[8];
  Point ptr = mesh->elements.ptr;
  uint i, j;
  int err = 0;
  for (i = 0; i < nelt && err == 0; i++) {
    for (j = 0; j < nv; j++)
      globalIds[j].id = ptr[i * nv + j].globalId + 1;

    sarray_sort(LongID, globalIds, nv, id, 1, bfr);

    for (j = 0; j < nv - 1; j++) {
      if (globalIds[j].id == globalIds[j + 1].id) {
        err = 1;
        break;
      }
    }
  }

  return err;
}
