#include "con-impl.h"
#include "parrsb-impl.h"
#include "sort.h"
#include <stdarg.h>

int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};
int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

void debug_print(struct comm *c, int verbose, const char *fmt, ...) {
  comm_barrier(c);
  va_list vargs;
  va_start(vargs, fmt);
  if (c->id == 0 && verbose > 0) {
    vprintf(fmt, vargs);
    fflush(stdout);
  }
  va_end(vargs);
}

double diff_sqr(double x, double y) { return (x - y) * (x - y); }

//==============================================================================
// Mesh struct
//
static struct mesh_t *mesh_init(int nelt, int ndim, double *coord,
                                long long *pinfo, int npinfo,
                                const struct comm *c) {
  struct mesh_t *m = tcalloc(struct mesh_t, 1);
  m->nelt = nelt, m->ndim = ndim, m->nnbrs = ndim;
  m->nv = (ndim == 2) ? 4 : 8;

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong start = out[0][0];
  m->nelgt = out[1][0];

  int nv = m->nv;
  array_init(struct point_t, &m->elements, nelt * nv);
  struct point_t p = {.origin = c->id};
  for (uint i = 0; i < nelt; i++) {
    for (uint k = 0; k < nv; k++) {
      uint j = PRE_TO_SYM_VERTEX[k];
      for (uint l = 0; l < ndim; l++)
        p.x[l] = coord[i * nv * ndim + j * ndim + l];
      p.elementId = start + i, p.sequenceId = nv * (start + i) + k;
      array_cat(struct point_t, &m->elements, &p, 1);
    }
  }

  array_init(struct boundary_t, &m->boundary, npinfo);
  struct boundary_t b;
  for (uint i = 0; i < npinfo; i++) {
    b.elementId = pinfo[4 * i + 0] - 1;
    b.faceId = PRE_TO_SYM_FACE[pinfo[4 * i + 1] - 1];
    b.bc[0] = pinfo[4 * i + 2] - 1;
    b.bc[1] = PRE_TO_SYM_FACE[pinfo[4 * i + 3] - 1];
    array_cat(struct boundary_t, &m->boundary, &b, 1);
  }

  return m;
}

static int mesh_free(struct mesh_t *m) {
  array_free(&m->elements), array_free(&m->boundary), free(m);
  return 0;
}

//==============================================================================
// Find the minimum distance between a vertex and its neighbors
//
static inline double distance_2d(struct point_t *a, struct point_t *b) {
  return diff_sqr(a->x[0], b->x[0]) + diff_sqr(a->x[1], b->x[1]);
}

static inline double distance_3d(struct point_t *a, struct point_t *b) {
  return distance_2d(a, b) + diff_sqr(a->x[2], b->x[2]);
}

int findMinNeighborDistance(Mesh mesh) {
  struct point_t *p = (struct point_t *)mesh->elements.ptr;
  int ndim = mesh->ndim;
  int nv = mesh->nv;

  uint i, j, k;
  int neighbor;
  scalar d;

  if (ndim == 3) {
    for (i = 0; i < mesh->elements.n; i += nv) {
      for (j = 0; j < nv; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nnbrs; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          d = distance_3d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  } else if (ndim == 2) {
    for (i = 0; i < mesh->elements.n; i += nv) {
      for (j = 0; j < nv; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nnbrs; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          d = distance_2d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  } else {
    return 1;
  }

  return 0;
}

//==============================================================================
// Global numbering
//
static int setGlobalID(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = (struct point_t *)mesh->elements.ptr;

  sint bin = (nPoints > 0);
  struct comm nonZeroRanks;
  comm_split(c, bin, c->id, &nonZeroRanks);

  sint rank = nonZeroRanks.id;
  sint size = nonZeroRanks.np;

  if (bin == 1) {
    slong count = 0;
    for (uint i = 0; i < nPoints; i++)
      if (points[i].ifSegment)
        count++;

    slong in = count, out[2][1], buf[2][1];
    comm_scan(out, &nonZeroRanks, gs_long, gs_add, &in, 1, buf);
    slong start = out[0][0];

    count = -1;
    for (uint i = 0; i < nPoints; i++) {
      if (points[i].ifSegment)
        count++;
      assert(start + count >= 0);
      points[i].globalId = start + count;
    }
  }

  comm_free(&nonZeroRanks);

  return 0;
}

int send_back(Mesh mesh, struct comm *c, buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct point_t, &mesh->elements, origin, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct point_t, mesh->elements.ptr, mesh->elements.n, sequenceId,
              1, bfr);

  return 0;
}

static int transferBoundaryFaces(Mesh mesh, struct comm *c) {
  uint size = c->np;

  struct array *boundary = &mesh->boundary;
  BoundaryFace ptr = (struct boundary_t *)boundary->ptr;
  int nFaces = boundary->n;

  slong nelgt = mesh->nelgt;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  sint i;
  slong eid;
  for (i = 0; i < nFaces; i++) {
    eid = ptr[i].elementId;
    if (eid < N)
      ptr[i].proc = eid / nelt;
    else
      ptr[i].proc = (eid - N) / (nelt + 1) + size - nrem;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct boundary_t, boundary, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

//==============================================================================
// C interface to find_conn
//
#define check_error(call, msg)                                                 \
  {                                                                            \
    sint err = (call);                                                         \
    sint buf;                                                                  \
    comm_allreduce(&c, gs_int, gs_max, &err, 1, &buf);                         \
    if (err) {                                                                 \
      buffer_free(&bfr), mesh_free(mesh), comm_free(&c);                       \
      return err;                                                              \
    }                                                                          \
  }

// Input:
//   nelt: Number of elements, nv: Number of vertices in an element
//   coord [nelt, nv, ndim]: Coordinates of elements vertices in preprocessor
//     ordering, nv = 8 if ndim == 3 (Hex) or nv = 4 if ndim = 2 (Quad).
// Output:
//   vtx[nelt, nv]: Global numbering of vertices of elements
int parrsb_conn_mesh(long long *vtx, double *coord, int nelt, int ndim,
                     long long *pinfo, int npinfo, double tol, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  buffer bfr;
  buffer_init(&bfr, 1024);

  int verbose = 0;
  {
    const char *val = getenv("PARRSB_VERBOSE_LEVEL");
    if (val != NULL)
      verbose = atoi(val);
  }

  debug_print(&c, verbose, "Running parCon ...\n");

  parrsb_barrier(&c);
  double tall = comm_time(), t;

  double duration[8] = {0};
  const char *name[8] = {"transferBoundaryFaces", "findMinNbrDistance   ",
                         "find_unique_vertices ", "setGlobalId          ",
                         "elementCheck         ", "faceCheck            ",
                         "matchPeriodicFaces   ", "copyOutput           "};

  // debug_print(&c, verbose, "\t%s ...");
  // parrsb_barrier(&c), t = comm_time();
  Mesh mesh = mesh_init(nelt, ndim, coord, pinfo, npinfo, &c);
  // duration[0] = comm_time() - t;
  // debug_print(&c, verbose, "done.\n");

  debug_print(&c, verbose, "\t%s ...", name[0]);
  parrsb_barrier(&c), t = comm_time();
  check_error(transferBoundaryFaces(mesh, &c), name[0]);
  duration[0] = comm_time() - t;
  debug_print(&c, verbose, "done.\n");

  debug_print(&c, verbose, "\t%s ...", name[1]);
  parrsb_barrier(&c), t = comm_time();
  check_error(findMinNeighborDistance(mesh), name[1]);
  duration[1] = comm_time() - t;
  debug_print(&c, verbose, "done.\n");

  debug_print(&c, verbose, "\t%s ...\n", name[2]);
  parrsb_barrier(&c), t = comm_time();
  check_error(find_unique_vertices(mesh, &c, tol, verbose, &bfr), name[2]);
  duration[2] = comm_time() - t;

  debug_print(&c, verbose, "\t%s ...", name[3]);
  parrsb_barrier(&c), t = comm_time();
  setGlobalID(mesh, &c);
  send_back(mesh, &c, &bfr);
  duration[3] = comm_time() - t;
  debug_print(&c, verbose, "done.\n");

  debug_print(&c, verbose, "\t%s ...", name[4]);
  parrsb_barrier(&c), t = comm_time();
  check_error(elementCheck(mesh, &c, &bfr), name[4]);
  duration[4] = comm_time() - t;
  debug_print(&c, verbose, "done.\n");

  debug_print(&c, verbose, "\t%s ...", name[5]);
  parrsb_barrier(&c), t = comm_time();
  check_error(faceCheck(mesh, &c, &bfr), name[5]);
  duration[5] = comm_time() - t;
  debug_print(&c, verbose, "done.\n");

  debug_print(&c, verbose, "\t%s ...", name[6]);
  parrsb_barrier(&c), t = comm_time();
  check_error(matchPeriodicFaces(mesh, &c, &bfr), name[6]);
  duration[6] = comm_time() - t;
  debug_print(&c, verbose, "done.\n");

  debug_print(&c, verbose, "\t%s ...", name[7]);
  parrsb_barrier(&c), t = comm_time();
  Point ptr = mesh->elements.ptr;
  for (uint i = 0; i < nelt; i++) {
    for (uint j = 0; j < mesh->nv; j++)
      vtx[i * mesh->nv + j] = ptr[i * mesh->nv + j].globalId + 1;
  }
  duration[7] = comm_time() - t;
  debug_print(&c, verbose, "done.\n");

  // Report timing info and finish
  double gmin[8], gmax[8], buf[8];
  for (unsigned i = 0; i < 8; i++)
    gmax[i] = gmin[i] = duration[i];
  comm_allreduce(&c, gs_double, gs_min, gmin, 8, buf);
  comm_allreduce(&c, gs_double, gs_max, gmax, 8, buf);

  if (c.id == 0 && verbose > 1) {
    for (unsigned i = 0; i < 7; i++)
      printf("%s: %e %e (min max)\n", name[i], gmin[i], gmax[i]);
    fflush(stdout);
  }

  parrsb_barrier(&c), tall = comm_time() - tall;
  if (c.id == 0) {
    printf("parCon (tol = %e) finished in %g s\n", tol, tall);
    fflush(stdout);
  }

  buffer_free(&bfr), mesh_free(mesh), comm_free(&c);

  return 0;
}
#undef check_error

//=============================================================================
// Fortran interface
//
void fparrsb_conn_mesh(long long *vtx, double *coord, int *nelt, int *ndim,
                       long long *pinfo, int *npinfo, double *tol,
                       MPI_Fint *fcomm, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*fcomm);
  *err = parrsb_conn_mesh(vtx, coord, *nelt, *ndim, pinfo, *npinfo, *tol, c);
}
